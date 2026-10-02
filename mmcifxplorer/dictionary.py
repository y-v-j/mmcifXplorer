"""The wwPDB PDBx/mmCIF dictionary (DDL2) as a queryable Python object.

The machine-readable dictionary ``mmcif_pdbx_v50.dic`` is parsed once with
gemmi and condensed into a JSON cache next to it.  Every explanation shown in
the dashboard for a data item comes from here: the definition text, data type
(with its regular expression), units, enumerations, permitted range, examples,
aliases to the core CIF dictionary and the parent/child links that tie
categories together.
"""
from __future__ import annotations

import gzip
import json
import re
import textwrap
import warnings
from functools import lru_cache
from pathlib import Path

import gemmi
import requests

from . import config

CACHE_SCHEMA = 4


# --------------------------------------------------------------------------- #
# Download / parse
# --------------------------------------------------------------------------- #
def ensure_dictionary_file(path: Path = config.DICT_FILE, url: str = config.DICT_URL) -> Path:
    if path.exists() and path.stat().st_size > 0:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    r = requests.get(url, timeout=config.HTTP_TIMEOUT)
    r.raise_for_status()
    data = r.content if url.endswith(".gz") else gzip.compress(r.content)
    tmp = path.with_suffix(".part")
    tmp.write_bytes(data)
    tmp.replace(path)
    return path


def merge_extension(base: dict, ext: dict, source: str) -> dict:
    """Add categories/items/groups/types/units of ``ext`` that ``base`` lacks."""
    for key in ("types", "units", "groups"):
        for k, v in ext[key].items():
            base[key].setdefault(k, v)
    for k, v in ext["categories"].items():
        if k not in base["categories"]:
            base["categories"][k] = dict(v, source=source)
            for g in v["groups"]:
                base["groups"].setdefault(g, {"id": g, "parent": "", "description": "", "categories": []})
                if v["id"] not in base["groups"][g]["categories"]:
                    base["groups"][g]["categories"].append(v["id"])
    for k, v in ext["items"].items():
        if k not in base["items"]:
            base["items"][k] = dict(v, source=source)
    base["extensions"] = base.get("extensions", []) + [f"{source} v{ext['version']}"]
    return base


def _s(value: str | None) -> str:
    """Unquote a raw CIF value; '?' and '.' become ''."""
    if value is None or value in ("?", "."):
        return ""
    return gemmi.cif.as_string(value)


def _text(value: str | None) -> str:
    """Clean a dictionary text field: unquote, dedent, trim blank edges."""
    s = _s(value)
    if not s:
        return ""
    lines = s.split("\n")
    # The first line of a text field is usually empty or flush-left while the
    # rest is indented by ~14 spaces; dedent the body separately.
    head, body = lines[0], "\n".join(lines[1:])
    body = textwrap.dedent(body)
    out = (head.strip() + "\n" + body) if head.strip() else body
    return out.strip("\n").rstrip()


def _loop_rows(frame, prefix: str, fields: list[str]) -> list[list[str]]:
    tab = frame.find(prefix, fields)
    rows = []
    for row in tab:
        rows.append([row[i] if row.has(i) else "?" for i in range(len(fields))])
    return rows


def _opt_rows(frame, prefix: str, required: list[str], optional: list[str]):
    """Like _loop_rows but tolerates missing optional columns."""
    tab = frame.find(prefix, required + ["?" + o for o in optional])
    rows = []
    for row in tab:
        vals = []
        for i in range(len(required) + len(optional)):
            vals.append(row[i] if row.has(i) else "?")
        rows.append(vals)
    return rows


def posix_to_python_regex(construct: str) -> str:
    """Translate a DDL2 (POSIX ERE) construct to Python ``re`` syntax.

    Inside POSIX bracket expressions a backslash is a literal character, while
    Python treats it as an escape, so backslashes inside brackets are doubled.
    """
    out, i, n = [], 0, len(construct)
    in_bracket = False
    while i < n:
        c = construct[i]
        if not in_bracket:
            if c == "[":
                in_bracket = True
                out.append(c)
                i += 1
                # A ']' (optionally after '^') right after '[' is a literal.
                if i < n and construct[i] == "^":
                    out.append("^")
                    i += 1
                if i < n and construct[i] == "]":
                    out.append(r"\]")
                    i += 1
                continue
            out.append(c)
        else:
            if c == "\\":
                out.append(r"\\")
            elif c == "[":
                out.append(r"\[")
            elif c == "]":
                in_bracket = False
                out.append(c)
            else:
                out.append(c)
        i += 1
    return "".join(out)


def build_cache(dic_path: Path) -> dict:
    doc = gemmi.cif.read(str(dic_path))
    block = doc[0]

    version = _s(block.find_value("_dictionary.version"))

    groups: dict[str, dict] = {}
    for row in block.find("_category_group_list.", ["id", "?parent_id", "?description"]):
        gid = _s(row[0])
        groups[gid] = {
            "id": gid,
            "parent": _s(row[1]) if row.has(1) else "",
            "description": _text(row[2]) if row.has(2) else "",
            "categories": [],
        }

    types: dict[str, dict] = {}
    for row in block.find("_item_type_list.", ["code", "primitive_code", "construct", "?detail"]):
        types[_s(row[0])] = {
            "primitive": _s(row[1]),
            "construct": _s(row[2]).strip("\n"),
            "detail": _text(row[3]) if row.has(3) else "",
        }

    units: dict[str, str] = {}
    for row in block.find("_item_units_list.", ["code", "detail"]):
        units[_s(row[0])] = _text(row[1])

    categories: dict[str, dict] = {}
    items: dict[str, dict] = {}
    children: dict[str, set] = {}
    parents: dict[str, set] = {}

    def link(child: str, parent: str):
        if not child or not parent:
            return
        parents.setdefault(child.lower(), set()).add(parent)
        children.setdefault(parent.lower(), set()).add(child)

    # Linked groups declared at block level (pdbx extension).
    for row in block.find("_pdbx_item_linked_group_list.", ["child_name", "parent_name"]):
        link(_s(row[0]), _s(row[1]))

    deferred_names: list[tuple[str, str, str, str]] = []  # (name, category, mandatory, defining frame)

    for it in block:
        frame = it.frame
        if frame is None:
            continue
        fname = frame.name
        if not fname.startswith("_"):
            # ---- category frame ----
            cid = _s(frame.find_value("_category.id")) or fname
            keys = [_s(v) for v in frame.find_values("_category_key.name")]
            grp = [_s(v) for v in frame.find_values("_category_group.id")]
            examples = []
            for case, detail in _opt_rows(frame, "_category_examples.", ["case"], ["detail"]):
                examples.append({"case": textwrap.dedent(_s(case)).strip("\n"), "detail": _text(detail)})
            ctx = _s(frame.find_value("_pdbx_category_context.type"))
            categories[cid.lower()] = {
                "id": cid,
                "description": _text(frame.find_value("_category.description")),
                "mandatory": _s(frame.find_value("_category.mandatory_code")),
                "keys": keys,
                "groups": grp,
                "examples": examples,
                "context": ctx,
                "items": [],
            }
            for g in grp:
                groups.setdefault(g, {"id": g, "parent": "", "description": "", "categories": []})
                groups[g]["categories"].append(cid)
            continue

        # ---- item frame ----
        names = []
        for row in frame.find("_item.", ["name", "?category_id", "?mandatory_code"]):
            names.append((
                _s(row[0]),
                _s(row[1]) if row.has(1) else "",
                _s(row[2]) if row.has(2) else "",
            ))
        if not names:
            names = [(fname, fname[1:].split(".")[0], "")]
        own = next((n for n in names if n[0].lower() == fname.lower()), names[0])
        name, cat, mandatory = own
        if not cat:
            cat = name[1:].split(".")[0]

        enum = [{"value": _s(v), "detail": _text(d)}
                for v, d in _opt_rows(frame, "_item_enumeration.", ["value"], ["detail"])]
        if not enum:
            enum = [{"value": _s(v), "detail": _text(d)}
                    for v, d in _opt_rows(frame, "_pdbx_item_enumeration.", ["value"], ["detail"])]
        ranges = [(_s(mn), _s(mx)) for mn, mx in _loop_rows(frame, "_item_range.", ["minimum", "maximum"])]
        pranges = [(_s(mn), _s(mx)) for mn, mx in _opt_rows(frame, "_pdbx_item_range.", ["minimum", "maximum"], [])]
        ex = [{"case": _s(c), "detail": _text(d)}
              for c, d in _opt_rows(frame, "_item_examples.", ["case"], ["detail"])]
        if not ex:
            ex = [{"case": _s(c), "detail": _text(d)}
                  for c, d in _opt_rows(frame, "_pdbx_item_examples.", ["case"], ["detail"])]
        aliases = [{"alias": _s(a), "dictionary": _s(d)}
                   for a, d in _opt_rows(frame, "_item_aliases.", ["alias_name"], ["dictionary"])]
        related = [{"name": _s(n), "function": _s(f)}
                   for n, f in _opt_rows(frame, "_item_related.", ["related_name", "function_code"], [])]
        for c, p in _opt_rows(frame, "_item_linked.", ["child_name", "parent_name"], []):
            link(_s(c), _s(p))
        type_code = _s(frame.find_value("_item_type.code"))
        unit = _s(frame.find_value("_item_units.code"))
        items[name.lower()] = {
            "name": name,
            "category": cat,
            "description": _text(frame.find_value("_item_description.description")),
            "pdbx_description": _text(frame.find_value("_pdbx_item_description.description")),
            "mandatory": mandatory or _s(frame.find_value("_pdbx_item.mandatory_code")),
            "type": type_code,
            "units": unit,
            "enum": enum,
            "range": ranges or pranges,
            "examples": ex,
            "aliases": aliases,
            "related": related,
            "default": _s(frame.find_value("_item_default.value")),
            "context": _s(frame.find_value("_pdbx_item_context.type")),
            "sub_category": [_s(v) for v in frame.find_values("_item_sub_category.id")],
            "defined_in": name,
        }
        for n, c, m in names:
            if n.lower() != name.lower():
                deferred_names.append((n, c or n[1:].split(".")[0], m, name))

    # Items that are only listed inside another (parent) frame get a stub that
    # points to the definition they share.
    for n, c, m, origin in deferred_names:
        if n.lower() in items:
            continue
        base = items[origin.lower()]
        items[n.lower()] = dict(base, name=n, category=c, mandatory=m or base["mandatory"],
                                defined_in=origin, enum=base["enum"], related=[], aliases=[])

    for key, item in items.items():
        item["parents"] = sorted(parents.get(key, ()))
        item["children"] = sorted(children.get(key, ()))
        cat = categories.get(item["category"].lower())
        if cat is not None:
            cat["items"].append(item["name"])
    for cat in categories.values():
        cat["items"].sort(key=str.lower)

    return {
        "schema": CACHE_SCHEMA,
        "version": version,
        "groups": groups,
        "types": types,
        "units": units,
        "categories": categories,
        "items": items,
    }


# --------------------------------------------------------------------------- #
# Query API
# --------------------------------------------------------------------------- #
class MmcifDictionary:
    def __init__(self, data: dict):
        self.data = data
        self.version: str = data["version"]
        self.items: dict[str, dict] = data["items"]
        self.categories: dict[str, dict] = data["categories"]
        self.groups: dict[str, dict] = data["groups"]
        self.types: dict[str, dict] = data["types"]
        self.units: dict[str, str] = data["units"]
        self.extensions: list[str] = data.get("extensions", [])
        self._regex: dict[str, re.Pattern | None] = {}

    # -- lookups ---------------------------------------------------------- #
    @staticmethod
    def norm_item(name: str) -> str:
        name = name.strip()
        if not name.startswith("_"):
            name = "_" + name
        return name.lower()

    def item(self, name: str) -> dict | None:
        return self.items.get(self.norm_item(name))

    def category(self, cat: str) -> dict | None:
        return self.categories.get(cat.strip().lstrip("_").lower())

    def category_group(self, cat: str) -> str:
        """The most specific group of a category (not the 'inclusive_group')."""
        c = self.category(cat)
        if not c:
            return "other"
        for g in c["groups"]:
            if g != "inclusive_group":
                return g
        return c["groups"][0] if c["groups"] else "other"

    def type_info(self, code: str) -> dict:
        return self.types.get(code, {})

    def type_regex(self, code: str) -> re.Pattern | None:
        if code in self._regex:
            return self._regex[code]
        construct = self.types.get(code, {}).get("construct", "")
        pat = None
        if construct:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    pat = re.compile(posix_to_python_regex(construct), re.DOTALL)
                except re.error:
                    pat = None
        self._regex[code] = pat
        return pat

    def check_value(self, item_name: str, value: str) -> dict:
        """Check a value against the item's type regex, enumeration and range."""
        out = {"type_ok": None, "enum_ok": None, "range_ok": None}
        it = self.item(item_name)
        if it is None or value in ("?", ".", "", None):
            return out
        pat = self.type_regex(it["type"])
        if pat is not None:
            out["type_ok"] = pat.fullmatch(value) is not None
        if it["enum"]:
            vals = {e["value"] for e in it["enum"]}
            out["enum_ok"] = value in vals or value.lower() in {v.lower() for v in vals}
        if it["range"] and self.type_info(it["type"]).get("primitive") == "numb":
            try:
                x = float(re.sub(r"\(\d+\)$", "", value))
            except ValueError:
                return out
            out["range_ok"] = self.in_range(x, it["range"])
        return out

    @staticmethod
    def in_range(x: float, ranges: list) -> bool:
        """DDL2 range semantics: rows with min == max are inclusive points,
        other rows are open intervals; '' means unbounded."""
        for mn, mx in ranges:
            lo = float(mn) if mn else float("-inf")
            hi = float(mx) if mx else float("inf")
            if mn and mx and lo == hi:
                if x == lo:
                    return True
            elif lo < x < hi:
                return True
        return False

    @staticmethod
    def range_text(ranges: list) -> str:
        if not ranges:
            return ""
        lows = [float(mn) for mn, _ in ranges if mn]
        highs = [float(mx) for _, mx in ranges if mx]
        lo = min(lows) if len(lows) == len(ranges) else None
        hi = max(highs) if len(highs) == len(ranges) else None
        incl_lo = lo is not None and any(mn and mx and float(mn) == float(mx) == lo for mn, mx in ranges)
        incl_hi = hi is not None and any(mn and mx and float(mn) == float(mx) == hi for mn, mx in ranges)
        left = ("[" if incl_lo else "(") + (f"{lo:g}" if lo is not None else "−∞")
        right = (f"{hi:g}" if hi is not None else "∞") + ("]" if incl_hi else ")")
        return f"{left}, {right}"

    def search(self, query: str, group: str | None = None, limit: int = 200,
               include_local: bool = False) -> list[dict]:
        """Rank items by where the query matches (name > category > text)."""
        q = query.strip().lower().lstrip("_")
        res = []
        for key, it in self.items.items():
            if not include_local and "WWPDB_LOCAL" in it["context"] and "EMDB" not in it["context"]:
                continue
            if group and self.category_group(it["category"]) != group:
                continue
            if not q:
                score = 3
            elif key.lstrip("_") == q:
                score = 0
            elif q in key:
                score = 1 if key.lstrip("_").startswith(q) else 2
            elif q in it["description"].lower():
                score = 4
            else:
                continue
            res.append((score, key, it))
        res.sort(key=lambda t: (t[0], t[1]))
        return [it for _, _, it in res[:limit]]


@lru_cache(maxsize=1)
def load_dictionary(dic_path: str | None = None, with_modelcif: bool = True) -> MmcifDictionary:
    path = Path(dic_path) if dic_path else ensure_dictionary_file()
    ext_path = None
    if with_modelcif and dic_path is None:
        try:
            ext_path = ensure_dictionary_file(config.EXT_DICT_FILE, config.EXT_DICT_URL)
        except requests.RequestException:
            ext_path = None          # offline: PDBx only
    cache = path.with_name(path.name.replace(".dic.gz", "").replace(".dic", "") + ".cache.json.gz")
    newest = max(path.stat().st_mtime, ext_path.stat().st_mtime if ext_path else 0)
    if cache.exists() and cache.stat().st_mtime >= newest:
        with gzip.open(cache, "rt", encoding="utf-8") as fh:
            data = json.load(fh)
        if data.get("schema") == CACHE_SCHEMA and bool(data.get("extensions")) == bool(ext_path):
            return MmcifDictionary(data)
    data = build_cache(path)
    if ext_path is not None:
        data = merge_extension(data, build_cache(ext_path), "ModelCIF")
    tmp = cache.with_suffix(".part")
    with gzip.open(tmp, "wt", encoding="utf-8") as fh:
        json.dump(data, fh)
    tmp.replace(cache)
    return MmcifDictionary(data)
