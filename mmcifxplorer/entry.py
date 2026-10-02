"""Fetching, caching and parsing one PDBx/mmCIF entry.

The ``_atom_site`` loop is kept as the single source of truth for atoms: the
3D viewer, the atom inspector and every "show this row in 3D" lookup refer to
row indices of that loop, so a picked atom always maps back to exactly one
line of the file.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import re
import threading
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
import requests

from . import config

PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
EXT_PDB_ID_RE = re.compile(r"^pdb_0000([0-9][A-Za-z0-9]{3})$", re.I)
UNIPROT_RE = re.compile(
    r"^(?:AF-)?([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-F1)?$", re.I)


class EntryError(Exception):
    """A user-facing problem (bad ID, download failed, unparsable file)."""


# --------------------------------------------------------------------------- #
# Identifier handling and download
# --------------------------------------------------------------------------- #
def normalize_id(text: str) -> tuple[str, str]:
    """Return (kind, id) where kind is 'pdb' or 'afdb'."""
    s = (text or "").strip()
    if PDB_ID_RE.match(s):
        return "pdb", s.upper()
    m = EXT_PDB_ID_RE.match(s)
    if m:
        return "pdb", m.group(1).upper()
    m = UNIPROT_RE.match(s.upper())
    if m:
        return "afdb", m.group(1).upper()
    raise EntryError(
        f"'{text}' is not a PDB ID (e.g. 1TUP, pdb_00001tup) or a UniProt accession (e.g. P04637).")


def _download(url: str) -> bytes:
    with requests.get(url, timeout=config.HTTP_TIMEOUT, stream=True) as r:
        if r.status_code == 404:
            raise EntryError(f"Not found: {url}")
        r.raise_for_status()
        buf = io.BytesIO()
        for chunk in r.iter_content(1 << 16):
            buf.write(chunk)
            if buf.tell() > config.MAX_DOWNLOAD_BYTES:
                raise EntryError("The file is larger than the configured download limit.")
        return buf.getvalue()


def fetch_entry_file(kind: str, ident: str) -> tuple[Path, str]:
    """Download (or reuse the cached copy of) an entry; return (path, source)."""
    config.ENTRY_DIR.mkdir(parents=True, exist_ok=True)
    if kind == "pdb":
        path = config.ENTRY_DIR / f"{ident}.cif.gz"
        if path.exists():
            return path, "RCSB PDB (cached)"
        try:
            data = _download(config.RCSB_URL.format(id=ident))
        except (requests.RequestException, EntryError) as exc:
            bundled = config.BUNDLED_ENTRIES.get(ident)
            if bundled and bundled.exists():
                return bundled, "bundled copy (offline)"
            raise EntryError(f"Could not download {ident} from RCSB: {exc}") from exc
        tmp = path.with_suffix(".part")
        tmp.write_bytes(data)
        tmp.replace(path)
        return path, "RCSB PDB"
    # AlphaFold DB model for a UniProt accession
    path = config.ENTRY_DIR / f"AF-{ident}-F1.cif.gz"
    if path.exists():
        return path, "AlphaFold DB (cached)"
    try:
        meta = requests.get(config.ALPHAFOLD_API.format(acc=ident), timeout=config.HTTP_TIMEOUT)
        if meta.status_code == 404:
            raise EntryError(f"No AlphaFold DB model for UniProt {ident}.")
        meta.raise_for_status()
        info = meta.json()
        cif_url = info[0]["cifUrl"]
        data = _download(cif_url)
    except (requests.RequestException, KeyError, IndexError, ValueError) as exc:
        raise EntryError(f"Could not fetch the AlphaFold DB model for {ident}: {exc}") from exc
    with open(path, "wb") as fh:
        fh.write(gzip.compress(data))
    return path, "AlphaFold DB"


# --------------------------------------------------------------------------- #
# Entry
# --------------------------------------------------------------------------- #
ATOM_COLUMNS = {
    # our name: (mmCIF item, required)
    "group": ("group_PDB", False),
    "id": ("id", True),
    "elem": ("type_symbol", True),
    "atom": ("label_atom_id", True),
    "alt": ("label_alt_id", False),
    "comp": ("label_comp_id", True),
    "asym": ("label_asym_id", True),
    "entity": ("label_entity_id", False),
    "seq": ("label_seq_id", False),
    "ins": ("pdbx_PDB_ins_code", False),
    "x": ("Cartn_x", True),
    "y": ("Cartn_y", True),
    "z": ("Cartn_z", True),
    "occ": ("occupancy", False),
    "b": ("B_iso_or_equiv", False),
    "auth_seq": ("auth_seq_id", False),
    "auth_comp": ("auth_comp_id", False),
    "auth_asym": ("auth_asym_id", False),
    "auth_atom": ("auth_atom_id", False),
    "model": ("pdbx_PDB_model_num", False),
}

METALS = {e for e in ("LI NA K RB CS BE MG CA SR BA AL GA V CR MN FE CO NI CU ZN MO W CD HG PT AU AG "
                      "PB TL SN SB BI RU RH PD OS IR LA CE PR ND SM EU GD TB DY HO ER TM YB LU Y ZR NB "
                      "HF TA RE SC TI U").split()}


def _str_values(raw: list[str]) -> list[str]:
    """Unquote raw CIF values but keep '?' and '.' visible (they are meaningful)."""
    out = []
    for v in raw:
        if v in ("?", "."):
            out.append(v)
        else:
            out.append(gemmi.cif.as_string(v))
    return out


@dataclass
class Entry:
    key: str
    entry_id: str
    source: str
    block: gemmi.cif.Block
    doc: gemmi.cif.Document
    atoms: pd.DataFrame = field(repr=False)
    _tables: dict = field(default_factory=dict, repr=False)
    _structure: gemmi.Structure | None = field(default=None, repr=False)

    # -- construction ----------------------------------------------------- #
    @classmethod
    def from_document(cls, doc: gemmi.cif.Document, key: str, source: str) -> "Entry":
        if len(doc) == 0:
            raise EntryError("The file contains no data block.")
        block = doc.sole_block() if len(doc) == 1 else doc[0]
        entry_id = gemmi.cif.as_string(block.find_value("_entry.id") or "") or block.name
        atoms = cls._atom_table(block)
        ent = cls(key=key, entry_id=entry_id, source=source, block=block, doc=doc, atoms=atoms)
        ent._annotate_atoms()
        return ent

    @staticmethod
    def _atom_table(block: gemmi.cif.Block) -> pd.DataFrame:
        cat = block.get_mmcif_category("_atom_site.", raw=True)
        if not cat:
            raise EntryError("This file has no _atom_site loop (no coordinates).")
        n = len(next(iter(cat.values())))
        cols = {}
        for ours, (tag, required) in ATOM_COLUMNS.items():
            if tag in cat:
                cols[ours] = _str_values(cat[tag])
            elif required:
                raise EntryError(f"_atom_site.{tag} is missing; cannot build the model.")
            else:
                cols[ours] = ["?"] * n
        df = pd.DataFrame(cols)
        for c in ("x", "y", "z"):
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
        df["occ"] = pd.to_numeric(df["occ"].replace({"?": "1", ".": "1"}), errors="coerce").fillna(1.0)
        df["b"] = pd.to_numeric(df["b"], errors="coerce").fillna(0.0)
        df["seq"] = pd.to_numeric(df["seq"].replace({"?": "-1", ".": "-1"}), errors="coerce").fillna(-1).astype(int)
        df["model"] = pd.to_numeric(df["model"].replace({"?": "1", ".": "1"}), errors="coerce").fillna(1).astype(int)
        for c in ("alt", "ins"):
            df[c] = df[c].replace({"?": "", ".": ""})
        for src, dst in (("auth_asym", "asym"), ("auth_seq", "seq"), ("auth_atom", "atom"), ("auth_comp", "comp")):
            missing = df[src] == "?"
            if missing.any():
                df.loc[missing, src] = df.loc[missing, dst].astype(str)
        df["elem"] = df["elem"].str.upper()
        return df

    def _annotate_atoms(self):
        """Add molecule kind and file-declared secondary structure per atom."""
        df = self.atoms
        etype = {r.get("id"): r.get("type", "") for r in self.rows("entity")}
        ptype = {r.get("entity_id"): r.get("type", "") for r in self.rows("entity_poly")}
        kind = np.full(len(df), "ligand", dtype=object)
        ent = df["entity"].to_numpy()
        comp = df["comp"].to_numpy()
        for eid in pd.unique(ent):
            m = ent == eid
            et = etype.get(eid, "")
            pt = ptype.get(eid, "")
            if et == "water":
                kind[m] = "water"
            elif et == "branched":
                kind[m] = "branched"
            elif et == "polymer" or pt:
                kind[m] = "nucleic" if "nucleotide" in pt else ("protein" if "peptide" in pt else "polymer")
        water = np.isin(comp, ("HOH", "DOD", "WAT"))
        kind[water] = "water"
        # single-atom metal residues are ions
        lig = kind == "ligand"
        if lig.any():
            metal = df["elem"].isin(METALS).to_numpy()
            res_key = df["asym"] + "|" + df["auth_seq"] + "|" + df["model"].astype(str)
            counts = res_key[lig].value_counts()
            single = res_key.map(counts).fillna(0).to_numpy() == 1
            kind[lig & metal & single] = "ion"
        df["kind"] = kind

        ss = np.full(len(df), "", dtype=object)
        asym = df["asym"].to_numpy()
        seq = df["seq"].to_numpy()
        for r in self.rows("struct_conf"):
            t = r.get("conf_type_id", "")
            code = "H" if t.upper().startswith("HELX") else ("T" if t.upper().startswith("TURN") else "")
            if not code:
                continue
            m = self._range_mask(asym, seq, r.get("beg_label_asym_id"), r.get("beg_label_seq_id"),
                                 r.get("end_label_seq_id"))
            ss[m] = code
        for r in self.rows("struct_sheet_range"):
            m = self._range_mask(asym, seq, r.get("beg_label_asym_id"), r.get("beg_label_seq_id"),
                                 r.get("end_label_seq_id"))
            ss[m] = "E"
        df["ss"] = ss

    @staticmethod
    def _range_mask(asym, seq, a, b, e):
        try:
            b, e = int(b), int(e)
        except (TypeError, ValueError):
            return np.zeros(len(asym), dtype=bool)
        return (asym == a) & (seq >= b) & (seq <= e)

    # -- generic category access ------------------------------------------ #
    def category_names(self) -> list[str]:
        return [c.rstrip(".").lstrip("_") for c in self.block.get_mmcif_category_names()]

    def category_info(self) -> list[dict]:
        info = []
        for c in self.block.get_mmcif_category_names():
            cat = self.block.get_mmcif_category(c, raw=True)
            nrows = len(next(iter(cat.values()))) if cat else 0
            is_loop = self.block.find_loop_item(c + next(iter(cat))) is not None if cat else False
            info.append({"name": c.rstrip(".").lstrip("_"), "rows": nrows, "items": len(cat), "loop": is_loop})
        return info

    def table(self, category: str) -> pd.DataFrame:
        """Category as a DataFrame of display strings ('?' / '.' kept)."""
        category = category.lstrip("_").rstrip(".")
        if category not in self._tables:
            cat = self.block.get_mmcif_category("_" + category + ".", raw=True)
            if not cat:
                self._tables[category] = pd.DataFrame()
            else:
                self._tables[category] = pd.DataFrame({k: _str_values(v) for k, v in cat.items()})
        return self._tables[category]

    def rows(self, category: str) -> list[dict]:
        df = self.table(category)
        return df.to_dict("records") if not df.empty else []

    def value(self, tag: str, default: str = "") -> str:
        """First value of a tag, unquoted; '?' and '.' give ``default``."""
        cat, _, item = tag.lstrip("_").partition(".")
        df = self.table(cat)
        if df.empty or item not in df.columns:
            return default
        v = df[item].iloc[0]
        return default if v in ("?", ".") else v

    def number(self, tag: str) -> float | None:
        v = self.value(tag)
        try:
            return float(re.sub(r"\(\d+\)$", "", v))
        except (TypeError, ValueError):
            return None

    def column(self, tag: str) -> list[str]:
        cat, _, item = tag.lstrip("_").partition(".")
        df = self.table(cat)
        return list(df[item]) if item in df.columns else []

    def raw_category_text(self, category: str, max_rows: int = 25) -> tuple[str, int]:
        """The category re-serialised in PDBx style (truncated to max_rows)."""
        category = category.lstrip("_").rstrip(".")
        cat = self.block.get_mmcif_category("_" + category + ".", raw=True)
        if not cat:
            return "", 0
        n = len(next(iter(cat.values())))
        d = gemmi.cif.Document()
        b = d.add_new_block(self.block.name)
        b.set_mmcif_category("_" + category + ".", {k: v[:max_rows] for k, v in cat.items()}, raw=True)
        opts = gemmi.cif.WriteOptions()
        opts.prefer_pairs = True
        opts.align_pairs = 33
        opts.align_loops = 30
        txt = d.as_string(opts)
        # drop the 'data_' header that the throw-away document adds
        txt = "\n".join(l for l in txt.splitlines() if not l.startswith("data_")).strip("\n")
        return txt, n

    # -- structure-level helpers ------------------------------------------ #
    @property
    def structure(self) -> gemmi.Structure:
        if self._structure is None:
            st = gemmi.make_structure_from_block(self.block)
            st.setup_entities()
            self._structure = st
        return self._structure

    @property
    def models(self) -> list[int]:
        return sorted(int(m) for m in pd.unique(self.atoms["model"]))

    @property
    def methods(self) -> list[str]:
        return [m for m in self.column("exptl.method") if m not in ("?", ".")]

    @property
    def is_computed_model(self) -> bool:
        return bool(self.table("ma_qa_metric_global").shape[0]) or self.source.startswith("AlphaFold")

    @property
    def title(self) -> str:
        return self.value("struct.title") or self.value("struct.pdbx_descriptor") or self.entry_id


# --------------------------------------------------------------------------- #
# In-memory LRU of parsed entries (single-user local app)
# --------------------------------------------------------------------------- #
_CACHE: "OrderedDict[str, Entry]" = OrderedDict()
_LOCK = threading.Lock()


def _remember(entry: Entry) -> Entry:
    with _LOCK:
        _CACHE[entry.key] = entry
        _CACHE.move_to_end(entry.key)
        while len(_CACHE) > config.ENTRY_CACHE_SIZE:
            _CACHE.popitem(last=False)
    return entry


def get_entry(key: str) -> Entry | None:
    """Return a parsed entry; re-parse from the on-disk cache if evicted."""
    with _LOCK:
        if key in _CACHE:
            _CACHE.move_to_end(key)
            return _CACHE[key]
    if key.startswith("upload:"):
        path = config.ENTRY_DIR / "uploads" / (key.split(":", 1)[1] + ".cif.gz")
        if path.exists():
            return _remember(Entry.from_document(gemmi.cif.read(str(path)), key, "uploaded file"))
        return None
    kind, _, ident = key.partition(":")
    try:
        return load_entry(ident if kind == "pdb" else ident)
    except EntryError:
        return None


def load_entry(text: str) -> Entry:
    kind, ident = normalize_id(text)
    key = f"{kind}:{ident}"
    with _LOCK:
        if key in _CACHE:
            _CACHE.move_to_end(key)
            return _CACHE[key]
    path, source = fetch_entry_file(kind, ident)
    try:
        doc = gemmi.cif.read(str(path))
    except (RuntimeError, ValueError) as exc:
        raise EntryError(f"Could not parse {path.name}: {exc}") from exc
    return _remember(Entry.from_document(doc, key, source))


def load_upload(data: bytes, filename: str = "upload.cif") -> Entry:
    if len(data) > config.MAX_UPLOAD_BYTES:
        raise EntryError("The uploaded file is larger than the configured limit.")
    if filename.endswith(".gz") or data[:2] == b"\x1f\x8b":
        try:
            data = gzip.decompress(data)
        except OSError as exc:
            raise EntryError(f"Could not decompress {filename}: {exc}") from exc
    digest = hashlib.sha1(data).hexdigest()[:16]
    key = f"upload:{digest}"
    try:
        doc = gemmi.cif.read_string(data.decode("utf-8", errors="replace"))
    except (RuntimeError, ValueError) as exc:
        raise EntryError(f"{filename} is not a readable mmCIF file: {exc}") from exc
    entry = Entry.from_document(doc, key, f"uploaded file ({filename})")
    up = config.ENTRY_DIR / "uploads"
    up.mkdir(parents=True, exist_ok=True)
    (up / f"{digest}.cif.gz").write_bytes(gzip.compress(data))
    return _remember(entry)
