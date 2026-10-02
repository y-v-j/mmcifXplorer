"""Resolve a row of any mmCIF category to the atoms it describes.

mmCIF categories point at atoms and residues through families of identifier
items (``label_asym_id``/``label_seq_id``/``label_atom_id`` and their
``auth_*`` twins), often prefixed per partner (``ptnr1_``, ``beg_``,
``range_2_`` ...) or suffixed (``_1``, ``_2``).  Instead of hard-coding every
category, the resolver discovers those families from the column names, so it
works for categories it has never seen.  A few categories that identify atoms
indirectly (entity, struct_asym, assembly operators, scheme tables) are handled
explicitly, and categories with no atom identifiers at all are resolved through
their dictionary children (e.g. ``struct_site`` -> ``struct_site_gen``).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .entry import Entry

ASYM_RE = re.compile(r"^(?P<pre>.*?)(?P<kind>label|auth)_asym_id(?P<suf>.*)$")
INVALID = {"", "?", "."}


@dataclass
class Partner:
    """One atom/residue reference found in a row (e.g. ptnr1 of a struct_conn)."""
    label: str
    mask: np.ndarray
    atom_level: bool


@dataclass
class Resolution:
    mask: np.ndarray                      # bool over entry.atoms rows
    level: str = "none"                   # atom | residue | range | entity | chain | none
    partners: list[Partner] = field(default_factory=list)
    note: str = ""
    measure: dict | None = None           # {'kind': 'distance'|'angle', 'reported': str, 'item': str}
    special: str = ""                     # e.g. 'cell' for unit-cell rows

    @property
    def n(self) -> int:
        return int(self.mask.sum())


def _valid(v) -> bool:
    return v is not None and str(v) not in INVALID


def _get(row: dict, pre: str, field_: str, suf: str):
    variants = [pre, "pdbx_" + pre]
    if pre.startswith("pdbx_"):
        variants.append(pre[len("pdbx_"):])
    for p in variants:
        v = row.get(f"{p}{field_}{suf}")
        if _valid(v):
            return str(v)
    return None


def _partner_mask(entry: Entry, row: dict, pre: str, suf: str) -> Partner | None:
    a = entry.atoms
    n = len(a)
    m = np.ones(n, dtype=bool)
    used = False
    la = _get(row, pre, "label_asym_id", suf)
    aa = _get(row, pre, "auth_asym_id", suf)
    if la:
        m &= a["asym"].to_numpy() == la
        used = True
    if aa:
        m &= a["auth_asym"].to_numpy() == aa
        used = True
    if not used:
        return None
    ls = _get(row, pre, "label_seq_id", suf)
    as_ = _get(row, pre, "auth_seq_id", suf)
    if ls:
        try:
            m &= a["seq"].to_numpy() == int(ls)
        except ValueError:
            pass
    if as_:
        m &= a["auth_seq"].to_numpy() == as_
    ins = _get(row, pre, "PDB_ins_code", suf)
    if ins:
        m &= a["ins"].to_numpy() == ins
    comp = _get(row, pre, "label_comp_id", suf) or _get(row, pre, "auth_comp_id", suf)
    if comp:
        m &= (a["comp"].to_numpy() == comp) | (a["auth_comp"].to_numpy() == comp)
    model = row.get("PDB_model_num") or row.get("pdbx_PDB_model_num")
    if _valid(model):
        try:
            m &= a["model"].to_numpy() == int(model)
        except ValueError:
            pass
    residue_level = bool(ls or as_)
    atom = _get(row, pre, "label_atom_id", suf) or _get(row, pre, "auth_atom_id", suf)
    label_bits = [x for x in (aa or la, comp, as_ or ls) if x]
    label = " ".join(label_bits)
    if atom and residue_level:
        am = m & (a["atom"].to_numpy() == atom)
        alt = _get(row, pre, "label_alt_id", suf)
        if alt:
            am &= np.isin(a["alt"].to_numpy(), (alt, ""))
        if am.any():
            return Partner(label=f"{label} {atom}".strip(), mask=am, atom_level=True)
    if not residue_level:
        # Only a chain was given (rare); treat as chain-level reference.
        return Partner(label=label, mask=m, atom_level=False)
    return Partner(label=label, mask=m, atom_level=False)


def _find_groups(row: dict) -> list[tuple[str, str]]:
    groups = []
    for col in row:
        mt = ASYM_RE.match(col)
        if not mt:
            continue
        pre, suf = mt.group("pre"), mt.group("suf")
        canon = pre[len("pdbx_"):] if pre.startswith("pdbx_") and f"{pre[5:]}label_asym_id{suf}" in row else pre
        key = (canon, suf)
        if key not in groups:
            groups.append(key)
    return groups


def _seq_range_mask(entry: Entry, asym: str, beg: str, end: str, auth: bool = False) -> np.ndarray:
    a = entry.atoms
    try:
        b, e = int(beg), int(end)
    except (TypeError, ValueError):
        return np.zeros(len(a), dtype=bool)
    if auth:
        s = np.array([int(re.sub(r"\D+$", "", x) or -10**9) if re.match(r"^-?\d", x) else -10**9
                      for x in a["auth_seq"]])
        return (a["auth_asym"].to_numpy() == asym) & (s >= b) & (s <= e)
    return (a["asym"].to_numpy() == asym) & (a["seq"].to_numpy() >= b) & (a["seq"].to_numpy() <= e)


# --------------------------------------------------------------------------- #
def resolve_row(entry: Entry, category: str, row: dict, dictionary=None, _depth: int = 0) -> Resolution:
    category = category.lstrip("_")
    a = entry.atoms
    none = np.zeros(len(a), dtype=bool)

    # ---- categories identified directly by atom serial -----------------------
    if category in ("atom_site", "atom_site_anisotrop") and _valid(row.get("id")):
        m = a["id"].to_numpy() == str(row["id"])
        return Resolution(m, "atom", [Partner(f"atom {row['id']}", m, True)])

    # ---- cell / symmetry: not atoms, but drawable ---------------------------
    if category in ("cell", "symmetry", "space_group", "space_group_symop", "atom_sites"):
        return Resolution(none, "none", special="cell",
                          note="Shown as the unit-cell box around the model.")

    # ---- structural ranges (helices, strands, sequence references) -----------
    if category in ("struct_conf", "struct_sheet_range", "pdbx_struct_sheet_range"):
        asym = row.get("beg_label_asym_id")
        m = _seq_range_mask(entry, asym, row.get("beg_label_seq_id"), row.get("end_label_seq_id"))
        if not m.any() and _valid(row.get("beg_auth_asym_id")):
            m = _seq_range_mask(entry, row["beg_auth_asym_id"], row.get("beg_auth_seq_id"),
                                row.get("end_auth_seq_id"), auth=True)
        lab = (f"{row.get('beg_auth_asym_id', asym)} {row.get('beg_auth_comp_id', '')}{row.get('beg_auth_seq_id', '')}"
               f" → {row.get('end_auth_comp_id', '')}{row.get('end_auth_seq_id', '')}")
        return Resolution(m, "range", [Partner(lab, m, False)])

    if category == "struct_ref_seq":
        strand = row.get("pdbx_strand_id")
        asyms = set(a.loc[a["auth_asym"] == strand, "asym"])
        m = none.copy()
        for asym in asyms:
            m |= _seq_range_mask(entry, asym, row.get("seq_align_beg"), row.get("seq_align_end"))
        return Resolution(m, "range", [Partner(f"chain {strand} {row.get('seq_align_beg')}–{row.get('seq_align_end')}", m, False)])

    if category == "struct_ref_seq_dif":
        strand, num = row.get("pdbx_pdb_strand_id"), row.get("seq_num")
        m = none.copy()
        if _valid(strand) and _valid(num):
            m = (a["auth_asym"].to_numpy() == strand) & (a["seq"].to_numpy() == int(num))
        elif _valid(strand) and _valid(row.get("pdbx_auth_seq_num")):
            m = (a["auth_asym"].to_numpy() == strand) & (a["auth_seq"].to_numpy() == row["pdbx_auth_seq_num"])
        note = "" if m.any() else "This difference refers to a residue that is not in the model (e.g. an expression tag)."
        return Resolution(m, "residue", [Partner(f"{strand} {row.get('mon_id', '')}{row.get('pdbx_auth_seq_num', '')}", m, False)], note=note)

    if category in ("pdbx_poly_seq_scheme", "pdbx_nonpoly_scheme", "pdbx_branch_scheme"):
        asym = row.get("asym_id")
        m = a["asym"].to_numpy() == asym
        if category == "pdbx_poly_seq_scheme" and _valid(row.get("seq_id")):
            m &= a["seq"].to_numpy() == int(row["seq_id"])
        elif _valid(row.get("pdb_seq_num")):
            m &= a["auth_seq"].to_numpy() == row["pdb_seq_num"]
        note = "" if m.any() else "Not observed: this residue is in the sequence but has no coordinates."
        return Resolution(m, "residue", [Partner(f"{row.get('pdb_strand_id', asym)} {row.get('mon_id', '')}"
                                                 f"{row.get('pdb_seq_num', row.get('seq_id', ''))}", m, False)], note=note)

    if category in ("pdbx_struct_assembly_gen",):
        ids = [s.strip() for s in str(row.get("asym_id_list", "")).split(",") if s.strip()]
        m = np.isin(a["asym"].to_numpy(), ids)
        return Resolution(m, "chain", [Partner("asym " + ",".join(ids), m, False)],
                          note=f"Operators applied to these chains: {row.get('oper_expression', '?')}")

    if category in ("struct_asym", "pdbx_struct_asym"):
        m = a["asym"].to_numpy() == row.get("id")
        return Resolution(m, "chain", [Partner(f"asym {row.get('id')}", m, False)])

    if category == "atom_type" and _valid(row.get("symbol")):
        m = a["elem"].to_numpy() == str(row["symbol"]).upper()
        return Resolution(m, "atom", [Partner(f"every {row['symbol']} atom", m, True)])

    if category == "struct_site" and _valid(row.get("id")):
        gen = entry.table("struct_site_gen")
        m = none.copy()
        if not gen.empty and "site_id" in gen.columns:
            for r in gen[gen["site_id"] == row["id"]].to_dict("records"):
                m |= resolve_row(entry, "struct_site_gen", r, dictionary, _depth=1).mask
        own = [_partner_mask(entry, row, pre, suf) for pre, suf in _find_groups(row)]
        for p in own:
            if p is not None:
                m |= p.mask
        return Resolution(m, "residue", [Partner(f"site {row['id']}", m, False)],
                          note="Site residues from _struct_site_gen plus the ligand the site is defined for.")

    if category == "chem_comp":
        m = a["comp"].to_numpy() == row.get("id")
        return Resolution(m, "entity", [Partner(f"all {row.get('id')}", m, False)])

    if category in ("chem_comp_atom", "chem_comp_bond", "chem_comp_angle"):
        comp = row.get("comp_id")
        names = [row.get(k) for k in ("atom_id", "atom_id_1", "atom_id_2", "atom_id_3") if _valid(row.get(k))]
        base = a["comp"].to_numpy() == comp
        m = base & np.isin(a["atom"].to_numpy(), names)
        return Resolution(m, "atom", [Partner(f"{comp} {'-'.join(names)} (every copy)", m, True)])

    # ---- generic: discover identifier families from the column names --------
    groups = _find_groups(row)
    partners = []
    for pre, suf in groups:
        p = _partner_mask(entry, row, pre, suf)
        if p is not None and p.mask.any():
            partners.append(p)
        elif p is not None:
            partners.append(p)

    # Range form: beg_* / end_* pairs
    begs = [g for g in groups if "beg" in g[0]]
    ends = [g for g in groups if "end" in g[0]]
    if begs and ends:
        bp, ep = begs[0], ends[0]
        asym = _get(row, bp[0], "label_asym_id", bp[1])
        b = _get(row, bp[0], "label_seq_id", bp[1])
        e = _get(row, ep[0], "label_seq_id", ep[1])
        m = _seq_range_mask(entry, asym, b, e) if (asym and b and e) else none
        if m.any():
            return Resolution(m, "range", [Partner(f"{asym} {b}→{e}", m, False)])

    if partners:
        mask = np.zeros(len(a), dtype=bool)
        for p in partners:
            mask |= p.mask
        level = "atom" if all(p.atom_level for p in partners) else "residue"
        res = Resolution(mask, level, partners)
        # Measurements reported in the row (distance / angle) can be checked.
        atom_partners = [p for p in partners if p.atom_level and p.mask.any()]
        sym = [row.get(k) for k in row if k.endswith("symmetry") and _valid(row.get(k))]
        if len(atom_partners) == 2:
            rep_item = next((k for k in ("pdbx_dist_value", "dist", "distance", "bond_value", "value")
                             if _valid(row.get(k))), None)
            res.measure = {"kind": "distance", "reported": row.get(rep_item) if rep_item else None,
                           "item": rep_item}
        elif len(atom_partners) == 3:
            rep_item = next((k for k in ("value", "angle_value", "angle") if _valid(row.get(k))), None)
            res.measure = {"kind": "angle", "reported": row.get(rep_item) if rep_item else None,
                           "item": rep_item}
        if any(s != "1_555" for s in sym):
            res.note = ("One partner is in a symmetry-related copy (" + ", ".join(sorted(set(sym))) +
                        "); the copy in the asymmetric unit is highlighted, so the drawn distance is "
                        "not the reported one.")
            if res.measure:
                res.measure["symmetry"] = True
        if not mask.any():
            res.note = res.note or "The referenced residues/atoms are not present in the model coordinates."
        return res

    # ---- entity level --------------------------------------------------------
    if category == "entity" and _valid(row.get("id")):
        m = a["entity"].to_numpy() == str(row["id"])
        return Resolution(m, "entity", [Partner(f"entity {row['id']}", m, False)])
    if _valid(row.get("entity_id")):
        m = a["entity"].to_numpy() == str(row["entity_id"])
        num = row.get("num")
        if category == "entity_poly_seq" and _valid(num):
            m &= a["seq"].to_numpy() == int(num)
            return Resolution(m, "residue", [Partner(f"entity {row['entity_id']} #{num} {row.get('mon_id', '')}", m, False)],
                              note="" if m.any() else "Residue in the sequence without coordinates.")
        return Resolution(m, "entity", [Partner(f"entity {row['entity_id']}", m, False)])
    if _valid(row.get("comp_id")) and category.startswith(("pdbx_entity_nonpoly", "pdbx_struct_mod")):
        m = a["comp"].to_numpy() == row["comp_id"]
        return Resolution(m, "entity", [Partner(row["comp_id"], m, False)])

    # ---- indirect: follow dictionary children one level down -----------------
    if dictionary is not None and _depth == 0:
        cat = dictionary.category(category)
        if cat:
            mask = none.copy()
            labels = []
            for key in cat["keys"]:
                kitem = dictionary.item(key)
                val = row.get(key.split(".", 1)[1])
                if not kitem or not _valid(val):
                    continue
                for child in kitem["children"]:
                    ccat, _, citem = child.lstrip("_").partition(".")
                    if ccat == category:
                        continue
                    df = entry.table(ccat)
                    if df.empty or citem not in df.columns:
                        continue
                    sub = df[df[citem] == val]
                    if sub.empty or len(sub) > 300:
                        continue
                    for r in sub.to_dict("records"):
                        rr = resolve_row(entry, ccat, r, dictionary, _depth=1)
                        mask |= rr.mask
                    if mask.any():
                        labels.append(f"{len(sub)} row(s) of _{ccat}")
            if mask.any():
                return Resolution(mask, "residue", [Partner(f"{category} {row.get(cat['keys'][0].split('.', 1)[1], '')}", mask, False)],
                                  note="Resolved through " + "; ".join(labels) + ".")
    return Resolution(none, "none", note="This row does not refer to atoms in the model.")


def can_map_category(entry: Entry, category: str, dictionary=None) -> bool:
    """True when sample rows of the category resolve to atoms (or the cell)."""
    df = entry.table(category)
    if df.empty:
        return False
    if category in ("cell", "symmetry", "atom_site", "atom_type", "entity", "struct_asym", "chem_comp"):
        return True
    if any(ASYM_RE.match(c) for c in df.columns) or "entity_id" in df.columns:
        # identifier columns present: check the first row only
        return bool(resolve_row(entry, category, df.iloc[0].to_dict(), dictionary).mask.any()) or \
            any(resolve_row(entry, category, r, dictionary).mask.any() for r in df.iloc[1:3].to_dict("records"))
    rows = df.head(1).to_dict("records")
    for r in rows:
        res = resolve_row(entry, category, r, dictionary)
        if res.special or res.mask.any():
            return True
    return False


def residue_references(entry: Entry, asym: str, auth_asym: str, seq: int, auth_seq: str, ins: str,
                       categories: list[str]) -> list[tuple[str, int, str]]:
    """Rows of ``categories`` that mention one residue, found by comparing the
    identifier columns directly (fast; no per-row atom masks)."""
    out = []
    seq_s = str(seq) if seq is not None and seq >= 0 else None
    for cat in categories:
        df = entry.table(cat)
        if df.empty:
            continue
        cols = set(df.columns)
        hit = np.zeros(len(df), dtype=bool)
        # explicit ranges (helices, strands)
        if {"beg_label_asym_id", "beg_label_seq_id", "end_label_seq_id"} <= cols and seq_s:
            b = pd.to_numeric(df["beg_label_seq_id"], errors="coerce").to_numpy()
            e = pd.to_numeric(df["end_label_seq_id"], errors="coerce").to_numpy()
            hit |= (df["beg_label_asym_id"].to_numpy() == asym) & (b <= seq) & (e >= seq)
        else:
            for pre, suf in _find_groups({c: None for c in df.columns}):
                def col(field_):
                    for p in (pre, "pdbx_" + pre, pre[5:] if pre.startswith("pdbx_") else None):
                        if p is not None and f"{p}{field_}{suf}" in cols:
                            return df[f"{p}{field_}{suf}"].to_numpy()
                    return None
                la, aa = col("label_asym_id"), col("auth_asym_id")
                ls, as_ = col("label_seq_id"), col("auth_seq_id")
                m = np.ones(len(df), dtype=bool)
                if la is not None:
                    m &= la == asym
                elif aa is not None:
                    m &= aa == auth_asym
                else:
                    continue
                if as_ is not None:
                    m &= as_ == auth_seq
                elif ls is not None and seq_s:
                    m &= ls == seq_s
                else:
                    continue
                hit |= m
        first = df.columns[0]
        for i in np.flatnonzero(hit)[:50]:
            r = df.iloc[i]
            label = str(r["id"]) if "id" in cols else str(r[first])
            out.append((cat, int(i), label))
    return out
