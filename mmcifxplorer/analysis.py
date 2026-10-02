"""Derived structural and crystallographic quantities for one entry.

Everything here is computed from the entry itself so the dashboard can put a
reported value next to an independent recomputation (cell volume, Matthews
coefficient, solvent content, bond distances, mean B ...).
"""
from __future__ import annotations

import math
import re

import gemmi
import numpy as np
import pandas as pd

from .entry import Entry

# Protein partial specific volume used by Matthews (1968): 1/1.230 Å^3/Da
MATTHEWS_K = 1.230


def _num(x):
    try:
        return float(re.sub(r"\(\d+\)$", "", str(x)))
    except (TypeError, ValueError):
        return None


def first_number(entry: Entry, *tags: str):
    for t in tags:
        v = entry.number(t)
        if v is not None:
            return v, t
    return None, None


# --------------------------------------------------------------------------- #
# Summary
# --------------------------------------------------------------------------- #
def summary(entry: Entry) -> dict:
    a = entry.atoms
    methods = entry.methods
    res, res_tag = first_number(entry, "refine.ls_d_res_high", "em_3d_reconstruction.resolution",
                                "reflns.d_resolution_high", "pdbx_refine.ls_d_res_high")
    rwork, rwork_tag = first_number(entry, "refine.ls_R_factor_R_work", "refine.ls_R_factor_obs",
                                    "refine.ls_R_factor_all")
    rfree, rfree_tag = first_number(entry, "refine.ls_R_factor_R_free")
    m1 = a[a["model"] == entry.models[0]]
    organisms = sorted({v for v in entry.column("entity_src_gen.pdbx_gene_src_scientific_name")
                        + entry.column("entity_src_nat.pdbx_organism_scientific_name")
                        + entry.column("pdbx_entity_src_syn.organism_scientific")
                        if v not in ("?", ".")})
    plddt = None
    if entry.is_computed_model:
        plddt = entry.number("ma_qa_metric_global.metric_value")
    return {
        "id": entry.entry_id,
        "title": entry.title,
        "methods": methods,
        "resolution": res, "resolution_tag": res_tag,
        "r_work": rwork, "r_work_tag": rwork_tag,
        "r_free": rfree, "r_free_tag": rfree_tag,
        "space_group": entry.value("symmetry.space_group_name_H-M"),
        "cell": [entry.number(f"cell.{k}") for k in ("length_a", "length_b", "length_c",
                                                     "angle_alpha", "angle_beta", "angle_gamma")],
        "n_models": len(entry.models),
        "n_atoms_model1": int(len(m1)),
        "n_atoms_total": int(len(a)),
        "n_entities": len(entry.rows("entity")),
        "n_chains": int(m1["auth_asym"].nunique()),
        "kinds": m1["kind"].value_counts().to_dict(),
        "mean_b": float(m1["b"].mean()) if len(m1) else None,
        "deposited": entry.value("pdbx_database_status.recvd_initial_deposition_date")
        or entry.value("pdbx_audit_revision_history.revision_date"),
        "organisms": organisms,
        "keywords": entry.value("struct_keywords.pdbx_keywords"),
        "source": entry.source,
        "plddt": plddt,
        "computed": entry.is_computed_model,
    }


# --------------------------------------------------------------------------- #
# Crystallography
# --------------------------------------------------------------------------- #
def unit_cell(entry: Entry) -> gemmi.UnitCell | None:
    vals = [entry.number(f"cell.{k}") for k in ("length_a", "length_b", "length_c",
                                               "angle_alpha", "angle_beta", "angle_gamma")]
    if any(v is None for v in vals):
        return None
    cell = gemmi.UnitCell(*vals)
    if cell.volume < 10 or (vals[0] == vals[1] == vals[2] == 1.0):
        return None
    return cell


def space_group(entry: Entry) -> gemmi.SpaceGroup | None:
    name = entry.value("symmetry.space_group_name_H-M")
    if not name:
        return None
    sg = gemmi.find_spacegroup_by_name(name)
    return sg


def asu_mass(entry: Entry) -> tuple[float | None, list[dict]]:
    """Mass (Da) of polymer + non-polymer (non-water) content of the ASU.

    Uses _entity.formula_weight × _entity.pdbx_number_of_molecules, i.e. the
    depositor-declared content, which is what Matthews' Vm is defined on.
    """
    rows = []
    total = 0.0
    for r in entry.rows("entity"):
        if r.get("type") == "water":
            continue
        fw, n = _num(r.get("formula_weight")), _num(r.get("pdbx_number_of_molecules"))
        if fw is None or n is None:
            continue
        rows.append({"entity": r.get("id"), "type": r.get("type"), "description": r.get("pdbx_description", ""),
                     "formula_weight": fw, "copies": int(n), "mass": fw * n})
        total += fw * n
    return (total or None), rows


def matthews(entry: Entry) -> dict | None:
    """Recompute Vm and solvent fraction and compare with _exptl_crystal."""
    cell, sg = unit_cell(entry), space_group(entry)
    if cell is None or sg is None:
        return None
    z = len(list(sg.operations()))
    mass, parts = asu_mass(entry)
    if not mass:
        return None
    vm = cell.volume / (z * mass)
    solvent = 1.0 - MATTHEWS_K / vm
    return {
        "volume": cell.volume, "z_sym": z, "mass": mass, "parts": parts,
        "vm": vm, "solvent": solvent * 100.0,
        "vm_reported": entry.number("exptl_crystal.density_Matthews"),
        "solvent_reported": entry.number("exptl_crystal.density_percent_sol"),
        "z_pdb": entry.number("cell.Z_PDB"),
    }


def crystal_packing(entry: Entry, max_points: int = 6000) -> dict | None:
    """Backbone trace of the ASU and its symmetry copies inside one unit cell."""
    cell, sg = unit_cell(entry), space_group(entry)
    if cell is None or sg is None:
        return None
    a = entry.atoms
    m = (a["model"] == entry.models[0]) & a["atom"].isin(["CA", "P"])
    sub = a[m]
    if sub.empty:
        return None
    if len(sub) * len(list(sg.operations())) > max_points:
        step = int(math.ceil(len(sub) * len(list(sg.operations())) / max_points))
        sub = sub.iloc[::step]
    xyz = sub[["x", "y", "z"]].to_numpy()
    frac_m = np.array(cell.frac.mat.tolist())
    frac_v = np.array(cell.frac.vec.tolist())
    orth_m = np.array(cell.orth.mat.tolist())
    orth_v = np.array(cell.orth.vec.tolist())
    frac = xyz @ frac_m.T + frac_v
    copies = []
    for op in sg.operations():
        rot = np.array(op.rot, dtype=float) / op.DEN
        tran = np.array(op.tran, dtype=float) / op.DEN
        f = frac @ rot.T + tran
        shift = -np.floor(f.mean(axis=0))      # bring the copy's centroid into [0,1)
        f = f + shift
        cart = f @ orth_m.T + orth_v
        copies.append({"op": op.triplet(), "xyz": cart, "chain": sub["auth_asym"].to_numpy()})
    return {"copies": copies, "cell": cell, "sg": sg.xhm()}


def shells(entry: Entry) -> dict[str, pd.DataFrame]:
    out = {}
    for cat in ("reflns_shell", "refine_ls_shell"):
        df = entry.table(cat)
        if df.empty:
            continue
        num = df.apply(lambda c: pd.to_numeric(c, errors="coerce"))
        out[cat] = num
    return out


# --------------------------------------------------------------------------- #
# Per-residue properties
# --------------------------------------------------------------------------- #
def residue_table(entry: Entry, model: int | None = None) -> pd.DataFrame:
    a = entry.atoms
    model = model or entry.models[0]
    a = a[(a["model"] == model) & a["kind"].isin(["protein", "nucleic", "polymer"])]
    if a.empty:
        return pd.DataFrame()
    g = a.groupby(["auth_asym", "asym", "seq", "auth_seq", "ins", "comp"], sort=False)
    df = g.agg(b_mean=("b", "mean"), b_max=("b", "max"), occ_min=("occ", "min"), n_atoms=("b", "size"),
               ss=("ss", "first"), kind=("kind", "first")).reset_index()
    return df


def _peptide_linked(prev: gemmi.Residue, cur: gemmi.Residue) -> bool:
    c = prev.find_atom("C", "*")
    n = cur.find_atom("N", "*")
    return bool(c and n and c.pos.dist(n.pos) < 2.0)


def ramachandran(entry: Entry, model_index: int = 0) -> pd.DataFrame:
    """phi/psi for every amino-acid residue with both angles defined.

    Residues whose neighbour is missing or not peptide-bonded (chain breaks)
    get no angle, which mirrors how validation software counts 'analysed'
    residues.
    """
    st = entry.structure
    if len(st) == 0:
        return pd.DataFrame()
    model = st[model_index]
    rows = []
    for ch in model:
        pol = ch.get_polymer()
        if pol.check_polymer_type() not in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD):
            continue
        n = len(pol)
        for i in range(n):
            res = pol[i]
            prev = pol[i - 1] if i > 0 and _peptide_linked(pol[i - 1], res) else None
            nxt = pol[i + 1] if i + 1 < n and _peptide_linked(res, pol[i + 1]) else None
            if prev is None or nxt is None:
                continue
            phi, psi = gemmi.calculate_phi_psi(prev, res, nxt)
            if math.isnan(phi) or math.isnan(psi):
                continue
            nxt_name = nxt.name
            cls = ("Glycine" if res.name == "GLY" else "Proline" if res.name == "PRO"
                   else "Pre-proline" if nxt_name == "PRO" else
                   "Ile/Val" if res.name in ("ILE", "VAL") else "General")
            rows.append({"chain": ch.name, "asym": res.subchain, "seq": res.label_seq if res.label_seq is not None else -1,
                         "auth_seq": f"{res.seqid.num}{res.seqid.icode.strip()}", "comp": res.name,
                         "phi": math.degrees(phi), "psi": math.degrees(psi), "class": cls})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    flagged = set()
    for r in entry.rows("pdbx_validate_torsion"):
        if str(r.get("PDB_model_num", "1")) in ("1", "?", ".") or model_index == 0:
            flagged.add((r.get("auth_asym_id"), str(r.get("auth_seq_id")) + (r.get("PDB_ins_code") or "").replace("?", "")))
    df["flagged"] = [(c, s) in flagged for c, s in zip(df["chain"], df["auth_seq"])]
    return df


def ensemble_rmsf(entry: Entry) -> pd.DataFrame:
    """Per-residue CA spread across NMR/ensemble models (no re-superposition;
    deposited ensembles are already superposed)."""
    a = entry.atoms
    if len(entry.models) < 2:
        return pd.DataFrame()
    ca = a[a["atom"].isin(["CA", "P"]) & (a["alt"].isin(["", "A"]))]
    rows = []
    for (chain, asym, seq, auth_seq, comp), g in ca.groupby(["auth_asym", "asym", "seq", "auth_seq", "comp"], sort=False):
        if g["model"].nunique() < 2:
            continue
        xyz = g[["x", "y", "z"]].to_numpy()
        rmsf = float(np.sqrt(((xyz - xyz.mean(axis=0)) ** 2).sum(axis=1).mean()))
        rows.append({"chain": chain, "asym": asym, "seq": seq, "auth_seq": auth_seq, "comp": comp,
                     "rmsf": rmsf, "n_models": int(g["model"].nunique())})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Composition
# --------------------------------------------------------------------------- #
def composition(entry: Entry) -> pd.DataFrame:
    a = entry.atoms
    m1 = a[a["model"] == entry.models[0]]
    poly = {r.get("entity_id"): r for r in entry.rows("entity_poly")}
    seq_len = {}
    for r in entry.rows("entity_poly_seq"):
        seq_len[r.get("entity_id")] = seq_len.get(r.get("entity_id"), 0) + 1
    rows = []
    for r in entry.rows("entity"):
        eid = r.get("id")
        sub = m1[m1["entity"] == eid]
        p = poly.get(eid, {})
        observed = sub[sub["seq"] >= 0].groupby(["asym", "seq"]).ngroups if len(sub) else 0
        rows.append({
            "entity": eid,
            "type": r.get("type", ""),
            "polymer_type": p.get("type", ""),
            "description": r.get("pdbx_description", ""),
            "formula_weight": _num(r.get("formula_weight")),
            "copies": _num(r.get("pdbx_number_of_molecules")),
            "auth_chains": ",".join(dict.fromkeys(sub["auth_asym"])),
            "label_asyms": ",".join(dict.fromkeys(sub["asym"])),
            "seq_length": seq_len.get(eid, 0),
            "observed_residues": observed,
            "atoms": int(len(sub)),
            "sequence": (p.get("pdbx_seq_one_letter_code_can") or "").replace("\n", ""),
        })
    return pd.DataFrame(rows)


def sequence_coverage(entry: Entry, entity_id: str) -> tuple[str, list[dict]]:
    """One-letter sequence of an entity and its modelled segments (label numbering)."""
    p = next((r for r in entry.rows("entity_poly") if r.get("entity_id") == entity_id), None)
    if not p:
        return "", []
    seq = (p.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "").replace(" ", "")
    a = entry.atoms
    sub = a[(a["entity"] == entity_id) & (a["model"] == entry.models[0]) & (a["seq"] >= 0)]
    if sub.empty:
        return seq, []
    first_asym = sub["asym"].iloc[0]
    observed = sorted(set(sub.loc[sub["asym"] == first_asym, "seq"]))
    segs = []
    start = prev = None
    for s in observed:
        if start is None:
            start = prev = s
        elif s == prev + 1:
            prev = s
        else:
            segs.append((start, prev))
            start = prev = s
    if start is not None:
        segs.append((start, prev))
    return seq, [{"start": s, "end": e, "asym": first_asym} for s, e in segs]
