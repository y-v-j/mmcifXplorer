"""Build the data the dash-bio Molecule3dViewer needs from an Entry.

The viewer works with a *view*: the atoms of one model, optionally without
water, and (for very large entries) reduced to a backbone trace plus ligands.
``view.index`` maps viewer atom serials back to rows of ``entry.atoms`` and
therefore to lines of ``_atom_site``.
"""
from __future__ import annotations

from dataclasses import dataclass

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from . import config
from .entry import METALS, Entry

TRACE_ATOMS = {"CA", "P", "C4'"}


@dataclass
class View:
    model_data: dict
    index: np.ndarray        # viewer serial -> entry.atoms row
    meta: dict
    trace_only: bool


def _covalent_radii(elems: np.ndarray) -> np.ndarray:
    cache = {}
    out = np.empty(len(elems))
    for i, e in enumerate(elems):
        r = cache.get(e)
        if r is None:
            el = gemmi.Element(e if e else "C")
            r = el.covalent_r if el.covalent_r > 0 else 0.77
            cache[e] = r
        out[i] = r
    return out


def compute_bonds(xyz: np.ndarray, elems: np.ndarray, alts: np.ndarray, kinds: np.ndarray,
                  tolerance: float = 0.45) -> np.ndarray:
    """Distance-based covalent bonds (pairs of view indices).

    Metal atoms are excluded on purpose: their coordination is shown through
    _struct_conn records rather than drawn as covalent sticks.
    """
    if len(xyz) < 2:
        return np.empty((0, 2), dtype=int)
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(r=2.6, output_type="ndarray")
    if len(pairs) == 0:
        return np.empty((0, 2), dtype=int)
    i, j = pairs[:, 0], pairs[:, 1]
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    rad = _covalent_radii(elems)
    ok = (d > 0.4) & (d <= rad[i] + rad[j] + tolerance)
    metal = np.isin(elems, list(METALS))
    ok &= ~(metal[i] | metal[j])
    ok &= ~((kinds[i] == "water") | (kinds[j] == "water"))
    # never bond two different alternate conformers
    ai, aj = alts[i], alts[j]
    ok &= ~((ai != "") & (aj != "") & (ai != aj))
    return pairs[ok]


def build_view(entry: Entry, model: int | None = None, show_water: bool = False,
               max_atoms: int = config.MAX_VIEW_ATOMS) -> View:
    a = entry.atoms
    models = entry.models
    if model is None or model not in models:
        model = models[0]
    m = a["model"].to_numpy() == model
    if not show_water:
        m &= a["kind"].to_numpy() != "water"
    idx = np.flatnonzero(m)
    trace_only = False
    if len(idx) > max_atoms:
        trace_only = True
        sub = a.iloc[idx]
        keep = sub["atom"].isin(TRACE_ATOMS).to_numpy() | ~sub["kind"].isin(["protein", "nucleic", "polymer"]).to_numpy()
        idx = idx[keep]
        if len(idx) > max_atoms:            # still too big: drop ligands/water beyond the polymer trace
            sub = a.iloc[idx]
            idx = idx[sub["atom"].isin(TRACE_ATOMS).to_numpy()][:max_atoms]

    v = a.iloc[idx]
    xyz = v[["x", "y", "z"]].to_numpy()
    elems = v["elem"].to_numpy()
    kinds = v["kind"].to_numpy()
    bonds = np.empty((0, 2), dtype=int) if trace_only else compute_bonds(xyz, elems, v["alt"].to_numpy(), kinds)

    chains = list(dict.fromkeys(v["auth_asym"]))
    entities = list(dict.fromkeys(v["entity"]))
    chain_i = {c: i for i, c in enumerate(chains)}
    ent_i = {e: i for i, e in enumerate(entities)}
    resi = np.where(v["seq"].to_numpy() >= 0, v["seq"].to_numpy(),
                    [int(s) if s.lstrip("-").isdigit() else 0 for s in v["auth_seq"]])
    elem_title = [e.capitalize() if len(e) > 1 else e for e in elems]
    weights = {}
    atoms = []
    for k, (row, x, y, z) in enumerate(zip(v.itertuples(index=False), xyz[:, 0], xyz[:, 1], xyz[:, 2])):
        el = elem_title[k]
        w = weights.get(el)
        if w is None:
            w = gemmi.Element(el).weight if el else 12.0
            weights[el] = w
        atoms.append({
            "serial": k,
            "name": row.atom,
            "elem": el,
            "chain": row.auth_asym,
            "residue_index": int(resi[k]),
            "residue_name": row.comp,
            "positions": [round(float(x), 3), round(float(y), 3), round(float(z), 3)],
            "mass_magnitude": w,
            # extra per-atom fields used by the clientside colouring code
            "b": round(float(row.b), 2),
            "q": round(float(row.occ), 2),
            "c": chain_i[row.auth_asym],
            "e": ent_i[row.entity],
            "k": row.kind,
            "s": row.ss,
            "alt": row.alt,
        })
    bond_list = [{"atom1_index": int(i), "atom2_index": int(j), "bond_order": 1} for i, j in bonds]

    poly = v["kind"].isin(["protein", "nucleic", "polymer"]).to_numpy()
    bsrc = v["b"].to_numpy()[poly] if poly.any() else v["b"].to_numpy()
    blo, bhi = (float(np.percentile(bsrc, 2)), float(np.percentile(bsrc, 98))) if len(bsrc) else (0.0, 1.0)
    if bhi - blo < 1e-6:
        bhi = blo + 1.0
    meta = {
        "model": model,
        "n_atoms": len(idx),
        "chains": chains,
        "entities": entities,
        "b_range": [round(blo, 2), round(bhi, 2)],
        "plddt": entry.is_computed_model,
        "trace_only": trace_only,
        "n_bonds": len(bond_list),
    }
    return View({"atoms": atoms, "bonds": bond_list}, idx, meta, trace_only)


def to_view_indices(view_index: np.ndarray, mask: np.ndarray) -> list[int]:
    """Viewer serials of the atoms selected by a mask over entry.atoms."""
    return np.flatnonzero(mask[view_index]).tolist()


def cell_edges(entry: Entry) -> list[tuple[list[float], list[float]]]:
    """The 12 edges of the unit cell in Cartesian coordinates (Å)."""
    try:
        cell = gemmi.UnitCell(*(entry.number(f"cell.{k}") or 0 for k in
                                ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")))
    except (TypeError, ValueError, RuntimeError):
        return []
    if cell.volume < 10 or (cell.a == 1 and cell.b == 1 and cell.c == 1):
        return []
    corners = {}
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                p = cell.orthogonalize(gemmi.Fractional(i, j, k))
                corners[(i, j, k)] = [p.x, p.y, p.z]
    edges = []
    for (i, j, k), p in corners.items():
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            q = (i + d[0], j + d[1], k + d[2])
            if q in corners:
                edges.append((p, corners[q]))
    return edges
