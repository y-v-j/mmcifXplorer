"""Offline tests for mmcifXplorer.

Oracles are independent of the code under test:
* the wwPDB validation report shipped with the project (docs/1tup_full_validation.pdf)
  for atom counts, chain composition and Ramachandran 'analysed' residues;
* values the depositor/annotator wrote into the file itself (struct_conn
  distances, pdbx_struct_conn_angle angles), recomputed from coordinates;
* textbook formulas (cell volume, Matthews coefficient).
"""
from __future__ import annotations

import math
import warnings
from pathlib import Path

import gemmi
import numpy as np
import pytest

warnings.filterwarnings("ignore", message="You may be importing Biopython")

from mmcifxplorer import analysis as an  # noqa: E402
from mmcifxplorer.dictionary import load_dictionary, posix_to_python_regex  # noqa: E402
from mmcifxplorer.entry import Entry, EntryError, normalize_id  # noqa: E402
from mmcifxplorer.knowledge import section_of  # noqa: E402
from mmcifxplorer.mapping import residue_references, resolve_row  # noqa: E402
from mmcifxplorer.viewer import build_view, to_view_indices  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="session")
def dic():
    return load_dictionary()


@pytest.fixture(scope="session")
def tup():
    return Entry.from_document(gemmi.cif.read(str(ROOT / "1tup.cif")), "test:1TUP", "bundled")


# --------------------------------------------------------------------------- dictionary
def test_dictionary_contents(dic):
    assert dic.version.startswith("5.")
    assert len(dic.items) > 6500 and len(dic.categories) > 600
    it = dic.item("_refine.ls_R_factor_R_free")
    assert it["type"] == "float" and dic.range_text(it["range"]) == "[0, 1]"
    assert [e["value"] for e in dic.item("_struct_conn.conn_type_id")["enum"]] == ["covale", "disulf", "metalc", "hydrog"]
    assert "_entity.id" in dic.item("_atom_site.label_entity_id")["parents"]
    assert "_atom_site.label_entity_id" in dic.item("_entity.id")["children"]
    assert dic.category("atom_site")["keys"] == ["_atom_site.id"]


def test_every_type_regex_compiles(dic):
    bad = [code for code in dic.types if dic.type_regex(code) is None]
    assert bad == []


def test_value_checks(dic):
    assert dic.check_value("_refine.ls_R_factor_R_free", "0.25")["range_ok"] is True
    assert dic.check_value("_refine.ls_R_factor_R_free", "1.5")["range_ok"] is False
    assert dic.check_value("_atom_site.label_asym_id", "A")["type_ok"] is True
    assert dic.check_value("_atom_site.label_asym_id", "A B")["type_ok"] is False   # code: no spaces
    assert dic.check_value("_struct_conn.conn_type_id", "metalc")["enum_ok"] is True
    assert dic.check_value("_struct_conn.conn_type_id", "salt")["enum_ok"] is False
    assert dic.check_value("_cell.angle_beta", "105.50(2)")["range_ok"] is True    # su in parentheses


def test_posix_bracket_backslash_is_literal():
    import re
    pat = re.compile(posix_to_python_regex(r"[\{}]*"))
    assert pat.fullmatch("\\{}")


def test_modelcif_extension_merged(dic):
    if not dic.extensions:
        pytest.skip("ModelCIF dictionary not available offline")
    assert dic.item("_ma_qa_metric_local.metric_value")["source"] == "ModelCIF"
    assert "source" not in dic.item("_atom_site.id")          # PDBx definitions win


def test_every_item_in_1tup_is_defined(dic, tup):
    missing = [f"_{c}.{col}" for c in tup.category_names() for col in tup.table(c).columns
               if dic.item(f"_{c}.{col}") is None]
    assert missing == []


# --------------------------------------------------------------------------- entry vs validation report
def test_composition_matches_validation_report(tup):
    a = tup.atoms
    assert len(a) == 5828                                    # report: "5828 atoms"
    poly = a[a["kind"].isin(["protein", "nucleic"])]
    assert poly.groupby("auth_asym").size().to_dict() == {"A": 1535, "B": 1522, "C": 1529, "E": 420, "F": 435}
    assert (a["kind"] == "ion").sum() == 3 and set(a.loc[a["kind"] == "ion", "comp"]) == {"ZN"}
    assert (a["kind"] == "water").sum() == 22 + 14 + 198 + 70 + 80
    res = poly.groupby(["auth_asym", "auth_seq"]).ngroups
    assert res == 21 + 21 + 196 + 194 + 195                 # residues per chain in the report


def test_ramachandran_analysed_counts_match_report(tup):
    rama = an.ramachandran(tup)
    assert rama.groupby("chain").size().to_dict() == {"A": 194, "B": 192, "C": 193}
    assert rama["phi"].between(-180, 180).all() and rama["psi"].between(-180, 180).all()


def test_summary_values(tup):
    s = an.summary(tup)
    assert s["resolution"] == pytest.approx(2.2) and s["r_work"] == pytest.approx(0.202)
    assert s["r_free"] is None and s["space_group"] == "C 1 2 1"


def test_label_and_auth_numbering_differ_for_p53(tup):
    a = tup.atoms
    first = a[(a["auth_asym"] == "A") & (a["kind"] == "protein")].iloc[0]
    assert first["asym"] == "C" and int(first["auth_seq"]) - first["seq"] == 93


# --------------------------------------------------------------------------- crystallography
def test_cell_volume_formula(tup):
    cell = an.unit_cell(tup)
    a, b, c, al, be, ga = (cell.a, cell.b, cell.c, *(math.radians(x) for x in (cell.alpha, cell.beta, cell.gamma)))
    v = a * b * c * math.sqrt(1 - math.cos(al) ** 2 - math.cos(be) ** 2 - math.cos(ga) ** 2
                              + 2 * math.cos(al) * math.cos(be) * math.cos(ga))
    assert cell.volume == pytest.approx(v, rel=1e-9)


def test_matthews_definition(tup):
    m = an.matthews(tup)
    assert m["z_sym"] == 4
    assert m["vm"] == pytest.approx(m["volume"] / (4 * m["mass"]))
    assert m["solvent"] == pytest.approx(100 * (1 - 1.230 / m["vm"]))
    assert m["vm_reported"] == pytest.approx(2.64)            # reported value is kept, not overwritten


def test_packing_identity_operator_reproduces_coordinates(tup):
    pack = an.crystal_packing(tup, max_points=10**6)
    assert [c["op"] for c in pack["copies"]][0] == "x,y,z" and len(pack["copies"]) == 4
    a = tup.atoms
    ca = a[a["atom"].isin(["CA", "P"])][["x", "y", "z"]].to_numpy()
    shift = pack["copies"][0]["xyz"] - ca
    # identity copy may be moved by whole lattice vectors only
    assert np.allclose(shift, shift[0], atol=1e-6)


# --------------------------------------------------------------------------- row -> atoms
def test_struct_conn_distances_match_file(dic, tup):
    a = tup.atoms
    n = 0
    for r in tup.rows("struct_conn"):
        if r["conn_type_id"] != "metalc":
            continue
        res = resolve_row(tup, "struct_conn", r, dic)
        p = [x for x in res.partners if x.atom_level]
        assert len(p) == 2 and all(x.mask.sum() == 1 for x in p)
        d = np.linalg.norm(a.loc[p[0].mask, ["x", "y", "z"]].to_numpy()[0] - a.loc[p[1].mask, ["x", "y", "z"]].to_numpy()[0])
        assert d == pytest.approx(float(r["pdbx_dist_value"]), abs=1e-3)
        n += 1
    assert n == 12


def test_conn_angles_match_file(dic, tup):
    a = tup.atoms
    for r in tup.rows("pdbx_struct_conn_angle"):
        res = resolve_row(tup, "pdbx_struct_conn_angle", r, dic)
        pts = [a.loc[p.mask, ["x", "y", "z"]].to_numpy()[0] for p in res.partners if p.atom_level]
        assert len(pts) == 3
        v1, v2 = pts[0] - pts[1], pts[2] - pts[1]
        ang = math.degrees(math.acos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
        assert ang == pytest.approx(float(r["value"]), abs=0.05)


def test_every_category_resolves_without_error(dic, tup):
    for c in tup.category_names():
        for r in tup.rows(c)[:25]:
            resolve_row(tup, c, r, dic)


def test_helix_ranges_cover_declared_residues(dic, tup):
    for r in tup.rows("struct_conf"):
        res = resolve_row(tup, "struct_conf", r, dic)
        n_res = tup.atoms[res.mask].groupby(["asym", "seq"]).ngroups
        assert n_res == int(r["pdbx_PDB_helix_length"])


def test_residue_references(tup):
    a = tup.atoms
    sg = a[(a["auth_asym"] == "A") & (a["auth_seq"] == "176") & (a["atom"] == "SG")].iloc[0]
    refs = residue_references(tup, sg["asym"], "A", int(sg["seq"]), "176", "", ["struct_conn", "struct_site_gen"])
    labels = {(c, lab) for c, _, lab in refs}
    assert ("struct_conn", "metalc1") in labels
    assert any(c == "struct_site_gen" for c, _ in labels)


# --------------------------------------------------------------------------- viewer
def test_view_index_maps_back_to_atom_site(dic, tup):
    v = build_view(tup, show_water=False)
    assert v.meta["n_atoms"] == 5828 - 384
    for k in (0, 100, len(v.index) - 1):
        atom = v.model_data["atoms"][k]
        row = tup.atoms.iloc[v.index[k]]
        assert atom["serial"] == k and atom["name"] == row["atom"] and atom["chain"] == row["auth_asym"]
    m = tup.atoms["comp"].to_numpy() == "ZN"
    assert len(to_view_indices(v.index, m)) == 3


def test_bonds_exclude_metals_and_are_reasonable(tup):
    v = build_view(tup)
    atoms = v.model_data["atoms"]
    for b in v.model_data["bonds"]:
        a1, a2 = atoms[b["atom1_index"]], atoms[b["atom2_index"]]
        assert "Zn" not in (a1["elem"], a2["elem"])
        assert np.linalg.norm(np.subtract(a1["positions"], a2["positions"])) < 2.3
    # a protein residue has ~ as many bonds as atoms; the whole model should too
    assert 0.9 < len(v.model_data["bonds"]) / len(atoms) < 1.15


# --------------------------------------------------------------------------- misc
@pytest.mark.parametrize("text,expected", [("1tup", ("pdb", "1TUP")), ("pdb_00001tup", ("pdb", "1TUP")),
                                           ("P04637", ("afdb", "P04637")), ("AF-P69905-F1", ("afdb", "P69905"))])
def test_normalize_id(text, expected):
    assert normalize_id(text) == expected


@pytest.mark.parametrize("bad", ["", "../etc", "12345", "1TU P"])
def test_normalize_id_rejects(bad):
    with pytest.raises(EntryError):
        normalize_id(bad)


def test_sections():
    assert section_of("exptl_crystal_grow") == "crystal"
    assert section_of("exptl") == "experiment"
    assert section_of("pdbx_struct_assembly_gen") == "assembly"
    assert section_of("struct_conn") == "struct"
    assert section_of("atom_sites") == "coords"
    assert section_of("refln") == "experiment"


def test_table_filter():
    import pandas as pd
    from mmcifxplorer.callbacks import apply_filter
    df = pd.DataFrame({"x": ["1.5", "2.5", "?"], "s": ["CYS", "HIS", "ZN"]})
    assert apply_filter(df, "{x} gt 2").index.tolist() == [1]
    assert apply_filter(df, '{s} contains "y"').index.tolist() == [0]
    assert apply_filter(df, "{s} = ZN").index.tolist() == [2]
