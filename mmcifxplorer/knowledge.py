"""Curated teaching material layered on top of the formal dictionary.

The dictionary says *what* an item is; the notes here say *why it matters*
and, where possible, interpret the value found in the loaded entry (e.g. a
B-factor turned into an RMS displacement, the Matthews coefficient
recomputed from the cell and the declared content).  Sources: International
Tables for Crystallography Vol. G ch. 3.6 and 4.5 (Fitzgerald et al.), the
PDB exchange dictionary appendix (Westbrook et al.), Westbrook et al. (2022)
J. Mol. Biol. 434:167599, and the wwPDB validation-report documentation.
"""
from __future__ import annotations

import math

# --------------------------------------------------------------------------- #
# Topical sections used to organise the categories of an entry
# --------------------------------------------------------------------------- #
SECTIONS = [
    ("coords", "Atoms & coordinates", "⚛", ["atom_site", "atom_type", "atom_sites"]),
    ("entry", "Entry, history & citation", "🗂", ["entry", "audit", "database", "pdbx_database", "pdbx_audit",
                                                  "citation", "software", "pdbx_deposit", "pdbx_contact",
                                                  "pdbx_initial_refinement", "computing", "journal"]),
    ("crystal", "Crystal, cell & symmetry", "💎", ["cell", "symmetry", "space_group", "exptl_crystal"]),
    ("experiment", "Experiment & data collection", "🔬", ["exptl", "diffrn", "reflns", "refln", "em_", "pdbx_nmr",
                                                         "pdbx_serial", "pdbx_diffrn", "pdbx_soln_scatter",
                                                         "pdbx_phasing", "phasing"]),
    ("refine", "Refinement", "📐", ["refine", "pdbx_refine", "pdbx_reflns"]),
    ("chem", "Molecules: entities & chemistry", "🧪", ["entity", "pdbx_entity", "chem_comp", "pdbx_chem_comp",
                                                       "pdbx_poly_seq_scheme", "pdbx_nonpoly_scheme",
                                                       "pdbx_branch_scheme", "struct_ref", "pdbx_struct_mod_residue",
                                                       "pdbx_molecule", "pdbx_entity_branch"]),
    ("struct", "Structure: chains, folds & contacts", "🧬", ["struct", "pdbx_struct", "ndb_", "geom"]),
    ("assembly", "Assemblies & symmetry operators", "🔁", ["pdbx_struct_assembly", "pdbx_struct_oper_list",
                                                          "struct_ncs", "struct_biol"]),
    ("validate", "Validation & completeness", "✅", ["pdbx_validate", "pdbx_unobs_or_zero_occ"]),
    ("model", "Computed-model quality", "🤖", ["ma_"]),
]
SECTION_ORDER = [s[0] for s in SECTIONS] + ["other"]
SECTION_LABEL = {s[0]: s[1] for s in SECTIONS} | {"other": "Other categories"}
SECTION_ICON = {s[0]: s[2] for s in SECTIONS} | {"other": "•"}


def section_of(category: str) -> str:
    """Longest matching prefix wins, so 'exptl_crystal' is crystal, not experiment."""
    best, best_len = "other", -1
    for sid, _, _, prefixes in SECTIONS:
        for p in prefixes:
            hit = category == p or category.startswith(p if p.endswith("_") else p + "_")
            if hit and len(p) > best_len:
                best, best_len = sid, len(p)
    return best


# --------------------------------------------------------------------------- #
# Curated notes: why an item/category matters
# --------------------------------------------------------------------------- #
CATEGORY_NOTES = {
    "atom_site": "The heart of the file: one row per atom (per model). Coordinates, occupancy and B-factor sit next "
                 "to **two parallel identifier systems** — `label_*` (archive-assigned, tied to the dictionary's "
                 "entity/struct_asym bookkeeping) and `auth_*` (what the authors and most papers use). "
                 "Click any atom in the 3D view to see its row.",
    "cell": "The repeating box of the crystal lattice. Three edge lengths (Å) and three angles (°) define it; "
            "its volume and the space group decide how many copies of the asymmetric unit pack into it.",
    "symmetry": "The space group: the set of rotations/screw axes/translations that generate the whole crystal "
                "from the asymmetric unit. Biological macromolecules are chiral, so only the 65 Sohncke groups "
                "(no mirrors or inversion centres) occur.",
    "exptl": "How the structure was determined (X-ray, NMR, cryo-EM, neutron, ...). An entry can list more than one "
             "method (e.g. joint X-ray/neutron refinement).",
    "exptl_crystal": "Properties of the crystal itself, including the Matthews coefficient and solvent content — "
                     "quick sanity checks that the declared content fits in the cell.",
    "exptl_crystal_grow": "Crystallisation conditions (method, pH, temperature, reservoir). Useful when reproducing "
                          "or improving crystals.",
    "reflns": "Statistics of the merged diffraction data (before model building): resolution limits, completeness, "
              "multiplicity, ⟨I/σ(I)⟩, R-merge, CC½, Wilson B.",
    "reflns_shell": "The same data-quality statistics split into resolution shells; the highest-resolution shell "
                    "shows where the data run out.",
    "refine": "Model-vs-data agreement after refinement: resolution used, R-work, R-free, number of reflections, "
              "overall B and coordinate-error estimates.",
    "refine_ls_shell": "R-factors and completeness of the refinement per resolution shell.",
    "refine_ls_restr": "Deviations from ideal stereochemistry (bonds, angles, planes ...) reported by the "
                       "refinement program, as RMS deviations.",
    "refine_hist": "Counts of protein / nucleic / ligand / solvent atoms used in the final refinement cycle.",
    "entity": "Each chemically distinct molecule is described **once**: a polymer sequence, a ligand, water. "
              "`pdbx_number_of_molecules` says how many copies are in the deposited model.",
    "entity_poly": "Sequence-level description of polymer entities, including the one-letter sequence (with "
                   "modified residues in parentheses) and the author chain IDs that carry it.",
    "entity_poly_seq": "The full sample sequence, residue by residue — including residues that were never seen "
                       "in the density. `num` here *is* `label_seq_id` in `_atom_site`.",
    "struct_asym": "One row per molecule *instance* in the model (each protein chain, each ligand copy, each "
                   "water group). Its `id` is what `label_asym_id` points to.",
    "pdbx_poly_seq_scheme": "The Rosetta stone between numbering systems: entity sequence position ↔ label_seq_id "
                            "↔ author residue number/chain, with '?' for residues not modelled.",
    "pdbx_nonpoly_scheme": "The same mapping for ligands, ions and waters.",
    "chem_comp": "Every chemical component (amino acid, nucleotide, ligand) used in the entry, with formula and "
                 "weight. Full definitions live in the wwPDB Chemical Component Dictionary.",
    "struct_conf": "Helices (and turns) as annotated in the entry, residue ranges in both numbering schemes.",
    "struct_sheet_range": "β-strands, grouped into sheets.",
    "struct_conn": "Explicit connections that are not ordinary polymer bonds: disulfides, covalent links to "
                   "ligands, metal coordination and (for nucleic acids) base-pair hydrogen bonds. Selecting a row "
                   "draws the link and checks the reported distance against the coordinates.",
    "struct_site": "Binding/active sites, often generated from ligand environments.",
    "struct_site_gen": "The residues that make up each site.",
    "pdbx_struct_assembly": "Biological assemblies: the oligomer believed to be functional, which may be bigger or "
                            "smaller than the asymmetric unit.",
    "pdbx_struct_assembly_gen": "Recipe for building an assembly: which chains (asym IDs) to transform with which "
                                "operators.",
    "pdbx_struct_oper_list": "The rotation matrices + translation vectors referenced by assembly recipes.",
    "struct_ncs_oper": "Non-crystallographic symmetry: operators relating copies inside the asymmetric unit.",
    "pdbx_validate_torsion": "Residues with unusual backbone (φ, ψ) found during annotation.",
    "pdbx_validate_close_contact": "Atom pairs closer than expected inside the asymmetric unit (possible clashes).",
    "pdbx_validate_symm_contact": "Close contacts between the model and a symmetry-related copy.",
    "pdbx_validate_rmsd_bond": "Covalent bonds whose length deviates strongly from the ideal value.",
    "pdbx_validate_rmsd_angle": "Bond angles that deviate strongly from the ideal value.",
    "pdbx_unobs_or_zero_occ_residues": "Residues present in the sample but missing from the model (disordered).",
    "pdbx_database_status": "Archive bookkeeping: deposition date, release status, processing site.",
    "pdbx_audit_revision_history": "Every public revision of the entry (remediation, corrections ...).",
    "citation": "The primary publication (id 'primary') and other references.",
    "software": "Programs used for data collection, processing, phasing and refinement.",
    "ma_qa_metric_local": "Per-residue model-confidence scores of a computed model (pLDDT for AlphaFold).",
    "em_3d_reconstruction": "Cryo-EM map reconstruction: resolution (usually by the FSC 0.143 criterion), "
                            "number of particles, symmetry imposed.",
    "pdbx_nmr_ensemble": "How the NMR ensemble was selected and how many conformers were deposited.",
}

ITEM_NOTES = {
    "_atom_site.group_pdb": "Legacy PDB record type: `ATOM` for standard polymer residues, `HETATM` for everything "
                            "else (ligands, water, modified residues).",
    "_atom_site.label_asym_id": "Archive chain ID → `_struct_asym.id`. Every ligand copy and every water group gets "
                                "its own label_asym_id, unlike author chains.",
    "_atom_site.auth_asym_id": "Author chain ID — the 'chain A' of papers and the legacy PDB format. Ligands and "
                               "waters are assigned to the author chain they sit next to.",
    "_atom_site.label_seq_id": "Position in `_entity_poly_seq` (1…N, no gaps, no insertion codes). '.' for "
                               "non-polymers. Use it to align coordinates to the sequence.",
    "_atom_site.auth_seq_id": "Author residue number. Can start anywhere, skip numbers, go negative, or repeat "
                              "with insertion codes (e.g. antibody Kabat numbering 52, 52A, 52B).",
    "_atom_site.pdbx_pdb_ins_code": "Insertion code that disambiguates repeated author residue numbers.",
    "_atom_site.label_atom_id": "Atom name, defined per chemical component in the Chemical Component Dictionary "
                                "(e.g. CA = alpha carbon, OG1 = threonine hydroxyl oxygen).",
    "_atom_site.label_alt_id": "Alternate-conformation flag. Atoms with 'A'/'B' are alternative positions; their "
                               "occupancies add up to 1.",
    "_atom_site.cartn_x": "Orthogonal coordinate in Å. Fractional coordinates are obtained with "
                          "`_atom_sites.fract_transf_matrix`.",
    "_atom_site.occupancy": "Fraction of unit cells (or molecules) in which this atom sits at this position.",
    "_atom_site.b_iso_or_equiv": r"Atomic displacement parameter, $B = 8\pi^2\langle u^2\rangle$. It mixes "
                                 "thermal vibration with static disorder; high B means poorly ordered.",
    "_atom_site.pdbx_pdb_model_num": "Model number. NMR ensembles and multi-state models have several.",
    "_refine.ls_d_res_high": "The resolution of the data used in refinement — the single most important "
                             "indicator of how much detail the model can be trusted to show.",
    "_refine.ls_r_factor_r_work": r"Agreement between observed and model structure-factor amplitudes for the "
                                  r"reflections used in refinement: $R=\frac{\sum\,||F_o|-|F_c||}{\sum |F_o|}$.",
    "_refine.ls_r_factor_r_free": "The same R-factor computed on a test set (typically 5 %) never used in "
                                  "refinement (Brünger, 1992). It guards against over-fitting.",
    "_refine.ls_r_factor_obs": "R-factor over all observed reflections (work + test).",
    "_refine.b_iso_mean": "Average B-factor of the model.",
    "_refine.overall_su_ml": "Maximum-likelihood estimate of the overall coordinate error (Å).",
    "_refine.pdbx_ls_cross_valid_method": "How the test set for R-free was chosen (THROUGHOUT = random across the "
                                          "whole resolution range).",
    "_reflns.d_resolution_high": "High-resolution limit of the merged data.",
    "_reflns.percent_possible_obs": "Completeness: the fraction of the theoretically measurable reflections that "
                                    "were actually measured.",
    "_reflns.pdbx_redundancy": "Multiplicity: average number of measurements per unique reflection.",
    "_reflns.pdbx_rmerge_i_obs": r"$R_{merge}=\frac{\sum_{hkl}\sum_i|I_i-\bar I|}{\sum_{hkl}\sum_i I_i}$. "
                                 "Grows with multiplicity, so CC½ and R-meas/R-pim are better modern measures.",
    "_reflns.pdbx_netI_over_sigmaI": "Mean signal-to-noise ⟨I/σ(I)⟩ of the merged data.",
    "_reflns.pdbx_cc_half": "Correlation between two random half data sets; the modern resolution-cutoff criterion.",
    "_reflns.b_iso_wilson_estimate": "Overall B estimated from how fast intensity falls off with resolution "
                                     "(Wilson plot), independent of the model.",
    "_cell.length_a": "Unit-cell edge a (Å).",
    "_cell.angle_beta": "Unit-cell angle β between edges a and c.",
    "_cell.z_pdb": "Number of polymer chains in the unit cell (= chains in the ASU × symmetry operators).",
    "_symmetry.space_group_name_h-m": "Hermann–Mauguin symbol of the space group, e.g. P 21 21 21 or C 1 2 1.",
    "_symmetry.int_tables_number": "Space-group number in International Tables Vol. A (1–230).",
    "_exptl_crystal.density_matthews": r"Matthews coefficient $V_M = V_{cell}/(Z\cdot M)$ in Å³/Da. Matthews "
                                       "(1968) found 1.68–3.53 for proteins.",
    "_exptl_crystal.density_percent_sol": r"Solvent fraction $1 - 1.230/V_M$ (×100), assuming the density of protein.",
    "_exptl.method": "Experimental technique used to determine the structure.",
    "_diffrn.ambient_temp": "Data-collection temperature in kelvin; ~100 K means a cryo-cooled crystal.",
    "_diffrn_radiation_wavelength.wavelength": r"X-ray wavelength λ (Å). Photon energy: $E\,[\mathrm{keV}] = 12.398/\lambda\,[\mathrm{Å}]$.",
    "_diffrn_source.source": "X-ray source: rotating anode / sealed tube (home source), synchrotron, XFEL.",
    "_entity.type": "polymer, non-polymer, branched (oligosaccharide), macrolide or water.",
    "_entity.pdbx_number_of_molecules": "Copies of this entity in the deposited model.",
    "_entity.formula_weight": "Molecular weight (Da) of one copy of the entity, from the sample sequence/formula.",
    "_entity_poly.pdbx_seq_one_letter_code": "Sample sequence; non-standard residues are written in parentheses, "
                                             "e.g. (MSE) for selenomethionine.",
    "_entity_poly.pdbx_seq_one_letter_code_can": "Canonical sequence: modified residues replaced by their parent.",
    "_entity_poly.pdbx_strand_id": "Author chain IDs that carry this polymer entity.",
    "_struct_conn.conn_type_id": "covale (covalent), disulf (disulfide), metalc (metal coordination), "
                                 "hydrog (hydrogen bond, mostly base pairs).",
    "_struct_conn.pdbx_dist_value": "Distance between the two partners in Å, as annotated.",
    "_struct_conn.ptnr1_symmetry": "Symmetry operator applied to partner 1, written n_klm: operator n of the "
                                   "space group plus lattice translations k-5, l-5, m-5. 1_555 = identity.",
    "_struct_conf.conf_type_id": "HELX_P = helix (protein), TURN_P = turn. The helix class is in "
                                 "pdbx_PDB_helix_class (1 = right-handed α, 5 = 3₁₀).",
    "_pdbx_struct_assembly.oligomeric_details": "Oligomeric state of the biological assembly (monomeric, dimeric ...).",
    "_pdbx_struct_assembly.method_details": "How the assembly was derived (author-defined, PISA, ...).",
    "_pdbx_struct_assembly_gen.oper_expression": "Operators to apply: '1' identity; '(1-60)' a range; "
                                                 "'(1,2)(3,4)' the product of two lists.",
    "_pdbx_database_status.recvd_initial_deposition_date": "Date the entry was first deposited.",
    "_ma_qa_metric_local.metric_value": "Per-residue confidence. For AlphaFold pLDDT: >90 very high, 70–90 "
                                        "confident, 50–70 low, <50 very low (often disordered).",
    "_em_3d_reconstruction.resolution": "Map resolution, normally where the half-map Fourier shell correlation "
                                        "drops below 0.143.",
}


# --------------------------------------------------------------------------- #
# Live interpretation of a value found in the entry
# --------------------------------------------------------------------------- #
def _f(v):
    try:
        return float(str(v).split("(")[0])
    except (TypeError, ValueError):
        return None


def _resolution_text(d: float) -> str:
    if d <= 1.2:
        s = "atomic resolution: individual atoms resolved, many hydrogens visible, anisotropic B refinable"
    elif d <= 1.8:
        s = "high resolution: side chains and waters well defined, alternate conformations visible"
    elif d <= 2.5:
        s = "medium resolution: side chains mostly defined, ordered waters visible"
    elif d <= 3.5:
        s = "low-to-medium resolution: backbone reliable, side-chain placement partly inferred"
    else:
        s = "low resolution: overall fold and domain placement; atomic details are model-dependent"
    return f"**{d:.2f} Å** — {s}."


def interpret(item: str, value: str, entry=None) -> str | None:
    """Markdown interpreting one value; None when nothing useful can be said."""
    key = item.lower()
    x = _f(value)
    if key in ("_refine.ls_d_res_high", "_reflns.d_resolution_high", "_em_3d_reconstruction.resolution") and x:
        return _resolution_text(x)
    if key == "_refine.ls_r_factor_r_free" and x is not None:
        out = f"R-free = **{x:.3f}**."
        if entry is not None:
            rw = entry.number("refine.ls_R_factor_R_work") or entry.number("refine.ls_R_factor_obs")
            d = entry.number("refine.ls_d_res_high")
            if rw is not None:
                gap = x - rw
                out += f" Gap to R-work ({rw:.3f}) = **{gap:.3f}**"
                out += " — normal (< 0.05)." if gap < 0.05 else (
                    " — slightly large." if gap < 0.07 else " — large; possible over-fitting.")
            if d:
                out += f" A common rule of thumb expects R-free ≈ resolution/10 ≈ {d / 10:.2f}."
        return out
    if key in ("_refine.ls_r_factor_r_work", "_refine.ls_r_factor_obs") and x is not None:
        out = f"**{x * 100:.1f} %** of the observed amplitude is not explained by the model."
        if entry is not None and entry.number("refine.ls_R_factor_R_free") is None:
            out += " This entry reports **no R-free** (common before ~1995, when cross-validation was introduced)."
        return out
    if key in ("_atom_site.b_iso_or_equiv", "_refine.b_iso_mean", "_reflns.b_iso_wilson_estimate") and x is not None and x > 0:
        u = math.sqrt(x / (8 * math.pi ** 2))
        out = f"B = {x:.2f} Å² ⇒ RMS displacement √(B/8π²) = **{u:.2f} Å**."
        if entry is not None and entry.is_computed_model and key == "_atom_site.b_iso_or_equiv":
            out = f"In this computed model the B column holds **pLDDT = {x:.1f}**, not a displacement parameter."
        return out
    if key == "_atom_site.occupancy" and x is not None:
        if x >= 0.999:
            return "Fully occupied position."
        if x == 0:
            return "**Zero occupancy**: the atom was placed for chemical completeness but is not supported by data."
        return f"**Partially occupied** ({x * 100:.0f} %): disorder or alternate conformations."
    if key == "_exptl_crystal.density_matthews" and entry is not None:
        from .analysis import matthews
        m = matthews(entry)
        if m:
            return (f"Recomputed from this entry: V = {m['volume']:,.0f} Å³, Z = {m['z_sym']} symmetry operators, "
                    f"ASU mass = {m['mass']:,.0f} Da (entity formula weights × copies) ⇒ "
                    f"**V_M = {m['vm']:.2f} Å³/Da**, solvent ≈ **{m['solvent']:.1f} %**. "
                    f"Reported: {value} Å³/Da. Differences usually come from the mass assumed "
                    f"(e.g. nucleic acid or ligands counted differently, a different construct).")
    if key == "_exptl_crystal.density_percent_sol" and x is not None:
        return f"About **{x:.0f} %** of the crystal volume is solvent channels."
    if key == "_diffrn_radiation_wavelength.wavelength" and x:
        e = 12.398 / x
        hint = " (Cu Kα home source)" if abs(x - 1.5418) < 0.003 else (
            " (Se K-edge region, typical for SeMet SAD/MAD)" if 0.97 < x < 0.985 else "")
        return f"λ = {x:.4f} Å ⇒ photon energy **{e:.2f} keV**{hint}."
    if key == "_diffrn.ambient_temp" and x:
        return ("**Cryo-cooled** data collection (reduces radiation damage)." if x < 150
                else "**Room-temperature** (or near) data collection.")
    if key in ("_cell.length_a", "_cell.length_b", "_cell.length_c", "_cell.angle_alpha", "_cell.angle_beta",
               "_cell.angle_gamma", "_cell.volume") and entry is not None:
        from .analysis import unit_cell
        c = unit_cell(entry)
        if c:
            return (f"Cell a={c.a:.2f}, b={c.b:.2f}, c={c.c:.2f} Å; α={c.alpha:.2f}, β={c.beta:.2f}, γ={c.gamma:.2f}° "
                    f"⇒ volume **{c.volume:,.0f} Å³**.")
    if key in ("_symmetry.space_group_name_h-m",) and value not in ("?", "."):
        import gemmi
        sg = gemmi.find_spacegroup_by_name(value)
        if sg:
            n_ops = len(list(sg.operations()))
            sohncke = sg.is_sohncke() if hasattr(sg, "is_sohncke") else None
            sys = sg.crystal_system_str()
            out = (f"Space group No. **{sg.number}**, {sys} crystal system, **{n_ops} symmetry operators** "
                   f"(copies of the asymmetric unit per cell).")
            if sohncke is not None:
                out += " Sohncke (chiral-compatible) group." if sohncke else " Not a Sohncke group — unusual for macromolecules."
            return out
    if key == "_reflns.percent_possible_obs" and x is not None:
        return f"**{x:.1f} %** complete" + (" — good." if x >= 95 else " — some data missing.")
    if key == "_reflns.pdbx_netI_over_sigmaI" and x is not None:
        return f"⟨I/σ⟩ = **{x:.1f}** over all data."
    if key == "_struct_conn.pdbx_dist_value" and x is not None:
        return ("Typical distances: S–S disulfide ≈ 2.05 Å; Zn–S(Cys) ≈ 2.3 Å; Zn–N(His) ≈ 2.0–2.1 Å; "
                "Mg–O ≈ 2.1 Å; hydrogen bonds 2.6–3.2 Å (donor–acceptor).")
    if key == "_ma_qa_metric_local.metric_value" and x is not None:
        band = ("very high" if x > 90 else "confident" if x > 70 else "low" if x > 50 else "very low")
        return f"pLDDT {x:.1f} — **{band}** confidence."
    return None


# --------------------------------------------------------------------------- #
# Lessons (Learn tab)
# --------------------------------------------------------------------------- #
LESSONS = [
    {
        "id": "anatomy",
        "title": "1 · Anatomy of an mmCIF file",
        "body": r"""
An mmCIF file is plain text made of **data blocks**, one per entry, starting with `data_1TUP`.
Inside, every value is labelled by a **data name** of the form `_category.item`:

```
_cell.length_a      117.500        <- single-row category: name/value pairs
_cell.angle_beta    105.50

loop_                              <- multi-row category: a table
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.pdbx_dist_value
metalc1 metalc 2.326
metalc2 metalc 1.992
```

* A **category** (`cell`, `struct_conn`) is a table; an **item** (`length_a`) is a column.
* `loop_` introduces a table with several rows; single-row categories are usually written as pairs.
* Two special values: **`?`** = unknown/missing, **`.`** = not applicable (e.g. `label_seq_id` of a ligand).
* Strings with spaces are quoted `'...'`; long text sits between lines starting with `;`.

Open any category in the **Fields** tab and press *Raw CIF* to see exactly this syntax for the loaded entry.
""",
    },
    {
        "id": "dictionary",
        "title": "2 · The dictionary behind the names (DDL2)",
        "body": r"""
Every data name is defined in the **PDBx/mmCIF dictionary**, itself written in CIF and governed by the
Dictionary Definition Language **DDL2**. For each item the dictionary records:

| attribute | example |
|---|---|
| definition | "Residual factor R for reflections ... used as the test reflections" |
| type (with a regular expression) | `float`, `code`, `line`, `yyyy-mm-dd` |
| units | `angstroms`, `degrees`, `8pi2_angstroms_squared` |
| enumeration | `covale`, `disulf`, `metalc`, `hydrog` |
| permitted range | R-free in [0, 1] |
| category key | `_atom_site.id` uniquely identifies an atom row |
| parent/child links | `_atom_site.label_entity_id` → `_entity.id` |

The parent/child links turn a flat file into a **relational database**: `_atom_site.label_asym_id` points to
`_struct_asym.id`, which points to `_entity.id`, which is described by `_entity_poly` and `_entity_poly_seq`.
Items added by the PDB beyond the IUCr standard carry the **`pdbx_`** prefix (Westbrook et al., ITC G A3.6.2).

The explanation panel shows all of this for any field you select, and checks the entry's value against the type,
enumeration and range.
""",
    },
    {
        "id": "identifiers",
        "title": "3 · Entity, asym, chain: who is who",
        "body": r"""
mmCIF separates **what** a molecule is from **where** its copies are:

* **entity** (`_entity`) — a chemically distinct molecule, described once (p53 core domain, Zn²⁺, water).
* **struct_asym** (`label_asym_id`) — one *instance* of an entity in the model. Three p53 chains and three Zn ions
  give six label_asym_ids. Waters near each chain form their own asym groups.
* **author chain** (`auth_asym_id`) — the familiar "chain A". A ligand bound to chain A gets auth chain A but its
  own label_asym_id.

Residues are numbered twice as well:

| | label | auth |
|---|---|---|
| chain | `label_asym_id` (archive) | `auth_asym_id` (author) |
| residue number | `label_seq_id` = position in `_entity_poly_seq`, 1…N | `auth_seq_id`, anything the authors chose |
| insertion codes | never | `pdbx_PDB_ins_code` |
| ligands/water | `label_seq_id` = `.` | numbered |

`_pdbx_poly_seq_scheme` maps the two systems residue by residue and shows unobserved residues with `?`.
Click an atom in the 3D view: the atom inspector shows both identities side by side.
""",
    },
    {
        "id": "crystal",
        "title": "4 · Crystallographic parameters",
        "body": r"""
**Unit cell and space group.** The crystal is a lattice of identical boxes (`_cell`, edges a, b, c and angles
α, β, γ). The **asymmetric unit** (ASU — what the file's coordinates describe) is copied into the cell by the
space-group operators (`_symmetry`). Volume of a general cell:

$$V = abc\sqrt{1-\cos^2\alpha-\cos^2\beta-\cos^2\gamma+2\cos\alpha\cos\beta\cos\gamma}$$

**Matthews coefficient** $V_M = V/(Z\,M)$ (Å³/Da, Z = symmetry operators, M = ASU mass) and solvent fraction
$1-1.230/V_M$: proteins crystallise with roughly 30–70 % solvent. The *Crystal & quality* tab recomputes both.

**Resolution** $d_{min}$: the finest spacing of lattice planes whose reflections were used. Smaller is better.

**R-factors.** $R=\sum||F_o|-|F_c||/\sum|F_o|$ measures model–data agreement; **R-free** uses ~5 % of
reflections left out of refinement and is the honest number. A gap R-free − R-work above ~0.05–0.07 suggests
over-fitting.

**B-factor** $B=8\pi^2\langle u^2\rangle$ — how much an atom is smeared out (vibration + disorder). Loops and
surface side chains have high B; the core is low. **Occupancy** < 1 marks alternate conformations or partial binding.

**Data quality**: completeness, multiplicity (redundancy), ⟨I/σ(I)⟩, R-merge and CC½ (`_reflns`, `_reflns_shell`).
""",
    },
    {
        "id": "quality",
        "title": "5 · Judging a model's quality",
        "body": r"""
Global indicators (Crystal & quality tab): resolution, R-work/R-free, mean B, Wilson B.

Local indicators recorded in the file during wwPDB annotation:

* `_pdbx_validate_torsion` — residues with unusual (φ, ψ); shown as red points on the Ramachandran plot.
* `_pdbx_validate_close_contact` / `_pdbx_validate_symm_contact` — clashes within the ASU / with symmetry mates.
* `_pdbx_validate_rmsd_bond`, `_pdbx_validate_rmsd_angle` — geometry far from ideal values.
* `_pdbx_unobs_or_zero_occ_residues` — parts of the sample that are not in the model.

The full wwPDB validation report adds MolProbity clashscore, Ramachandran and rotamer outliers and fit to
electron density (RSRZ), with percentile ranks against the archive. For 1TUP the report counts 194/192/193 analysed
residues in chains A/B/C — the Ramachandran plot here reproduces exactly those counts.

Rule of thumb: never trust a ligand, a loop or a side chain more than its B-factors, occupancy and the
resolution allow.
""",
    },
    {
        "id": "methods",
        "title": "6 · Beyond X-ray: NMR, cryo-EM and computed models",
        "body": r"""
PDBx/mmCIF covers every method in the archive (Westbrook et al., 2022):

* **NMR** — several models (`pdbx_PDB_model_num`), restraints and ensemble details in `_pdbx_nmr_*`. No unit cell
  (a dummy 1 Å cube is written). The *Crystal & quality* tab shows per-residue spread across the ensemble instead.
* **Cryo-EM** — `_em_*` categories: specimen, imaging, particle counts, reconstruction resolution (FSC 0.143).
* **Integrative/hybrid** models (PDB-IHM) and **computed models** (ModelArchive, AlphaFold DB) extend the dictionary
  with `ihm_` and `ma_` (ModelCIF) categories. In AlphaFold files the B-factor column carries **pLDDT**.

Since 2014 mmCIF has been the PDB's master format: legacy PDB files cannot hold more than 62 chains or 99 999
atoms, so large structures exist *only* as mmCIF.
""",
    },
]


# --------------------------------------------------------------------------- #
# Guided tours: each step selects a category/item (and optionally a row)
# --------------------------------------------------------------------------- #
TOURS = {
    "xray": {
        "title": "Anatomy of an X-ray entry",
        "steps": [
            ("entry", "_entry.id", 0, "Every data block names the entry it describes."),
            ("exptl", "_exptl.method", 0, "First question for any structure: how was it determined?"),
            ("cell", "_cell.length_a", 0, "The unit cell — drawn as a box around the model in 3D."),
            ("symmetry", "_symmetry.space_group_name_H-M", 0, "The space group copies the ASU through the cell."),
            ("exptl_crystal", "_exptl_crystal.density_Matthews", 0, "Does the declared content fit the cell?"),
            ("reflns", "_reflns.d_resolution_high", 0, "The merged diffraction data."),
            ("refine", "_refine.ls_d_res_high", 0, "Resolution used in refinement."),
            ("refine", "_refine.ls_R_factor_R_work", 0, "How well the model explains the data."),
            ("entity", "_entity.type", 0, "Which distinct molecules are present?"),
            ("struct_asym", "_struct_asym.id", 0, "...and how many copies of each."),
            ("atom_site", "_atom_site.B_iso_or_equiv", 0, "Finally the atoms — with their B-factors."),
        ],
    },
    "ids": {
        "title": "From sequence to atoms (identifiers)",
        "steps": [
            ("entity_poly", "_entity_poly.pdbx_seq_one_letter_code", 0, "The sample sequence of a polymer entity."),
            ("entity_poly_seq", "_entity_poly_seq.num", 0, "Each position gets a number: this is label_seq_id."),
            ("struct_asym", "_struct_asym.entity_id", 0, "Instances of entities in the model."),
            ("pdbx_poly_seq_scheme", "_pdbx_poly_seq_scheme.auth_seq_num", 0,
             "The mapping between label and author numbering ('?' = not modelled)."),
            ("atom_site", "_atom_site.label_seq_id", 0, "Atoms carry both numbering systems."),
            ("atom_site", "_atom_site.auth_seq_id", 0, "...the author one is what papers cite."),
        ],
    },
    "structure": {
        "title": "Structure annotation & contacts",
        "steps": [
            ("struct_conf", "_struct_conf.conf_type_id", 0, "Helices annotated in the file."),
            ("struct_sheet_range", "_struct_sheet_range.sheet_id", 0, "β-strands grouped into sheets."),
            ("struct_conn", "_struct_conn.conn_type_id", 0, "Explicit links: metals, disulfides, covalent bonds."),
            ("struct_conn", "_struct_conn.pdbx_dist_value", 0, "Reported distance — checked against coordinates."),
            ("struct_site", "_struct_site.details", 0, "Binding sites."),
            ("pdbx_struct_assembly", "_pdbx_struct_assembly.oligomeric_details", 0, "The biological assembly."),
            ("pdbx_validate_torsion", "_pdbx_validate_torsion.phi", 0, "Residues flagged for unusual backbone torsions."),
        ],
    },
}


def cif_markup_to_markdown(text: str) -> str:
    """Convert CIF typographic markup (B~eq~ subscripts, x^2^ superscripts) to HTML tags."""
    import re
    text = re.sub(r"~([^~\n]{1,20})~", r"<sub>\1</sub>", text)
    text = re.sub(r"\^([^\^\n]{1,20})\^", r"<sup>\1</sup>", text)
    return text
