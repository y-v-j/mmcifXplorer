---
title: mmcifXplorer — User Manual
aliases: [mmcifXplorer manual, mmcifXplorer]
tags: [mmcif, pdb, structural-biology, dash, user-manual]
created: 2026-10-01
---

# mmcifXplorer — User Manual

> [!abstract] In one sentence
> mmcifXplorer loads a structure from the Protein Data Bank, shows it in 3D, and lists every PDBx/mmCIF data field of
> that entry. Click a field to read its official definition and see which atoms it describes.

**Contents:** [[#What the application does]] · [[#Installation]] · [[#Starting the app]] · [[#Using the app]] ·
[[#Tips]] · [[#Troubleshooting]]

---

## What the application does

PDBx/mmCIF is the standard file format of the Protein Data Bank. Each value in a file is labelled
`_category.item`, for example `_refine.ls_d_res_high` (resolution) or `_atom_site.B_iso_or_equiv` (B-factor).
mmcifXplorer helps you learn the format and the structure at the same time:

- **3D structure** of any PDB entry, or an AlphaFold model from a UniProt accession.
- **Every category and field in the entry**, with its definition from the wwPDB dictionary (PDBx v5.417,
  plus ModelCIF for computed models).
- **Value checks**: is the value the right type, an allowed option, and within the permitted range?
- **Plain-language interpretation** of key numbers: resolution, R-factors, B-factors, cell, space group,
  Matthews coefficient, wavelength.
- **Fields linked to atoms**: selecting a table row highlights its atoms in 3D. Reported distances and angles are
  recalculated from the coordinates.
- **Crystallographic and quality plots**, **composition** views, a **dictionary search**, **lessons** and
  **guided tours**.

---

## Installation

> [!info] Requirements
> Python ≥ 3.11 and internet access the first time (to download the dictionaries and entries).
> The bundled example **1TUP** also works offline.

**Option A — use the existing `biopyc` conda environment** (it already has every package):

```bash
conda activate biopyc
```

**Option B — create a new environment** from the file in the project folder:

```bash
cd mmcifXplorer
conda env create -f environment.yml
conda activate mmcifxplorer
```

Packages used: `dash`, `dash-bio`, `dash-bootstrap-components`, `plotly`, `gemmi`, `numpy`, `pandas`, `scipy`,
`requests` (and `pytest` for the tests).

---

## Starting the app

```bash
cd mmcifXplorer
python app.py                 # then open http://127.0.0.1:8050 in a browser
```

| Option | Meaning |
|---|---|
| `--port 8060` | use another port |
| `--host 0.0.0.0` | allow other computers on the network to connect |
| `--debug` | developer mode (auto-reload, error pop-ups) |

> [!note] First start
> The first start downloads the two dictionaries (~1.2 MB) into `datasets/dictionary/` and builds a cache,
> which takes a few seconds. Entries you load are cached in `datasets/entries/`, so they open instantly next time.

Stop the app with `Ctrl+C` in the terminal.

---

## Using the app

### 1. Load a structure

In the header:

- Type an ID and press **Enter** or **Load**. Accepted forms:
  - PDB ID: `1TUP`
  - extended PDB ID: `pdb_00001tup`
  - UniProt accession, which loads the AlphaFold model: `P69905`
- Or pick one of the **Examples** (X-ray, NMR ensemble, cryo-EM, DNA, AlphaFold).
- Or click **open local .cif / .cif.gz** to load your own file.

The **summary tiles** below the header show method, resolution, R-factors, cell, atom count and more.
Click a tile to read the definition of the field behind it.

### 2. The 3D viewer (left)

| Control | What it does |
|---|---|
| **Cartoon / Sticks / Spheres** | representation (cartoon also shows ligands as sticks, ions as spheres) |
| **Colour menu** | molecule kind, chain, entity, secondary structure (from the file), B-factor/pLDDT, occupancy, element, residue property |
| **model** | choose a model of an NMR ensemble |
| **water / unit cell** | show water molecules; draw the crystal unit cell |
| **Clear / Reset view** | remove the highlight; re-centre the molecule |

- **Rotate** by dragging, **zoom** with the mouse wheel.
- **Click an atom** to open the **Picked atom** tab. It shows the atom's full `_atom_site` row, its archive (`label_*`)
  and author (`auth_*`) identifiers side by side, and the other rows in the file that mention this residue.
- **Click 2 atoms** to measure a distance, or **3 atoms** for an angle.

### 3. Fields of this entry (right, first tab)

1. Pick a **category** from the list on the left of the tab. Categories are grouped by topic, and **◆** marks those
   whose rows can be shown in 3D. Use the filter box to find a category or item.
2. The panel shows the category description, its key, the **item chips** and the data **table**:
   - Click any **cell**: the field is explained in the *Field explanation* panel (under the viewer), and the row's
     atoms are highlighted in 3D.
   - Hover a column header to see a short definition.
   - Type in the grey row under the headers to filter (`A`, `>2.5`, `contains HOH`). Click a header to sort.
   - Single-row categories (e.g. `_refine`, `_cell`) are shown as an *item / value / definition* list.
3. Open **Raw CIF text** to see how the category is written in the file, or **Dictionary examples** for
   reference examples.

> [!example] Try this
> Load **1TUP** → category `_struct_conn` → click the first value of `pdbx_dist_value`.
> The zinc–cysteine bond is drawn in 3D and its distance is recalculated: *2.326 Å ✓ agrees*.

### 4. Field explanation (under the viewer)

For the selected field you see:
- the definition, data type, units and permitted range;
- the value in this entry, with ✓/✗ checks against the type, allowed values and range;
- an interpretation of the value, where one applies;
- the allowed values (the current one highlighted);
- the fields it **points to** and the fields that **refer to it**. These are clickable, so you can follow links
  through the dictionary.

### 5. Other tabs (right)

| Tab | Contents |
|---|---|
| **Crystal & quality** | key numbers; Matthews coefficient reported vs recalculated; crystal packing in 3D; per-residue B-factor; Ramachandran plot; NMR ensemble spread; per-shell statistics. Click a residue or point to show it in 3D. Use the **Chain** and **Ramachandran class** menus to filter. |
| **Composition** | molecule sunburst and entity table (click to highlight); sequence viewer with modelled residues underlined (select a stretch to highlight it); the largest categories (click to open one). |
| **Dictionary** | search all dictionary fields by name or definition; filter by category group or *in this entry*; click a result to read it. |
| **Learn & tours** | six short lessons on the format and crystallography, and three **guided tours**. A tour bar appears at the top: use **Next ▶ / ◀ Prev / Exit**. |

---

## Tips

> [!tip]
> - `?` in a value means *unknown*; `.` means *not applicable* (e.g. no sequence number for a ligand).
> - `label_*` identifiers are assigned by the archive; `auth_*` are the authors' (the "chain A, residue 176" of papers).
> - For AlphaFold models the B-factor column holds **pLDDT** (confidence), and the viewer colours by it automatically.
> - Entries with more than 40,000 atoms are drawn as a backbone trace; tables and plots still use every atom.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| "Could not download …" | Check the ID and your internet connection. 1TUP works offline. |
| Nothing highlighted after clicking a row | The row may describe residues that are not in the model (e.g. disordered loops), or atoms hidden by the water/model settings. Read the note under the viewer. |
| Viewer has the wrong size after resizing the window | Resize once more or reload the page. |
| Port already in use | Start with `python app.py --port 8060`. |

---

*See also the project `README.md` for technical details and verification results.*
