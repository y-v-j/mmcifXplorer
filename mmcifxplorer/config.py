"""Paths and limits shared by the whole application."""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Every downloaded file lives under datasets/ (dictionary + entry cache).
DATASETS = Path(os.environ.get("MMCIFXPLORER_DATASETS", ROOT / "datasets"))
DICT_DIR = DATASETS / "dictionary"
ENTRY_DIR = DATASETS / "entries"

DICT_URL = "https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic.gz"
DICT_FILE = DICT_DIR / "mmcif_pdbx_v50.dic.gz"
# ModelCIF (computed models, ma_* categories) — merged in for items PDBx lacks.
EXT_DICT_URL = "https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_ma.dic"
EXT_DICT_FILE = DICT_DIR / "mmcif_ma.dic.gz"

RCSB_URL = "https://files.rcsb.org/download/{id}.cif.gz"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

# Offline fallback shipped with the project.
BUNDLED_ENTRIES = {"1TUP": ROOT / "1tup.cif"}

HTTP_TIMEOUT = 30            # seconds
MAX_DOWNLOAD_BYTES = 150e6   # refuse anything larger (compressed size)
MAX_UPLOAD_BYTES = 60e6

# The 3D viewer receives one JSON record per atom.  Beyond this many atoms
# (after removing water) only the backbone trace (CA / P) and ligands are sent.
MAX_VIEW_ATOMS = 40_000
ENTRY_CACHE_SIZE = 4         # parsed entries kept in memory (LRU)
TABLE_PAGE_SIZE = 12
