"""Static page layout.  Dynamic content is filled by callbacks."""
from __future__ import annotations

import dash_bio
import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html

from . import config
from .figures import INK, INK_2, MUTED, SURFACE
from .knowledge import LESSONS, TOURS

EXAMPLES = [
    ("1TUP", "p53 core domain bound to DNA · X-ray 2.2 Å (bundled offline)"),
    ("4HHB", "Haemoglobin · X-ray 1.74 Å, haem ligands"),
    ("1CRN", "Crambin · X-ray 1.5 Å, disulfides, alternate conformers"),
    ("2LYZ", "Hen lysozyme · classic X-ray"),
    ("1BNA", "B-DNA dodecamer (Dickerson) · nucleic acid"),
    ("6LU7", "SARS-CoV-2 main protease + covalent inhibitor"),
    ("1D3Z", "Ubiquitin · NMR ensemble (10 models)"),
    ("5A63", "γ-secretase · cryo-EM"),
    ("P69905", "Haemoglobin α · AlphaFold DB model (UniProt)"),
]

COLOR_MODES = [
    {"label": "Molecule kind", "value": "kind"},
    {"label": "Chain (auth_asym_id)", "value": "chain"},
    {"label": "Entity (label_entity_id)", "value": "entity"},
    {"label": "Secondary structure (file)", "value": "ss"},
    {"label": "B-factor / pLDDT", "value": "bfactor"},
    {"label": "Occupancy", "value": "occupancy"},
    {"label": "Element", "value": "element"},
    {"label": "Residue property", "value": "residue"},
]

TABLE_STYLE = dict(
    style_table={"overflowX": "auto"},
    style_header={"backgroundColor": "#232321", "color": INK, "fontWeight": "600", "border": "1px solid #383835",
                  "fontFamily": "ui-monospace, Menlo, monospace", "fontSize": "11.5px", "cursor": "pointer"},
    style_cell={"backgroundColor": SURFACE, "color": INK_2, "border": "1px solid #2c2c2a",
                "fontFamily": "ui-monospace, Menlo, monospace", "fontSize": "11.5px", "textAlign": "left",
                "maxWidth": "260px", "overflow": "hidden", "textOverflow": "ellipsis", "padding": "3px 8px"},
    style_filter={"backgroundColor": "#0f0f0e", "color": INK},
    style_data_conditional=[
        {"if": {"state": "active"}, "backgroundColor": "#14243a", "border": "1px solid #3987e5", "color": INK},
        {"if": {"filter_query": '{_value} = "?"'}, "color": MUTED},
    ],
    css=[{"selector": ".dash-spreadsheet-menu", "rule": "display:none"},
         {"selector": ".previous-next-container", "rule": "float:right"}],
)


def header():
    return html.Div(className="app-header", children=[
        html.Div(className="brand", children=[
            html.H1("mmcifXplorer"),
            html.Div(id="brand-sub", className="sub", children="PDBx/mmCIF dictionary ↔ live 3D structure"),
        ]),
        html.Div(className="load-group", children=[
            dbc.Input(id="pdb-input", value="1TUP", placeholder="PDB ID or UniProt", debounce=False,
                      n_submit=0, size="sm", maxLength=16, autoComplete="off"),
            dbc.Button("Load", id="btn-load", color="primary", size="sm"),
            dbc.DropdownMenu(label="Examples", size="sm", color="secondary", children=[
                dbc.DropdownMenuItem([html.Span(pid, className="mono me-2", style={"color": INK}),
                                      html.Small(desc, style={"color": MUTED})],
                                     id={"type": "example", "id": pid}) for pid, desc in EXAMPLES
            ]),
            dcc.Upload(id="upload", className="upload-box", multiple=False,
                       children=html.Span(["⤒ open local .cif / .cif.gz"])),
        ]),
        html.Div(id="status", className="status", children=[html.Span(className="dot"), "ready"]),
    ])


def viewer_panel():
    return html.Div(className="panel", children=[
        html.Div(className="panel-head", children=[
            dbc.RadioItems(id="rep", className="btn-group", inputClassName="btn-check",
                           labelClassName="btn btn-outline-secondary btn-sm", labelCheckedClassName="active",
                           value="cartoon", options=[{"label": "Cartoon", "value": "cartoon"},
                                                     {"label": "Sticks", "value": "stick"},
                                                     {"label": "Spheres", "value": "sphere"}]),
            dbc.Select(id="color-mode", options=COLOR_MODES, value="kind", size="sm", style={"width": "190px"}),
            dbc.Select(id="model-select", options=[{"label": "model 1", "value": "1"}], value="1", size="sm",
                       style={"width": "105px"}),
            dbc.Switch(id="show-water", label="water", value=False, className="mb-0 small"),
            dbc.Switch(id="show-cell", label="unit cell", value=False, className="mb-0 small"),
            dbc.Button("Clear", id="btn-clear", size="sm", className="btn-ghost ms-auto"),
            dbc.Button("Reset view", id="btn-reset", size="sm", className="btn-ghost"),
        ]),
        html.Div(id="viewer-wrap", children=[
            dcc.Loading(type="circle", color="#3987e5", children=dash_bio.Molecule3dViewer(
                id="mol3d", modelData={"atoms": [], "bonds": []}, styles=[], backgroundColor=SURFACE,
                backgroundOpacity=1, selectionType="atom", selectedAtomIds=[], width=600, height=480,
                zoom={"factor": 0.95, "animationDuration": 0, "fixedPath": False},
                zoomTo={"sel": {}, "animationDuration": 0, "fixedPath": False}, labels=[], shapes=[])),
        ]),
        html.Div(id="viewer-legend", className="viewer-legend"),
        html.Div(id="hl-banner", className="hl-banner", style={"display": "none"}),
        html.Div("Click an atom to inspect its _atom_site row (click two or three atoms to measure a distance "
                 "or angle). Select a table row or field to see it in 3D.", className="viewer-hint"),
    ])


def left_tabs():
    return dbc.Tabs(id="left-tabs", active_tab="explain", className="mt-2", children=[
        dbc.Tab(label="Field explanation", tab_id="explain", children=html.Div(
            id="explain-panel", className="panel panel-body explain mt-2",
            children=html.Div("Select a category, a table column or a field chip to read its definition.",
                              style={"color": MUTED}))),
        dbc.Tab(label="Picked atom", tab_id="atom", children=html.Div(
            id="atom-panel", className="panel panel-body mt-2",
            children=html.Div("Click an atom in the 3D view.", style={"color": MUTED}))),
    ])


def fields_tab():
    return dbc.Row(className="g-2 mt-1", children=[
        dbc.Col(md=4, children=html.Div(className="panel panel-body", children=[
            dbc.Input(id="cat-search", placeholder="filter categories / items…", size="sm", debounce=True,
                      className="mb-2", style={"background": "#232321", "color": INK, "borderColor": "#383835"}),
            html.Div("◆ = can be shown in 3D", className="lbl-sm mb-1"),
            html.Div(id="cat-nav", className="cat-nav"),
        ])),
        dbc.Col(md=8, children=html.Div(className="panel panel-body", children=[
            html.Div(id="cat-header", className="cat-header"),
            html.Div(id="cat-chips", className="chips"),
            dash_table.DataTable(id="cat-table", columns=[], data=[], page_action="custom", page_current=0,
                                 page_size=config.TABLE_PAGE_SIZE, filter_action="custom", filter_query="",
                                 sort_action="custom", sort_mode="single", sort_by=[], cell_selectable=True,
                                 tooltip_delay=400, tooltip_duration=None, **TABLE_STYLE),
            html.Div(id="cat-table-note", className="lbl-sm mt-1"),
            dbc.Accordion(start_collapsed=True, always_open=True, className="mt-2", children=[
                dbc.AccordionItem(title="Raw CIF text of this category", item_id="raw",
                                  children=html.Pre(id="cat-raw", className="raw-cif")),
                dbc.AccordionItem(title="Dictionary examples for this category", item_id="ex",
                                  children=html.Div(id="cat-examples")),
            ]),
        ])),
    ])


def crystal_tab():
    return html.Div(className="mt-2", children=[
        html.Div(id="kpi-grid", className="kpi-grid mb-2"),
        dbc.Row(className="g-2", children=[
            dbc.Col(md=6, children=html.Div(className="panel panel-body", children=[
                html.Div("Matthews analysis — reported vs recomputed", className="lbl-sm mb-1"),
                html.Div(id="matthews-box"),
            ])),
            dbc.Col(md=6, children=html.Div(className="panel", children=dcc.Graph(
                id="fig-packing", config={"displaylogo": False}))),
        ]),
        html.Div(className="panel mt-2", children=[
            html.Div(className="panel-head", children=[
                html.Span("Chain", className="lbl-sm"),
                dbc.Select(id="chain-select", size="sm", style={"width": "110px"}, options=[], value=None),
                html.Span("Ramachandran class", className="lbl-sm ms-3"),
                dbc.Select(id="rama-class", size="sm", style={"width": "150px"}, value="__all__", options=[
                    {"label": "all residues", "value": "__all__"}] + [
                    {"label": c, "value": c} for c in ("General", "Glycine", "Proline", "Pre-proline", "Ile/Val")]),
                html.Span("click a point or residue to show it in 3D", className="lbl-sm ms-auto"),
            ]),
            dbc.Row(className="g-0", children=[
                dbc.Col(md=7, children=[dcc.Graph(id="fig-bprofile", config={"displaylogo": False}),
                                        dcc.Graph(id="fig-rmsf", config={"displaylogo": False},
                                                  style={"display": "none"}),
                                        dcc.Graph(id="fig-bhist", config={"displaylogo": False})]),
                dbc.Col(md=5, children=dcc.Graph(id="fig-rama", config={"displaylogo": False})),
            ]),
        ]),
        html.Div(id="shell-figs", className="mt-2"),
    ])


def composition_tab():
    return html.Div(className="mt-2", children=[
        dbc.Row(className="g-2", children=[
            dbc.Col(md=5, children=html.Div(className="panel", children=dcc.Graph(
                id="fig-sunburst", config={"displaylogo": False}))),
            dbc.Col(md=7, children=html.Div(className="panel panel-body", children=[
                html.Div("Entities (_entity) — click a row to show it in 3D", className="lbl-sm mb-1"),
                dash_table.DataTable(id="entity-table", columns=[], data=[], cell_selectable=True,
                                     page_size=10, **TABLE_STYLE),
            ])),
        ]),
        html.Div(className="panel panel-body mt-2", children=[
            html.Div(className="d-flex gap-2 align-items-center mb-2", children=[
                html.Span("Sequence of entity", className="lbl-sm"),
                dbc.Select(id="seq-entity", size="sm", style={"width": "300px"}, options=[], value=None),
                html.Span("blue underline = modelled residues · select a stretch to show it in 3D",
                          className="lbl-sm ms-auto"),
            ]),
            html.Div(id="seq-wrap", children=dash_bio.SequenceViewer(
                id="seq-viewer", sequence="-", showLineNumbers=True, charsPerLine=50, toolbar=False, search=True,
                badge=False, wrapAminoAcids=True)),
        ]),
        html.Div(className="panel mt-2", children=dcc.Graph(id="fig-catbar", config={"displaylogo": False})),
    ])


def dictionary_tab():
    return html.Div(className="mt-2 panel panel-body", children=[
        html.Div(id="dict-version", className="lbl-sm mb-2"),
        dbc.Row(className="g-2 mb-2", children=[
            dbc.Col(md=6, children=dbc.Input(id="dict-search", placeholder="search 6,800 data items, e.g. "
                                             "'resolution', 'B_iso', 'twin', 'helix'", debounce=True, size="sm",
                                             style={"background": "#232321", "color": INK, "borderColor": "#383835"})),
            dbc.Col(md=4, children=dbc.Select(id="dict-group", size="sm", value="", options=[])),
            dbc.Col(md=2, children=dbc.Switch(id="dict-only-entry", label="in this entry", value=False)),
        ]),
        html.Div(id="dict-results", style={"maxHeight": "620px", "overflowY": "auto"}),
    ])


def learn_tab():
    return html.Div(className="mt-2", children=[
        html.Div(className="panel panel-body mb-2", children=[
            html.Div("Guided tours — each step opens a category, explains a field and shows it in 3D",
                     className="lbl-sm mb-2"),
            html.Div(className="d-flex gap-2 flex-wrap", children=[
                dbc.Button(t["title"], id={"type": "tour-start", "id": tid}, size="sm", className="btn-ghost")
                for tid, t in TOURS.items()
            ]),
        ]),
        dbc.Accordion(id="lessons", start_collapsed=False, active_item=LESSONS[0]["id"], children=[
            dbc.AccordionItem(title=l["title"], item_id=l["id"],
                              children=dcc.Markdown(l["body"], mathjax=True, className="lesson"))
            for l in LESSONS
        ]),
    ])


def serve_layout():
    return html.Div(**{"data-bs-theme": "dark"}, children=[
        dcc.Store(id="store-entry"),
        dcc.Store(id="store-cat"),
        dcc.Store(id="store-item"),
        dcc.Store(id="store-hl"),
        dcc.Store(id="store-hl-view"),
        dcc.Store(id="store-view-meta"),
        dcc.Store(id="store-tour", data={"tour": None, "step": 0}),
        dcc.Store(id="viewer-resize"),
        header(),
        html.Div(id="summary", className="summary"),
        html.Div(id="tour-bar", className="tour-bar", style={"display": "none"}, children=[
            html.Strong(id="tour-title"),
            html.Span(id="tour-step", className="step"),
            html.Span(id="tour-text", style={"flex": "1"}),
            dbc.Button("◀ Prev", id="tour-prev", size="sm", className="btn-ghost"),
            dbc.Button("Next ▶", id="tour-next", size="sm", color="primary"),
            dbc.Button("Exit", id="tour-exit", size="sm", className="btn-ghost"),
        ]),
        dbc.Row(className="main g-3", children=[
            dbc.Col(lg=5, children=[viewer_panel(), left_tabs()]),
            dbc.Col(lg=7, children=dbc.Tabs(id="right-tabs", active_tab="fields", children=[
                dbc.Tab(label="Fields of this entry", tab_id="fields", children=fields_tab()),
                dbc.Tab(label="Crystal & quality", tab_id="crystal", children=crystal_tab()),
                dbc.Tab(label="Composition", tab_id="composition", children=composition_tab()),
                dbc.Tab(label="Dictionary", tab_id="dictionary", children=dictionary_tab()),
                dbc.Tab(label="Learn & tours", tab_id="learn", children=learn_tab()),
            ])),
        ]),
        html.Div(className="footer", children=[
            "Definitions: wwPDB PDBx/mmCIF dictionary (mmcif_pdbx_v50). Background: International Tables for "
            "Crystallography Vol. G ch. 3.6 & 4.5; Westbrook et al. (2022) J. Mol. Biol. 434:167599. "
            "Structures: RCSB PDB / AlphaFold DB. Viewer: dash-bio Molecule3dViewer (3Dmol.js).",
        ]),
    ])
