"""All Dash callbacks.  Server state lives in the entry LRU (entry.py) and a
small view cache here; the browser only holds keys and small selection specs."""
from __future__ import annotations

import base64
import math
import threading
import time
from collections import OrderedDict

import dash_bootstrap_components as dbc
import gemmi
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from dash import ALL, Input, Output, State, ctx, dcc, html, no_update

from . import analysis as an
from . import components as comp
from . import config
from . import figures as figs
from . import knowledge as kb
from .dictionary import load_dictionary
from .entry import Entry, EntryError, get_entry, load_entry, load_upload
from .layout import TABLE_STYLE
from .mapping import Resolution, can_map_category, residue_references, resolve_row
from .viewer import build_view, cell_edges, to_view_indices

DIC = load_dictionary()

# --------------------------------------------------------------------------- #
# server-side caches
# --------------------------------------------------------------------------- #
_VIEWS: "OrderedDict[tuple, object]" = OrderedDict()
_VLOCK = threading.Lock()


def get_view(entry: Entry, model: int, water: bool):
    key = (entry.key, model, bool(water))
    with _VLOCK:
        if key in _VIEWS:
            _VIEWS.move_to_end(key)
            return _VIEWS[key]
    v = build_view(entry, model, water)
    with _VLOCK:
        _VIEWS[key] = v
        while len(_VIEWS) > 6:
            _VIEWS.popitem(last=False)
    return v


def derived(entry: Entry, name: str, fn):
    cache = entry.__dict__.setdefault("_derived", {})
    if name not in cache:
        cache[name] = fn()
    return cache[name]


def _entry(store) -> Entry | None:
    if not store or not store.get("key"):
        return None
    return get_entry(store["key"])


def _status(ok: bool, text: str):
    return [html.Span(className="dot " + ("ok" if ok else "err")), text]


XREF_CATS = ["struct_conn", "struct_conf", "struct_sheet_range", "struct_site_gen", "pdbx_struct_conn_angle",
             "pdbx_validate_close_contact", "pdbx_validate_symm_contact", "pdbx_validate_rmsd_bond",
             "pdbx_validate_rmsd_angle", "pdbx_validate_torsion", "pdbx_validate_chiral", "pdbx_validate_planes",
             "pdbx_struct_mod_residue", "struct_mon_prot_cis", "pdbx_struct_sheet_hbond",
             "ndb_struct_na_base_pair", "struct_ref_seq_dif", "pdbx_validate_peptide_omega"]


def mappable(entry: Entry) -> dict:
    return derived(entry, "mappable", lambda: {c: can_map_category(entry, c, DIC) for c in entry.category_names()})


def default_category(entry: Entry) -> str:
    for c in ("refine", "em_3d_reconstruction", "pdbx_nmr_ensemble", "ma_qa_metric_global", "struct", "entity"):
        if not entry.table(c).empty:
            return c
    return entry.category_names()[0]


def _hl_from_rows(rows, label: str, level: str = "residue", note: str = "") -> dict:
    return {"source": "rows", "rows": [int(x) for x in rows], "label": label, "level": level, "note": note}


def _residue_rows(entry: Entry, asym: str, seq=None, auth_seq=None) -> np.ndarray:
    a = entry.atoms
    m = a["asym"].to_numpy() == asym
    if seq is not None and int(seq) >= 0:
        m &= a["seq"].to_numpy() == int(seq)
    elif auth_seq is not None:
        m &= a["auth_seq"].to_numpy() == str(auth_seq)
    return np.flatnonzero(m)


# --------------------------------------------------------------------------- #
# table filtering (DataTable custom filter syntax)
# --------------------------------------------------------------------------- #
OPERATORS = [["ge ", ">="], ["le ", "<="], ["lt ", "<"], ["gt ", ">"], ["ne ", "!="], ["eq ", "="],
             ["contains "], ["datestartswith "]]


def split_filter_part(filter_part):
    for operator_type in OPERATORS:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                name = name_part[name_part.find("{") + 1: name_part.rfind("}")]
                value_part = value_part.strip()
                v0 = value_part[0] if value_part else ""
                if v0 == value_part[-1:] and v0 in ("'", '"', "`"):
                    value = value_part[1:-1].replace("\\" + v0, v0)
                else:
                    value = value_part
                return name, operator_type[0].strip(), value
    return [None] * 3


def apply_filter(df: pd.DataFrame, query: str) -> pd.DataFrame:
    if not query:
        return df
    for part in query.split(" && "):
        col, op, val = split_filter_part(part)
        if col not in df.columns:
            continue
        s = df[col]
        if op in ("eq", "ne", "lt", "le", "gt", "ge"):
            num = pd.to_numeric(s, errors="coerce")
            try:
                fv = float(val)
                cmp = {"eq": num == fv, "ne": num != fv, "lt": num < fv, "le": num <= fv,
                       "gt": num > fv, "ge": num >= fv}[op]
                if op in ("eq", "ne") and num.isna().all():
                    raise ValueError
            except (TypeError, ValueError):
                cmp = (s == str(val)) if op == "eq" else (s != str(val)) if op == "ne" else s.str.contains(str(val))
            df = df[cmp]
        elif op == "contains":
            df = df[s.astype(str).str.contains(str(val), case=False, regex=False)]
        elif op == "datestartswith":
            df = df[s.astype(str).str.startswith(str(val))]
    return df


# --------------------------------------------------------------------------- #
def register(app):

    # ---------------- load an entry --------------------------------------- #
    @app.callback(
        Output("store-entry", "data"), Output("status", "children"), Output("pdb-input", "value"),
        Input("btn-load", "n_clicks"), Input("pdb-input", "n_submit"),
        Input({"type": "example", "id": ALL}, "n_clicks"), Input("upload", "contents"),
        State("pdb-input", "value"), State("upload", "filename"),
    )
    def load(_n, _s, _ex, contents, text, filename):
        trig = ctx.triggered_id
        t0 = time.time()
        try:
            if trig == "upload" and contents:
                data = base64.b64decode(contents.split(",", 1)[1])
                e = load_upload(data, filename or "upload.cif")
                text_out = e.entry_id
            else:
                if isinstance(trig, dict) and trig.get("type") == "example":
                    if not any(_ex):
                        return no_update, no_update, no_update
                    text = trig["id"]
                e = load_entry(text or "1TUP")
                text_out = text
        except EntryError as exc:
            return no_update, _status(False, str(exc)), no_update
        except Exception as exc:  # noqa: BLE001 - shown to the user, app keeps running
            return no_update, _status(False, f"Failed: {type(exc).__name__}: {exc}"), no_update
        n = len(e.atoms)
        return ({"key": e.key, "nonce": time.time()},
                _status(True, f"{e.entry_id} · {e.source} · {n:,} atoms · {len(e.category_names())} categories · "
                              f"{time.time() - t0:.1f}s"),
                text_out)

    # ---------------- entry-level outputs --------------------------------- #
    @app.callback(
        Output("summary", "children"), Output("model-select", "options"), Output("model-select", "value"),
        Output("chain-select", "options"), Output("chain-select", "value"),
        Output("seq-entity", "options"), Output("seq-entity", "value"),
        Output("store-cat", "data"), Output("store-item", "data"), Output("store-hl", "data"),
        Output("brand-sub", "children"), Output("mol3d", "selectedAtomIds"), Output("color-mode", "value"),
        Output("cat-search", "value"),
        Input("store-entry", "data"),
    )
    def on_entry(store):
        e = _entry(store)
        if e is None:
            return [no_update] * 14
        s = an.summary(e)
        models = [{"label": f"model {m}", "value": str(m)} for m in e.models]
        res = derived(e, "residues", lambda: an.residue_table(e))
        chains = list(dict.fromkeys(res["auth_asym"])) if not res.empty else []
        prot = list(dict.fromkeys(res.loc[res["kind"] == "protein", "auth_asym"])) if not res.empty else []
        chains = prot + [c for c in chains if c not in prot]
        comp_df = derived(e, "composition", lambda: an.composition(e))
        poly = comp_df[comp_df["type"] == "polymer"] if not comp_df.empty else comp_df
        if not poly.empty:  # proteins first
            poly = poly.iloc[np.argsort(~poly["polymer_type"].str.contains("peptide").to_numpy(), kind="stable")]
        seq_opts = [{"label": f"{r.entity}: {str(r.description)[:50]} ({r.polymer_type})", "value": str(r.entity)}
                    for r in poly.itertuples()]
        cat = default_category(e)
        return (comp.summary_tiles(s), models, str(e.models[0]),
                [{"label": c, "value": c} for c in chains], chains[0] if chains else None,
                seq_opts, seq_opts[0]["value"] if seq_opts else None,
                cat, {"category": cat}, None,
                f"dictionary v{DIC.version} · {e.entry_id} · {len(e.category_names())} categories present", [],
                "bfactor" if e.is_computed_model else "kind", "")

    # ---------------- 3D model -------------------------------------------- #
    @app.callback(
        Output("mol3d", "modelData"), Output("store-view-meta", "data"),
        Output("mol3d", "zoomTo", allow_duplicate=True), Output("mol3d", "shapes", allow_duplicate=True),
        Output("mol3d", "labels", allow_duplicate=True), Output("mol3d", "zoom", allow_duplicate=True),
        Input("store-entry", "data"), Input("model-select", "value"), Input("show-water", "value"),
        prevent_initial_call=True,
    )
    def model_data(store, model, water):
        e = _entry(store)
        if e is None:
            return [no_update] * 6
        try:
            model = int(model)
        except (TypeError, ValueError):
            model = e.models[0]
        if model not in e.models:
            model = e.models[0]
        v = get_view(e, model, water)
        meta = dict(v.meta, key=e.key, water=bool(water))
        # drop shapes/labels of the previous entry first, otherwise they enter the zoom bounds
        return (v.model_data, meta, {"sel": {}, "animationDuration": 0, "fixedPath": False}, [], [],
                {"factor": 0.95, "animationDuration": 0, "fixedPath": False})

    app.clientside_callback(
        "function(md, rep, mode, hl, meta){ return window.dash_clientside.mmcif.styles(md, rep, mode, hl, meta); }",
        Output("mol3d", "styles"),
        Input("mol3d", "modelData"), Input("rep", "value"), Input("color-mode", "value"), Input("store-hl-view", "data"),
        State("store-view-meta", "data"),
    )
    app.clientside_callback(
        "function(n, meta){ return window.dash_clientside.mmcif.viewerSize(n, meta); }",
        Output("mol3d", "width"), Output("mol3d", "height"),
        Input("viewer-resize", "data"), Input("store-view-meta", "data"),
    )

    @app.callback(Output("viewer-legend", "children"), Input("color-mode", "value"), Input("store-view-meta", "data"))
    def legend(mode, meta):
        if not meta:
            return no_update
        def sw(color, label):
            return html.Span([html.I(style={"background": color}), label])
        items = []
        if meta.get("trace_only"):
            items.append(html.Span("large entry: backbone trace (CA/P) + ligands shown", style={"color": figs.WARNING}))
        if mode == "kind":
            items += [sw(c, k) for k, c in figs.KIND_COLORS.items() if k != "polymer"]
        elif mode in ("chain", "entity"):
            vals = meta["chains"] if mode == "chain" else meta["entities"]
            items += [sw(figs.SERIES[i % 8], ("chain " if mode == "chain" else "entity ") + str(v))
                      for i, v in enumerate(vals[:16])]
            if len(vals) > 8:
                items.append(html.Span("(colours repeat after 8; click an atom to identify it)", style={"color": figs.MUTED}))
        elif mode == "ss":
            items += [sw(figs.SS_COLORS[k], figs.SS_NAMES[k]) for k in ("H", "E", "")]
        elif mode == "bfactor":
            if meta.get("plddt"):
                items += [sw("#0053d6", "pLDDT > 90"), sw("#65cbf3", "70–90"), sw("#ffdb13", "50–70"),
                          sw("#ff7d45", "< 50")]
            else:
                lo, hi = meta["b_range"]
                items += [html.Span([f"B {lo:.0f}", html.Span(className="ramp", style={
                    "background": "linear-gradient(90deg," + ",".join(figs.SEQ_BLUE) + ")"}), f"{hi:.0f} Å²",
                    html.Span("  (2nd–98th percentile of polymer atoms; brighter = more mobile)",
                              style={"color": figs.MUTED})])]
        elif mode == "occupancy":
            items += [sw("#5c6f8a", "occupancy 1"), sw("#fab219", "partial (alt. conformers)"), sw("#d03b3b", "zero")]
        elif mode == "element":
            items += [sw(c, el) for el, c in (("C", "#909090"), ("N", "#3050f8"), ("O", "#ff0d0d"), ("S", "#ffff30"),
                                              ("P", "#ff8000"), ("metal", "#7d80b0"))]
        elif mode == "residue":
            items += [sw(figs.SERIES[0], "hydrophobic"), sw(figs.SERIES[1], "positive"), sw(figs.SERIES[2], "polar"),
                      sw(figs.SERIES[4], "negative"), sw(figs.SERIES[3], "G / P / C")]
        items.append(html.Span(f"{meta['n_atoms']:,} atoms · {meta['n_bonds']:,} bonds · model {meta['model']}",
                               style={"color": figs.MUTED, "marginLeft": "auto"}))
        return items

    # ---------------- highlight → viewer ---------------------------------- #
    @app.callback(
        Output("store-hl-view", "data"), Output("mol3d", "shapes"), Output("mol3d", "labels"),
        Output("mol3d", "zoomTo", allow_duplicate=True), Output("hl-banner", "children"), Output("hl-banner", "style"),
        Output("mol3d", "zoom", allow_duplicate=True),
        Input("store-hl", "data"), Input("store-view-meta", "data"), Input("show-cell", "value"),
        prevent_initial_call=True,
    )
    def highlight(hl, meta, show_cell):
        e = get_entry(meta["key"]) if meta else None
        if e is None:
            return [no_update] * 7
        v = get_view(e, meta["model"], meta.get("water", False))
        shapes, labels, banner = [], [], []
        zoom_factor = no_update
        hidden = {"display": "none"}
        res = None
        want_cell = bool(show_cell)
        if hl and hl.get("source") == "row":
            df = e.table(hl["cat"])
            if not df.empty and 0 <= hl["row"] < len(df):
                row = df.iloc[hl["row"]].to_dict()
                res = resolve_row(e, hl["cat"], row, DIC)
                title = f"_{hl['cat']} row {hl['row'] + 1}"
                rid = row.get("id")
                if rid not in (None, "?", "."):
                    title += f" (id {rid})"
        elif hl and hl.get("source") == "rows":
            mask = np.zeros(len(e.atoms), dtype=bool)
            mask[np.asarray(hl["rows"], dtype=int)] = True
            res = Resolution(mask, hl.get("level", "residue"), note=hl.get("note", ""))
            title = hl.get("label", "selection")

        out_view = None
        zoom = no_update
        if res is not None:
            if res.special == "cell":
                want_cell = True
            serials = to_view_indices(v.index, res.mask)
            focus, sticks = [], []
            atoms = e.atoms
            if serials:
                vrows = v.index[serials]
                if res.level == "atom" and len(serials) <= 80:
                    focus = serials
                small = atoms.iloc[vrows][["asym", "auth_seq"]].drop_duplicates()
                if res.level in ("atom", "residue") and len(small) <= 30:
                    vis = atoms.iloc[v.index]
                    keys = set(zip(small["asym"], small["auth_seq"]))
                    side = ~vis["atom"].isin(["N", "C", "O", "CA", "P", "OP1", "OP2", "O5'", "C5'"]).to_numpy()
                    inres = np.array([(a_, s_) in keys for a_, s_ in zip(vis["asym"], vis["auth_seq"])])
                    sticks = np.flatnonzero(inres & side & vis["kind"].isin(["protein", "nucleic", "polymer"]).to_numpy()).tolist()
                zsel = serials
                if len(serials) <= 150:
                    # small selection: frame its 10 Å neighbourhood so the context stays visible
                    xyz = e.atoms.iloc[v.index][["x", "y", "z"]].to_numpy()
                    tree = cKDTree(xyz)
                    near = set()
                    for hits in tree.query_ball_point(xyz[serials], r=10.0):
                        near.update(hits)
                    zsel = sorted(near)
                zoom = {"sel": {"serial": zsel[:20000]}, "animationDuration": 500, "fixedPath": False}
                zoom_factor = {"factor": 0.95, "animationDuration": 500, "fixedPath": False}
            context = []
            if res.level == "atom" and serials and len(serials) <= 80:
                vis = atoms.iloc[v.index]
                keys = set(zip(atoms.iloc[v.index[serials]]["asym"], atoms.iloc[v.index[serials]]["auth_seq"]))
                context = np.flatnonzero([(a_, s_) in keys for a_, s_ in zip(vis["asym"], vis["auth_seq"])]).tolist()
            out_view = {"serials": serials, "focus": focus, "sticks": sticks, "context": context}

            # labels for partners
            for p in res.partners[:6]:
                pm = to_view_indices(v.index, p.mask)
                if not pm or not p.label:
                    continue
                pos = v.model_data["atoms"][pm[0]]["positions"]
                labels.append({"text": p.label, "position": {"x": pos[0], "y": pos[1], "z": pos[2]},
                               "backgroundColor": "#232321", "fontColor": "white", "fontSize": 11,
                               "backgroundOpacity": 0.85, "borderThickness": 0, "inFront": True})
            # measurement check
            measure_txt = None
            ap = [p for p in res.partners if p.atom_level and to_view_indices(v.index, p.mask)]
            if res.measure and len(ap) in (2, 3):
                pts = [np.array(v.model_data["atoms"][to_view_indices(v.index, p.mask)[0]]["positions"]) for p in ap]
                for i in range(len(pts) - 1):
                    shapes.append({"type": "Cylinder", "start": dict(zip("xyz", pts[i].tolist())),
                                   "end": dict(zip("xyz", pts[i + 1].tolist())), "radius": 0.07,
                                   "color": "#ffffff", "dashed": True, "alpha": 0.9})
                if len(pts) == 2:
                    val, unit = float(np.linalg.norm(pts[0] - pts[1])), "Å"
                else:
                    v1, v2 = pts[0] - pts[1], pts[2] - pts[1]
                    val = math.degrees(math.acos(max(-1, min(1, float(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2))))))
                    unit = "°"
                mid = (pts[0] + pts[-1]) / 2
                labels.append({"text": f"{val:.2f} {unit}", "position": dict(zip("xyz", mid.tolist())),
                               "backgroundColor": "#3987e5", "fontColor": "white", "fontSize": 12,
                               "backgroundOpacity": 0.95, "borderThickness": 0, "inFront": True})
                rep = res.measure.get("reported")
                measure_txt = f"Computed from coordinates: {val:.3f} {unit}"
                if rep not in (None, "?", "."):
                    try:
                        diff = abs(float(rep) - val)
                        ok = diff < (0.01 if unit == "Å" else 0.5)
                        measure_txt += f" · reported _{hl.get('cat', '')}.{res.measure['item']} = {rep}"
                        banner_cls = "ok" if ok else "bad"
                        if res.measure.get("symmetry"):
                            banner_cls = ""
                        measure_txt = html.Span(measure_txt + (" ✓ agrees" if ok else f" ✗ differs by {diff:.3f}"),
                                                className=banner_cls)
                    except ValueError:
                        pass
            n = len(serials)
            banner = [html.Strong(title), f" — {n:,} atom{'s' if n != 1 else ''} in view · {res.level}",
                      html.Span(res.note, className="note") if res.note else None,
                      html.Span(measure_txt, className="note") if measure_txt else None]
            if n == 0 and not res.special:
                banner.append(html.Span("Nothing to show in the current model/view (try toggling water or another model).",
                                        className="note"))

        if want_cell:
            edges = cell_edges(e)
            for s, t in edges:
                shapes.append({"type": "Cylinder", "start": dict(zip("xyz", s)), "end": dict(zip("xyz", t)),
                               "radius": 0.25, "color": "#c3c2b7", "alpha": 0.9})
            if edges:
                cell = an.unit_cell(e)
                for lab, f in (("a", (1, 0, 0)), ("b", (0, 1, 0)), ("c", (0, 0, 1)), ("O", (0, 0, 0))):
                    p = cell.orthogonalize(gemmi.Fractional(*f))
                    labels.append({"text": lab, "position": {"x": p.x, "y": p.y, "z": p.z}, "fontSize": 14,
                                   "backgroundColor": "#383835", "fontColor": "white", "backgroundOpacity": 0.9,
                                   "borderThickness": 0})
                if res is not None and res.special == "cell":
                    zoom = {"sel": {}, "animationDuration": 500, "fixedPath": False}
                    banner = [html.Strong(f"Unit cell of {e.entry_id}"),
                              f" — {an.summary(e)['space_group']}; origin O and axes a, b, c labelled. "
                              "The model is the asymmetric unit; symmetry copies fill the rest of the cell "
                              "(see Crystal & quality → packing)."]
            elif res is not None and res.special == "cell":
                banner = [html.Strong("No unit cell"), " — this entry has no crystal lattice (NMR, EM or computed model)."]
        style = {} if banner else hidden
        return out_view, shapes, labels, zoom, banner, style, zoom_factor

    @app.callback(Output("store-hl", "data", allow_duplicate=True), Input("btn-clear", "n_clicks"),
                  prevent_initial_call=True)
    def clear(_n):
        return None

    @app.callback(Output("mol3d", "zoomTo", allow_duplicate=True), Output("mol3d", "zoom", allow_duplicate=True),
                  Input("btn-reset", "n_clicks"), prevent_initial_call=True)
    def reset(_n):
        return ({"sel": {}, "animationDuration": 400, "fixedPath": False},
                {"factor": 0.95, "animationDuration": 400, "fixedPath": False})

    # ---------------- category navigator ---------------------------------- #
    @app.callback(Output("cat-nav", "children"), Input("store-entry", "data"), Input("cat-search", "value"),
                  Input("store-cat", "data"))
    def nav(store, query, current):
        e = _entry(store)
        if e is None:
            return no_update
        q = (query or "").strip().lower().lstrip("_")
        maps = mappable(e)
        groups: dict[str, list] = {}
        for info in e.category_info():
            name = info["name"]
            if q and q not in name.lower() and not any(q in c.lower() for c in e.table(name).columns):
                continue
            groups.setdefault(kb.section_of(name), []).append(info)
        out = []
        for sid in kb.SECTION_ORDER:
            if sid not in groups:
                continue
            out.append(html.Div(className="cat-section", children=[
                html.Div(f"{kb.SECTION_ICON[sid]} {kb.SECTION_LABEL[sid]}", className="cat-section-title"),
                *[html.Button(
                    [html.Span(["_" + i["name"], html.Span("◆", className="in3d") if maps.get(i["name"]) else None]),
                     html.Span(f"{i['rows']}×{i['items']}" if i["loop"] else f"{i['items']}", className="cnt")],
                    id={"type": "cat-link", "cat": i["name"]}, n_clicks=0,
                    className="cat-link" + (" active" if i["name"] == current else ""),
                    title=comp.first_sentence((DIC.category(i["name"]) or {}).get("description", ""), 200))
                  for i in groups[sid]],
            ]))
        return out or html.Div("No category matches.", className="lbl-sm")

    @app.callback(Output("store-cat", "data", allow_duplicate=True), Output("store-item", "data", allow_duplicate=True),
                  Output("right-tabs", "active_tab", allow_duplicate=True),
                  Input({"type": "cat-link", "cat": ALL}, "n_clicks"), Input({"type": "open-cat", "cat": ALL}, "n_clicks"),
                  Input("fig-catbar", "clickData"), prevent_initial_call=True)
    def pick_cat(_a, _b, click):
        trig = ctx.triggered_id
        if trig == "fig-catbar":
            if not click:
                return no_update, no_update, no_update
            cat = click["points"][0]["customdata"]
        elif isinstance(trig, dict):
            if not ctx.triggered[0]["value"]:
                return no_update, no_update, no_update
            cat = trig["cat"]
        else:
            return no_update, no_update, no_update
        return cat, {"category": cat}, "fields"

    # ---------------- category panel -------------------------------------- #
    @app.callback(
        Output("cat-header", "children"), Output("cat-table", "columns"), Output("cat-table", "page_current"),
        Output("cat-table", "filter_query"), Output("cat-table", "sort_by"), Output("cat-raw", "children"),
        Output("cat-examples", "children"), Output("cat-table", "tooltip_header"),
        Input("store-cat", "data"), State("store-entry", "data"), State("store-item", "data"),
    )
    def cat_panel(cat, store, item):
        e = _entry(store)
        if e is None or not cat:
            return [no_update] * 8
        df = e.table(cat)
        c = DIC.category(cat) or {"description": "", "keys": [], "examples": []}
        note = kb.CATEGORY_NOTES.get(cat)
        maps = mappable(e).get(cat)
        header = [
            html.Div(className="d-flex align-items-center gap-2 flex-wrap", children=[
                html.H3("_" + cat),
                dbc.Badge(f"{len(df)} row{'s' if len(df) != 1 else ''} × {len(df.columns)} items", color="secondary"),
                dbc.Badge("loop_" if len(df) > 1 else "key-value pairs", color="dark"),
                dbc.Badge("◆ rows map to atoms — click a cell" if maps else "no atom mapping",
                          color="success" if maps else "dark"),
                dbc.Badge(kb.SECTION_LABEL[kb.section_of(cat)], color="info"),
            ]),
            html.Div(comp.rich_text(c["description"]), className="cat-desc"),
            html.Div(dcc.Markdown(note, mathjax=True), className="note-box") if note else None,
            html.Div(["Key: "] + [html.Span(k + " ", className="mono", style={"color": figs.INK}) for k in c["keys"]],
                     className="lbl-sm"),
        ]
        if len(df) == 1:
            cols = [{"name": "item", "id": "_item"}, {"name": "value", "id": "_value"},
                    {"name": "definition (first sentence)", "id": "_def"}]
        else:
            cols = [{"name": "#", "id": "_n"}] + [{"name": col, "id": col} for col in df.columns]
        page = 0
        if item and item.get("cat") == cat and item.get("row") is not None:
            page = int(item["row"]) // config.TABLE_PAGE_SIZE
        raw, n = e.raw_category_text(cat, 40)
        if n > 40:
            raw += f"\n# … {n - 40} more rows not shown"
        examples = [html.Div([html.Div(comp.rich_text(x["detail"]), className="lbl-sm mt-2"),
                              html.Pre(x["case"], className="raw-cif")]) for x in c.get("examples", [])[:4]] \
            or html.Div("No examples in the dictionary.", className="lbl-sm")
        tips = {}
        for col in df.columns:
            it = DIC.item(f"_{cat}.{col}")
            if it:
                tips[col] = {"value": f"**_{cat}.{col}** ({it['type']}{', ' + it['units'] if it['units'] else ''})\n\n"
                                      + comp.first_sentence(it["description"], 220), "type": "markdown"}
        return header, cols, page, "", [], raw, examples, tips

    @app.callback(
        Output("cat-table", "data"), Output("cat-table", "page_count"), Output("cat-table-note", "children"),
        Input("cat-table", "page_current"), Input("cat-table", "page_size"), Input("cat-table", "filter_query"),
        Input("cat-table", "sort_by"), Input("cat-table", "columns"),
        State("store-cat", "data"), State("store-entry", "data"),
    )
    def table_data(page, size, query, sort_by, _cols, cat, store):
        e = _entry(store)
        if e is None or not cat:
            return [], 1, ""
        df = e.table(cat)
        if df.empty:
            return [], 1, ""
        if len(df) == 1:
            recs = []
            for col in df.columns:
                it = DIC.item(f"_{cat}.{col}")
                recs.append({"_item": col, "_value": df[col].iloc[0], "_row": 0,
                             "_def": comp.first_sentence(it["description"], 120) if it else "(not in dictionary)"})
            return recs, 1, (f"single-row category shown transposed: {len(recs)} items · click an item or value "
                             "for its full definition")
        df = df.copy()
        df["_row"] = np.arange(len(df))
        df = apply_filter(df, query or "")
        if sort_by:
            col = sort_by[0]["column_id"]
            if col in df.columns:
                num = pd.to_numeric(df[col], errors="coerce")
                key = num if num.notna().mean() > 0.8 else df[col]
                df = df.iloc[np.argsort(key.to_numpy(), kind="stable")]
                if sort_by[0]["direction"] == "desc":
                    df = df.iloc[::-1]
        page = page or 0
        size = size or config.TABLE_PAGE_SIZE
        part = df.iloc[page * size:(page + 1) * size].copy()
        part["_n"] = part["_row"] + 1
        n_pages = max(1, math.ceil(len(df) / size))
        note = f"{len(df):,} of {len(e.table(cat)):,} rows" + (" (filtered)" if query else "") + \
            " · type in the header row to filter (e.g. A, >2.5, contains HOH) · hover a header for its definition"
        return part.to_dict("records"), n_pages, note

    @app.callback(
        Output("store-item", "data", allow_duplicate=True), Output("store-hl", "data", allow_duplicate=True),
        Input("cat-table", "active_cell"), State("cat-table", "data"), State("store-cat", "data"),
        prevent_initial_call=True,
    )
    def table_click(cell, data, cat):
        if not cell or not data or cell["row"] >= len(data):
            return no_update, no_update
        rec = data[cell["row"]]
        row = int(rec["_row"])
        col = cell["column_id"]
        if "_item" in rec:
            col = rec["_item"]
        item = {"item": f"_{cat}.{col}" if col not in ("_n", "_row") else None, "row": row, "cat": cat}
        if item["item"] is None:
            item = {"category": cat, "row": row, "cat": cat}
        return item, {"source": "row", "cat": cat, "row": row}

    @app.callback(Output("cat-table", "style_data_conditional"), Input("store-item", "data"), State("store-cat", "data"))
    def row_style(item, cat):
        base = list(TABLE_STYLE["style_data_conditional"])
        if item and item.get("cat") == cat and item.get("item") and item.get("row") == 0:
            base.append({"if": {"filter_query": f'{{_item}} = "{item["item"].split(".", 1)[1]}"'},
                         "backgroundColor": "#14243a", "color": figs.INK})
        if item and item.get("cat") == cat and item.get("row") is not None:
            base.append({"if": {"filter_query": f"{{_row}} = {int(item['row'])}"},
                         "backgroundColor": "#14243a", "color": figs.INK})
        return base

    @app.callback(Output("cat-chips", "children"), Input("store-cat", "data"), Input("store-item", "data"),
                  State("store-entry", "data"))
    def chips(cat, item, store):
        e = _entry(store)
        if e is None or not cat:
            return no_update
        df = e.table(cat)
        c = DIC.category(cat) or {"keys": []}
        keys = {k.lower() for k in c["keys"]}
        active = (item or {}).get("item") or ""
        out = []
        for col in df.columns:
            name = f"_{cat}.{col}"
            cls = "chip"
            if name.lower() in keys:
                cls += " key"
            if name.lower() in kb.ITEM_NOTES:
                cls += " curated"
            if name.lower() == active.lower():
                cls += " active"
            out.append(html.Button(col, id={"type": "item-chip", "item": name}, n_clicks=0, className=cls,
                                   title=comp.first_sentence((DIC.item(name) or {}).get("description", ""), 160)))
        return [html.Span("items:", className="lbl-sm me-1")] + out + [
            html.Span("  (orange = category key, ★ = curated notes)", className="lbl-sm ms-1")]

    # ---------------- item selection from many places ---------------------- #
    @app.callback(
        Output("store-item", "data", allow_duplicate=True),
        Input({"type": "item-chip", "item": ALL}, "n_clicks"), Input({"type": "item-link", "item": ALL}, "n_clicks"),
        Input({"type": "cat-explain", "cat": ALL}, "n_clicks"), Input({"type": "tile", "item": ALL}, "n_clicks"),
        Input({"type": "kpi", "item": ALL}, "n_clicks"),
        State("store-item", "data"), State("store-cat", "data"),
        prevent_initial_call=True,
    )
    def pick_item(_a, _b, _c, _d, _e, current, cat):
        trig = ctx.triggered_id
        if not isinstance(trig, dict) or not ctx.triggered[0]["value"]:
            return no_update
        if trig["type"] == "cat-explain":
            return {"category": trig["cat"], "cat": cat}
        name = trig["item"]
        if not name:
            return no_update
        row = None
        icat = name.lstrip("_").split(".")[0]
        if current and current.get("cat") == icat:
            row = current.get("row")
        return {"item": name, "row": row, "cat": icat}

    @app.callback(
        Output("store-item", "data", allow_duplicate=True), Output("store-cat", "data", allow_duplicate=True),
        Input({"type": "atom-item", "item": ALL, "row": ALL}, "n_clicks"), prevent_initial_call=True,
    )
    def pick_atom_item(_n):
        trig = ctx.triggered_id
        if not isinstance(trig, dict) or not ctx.triggered[0]["value"]:
            return no_update, no_update
        return {"item": trig["item"], "row": trig["row"], "cat": "atom_site"}, "atom_site"

    @app.callback(
        Output("store-cat", "data", allow_duplicate=True), Output("store-item", "data", allow_duplicate=True),
        Output("store-hl", "data", allow_duplicate=True), Output("right-tabs", "active_tab", allow_duplicate=True),
        Input({"type": "xref", "cat": ALL, "row": ALL}, "n_clicks"), State("store-entry", "data"),
        prevent_initial_call=True,
    )
    def pick_xref(_n, store):
        trig = ctx.triggered_id
        if not isinstance(trig, dict) or not ctx.triggered[0]["value"]:
            return no_update, no_update, no_update, no_update
        cat, row = trig["cat"], trig["row"]
        e = _entry(store)
        first = e.table(cat).columns[0] if e is not None and not e.table(cat).empty else "id"
        key = (DIC.category(cat) or {}).get("keys") or [f"_{cat}.{first}"]
        return cat, {"item": key[0], "row": row, "cat": cat}, {"source": "row", "cat": cat, "row": row}, "fields"

    @app.callback(Output("explain-panel", "children"), Output("left-tabs", "active_tab", allow_duplicate=True),
                  Input("store-item", "data"), State("store-entry", "data"), prevent_initial_call=True)
    def explain(item, store):
        if not item:
            return no_update, no_update
        e = _entry(store)
        if item.get("item"):
            return comp.item_explanation(item["item"], DIC, e, item.get("row")), "explain"
        if item.get("category"):
            return comp.category_explanation(item["category"], DIC, e), "explain"
        return no_update, no_update

    # ---------------- picked atoms ---------------------------------------- #
    @app.callback(Output("atom-panel", "children"), Output("left-tabs", "active_tab", allow_duplicate=True),
                  Output("mol3d", "selectedAtomIds", allow_duplicate=True),
                  Input("mol3d", "selectedAtomIds"), State("store-view-meta", "data"), prevent_initial_call=True)
    def picked(ids, meta):
        if not ids or not meta:
            return no_update, no_update, no_update
        if len(ids) > 3:
            return no_update, no_update, ids[-1:]
        e = get_entry(meta["key"])
        if e is None:
            return no_update, no_update, no_update
        v = get_view(e, meta["model"], meta.get("water", False))
        rows = [int(v.index[i]) for i in ids if 0 <= i < len(v.index)]
        if not rows:
            return no_update, no_update, no_update
        at = e.atoms.iloc[rows[-1]]
        xr = residue_references(e, at["asym"], at["auth_asym"], int(at["seq"]), at["auth_seq"], at["ins"], XREF_CATS)
        return comp.atom_panel(e, DIC, rows, xr), "atom", no_update

    # ---------------- crystal & quality tab ------------------------------- #
    @app.callback(
        Output("kpi-grid", "children"), Output("matthews-box", "children"), Output("fig-packing", "figure"),
        Output("shell-figs", "children"), Output("fig-bhist", "figure"),
        Input("store-entry", "data"),
    )
    def crystal(store):
        e = _entry(store)
        if e is None:
            return [no_update] * 5
        s = an.summary(e)

        def kpi(label, val, item, sub=""):
            return html.Div(className="kpi", id={"type": "kpi", "item": item}, n_clicks=0,
                            children=[html.Div(label, className="lbl"), html.Div(val, className="val"),
                                      html.Div(sub, className="sub")])
        k = []
        if s["resolution"]:
            k.append(kpi("Resolution", f"{s['resolution']:.2f} Å", "_" + s["resolution_tag"], s["resolution_tag"]))
        if s["r_work"] is not None:
            k.append(kpi("R-work", f"{s['r_work']:.3f}", "_" + s["r_work_tag"], s["r_work_tag"]))
        k.append(kpi("R-free", f"{s['r_free']:.3f}" if s["r_free"] is not None else "not reported",
                     "_refine.ls_R_factor_R_free",
                     f"gap {s['r_free'] - s['r_work']:.3f}" if s["r_free"] is not None and s["r_work"] is not None else ""))
        for label, tag, unit in (("Wilson B", "reflns.B_iso_Wilson_estimate", " Å²"),
                                 ("Completeness", "reflns.percent_possible_obs", " %"),
                                 ("Multiplicity", "reflns.pdbx_redundancy", ""),
                                 ("⟨I/σ(I)⟩", "reflns.pdbx_netI_over_sigmaI", ""),
                                 ("R-merge", "reflns.pdbx_Rmerge_I_obs", ""),
                                 ("Temperature", "diffrn.ambient_temp", " K"),
                                 ("Wavelength", "diffrn_radiation_wavelength.wavelength", " Å"),
                                 ("Crystal pH", "exptl_crystal_grow.pH", ""),
                                 ("EM resolution", "em_3d_reconstruction.resolution", " Å"),
                                 ("NMR conformers", "pdbx_nmr_ensemble.conformers_submitted_total_number", "")):
            x = e.value(tag)
            if x:
                k.append(kpi(label, x + unit, "_" + tag, tag))
        if s["mean_b"] is not None:
            k.append(kpi("Mean B (computed)" if not s["computed"] else "Mean pLDDT (computed)", f"{s['mean_b']:.1f}",
                         "_atom_site.B_iso_or_equiv", "mean over all atoms of model 1"))

        mt = an.matthews(e)
        if mt:
            def row(label, rep, calc, fmt="{:.2f}"):
                return html.Tr([html.Td(label), html.Td("—" if rep is None else fmt.format(rep), className="v"),
                                html.Td(fmt.format(calc), className="v")])
            box = [
                dbc.Table([html.Thead(html.Tr([html.Th(""), html.Th("reported"), html.Th("recomputed")])),
                           html.Tbody([
                               row("cell volume (Å³)", None, mt["volume"], "{:,.0f}"),
                               row("symmetry operators", None, mt["z_sym"], "{:d}"),
                               row("ASU mass (Da)", None, mt["mass"], "{:,.0f}"),
                               row("V_M (Å³/Da)", mt["vm_reported"], mt["vm"]),
                               row("solvent (%)", mt["solvent_reported"], mt["solvent"], "{:.1f}"),
                           ])], size="sm", bordered=True, className="calc-table"),
                dcc.Markdown(
                    r"$V_M = V_{cell}/(Z\cdot M_{ASU})$, solvent $= 1-1.230/V_M$. "
                    "ASU mass = Σ `_entity.formula_weight` × `pdbx_number_of_molecules` (water excluded). "
                    + (f"_cell.Z_PDB = {mt['z_pdb']:.0f} chains per cell. " if mt["z_pdb"] else "")
                    + "A reported value that disagrees usually reflects a different mass assumption.",
                    mathjax=True, className="lbl-sm", style={"textTransform": "none", "letterSpacing": "0"}),
                html.Div("Asymmetric-unit content used:", className="lbl-sm mt-2"),
                html.Ul([html.Li(f"entity {p['entity']} ({p['type']}): {p['copies']} × {p['formula_weight']:,.0f} Da — "
                                 f"{str(p['description'])[:60]}", style={"fontSize": "12px"}) for p in mt["parts"]]),
            ]
        else:
            box = html.Div("No unit cell / space group in this entry (not a crystal structure), so no Matthews analysis.",
                           className="lbl-sm")
        pack = derived(e, "packing", lambda: an.crystal_packing(e))
        sh = an.shells(e)
        shell_figs = [html.Div(className="panel mb-2", children=dcc.Graph(
            figure=f, config={"displaylogo": False})) for f in figs.shells_figures(sh)]
        a = e.atoms
        m1 = a[a["model"] == e.models[0]]
        return k, box, figs.packing_figure(pack), shell_figs, figs.b_histogram(m1, s["computed"])

    @app.callback(
        Output("fig-bprofile", "figure"), Output("fig-rama", "figure"), Output("fig-rmsf", "figure"),
        Output("fig-rmsf", "style"),
        Input("chain-select", "value"), Input("rama-class", "value"), State("store-entry", "data"),
    )
    def chain_figs(chain, cls, store):
        e = _entry(store)
        if e is None:
            return [no_update] * 4
        res = derived(e, "residues", lambda: an.residue_table(e))
        rama = derived(e, "rama", lambda: an.ramachandran(e))
        rmsf = derived(e, "rmsf", lambda: an.ensemble_rmsf(e))
        chain = chain or (res["auth_asym"].iloc[0] if not res.empty else None)
        bp = figs.b_profile(res, chain, plddt=e.is_computed_model) if chain else figs.empty("No polymer chains.")
        rf = figs.ramachandran(rama, chain, cls)
        if not rmsf.empty and chain:
            return bp, rf, figs.rmsf_profile(rmsf, chain), {}
        return bp, rf, figs.empty(""), {"display": "none"}

    @app.callback(
        Output("store-hl", "data", allow_duplicate=True),
        Input("fig-bprofile", "clickData"), Input("fig-rama", "clickData"), Input("fig-rmsf", "clickData"),
        State("store-entry", "data"), prevent_initial_call=True,
    )
    def residue_click(c1, c2, c3, store):
        e = _entry(store)
        click = {"fig-bprofile": c1, "fig-rama": c2, "fig-rmsf": c3}.get(ctx.triggered_id)
        if e is None or not click:
            return no_update
        cd = click["points"][0].get("customdata")
        if cd is None:
            return no_update
        if ctx.triggered_id == "fig-rama":
            chain, auth_seq, compn, cls, asym, seq = cd
            rows = _residue_rows(e, asym, seq if int(seq) >= 0 else None, auth_seq)
            label = f"{compn} {chain}{auth_seq} (φ {click['points'][0]['x']:.0f}°, ψ {click['points'][0]['y']:.0f}°, {cls})"
        else:
            asym, seq, auth_seq, compn = cd[0], cd[1], cd[2], cd[3]
            rows = _residue_rows(e, asym, seq, auth_seq)
            label = f"{compn} {auth_seq} (label {asym}:{seq})"
        return _hl_from_rows(rows, label)

    # ---------------- composition tab ------------------------------------- #
    @app.callback(
        Output("fig-sunburst", "figure"), Output("entity-table", "columns"), Output("entity-table", "data"),
        Output("fig-catbar", "figure"),
        Input("store-entry", "data"),
    )
    def composition(store):
        e = _entry(store)
        if e is None:
            return [no_update] * 4
        c = derived(e, "composition", lambda: an.composition(e))
        a = e.atoms
        m1 = a[a["model"] == e.models[0]]
        show = c.drop(columns=["sequence"]) if not c.empty else c
        cols = [{"name": n, "id": n} for n in show.columns]
        return figs.composition_sunburst(c, m1), cols, show.astype(str).to_dict("records"), figs.category_bar(e.category_info())

    @app.callback(
        Output("store-hl", "data", allow_duplicate=True), Output("store-item", "data", allow_duplicate=True),
        Input("entity-table", "active_cell"), Input("fig-sunburst", "clickData"),
        State("entity-table", "data"), State("store-entry", "data"), prevent_initial_call=True,
    )
    def entity_pick(cell, sun, data, store):
        e = _entry(store)
        if e is None:
            return no_update, no_update
        a = e.atoms
        if ctx.triggered_id == "entity-table" and cell and data:
            eid = data[cell["row"]]["entity"]
            rows = np.flatnonzero(a["entity"].to_numpy() == eid)
            df = e.table("entity")
            ridx = int(np.flatnonzero(df["id"].to_numpy() == eid)[0]) if "id" in df else None
            col = cell["column_id"]
            item_map = {"entity": "_entity.id", "type": "_entity.type", "description": "_entity.pdbx_description",
                        "formula_weight": "_entity.formula_weight", "copies": "_entity.pdbx_number_of_molecules",
                        "polymer_type": "_entity_poly.type", "auth_chains": "_entity_poly.pdbx_strand_id",
                        "label_asyms": "_struct_asym.id", "seq_length": "_entity_poly_seq.num",
                        "observed_residues": "_pdbx_poly_seq_scheme.auth_seq_num", "atoms": "_atom_site.label_entity_id"}
            item = item_map.get(col, "_entity.id")
            icat = item.lstrip("_").split(".")[0]
            return (_hl_from_rows(rows, f"entity {eid}", "entity"),
                    {"item": item, "row": ridx if icat == "entity" else None, "cat": icat})
        if ctx.triggered_id == "fig-sunburst" and sun:
            pid = sun["points"][0].get("id", "")
            kind, _, val = pid.partition(":")
            col = {"k": "kind", "e": "entity", "a": "asym"}.get(kind)
            if not col:
                return no_update, no_update
            rows = np.flatnonzero(a[col].to_numpy() == val)
            return _hl_from_rows(rows, f"{ {'k': 'kind', 'e': 'entity', 'a': 'label_asym_id'}[kind]} {val}", "entity"), no_update
        return no_update, no_update

    @app.callback(Output("seq-viewer", "sequence"), Output("seq-viewer", "coverage"), Output("seq-viewer", "title"),
                  Input("seq-entity", "value"), State("store-entry", "data"))
    def seq(eid, store):
        e = _entry(store)
        if e is None or not eid:
            return "-", [], ""
        s, segs = an.sequence_coverage(e, str(eid))
        if not s:
            return "-", [], "no sequence"
        cov = [{"start": seg["start"] - 1, "end": seg["end"], "color": "#ffffff", "bgcolor": "#1c5cab",
                "underscore": True, "tooltip": f"modelled {seg['start']}–{seg['end']} (label_seq_id)"} for seg in segs]
        n_obs = sum(seg["end"] - seg["start"] + 1 for seg in segs)
        title = (f"entity {eid} · {len(s)} residues in _entity_poly_seq · {n_obs} modelled in "
                 f"{segs[0]['asym'] if segs else '-'} (first copy)")
        return s, cov, title

    @app.callback(Output("store-hl", "data", allow_duplicate=True), Input("seq-viewer", "mouseSelection"),
                  State("seq-entity", "value"), State("store-entry", "data"), prevent_initial_call=True)
    def seq_select(sel, eid, store):
        e = _entry(store)
        if e is None or not sel or not eid:
            return no_update
        start, end = int(sel.get("start", 0)), int(sel.get("end", 0))
        if end < start:
            start, end = end, start
        a = e.atoms
        m = (a["entity"].to_numpy() == str(eid)) & (a["seq"].to_numpy() >= start) & (a["seq"].to_numpy() <= end)
        return _hl_from_rows(np.flatnonzero(m), f"entity {eid} sequence {start}–{end} (label_seq_id, every copy)", "range")

    # ---------------- dictionary tab -------------------------------------- #
    @app.callback(Output("dict-group", "options"), Output("dict-version", "children"), Input("store-entry", "data"))
    def dict_opts(_store):
        groups = sorted((g for g, v in DIC.groups.items() if v["categories"] and g != "inclusive_group"))
        opts = [{"label": "all category groups", "value": ""}] + [
            {"label": f"{g.replace('_group', '').replace('_', ' ')} ({len(DIC.groups[g]['categories'])})", "value": g}
            for g in groups]
        txt = (f"PDBx/mmCIF dictionary v{DIC.version}" + (" + " + ", ".join(DIC.extensions) if DIC.extensions else "") +
               f" · {len(DIC.categories)} categories · {len(DIC.items):,} items · "
               f"{len(DIC.types)} data types · {len(DIC.units)} units")
        return opts, txt

    @app.callback(Output("dict-results", "children"),
                  Input("dict-search", "value"), Input("dict-group", "value"), Input("dict-only-entry", "value"),
                  State("store-entry", "data"))
    def dict_search(q, group, only, store):
        e = _entry(store)
        present = set()
        if e is not None:
            for c in e.category_names():
                present.update(f"_{c}.{col}".lower() for col in e.table(c).columns)
        q = (q or "").strip()
        out = []
        if q or group:
            cats = [c for c in DIC.categories.values()
                    if (not q or q.lower().lstrip("_") in c["id"].lower()) and (not group or group in c["groups"])][:15]
            for c in cats:
                inent = e is not None and not e.table(c["id"]).empty
                if only and not inent:
                    continue
                out.append(html.Div(className="dict-result", id={"type": "cat-explain", "cat": c["id"]}, n_clicks=0,
                                    children=[html.Span("category  ", className="lbl-sm"),
                                              html.Span("_" + c["id"], className="n"),
                                              dbc.Badge("in entry", color="success", className="ms-2") if inent else None,
                                              html.Div(comp.first_sentence(c["description"], 180), className="d")]))
        hits = DIC.search(q, group or None, limit=600 if only else 150)
        if only:
            hits = [h for h in hits if h["name"].lower() in present][:150]
        for it in hits:
            inent = it["name"].lower() in present
            out.append(html.Div(className="dict-result", id={"type": "item-link", "item": it["name"]}, n_clicks=0,
                                children=[html.Span(it["name"], className="n"),
                                          html.Span(f"  {it['type']}" + (f" · {it['units']}" if it["units"] else ""),
                                                    className="lbl-sm"),
                                          dbc.Badge("in entry", color="success", className="ms-2") if inent else None,
                                          html.Div(comp.first_sentence(it["description"], 180), className="d")]))
        if not out:
            return html.Div("Type a word to search names and definitions, or pick a category group.", className="lbl-sm")
        return out

    # ---------------- guided tours ---------------------------------------- #
    @app.callback(Output("store-tour", "data"),
                  Input({"type": "tour-start", "id": ALL}, "n_clicks"), Input("tour-prev", "n_clicks"),
                  Input("tour-next", "n_clicks"), Input("tour-exit", "n_clicks"), State("store-tour", "data"),
                  prevent_initial_call=True)
    def tour_nav(_s, _p, _n, _x, tour):
        trig = ctx.triggered_id
        tour = tour or {"tour": None, "step": 0}
        if isinstance(trig, dict):
            if not ctx.triggered[0]["value"]:
                return no_update
            return {"tour": trig["id"], "step": 0}
        if trig == "tour-exit":
            return {"tour": None, "step": 0}
        if not tour.get("tour"):
            return no_update
        n = len(kb.TOURS[tour["tour"]]["steps"])
        step = tour["step"] + (1 if trig == "tour-next" else -1)
        return {"tour": tour["tour"], "step": max(0, min(n - 1, step))}

    @app.callback(
        Output("tour-bar", "style"), Output("tour-title", "children"), Output("tour-step", "children"),
        Output("tour-text", "children"),
        Output("store-cat", "data", allow_duplicate=True), Output("store-item", "data", allow_duplicate=True),
        Output("store-hl", "data", allow_duplicate=True), Output("right-tabs", "active_tab", allow_duplicate=True),
        Input("store-tour", "data"), State("store-entry", "data"), prevent_initial_call=True,
    )
    def tour_step(tour, store):
        if not tour or not tour.get("tour"):
            return {"display": "none"}, "", "", "", no_update, no_update, no_update, no_update
        t = kb.TOURS[tour["tour"]]
        cat, item, row, text = t["steps"][tour["step"]]
        e = _entry(store)
        present = e is not None and not e.table(cat).empty
        if not present:
            text += f"  (_{cat} is not in this entry — showing the dictionary definition only.)"
        hl = {"source": "row", "cat": cat, "row": row} if present else None
        return ({}, f"Tour: {t['title']}", f"step {tour['step'] + 1}/{len(t['steps'])}", text,
                cat if present else no_update, {"item": item, "row": row if present else None, "cat": cat},
                hl, "fields" if present else no_update)
