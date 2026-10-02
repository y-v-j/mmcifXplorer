"""Builders for the dynamic HTML fragments (explanation card, category header,
atom inspector, summary tiles)."""
from __future__ import annotations

import re

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
from dash import dcc, html

from . import knowledge as kb
from .dictionary import MmcifDictionary
from .entry import Entry

SUBSUP_RE = re.compile(r"(~[^~\n]{1,24}~|\^[^\^\n]{1,24}\^)")


def rich_text(text: str):
    """Dictionary text with CIF markup (B~eq~, Å^2^) turned into <sub>/<sup>."""
    if not text:
        return []
    out = []
    for part in SUBSUP_RE.split(text):
        if not part:
            continue
        if part.startswith("~") and part.endswith("~") and len(part) > 2:
            out.append(html.Sub(part[1:-1]))
        elif part.startswith("^") and part.endswith("^") and len(part) > 2:
            out.append(html.Sup(part[1:-1]))
        else:
            out.append(part)
    return out


def first_sentence(text: str, limit: int = 110) -> str:
    t = " ".join((text or "").split())
    m = re.match(r"(.+?[.;])(\s|$)", t)
    s = m.group(1) if m else t
    return s if len(s) <= limit else s[: limit - 1] + "…"


def item_link(name: str, label: str | None = None, present: bool = True, kind: str = "item-link", **extra):
    return html.Button(label or name, id={"type": kind, "item": name, **extra}, n_clicks=0,
                       className="chip" + ("" if present else " absent"),
                       title="open definition" + ("" if present else " (not in this entry)"))


def _badge(text, color="secondary"):
    return dbc.Badge(text, color=color, className="me-1", pill=True)


def check_line(checks: dict, it: dict, dic: MmcifDictionary):
    parts = []
    names = {"type_ok": f"matches type '{it['type']}'", "enum_ok": "allowed enumeration value",
             "range_ok": f"inside permitted range {dic.range_text(it['range'])}"}
    for k, label in names.items():
        v = checks.get(k)
        if v is None:
            continue
        parts.append(html.Span(("✓ " if v else "✗ ") + label, className="check " + ("ok" if v else "bad")))
    return parts


# --------------------------------------------------------------------------- #
def category_explanation(cat_name: str, dic: MmcifDictionary, entry: Entry | None):
    c = dic.category(cat_name)
    if c is None:
        return html.Div(f"_{cat_name} is not defined in dictionary v{dic.version}.")
    note = kb.CATEGORY_NOTES.get(c["id"])
    present = entry is not None and not entry.table(c["id"]).empty
    return html.Div([
        html.H4("_" + c["id"]),
        html.Div(className="badges mb-2", children=[
            _badge("category"), _badge(kb.SECTION_LABEL[kb.section_of(c["id"])], "info"),
            _badge(dic.category_group(c["id"]).replace("_group", ""), "dark"),
            _badge("present in entry" if present else "not in this entry", "success" if present else "dark"),
        ]),
        html.Div(rich_text(c["description"]), className="def"),
        html.Div(dcc.Markdown(note, mathjax=True), className="note-box") if note else None,
        html.Div("Category key", className="section-title"),
        html.Div([item_link(k, present=present) for k in c["keys"]], className="chips"),
        html.Div(f"{len(c['items'])} items defined", className="section-title"),
        html.Div([item_link(i, i.split(".", 1)[1],
                            present=bool(entry is not None and i.split(".", 1)[1] in entry.table(c["id"]).columns))
                  for i in c["items"][:400]], className="chips"),
        dbc.Button("Open in Fields tab", id={"type": "open-cat", "cat": c["id"]}, size="sm",
                   className="btn-ghost mt-2") if present else None,
    ])


def item_explanation(item_name: str, dic: MmcifDictionary, entry: Entry | None, row: int | None = None):
    it = dic.item(item_name)
    if it is None:
        return html.Div([html.H4(item_name), html.Div("Not defined in the dictionary.", className="def")])
    name = it["name"]
    cat, _, col = name.lstrip("_").partition(".")
    df = entry.table(cat) if entry is not None else None
    present = df is not None and not df.empty and col in df.columns

    # ---- value(s) in this entry
    value_block = []
    value = None
    if present:
        series = df[col]
        if row is not None and 0 <= row < len(series):
            value = str(series.iloc[row])
            where = f"row {row + 1} of {len(series)}"
        else:
            value = str(series.iloc[0])
            where = "single value" if len(series) == 1 else f"row 1 of {len(series)}"
        value_block.append(html.Div(f"Value in {entry.entry_id} ({where})", className="section-title"))
        shown = value if len(value) < 1500 else value[:1500] + " …"
        value_block.append(html.Div(shown, className="value-box"))
        if value == "?":
            value_block.append(html.Small("'?' means the value is unknown or was not provided.", className="check na"))
        elif value == ".":
            value_block.append(html.Small("'.' means the item is not applicable here (e.g. no sequence number "
                                          "for a ligand).", className="check na"))
        else:
            value_block.append(html.Div(check_line(dic.check_value(name, value), it, dic), className="mt-1"))
        if len(series) > 1:
            num = pd.to_numeric(series.replace({"?": None, ".": None}), errors="coerce")
            if num.notna().sum() >= max(2, 0.8 * series.ne("?").sum()):
                value_block.append(html.Small(
                    f"Across {len(series)} rows: min {num.min():g}, mean {num.mean():.4g}, max {num.max():g}",
                    className="check na"))
            else:
                vc = series.value_counts().head(6)
                value_block.append(html.Small(
                    f"{series.nunique()} distinct values; most common: " +
                    ", ".join(f"{k} ({v})" for k, v in vc.items()), className="check na"))
        interp = kb.interpret(name, value, entry) if value not in (None, "?", ".") else None
        if interp:
            value_block.append(html.Div(dcc.Markdown(interp, mathjax=True), className="note-box"))
    elif entry is not None:
        value_block.append(html.Div(f"Not present in {entry.entry_id}.", className="check na mt-2"))

    # ---- enumeration
    enum_block = []
    if it["enum"]:
        rows = []
        for e in it["enum"][:80]:
            cur = value is not None and e["value"].lower() == str(value).lower()
            rows.append(html.Tr([html.Td(e["value"], className="mono"), html.Td(first_sentence(e["detail"], 160))],
                                className="current" if cur else ""))
        enum_block = [html.Div(f"Allowed values ({len(it['enum'])})", className="section-title"),
                      dbc.Table(html.Tbody(rows), size="sm", className="enum-table", bordered=True)]

    tinfo = dic.type_info(it["type"])
    facts = []
    if it["units"]:
        facts.append(html.Li(["Units: ", html.Span(it["units"], className="mono"), " — ",
                              dic.units.get(it["units"], "")]))
    if it["range"]:
        facts.append(html.Li(["Permitted range: ", html.Span(dic.range_text(it["range"]), className="mono")]))
    if it["default"]:
        facts.append(html.Li(["Default: ", html.Span(it["default"], className="mono")]))
    if it["mandatory"]:
        facts.append(html.Li(f"Mandatory: {it['mandatory']}"))
    if it["aliases"]:
        facts.append(html.Li(["Core-CIF alias: "] + [html.Span(a["alias"] + " ", className="mono") for a in it["aliases"]]))
    if it["defined_in"].lower() != name.lower():
        facts.append(html.Li(["Definition shared with parent item ", html.Span(it["defined_in"], className="mono")]))
    if "DEPRECATED" in it["context"]:
        facts.append(html.Li("Deprecated by the wwPDB.", style={"color": "#fab219"}))

    def _present(n):
        c, _, i = n.lstrip("_").partition(".")
        return entry is not None and i in entry.table(c).columns

    links = []
    if it["parents"]:
        links += [html.Div("Points to (parent items)", className="section-title"),
                  html.Div([item_link(p, present=_present(p)) for p in it["parents"]], className="chips")]
    if it["children"]:
        kids = it["children"]
        links += [html.Div(f"Referenced by (child items, {len(kids)})", className="section-title"),
                  html.Div([item_link(c, present=_present(c)) for c in kids[:60]], className="chips")]
    if it["related"]:
        links += [html.Div("Related items", className="section-title"),
                  html.Div([item_link(r["name"], f"{r['name']} · {r['function']}", present=_present(r["name"]))
                            for r in it["related"]], className="chips")]

    examples = []
    if it["examples"]:
        examples = [html.Div("Dictionary examples", className="section-title")] + [
            html.Div([html.Span(e["case"], className="mono", style={"color": "#ffffff"}),
                      html.Span("  " + first_sentence(e["detail"], 140), style={"color": "#898781"})])
            for e in it["examples"][:6]]

    note = kb.ITEM_NOTES.get(name.lower())
    cat_note = kb.CATEGORY_NOTES.get(cat)
    badges = [_badge(it["type"] or "no type", "primary"), _badge(it["units"], "dark") if it["units"] else None,
              _badge("pdbx extension", "warning") if ("pdbx" in col.lower() or cat.startswith("pdbx")) else None,
              _badge("category key", "danger") if dic.category(cat) and name in dic.category(cat)["keys"] else None,
              _badge("mandatory", "secondary") if it["mandatory"] == "yes" else None,
              _badge(f"from {it['source']} dictionary", "info") if it.get("source") else None]
    return html.Div([
        html.Div(className="d-flex align-items-start gap-2", children=[
            html.Div([html.H4(name),
                      html.Div(className="badges mb-2", children=[b for b in badges if b is not None])],
                     style={"flex": "1"}),
            dbc.Button("_" + cat, id={"type": "cat-explain", "cat": cat}, size="sm", className="btn-ghost mono",
                       title="about this category"),
        ]),
        html.Div(rich_text(it["description"]), className="def"),
        html.Div(dcc.Markdown(note, mathjax=True), className="note-box") if note else None,
        html.Div([html.Span("Deposition help text: ", style={"color": "#898781"}), it["pdbx_description"]],
                 className="mt-1", style={"fontSize": "12.5px"}) if it["pdbx_description"] else None,
        *value_block,
        *enum_block,
        html.Div("Definition facts", className="section-title") if facts else None,
        html.Ul(facts, style={"fontSize": "12.5px", "paddingLeft": "18px"}) if facts else None,
        *links,
        *examples,
        dbc.Accordion(start_collapsed=True, className="mt-2", children=[
            dbc.AccordionItem(title=f"Data type '{it['type']}' ({tinfo.get('primitive', '?')})", children=[
                html.Div(tinfo.get("detail", ""), className="def"),
                html.Div("Regular expression (POSIX) that every value must match:", className="section-title"),
                html.Pre(tinfo.get("construct", "")),
            ]),
            dbc.AccordionItem(title=f"About category _{cat}", children=html.Div(
                [html.Div(rich_text(dic.category(cat)["description"]) if dic.category(cat) else "", className="def"),
                 html.Div(dcc.Markdown(cat_note), className="note-box") if cat_note else None])),
        ]) if it["type"] else None,
    ])


# --------------------------------------------------------------------------- #
def summary_tiles(s: dict):
    def tile(label, val, item, sub=None):
        return html.Div(className="tile", id={"type": "tile", "item": item}, n_clicks=0, children=[
            html.Div(label, className="lbl"),
            html.Div([val, html.Small(" " + sub) if sub else None], className="val"),
        ])
    tiles = [html.Div(className="title", children=[s["title"], html.Span(f"{s['id']} · {s['source']}", className="src")])]
    method = ", ".join(m.title() for m in s["methods"]) or ("Computed model" if s["computed"] else "—")
    tiles.append(tile("Method", method, "_exptl.method" if s["methods"] else "_ma_model_list.model_type"))
    if s["resolution"]:
        tiles.append(tile("Resolution", f"{s['resolution']:.2f}", s["resolution_tag"] and "_" + s["resolution_tag"], "Å"))
    if s["r_work"] is not None:
        tiles.append(tile("R-work / R-free", f"{s['r_work']:.3f} / " + (f"{s['r_free']:.3f}" if s["r_free"] else "—"),
                          "_refine.ls_R_factor_R_free" if s["r_free"] else "_" + s["r_work_tag"]))
    if s["space_group"] and s["cell"][0] and s["cell"][0] != 1.0:
        tiles.append(tile("Space group", s["space_group"], "_symmetry.space_group_name_H-M"))
        a, b, c = s["cell"][:3]
        tiles.append(tile("Cell a b c", f"{a:.1f} {b:.1f} {c:.1f}", "_cell.length_a", "Å"))
    if s["computed"] and s["plddt"] is not None:
        tiles.append(tile("Mean pLDDT", f"{s['plddt']:.1f}", "_ma_qa_metric_global.metric_value"))
    tiles.append(tile("Atoms (model 1)", f"{s['n_atoms_model1']:,}", "_atom_site.id",
                      f"× {s['n_models']} models" if s["n_models"] > 1 else None))
    tiles.append(tile("Entities / chains", f"{s['n_entities']} / {s['n_chains']}", "_entity.id"))
    if s["mean_b"] is not None:
        tiles.append(tile("Mean pLDDT (B col.)" if s["computed"] else "Mean B (all atoms)", f"{s['mean_b']:.1f}",
                          "_atom_site.B_iso_or_equiv", "" if s["computed"] else "Å²"))
    if s["deposited"]:
        tiles.append(tile("Deposited", s["deposited"], "_pdbx_database_status.recvd_initial_deposition_date"))
    if s["organisms"]:
        tiles.append(tile("Organism", ", ".join(s["organisms"][:2]), "_entity_src_gen.pdbx_gene_src_scientific_name"))
    return tiles


# --------------------------------------------------------------------------- #
def atom_panel(entry: Entry, dic: MmcifDictionary, rows: list[int], xrefs: list[tuple[str, int, str]]):
    a = entry.atoms
    r = rows[-1]
    at = a.iloc[r]
    raw = entry.table("atom_site").iloc[r]
    head = f"{at['auth_asym']} {at['comp']} {at['auth_seq']}{at['ins']} · {at['atom']}" + (f" (alt {at['alt']})" if at["alt"] else "")
    measure = None
    if len(rows) >= 2:
        p = a.iloc[rows[-3:] if len(rows) >= 3 else rows[-2:]][["x", "y", "z"]].to_numpy()
        names = [f"{x.auth_asym}{x.auth_seq}:{x.atom}" for x in a.iloc[rows[-3:] if len(rows) >= 3 else rows[-2:]].itertuples()]
        if len(p) == 2:
            measure = f"Distance {names[0]} – {names[1]}: {np.linalg.norm(p[0] - p[1]):.3f} Å"
        else:
            v1, v2 = p[0] - p[1], p[2] - p[1]
            ang = np.degrees(np.arccos(np.clip(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2), -1, 1)))
            measure = f"Angle {names[0]} – {names[1]} – {names[2]}: {ang:.1f}°"
    table_rows = []
    for col, val in raw.items():
        it = dic.item(f"_atom_site.{col}")
        table_rows.append(html.Tr([
            html.Td(html.Button(col, id={"type": "atom-item", "item": f"_atom_site.{col}", "row": int(r)},
                                n_clicks=0, className="linkish")),
            html.Td(str(val), className="v"),
            html.Td(first_sentence(it["description"], 80) if it else "", style={"color": "#898781"}),
        ]))
    ident = dbc.Table([
        html.Thead(html.Tr([html.Th(""), html.Th("label_* (archive)"), html.Th("auth_* (author)")])),
        html.Tbody([
            html.Tr([html.Td("chain"), html.Td(at["asym"], className="v"), html.Td(at["auth_asym"], className="v")]),
            html.Tr([html.Td("residue no."), html.Td("." if at["seq"] < 0 else str(at["seq"]), className="v"),
                     html.Td(f"{at['auth_seq']}{at['ins']}", className="v")]),
            html.Tr([html.Td("residue"), html.Td(at["comp"], className="v"), html.Td(at["auth_comp"], className="v")]),
            html.Tr([html.Td("atom"), html.Td(at["atom"], className="v"), html.Td(at["auth_atom"], className="v")]),
        ]),
    ], size="sm", className="atom-table", bordered=True)
    xref_block = []
    if xrefs:
        xref_block = [html.Div("This residue is referenced by", className="lbl-sm mt-2"),
                      html.Div([html.Button(f"_{c} · {label}", id={"type": "xref", "cat": c, "row": int(i)},
                                            n_clicks=0, className="chip") for c, i, label in xrefs[:40]],
                               className="chips")]
    return html.Div([
        html.H5(head, className="mono", style={"color": "#ffffff"}),
        html.Div(f"_atom_site row {r + 1} of {len(a)} · {at['kind']} · B {at['b']:.2f} Å² · occupancy {at['occ']:.2f}",
                 className="lbl-sm mb-2"),
        html.Div(measure, className="note-box mono") if measure else None,
        ident,
        *xref_block,
        html.Div("Full _atom_site row — click a column name for its definition", className="lbl-sm mt-2"),
        dbc.Table(html.Tbody(table_rows), size="sm", className="atom-table", bordered=True),
    ])
