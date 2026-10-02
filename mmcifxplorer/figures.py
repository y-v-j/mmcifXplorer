"""Plotly figures.  Colours follow the validated dark-mode reference palette
(dataviz skill): categorical slots in fixed order, a one-hue sequential ramp,
status red reserved for flagged items (always with a symbol + label)."""
from __future__ import annotations

import math

import gemmi
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# ---- palette (dark mode) ---------------------------------------------------
SURFACE = "#1a1a19"
PAGE = "#0d0d0d"
INK = "#ffffff"
INK_2 = "#c3c2b7"
MUTED = "#898781"
GRID = "#2c2c2a"
AXIS = "#383835"
SERIES = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"]
SEQ_BLUE = ["#184f95", "#1c5cab", "#256abf", "#2a78d6", "#3987e5", "#5598e7", "#6da7ec", "#86b6ef",
            "#9ec5f4", "#b7d3f6", "#cde2fb"]
CRITICAL = "#d03b3b"
GOOD = "#0ca30c"
WARNING = "#fab219"

KIND_COLORS = {"protein": SERIES[0], "nucleic": SERIES[1], "ligand": SERIES[2], "ion": SERIES[3],
               "branched": SERIES[4], "polymer": SERIES[6], "water": MUTED}
SS_COLORS = {"H": SERIES[0], "E": SERIES[1], "T": SERIES[2], "": MUTED}
SS_NAMES = {"H": "helix (_struct_conf)", "E": "strand (_struct_sheet_range)", "T": "turn", "": "other"}

FONT = 'system-ui, -apple-system, "Segoe UI", sans-serif'


def base_layout(fig: go.Figure, title: str | None = None, height: int = 300, **kw) -> go.Figure:
    layout = dict(
        template="none",
        height=height,
        margin=dict(l=52, r=16, t=40 if title else 12, b=44),
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        font=dict(family=FONT, size=12, color=INK_2),
        title=dict(text=title, x=0.01, xanchor="left", font=dict(size=13, color=INK)) if title else None,
        hoverlabel=dict(bgcolor="#262624", bordercolor=AXIS, font=dict(color=INK, family=FONT)),
        legend=dict(orientation="h", y=1.02, yanchor="bottom", x=1, xanchor="right",
                    font=dict(color=INK_2), bgcolor="rgba(0,0,0,0)"),
    )
    layout.update(kw)
    fig.update_layout(**layout)
    fig.update_xaxes(gridcolor=GRID, zerolinecolor=AXIS, linecolor=AXIS, tickfont=dict(color=MUTED),
                     title_font=dict(color=INK_2))
    fig.update_yaxes(gridcolor=GRID, zerolinecolor=AXIS, linecolor=AXIS, tickfont=dict(color=MUTED),
                     title_font=dict(color=INK_2))
    return fig


def empty(message: str, height: int = 240) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=message, showarrow=False, font=dict(color=MUTED, size=13), x=0.5, y=0.5,
                       xref="paper", yref="paper")
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    return base_layout(fig, height=height)


# --------------------------------------------------------------------------- #
def b_profile(res: pd.DataFrame, chain: str, plddt: bool = False) -> go.Figure:
    sub = res[res["auth_asym"] == chain] if not res.empty else res
    if sub.empty:
        return empty("No polymer residues in this chain.")
    x = np.arange(len(sub))
    label = [f"{c}{s}{i}" for c, s, i in zip(sub["comp"], sub["auth_seq"], sub["ins"])]
    ylab = "pLDDT" if plddt else "mean B (Å²)"
    fig = go.Figure()
    # secondary-structure strip from the file's own annotation
    for code in ("H", "E"):
        m = (sub["ss"] == code).to_numpy()
        if not m.any():
            continue
        starts = np.flatnonzero(m & ~np.r_[False, m[:-1]])
        ends = np.flatnonzero(m & ~np.r_[m[1:], False])
        for s, e in zip(starts, ends):
            fig.add_vrect(x0=s - 0.5, x1=e + 0.5, fillcolor=SS_COLORS[code], opacity=0.14, line_width=0, layer="below")
    fig.add_trace(go.Scatter(
        x=x, y=sub["b_mean"], mode="lines", line=dict(color=SERIES[0], width=2), name=ylab,
        customdata=np.stack([sub["asym"], sub["seq"], sub["auth_seq"], sub["comp"], sub["n_atoms"]], axis=1),
        text=label,
        hovertemplate="<b>%{text}</b> (label %{customdata[0]}:%{customdata[1]})<br>" + ylab +
                      " %{y:.1f}<br>%{customdata[4]} atoms<extra></extra>",
    ))
    tick_step = max(1, int(math.ceil(len(sub) / 8)))
    fig.update_xaxes(tickmode="array", tickvals=x[::tick_step], ticktext=[str(s) for s in sub["auth_seq"].iloc[::tick_step]],
                     tickangle=0, title="author residue number (shaded: blue helix, orange strand from the file)")
    fig.update_yaxes(title=ylab)
    return base_layout(fig, f"Per-residue {'confidence' if plddt else 'B-factor'} · chain {chain}", height=270,
                       showlegend=False)


def rmsf_profile(rmsf: pd.DataFrame, chain: str) -> go.Figure:
    sub = rmsf[rmsf["chain"] == chain]
    if sub.empty:
        return empty("Single-model entry: no ensemble spread.")
    fig = go.Figure(go.Scatter(
        x=np.arange(len(sub)), y=sub["rmsf"], mode="lines", line=dict(color=SERIES[1], width=2),
        customdata=np.stack([sub["asym"], sub["seq"], sub["auth_seq"], sub["comp"], sub["n_models"]], axis=1),
        hovertemplate="<b>%{customdata[3]}%{customdata[2]}</b><br>CA spread %{y:.2f} Å over "
                      "%{customdata[4]} models<extra></extra>"))
    step = max(1, len(sub) // 12)
    fig.update_xaxes(tickmode="array", tickvals=np.arange(len(sub))[::step],
                     ticktext=[str(s) for s in sub["auth_seq"].iloc[::step]], title="author residue number")
    fig.update_yaxes(title="RMS deviation from mean (Å)")
    return base_layout(fig, f"Ensemble spread (CA/P) · chain {chain}", height=250, showlegend=False)


def ramachandran(df: pd.DataFrame, chain: str | None = None, cls: str | None = None) -> go.Figure:
    if df.empty:
        return empty("No protein backbone to analyse.")
    sub = df
    if chain and chain != "__all__":
        sub = sub[sub["chain"] == chain]
    if cls and cls != "__all__":
        sub = sub[sub["class"] == cls]
    if sub.empty:
        return empty("No protein residues with both φ and ψ in this selection.", height=380)
    fig = go.Figure()
    # Reference regions (approximate centres of the populated basins)
    for name, (x0, x1, y0, y1) in {"β": (-180, -45, 90, 180), "αR": (-160, -20, -80, 10),
                                   "αL": (40, 90, 10, 90)}.items():
        fig.add_shape(type="rect", x0=x0, x1=x1, y0=y0, y1=y1, line=dict(color=AXIS, width=1, dash="dot"),
                      fillcolor="rgba(57,135,229,0.06)", layer="below")
        fig.add_annotation(x=x0 + 6, y=y1 - 8, text=name, showarrow=False, xanchor="left",
                           font=dict(color=MUTED, size=12))
    for flagged, color, sym, name in ((False, SERIES[0], "circle", "residues"),
                                      (True, CRITICAL, "x", "✕ listed in _pdbx_validate_torsion")):
        s = sub[sub["flagged"] == flagged]
        if s.empty:
            continue
        fig.add_trace(go.Scattergl(
            x=s["phi"], y=s["psi"], mode="markers", name=name,
            marker=dict(color=color, size=8 if not flagged else 10, symbol=sym,
                        line=dict(color=SURFACE, width=1)),
            customdata=np.stack([s["chain"], s["auth_seq"], s["comp"], s["class"], s["asym"], s["seq"]], axis=1),
            hovertemplate="<b>%{customdata[2]} %{customdata[0]}%{customdata[1]}</b> (%{customdata[3]})"
                          "<br>φ %{x:.1f}°, ψ %{y:.1f}°<extra></extra>"))
    fig.update_xaxes(range=[-180, 180], dtick=60, title="φ (°)", zeroline=True)
    fig.update_xaxes(autorange=False, constrain="domain")
    fig.update_yaxes(range=[-180, 180], dtick=60, title="ψ (°)", zeroline=True, autorange=False)
    n = len(sub)
    return base_layout(fig, f"Ramachandran plot · {n} residues with φ and ψ", height=400,
                       legend=dict(orientation="h", y=-0.22, yanchor="top", x=0, xanchor="left",
                                   font=dict(color=INK_2), bgcolor="rgba(0,0,0,0)"),
                       margin=dict(l=52, r=16, t=40, b=90))


def b_histogram(atoms: pd.DataFrame, plddt: bool = False) -> go.Figure:
    if atoms.empty:
        return empty("No atoms.")
    fig = go.Figure(go.Histogram(x=atoms["b"], nbinsx=40, marker=dict(color=SERIES[0], line=dict(color=SURFACE, width=2)),
                                 hovertemplate=("pLDDT" if plddt else "B") + " %{x}<br>%{y} atoms<extra></extra>"))
    mean = atoms["b"].mean()
    fig.add_vline(x=mean, line=dict(color=INK_2, width=1, dash="dash"))
    fig.add_annotation(x=mean, y=1, yref="paper", text=f"mean {mean:.1f}", showarrow=False, yanchor="bottom",
                       font=dict(color=INK_2))
    fig.update_xaxes(title="pLDDT" if plddt else "B (Å²)")
    fig.update_yaxes(title="atoms")
    return base_layout(fig, "Distribution over atoms (model 1)", height=240, showlegend=False, bargap=0.02)


def shells_figures(shells: dict) -> list[go.Figure]:
    figs = []
    rs = shells.get("reflns_shell")
    ls = shells.get("refine_ls_shell")

    if ls is not None and "d_res_high" in ls:
        fig = go.Figure()
        for i, (col, name) in enumerate((("R_factor_R_work", "R-work"), ("R_factor_R_free", "R-free"))):
            if col in ls and ls[col].notna().any():
                fig.add_trace(go.Scatter(x=ls["d_res_high"], y=ls[col], mode="lines+markers", name=name,
                                         line=dict(color=SERIES[i], width=2), marker=dict(size=8),
                                         hovertemplate=name + " %{y:.3f} at %{x:.2f} Å<extra></extra>"))
        if fig.data:
            fig.update_xaxes(title="shell high-resolution limit (Å)", autorange="reversed")
            fig.update_yaxes(title="R")
            figs.append(base_layout(fig, "Refinement R per resolution shell", height=250))
    if rs is not None and "d_res_high" in rs:
        for col, name, color in (("percent_possible_all", "completeness (%)", SERIES[0]),
                                 ("meanI_over_sigI_obs", "⟨I/σ(I)⟩", SERIES[0]),
                                 ("Rmerge_I_obs", "R-merge", SERIES[0]),
                                 ("pdbx_CC_half", "CC½", SERIES[0])):
            if col in rs and rs[col].notna().any():
                fig = go.Figure(go.Bar(x=rs["d_res_high"].astype(str), y=rs[col],
                                       marker=dict(color=color, line=dict(color=SURFACE, width=2)),
                                       hovertemplate=name + " %{y} at %{x} Å<extra></extra>"))
                fig.update_xaxes(title="shell high-resolution limit (Å)", type="category")
                fig.update_yaxes(title=name)
                figs.append(base_layout(fig, f"Data: {name} per shell", height=230, showlegend=False))
    return figs


def packing_figure(pack: dict | None) -> go.Figure:
    if not pack:
        return empty("No crystal lattice: NMR, cryo-EM and computed models have no unit cell.", height=380)
    cell = pack["cell"]
    fig = go.Figure()
    corners = {}
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                p = cell.orthogonalize(gemmi.Fractional(i, j, k))
                corners[(i, j, k)] = (p.x, p.y, p.z)
    ex, ey, ez = [], [], []
    for (i, j, k), p in corners.items():
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            q = (i + d[0], j + d[1], k + d[2])
            if q in corners:
                ex += [p[0], corners[q][0], None]
                ey += [p[1], corners[q][1], None]
                ez += [p[2], corners[q][2], None]
    fig.add_trace(go.Scatter3d(x=ex, y=ey, z=ez, mode="lines", line=dict(color=INK_2, width=3),
                               name="unit cell", hoverinfo="skip"))
    for name, idx in (("a", (1, 0, 0)), ("b", (0, 1, 0)), ("c", (0, 0, 1))):
        p = corners[idx]
        fig.add_trace(go.Scatter3d(x=[p[0]], y=[p[1]], z=[p[2]], mode="text", text=[name], showlegend=False,
                                   textfont=dict(color=INK, size=14), hoverinfo="skip"))
    for n, cp in enumerate(pack["copies"]):
        xyz, ch = cp["xyz"], cp["chain"]
        xs, ys, zs = [], [], []
        for c in dict.fromkeys(ch):
            m = ch == c
            xs += list(xyz[m, 0]) + [None]
            ys += list(xyz[m, 1]) + [None]
            zs += list(xyz[m, 2]) + [None]
        ident = n == 0
        fig.add_trace(go.Scatter3d(
            x=xs, y=ys, z=zs, mode="lines",
            line=dict(color=SERIES[0] if ident else "#6b6a66", width=4 if ident else 2),
            name=("asymmetric unit (x,y,z)" if ident else f"copy {cp['op']}"),
            showlegend=ident or n == 1,
            legendgroup="asu" if ident else "copies",
            hovertemplate=("ASU" if ident else cp["op"]) + "<extra></extra>"))
    if len(pack["copies"]) > 1:
        fig.data[-1].name = "symmetry copies (grey)"
    scene_axis = dict(backgroundcolor=SURFACE, gridcolor=GRID, zerolinecolor=AXIS, color=MUTED, showbackground=True)
    base_layout(fig, f"Crystal packing · {pack['sg']} · {len(pack['copies'])} copies per cell", height=420)
    fig.update_layout(scene=dict(xaxis=dict(title="x (Å)", **scene_axis), yaxis=dict(title="y (Å)", **scene_axis),
                                 zaxis=dict(title="z (Å)", **scene_axis), aspectmode="data",
                                 bgcolor=SURFACE),
                      margin=dict(l=0, r=0, t=40, b=0),
                      legend=dict(orientation="h", y=0, yanchor="bottom", x=0, xanchor="left",
                                  font=dict(color=INK_2), bgcolor="rgba(0,0,0,0)"))
    return fig


def composition_sunburst(comp: pd.DataFrame, atoms: pd.DataFrame) -> go.Figure:
    if comp.empty:
        return empty("No entities.")
    ids, labels, parents, values, colors, hover = [], [], [], [], [], []
    a = atoms
    kinds = {}
    for _, r in comp.iterrows():
        sub = a[a["entity"] == r["entity"]]
        kind = sub["kind"].mode().iloc[0] if len(sub) else "ligand"
        kinds.setdefault(kind, 0)
        kinds[kind] += len(sub)
    for k, n in kinds.items():
        ids.append(f"k:{k}"); labels.append(k); parents.append(""); values.append(n)
        colors.append(KIND_COLORS.get(k, MUTED)); hover.append(f"{k}: {n} atoms")
    for _, r in comp.iterrows():
        sub = a[a["entity"] == r["entity"]]
        if sub.empty:
            continue
        kind = sub["kind"].mode().iloc[0]
        eid = f"e:{r['entity']}"
        desc = str(r["description"])[:40]
        ids.append(eid); labels.append(f"entity {r['entity']}"); parents.append(f"k:{kind}")
        values.append(len(sub)); colors.append(KIND_COLORS.get(kind, MUTED))
        hover.append(f"entity {r['entity']}: {desc}<br>{len(sub)} atoms")
        for asym, g in sub.groupby("asym", sort=False):
            ids.append(f"a:{asym}"); labels.append(f"{g['auth_asym'].iloc[0]} ({asym})"); parents.append(eid)
            values.append(len(g)); colors.append(KIND_COLORS.get(kind, MUTED))
            hover.append(f"label_asym_id {asym} · auth chain {g['auth_asym'].iloc[0]}<br>{len(g)} atoms")
    fig = go.Figure(go.Sunburst(ids=ids, labels=labels, parents=parents, values=values, branchvalues="total",
                                marker=dict(colors=colors, line=dict(color=SURFACE, width=2)),
                                hovertext=hover, hovertemplate="%{hovertext}<extra></extra>",
                                insidetextfont=dict(color=INK), maxdepth=3))
    return base_layout(fig, "Kind → entity → instance (label_asym_id), sized by atoms", height=360)


def category_bar(info: list[dict], top: int = 15) -> go.Figure:
    df = pd.DataFrame(info).sort_values("rows", ascending=False).head(top).iloc[::-1]
    if df.empty:
        return empty("No categories.")
    fig = go.Figure(go.Bar(
        y=["_" + n for n in df["name"]], x=df["rows"], orientation="h",
        marker=dict(color=SERIES[0], line=dict(color=SURFACE, width=2)),
        customdata=df["name"], text=df["rows"], textposition="outside", textfont=dict(color=INK_2),
        hovertemplate="<b>%{y}</b><br>%{x} rows · click to open<extra></extra>"))
    fig.update_xaxes(type="log", title="rows (log scale)")
    fig.update_yaxes(tickfont=dict(family="ui-monospace, Menlo, monospace", size=11, color=INK_2))
    return base_layout(fig, f"Largest categories in this entry (top {len(df)})", height=40 + 22 * len(df) + 60,
                       showlegend=False, margin=dict(l=260, r=50, t=40, b=44))
