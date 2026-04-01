#!/usr/bin/env python3
"""
Interactive VFS explorer with sliders.
Computes VFS directly from the formula.

Usage: python3 scripts/vfs_explorer.py
Opens in browser at http://localhost:8050
"""

import math
import numpy as np

try:
    import dash
    from dash import dcc, html
    from dash.dependencies import Input, Output
    import plotly.graph_objects as go
except ImportError:
    print("pip install dash plotly")
    exit(1)

# ---- VFS formula (from vfs.py + predictor_metrics.py) ----
IPCcbp0 = 8
CPIcbp0 = 0.0315
EPIcbp0 = 1000
ALPHA = 1.625
BETA = 4 * ALPHA / (ALPHA - 1) ** 2
GAMMA = 2 / (ALPHA - 1)
cbp_energy_ratio = 0.05
MISPREDICTION_PENALTY = 8

def compute_vfs(mpki, p2_ceil, epi, ipc_base=4.0):
    """
    Simplified VFS from MPKI, P2 ceil, EPI.
    ipc_base: base IPC (instructions per prediction block / p2_ceil).
    """
    mpi = mpki / 1000.0
    ipc = ipc_base / max(1, p2_ceil)
    cpi = mpi * (MISPREDICTION_PENALTY + p2_ceil)

    WPI0 = IPCcbp0 * CPIcbp0
    WPI = ipc * cpi
    speedup = (ipc / IPCcbp0) * (1 + WPI0) / (1 + WPI)
    LAMBDA = 1 / (1 + WPI0 / 2) - cbp_energy_ratio
    normalizedEPI = ((epi / EPIcbp0) * cbp_energy_ratio + LAMBDA * speedup ** GAMMA) * (1 + WPI / 2)
    vfs = speedup * ALPHA * (1 - 2 / (1 + math.sqrt(1 + BETA / (speedup * normalizedEPI))))
    return vfs

app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("VFS Explorer"),

    html.Div([
        html.Div([
            html.H3("Fixed Parameters"),
            html.Label("IPC base (instr per block cycle):"),
            dcc.Slider(id="ipc-base", min=1, max=10, step=0.5, value=4.0,
                      marks={i: str(i) for i in range(1, 11)}),
            html.Label("X axis:"),
            dcc.Dropdown(id="x-axis", options=[
                {"label": "MPKI", "value": "mpki"},
                {"label": "EPI (fJ)", "value": "epi"},
                {"label": "P2 latency (ceil cycles)", "value": "p2"},
            ], value="mpki"),
            html.Label("Y axis:"),
            dcc.Dropdown(id="y-axis", options=[
                {"label": "MPKI", "value": "mpki"},
                {"label": "EPI (fJ)", "value": "epi"},
                {"label": "P2 latency (ceil cycles)", "value": "p2"},
            ], value="epi"),
        ], style={"width": "30%", "display": "inline-block", "verticalAlign": "top", "padding": "20px"}),

        html.Div([
            html.H3("Sliders for fixed axis"),
            html.Label("MPKI (when not an axis):"),
            dcc.Slider(id="mpki-slider", min=1, max=20, step=0.5, value=5.85,
                      marks={i: str(i) for i in range(1, 21)}),
            html.Label("EPI (when not an axis):"),
            dcc.Slider(id="epi-slider", min=200, max=6000, step=100, value=1245,
                      marks={i: str(i) for i in range(0, 6001, 1000)}),
            html.Label("P2 ceil (when not an axis):"),
            dcc.Slider(id="p2-slider", min=1, max=5, step=1, value=2,
                      marks={i: str(i) for i in range(1, 6)}),
        ], style={"width": "30%", "display": "inline-block", "verticalAlign": "top", "padding": "20px"}),
    ]),

    html.Div([
        dcc.Graph(id="vfs-heatmap", style={"height": "700px"}),
    ]),

    html.Div([
        html.H3("Reference Points"),
        html.Table([
            html.Tr([html.Th("Config"), html.Th("MPKI"), html.Th("P2 ceil"), html.Th("EPI"), html.Th("VFS")]),
            html.Tr([html.Td("Direct Tage (best)"), html.Td("5.85"), html.Td("2"), html.Td("1245"), html.Td(id="ref-direct")]),
            html.Tr([html.Td("TageAhead (current)"), html.Td("13.57"), html.Td("2"), html.Td("3049"), html.Td(id="ref-ahead")]),
        ], style={"margin": "20px"}),
    ]),
])

@app.callback(
    [Output("vfs-heatmap", "figure"),
     Output("ref-direct", "children"),
     Output("ref-ahead", "children")],
    [Input("x-axis", "value"),
     Input("y-axis", "value"),
     Input("mpki-slider", "value"),
     Input("epi-slider", "value"),
     Input("p2-slider", "value"),
     Input("ipc-base", "value")]
)
def update_plot(x_axis, y_axis, mpki_fixed, epi_fixed, p2_fixed, ipc_base):
    ranges = {
        "mpki": np.linspace(1, 20, 50),
        "epi": np.linspace(200, 6000, 50),
        "p2": np.array([1, 2, 3, 4, 5]),
    }
    labels = {
        "mpki": "MPKI",
        "epi": "EPI (fJ)",
        "p2": "P2 latency (ceil cycles)",
    }

    if x_axis == y_axis:
        y_axis = [k for k in ["mpki", "epi", "p2"] if k != x_axis][0]

    x_vals = ranges[x_axis]
    y_vals = ranges[y_axis]
    Z = np.zeros((len(y_vals), len(x_vals)))

    for i, yv in enumerate(y_vals):
        for j, xv in enumerate(x_vals):
            params = {"mpki": mpki_fixed, "epi": epi_fixed, "p2": p2_fixed}
            params[x_axis] = xv
            params[y_axis] = yv
            Z[i, j] = compute_vfs(params["mpki"], params["p2"], params["epi"], ipc_base)

    fig = go.Figure(data=go.Heatmap(
        x=x_vals, y=y_vals, z=Z,
        colorscale="Viridis", zmin=0.4, zmax=1.0,
        colorbar=dict(title="VFS"),
        hovertemplate=f"{labels[x_axis]}: %{{x:.1f}}<br>{labels[y_axis]}: %{{y:.1f}}<br>VFS: %{{z:.4f}}<extra></extra>",
    ))

    # Add contour lines
    fig.add_trace(go.Contour(
        x=x_vals, y=y_vals, z=Z,
        contours=dict(start=0.5, end=1.0, size=0.05, showlabels=True),
        line=dict(width=1, color="white"),
        showscale=False,
        opacity=0.5,
        hoverinfo="skip",
    ))

    # Mark reference points if they're on the axes
    refs = {
        "Direct Tage": {"mpki": 5.85, "epi": 1245, "p2": 2},
        "TageAhead": {"mpki": 13.57, "epi": 3049, "p2": 2},
    }
    markers = {"Direct Tage": "star", "TageAhead": "x"}
    colors = {"Direct Tage": "white", "TageAhead": "red"}

    for name, vals in refs.items():
        fig.add_trace(go.Scatter(
            x=[vals[x_axis]], y=[vals[y_axis]],
            mode="markers+text",
            marker=dict(size=15, color=colors[name], symbol=markers[name],
                       line=dict(width=2, color="black")),
            text=[name], textposition="top center",
            textfont=dict(color="white", size=12),
            name=name,
        ))

    fig.update_layout(
        title=f"VFS: {labels[x_axis]} vs {labels[y_axis]} (fixed: MPKI={mpki_fixed}, EPI={epi_fixed}, P2={p2_fixed})",
        xaxis_title=labels[x_axis],
        yaxis_title=labels[y_axis],
        width=900, height=700,
    )

    vfs_direct = compute_vfs(5.85, 2, 1245, ipc_base)
    vfs_ahead = compute_vfs(13.57, 2, 3049, ipc_base)

    return fig, f"{vfs_direct:.4f}", f"{vfs_ahead:.4f}"

if __name__ == "__main__":
    print("Opening VFS Explorer at http://localhost:8050")
    app.run(debug=False, port=8050)
