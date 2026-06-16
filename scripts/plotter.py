import numpy as np
import plotly.graph_objects as go
import plotly.subplots as sp
import argparse

# =========================
# MODE FLAG
# =========================
parser = argparse.ArgumentParser()
parser.add_argument("--mode", type=str, default="3d", choices=["3d", "2d"])
args = parser.parse_args()

MODE = args.mode


# =========================
# VFS MODEL
# =========================
alpha = 1.625
beta = 16.64
gamma = 3.2

IPC_cbp_star = 0.6
CPI_cbp_star = 0.05
mu = 0.05

WPI_star = IPC_cbp_star * CPI_cbp_star
lambda_val = 1.0 / (1.0 + WPI_star / 2.0) - mu


def compute_fn(x, y, z):
    IPC_ratio = (1.0 + WPI_star) / (1.0/x + WPI_star * y)
    WPI = x * y * WPI_star
    EPI_ratio = (lambda_val * (IPC_ratio ** gamma) + mu * z) * (1.0 + WPI / 2.0)

    sqrt_term = np.sqrt(1.0 + beta / EPI_ratio)
    VFS = IPC_ratio * alpha * (1.0 - 2.0 / (1.0 + sqrt_term / IPC_ratio))
    return VFS


# =========================
# DOMAIN
# =========================
x_range = (0.4, 1.6)
y_range = (0.4, 1.6)
z_range = (0.4, 1.6)

n = 45

x = np.linspace(*x_range, n)
y = np.linspace(*y_range, n)
z = np.linspace(*z_range, n)


# ============================================================
# 3D MODE (TRUE ISOSURFACES)
# ============================================================
if MODE == "3d":

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    F = compute_fn(X, Y, Z)

    levels = np.percentile(F, [90, 95, 98, 99, 100])

    fig = go.Figure()

    fig.add_trace(go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=F.flatten(),

        isomin=levels[0],
        isomax=levels[-1],
        surface_count=5,

        opacity=0.25,
        colorscale="Viridis",
        caps=dict(x_show=False, y_show=False, z_show=False),
    ))

    # -------- optimization path --------
    def grad(p, eps=1e-3):
        x0, y0, z0 = p

        return np.array([
            compute_fn(x0+eps, y0, z0) - compute_fn(x0-eps, y0, z0),
            compute_fn(x0, y0+eps, z0) - compute_fn(x0, y0-eps, z0),
            compute_fn(x0, y0, z0+eps) - compute_fn(x0, y0, z0-eps),
        ]) / (2 * eps + 1e-12)

    p = np.array([
        np.mean(x_range),
        np.mean(y_range),
        np.mean(z_range)
    ])

    path = [p.copy()]
    momentum = np.zeros(3)

    lr = 0.03
    beta_m = 0.85

    for _ in range(80):
        g = grad(p)
        momentum = beta_m * momentum + (1 - beta_m) * g
        step = np.clip(lr * momentum, -0.1, 0.1)
        p = p + step
        path.append(p.copy())

    path = np.array(path)

    fig.add_trace(go.Scatter3d(
        x=path[:, 0],
        y=path[:, 1],
        z=path[:, 2],
        mode="lines+markers",
        line=dict(width=6),
        marker=dict(size=3)
    ))

    fig.update_layout(
        title="VFS Optimization Landscape (3D Isosurfaces)",
        scene=dict(
            xaxis_title="IPC",
            yaxis_title="CPI",
            zaxis_title="EPI"
        ),
        margin=dict(l=0, r=0, t=40, b=0)
    )

    fig.show()


# ============================================================
# 2D MODE (CLEAN + NO TRANSPOSE BUGS)
# ============================================================
elif MODE == "2d":

    percentiles = [90, 95, 98, 99, 100]

    # ---- slices (fixed-axis views) ----
    def slice_xy(z0):
        X, Y = np.meshgrid(x, y, indexing="xy")
        Z = np.full_like(X, z0)
        F2 = compute_fn(X, Y, Z)
        return X, Y, F2, "IPC_ratio", "CPI_ratio"

    def slice_xz(y0):
        X, Z = np.meshgrid(x, z, indexing="xy")
        Y = np.full_like(X, y0)
        F2 = compute_fn(X, Y, Z)
        return X, Z, F2, "IPC_ratio", "EPI_ratio"

    def slice_yz(x0):
        Y, Z = np.meshgrid(y, z, indexing="xy")
        X = np.full_like(Y, x0)
        F2 = compute_fn(X, Y, Z)
        return Y, Z, F2, "CPI_ratio", "EPI_ratio"


    x0 = np.mean(x_range)
    y0 = np.mean(y_range)
    z0 = np.mean(z_range)

    fig = sp.make_subplots(
        rows=1, cols=3,
        subplot_titles=[
            "(EPI_ratio fixed to 1)",
            "(CPI_ratio fixed to 1)",
            "(IPC_ratio fixed to 1)"
        ]
    )

    slices = [
        slice_xy(z0),
        slice_xz(y0),
        slice_yz(x0)
    ]

    for i, (X2, Y2, F2, xl, yl) in enumerate(slices):

        levels = np.percentile(F2, percentiles)

        fig.add_trace(
            go.Contour(
                x=x,          # FIXED: true axis vectors
                y=y if i == 0 else (z if i == 1 else z),
                z=F2,

                colorscale="Viridis",
                contours=dict(
                    coloring="heatmap",
                    showlabels=True
                ),

                showscale=(i == 2),
                colorbar=dict(title="VFS")
            ),
            row=1, col=i+1
        )

        fig.update_xaxes(title_text=xl, row=1, col=i+1)
        fig.update_yaxes(title_text=yl, row=1, col=i+1)


    fig.update_layout(
        title="",
        height=450,
        margin=dict(l=10, r=10, t=50, b=10)
    )

    fig.show()

# fig.write_image("my_vfs_plot.png")   # Raster (recommended for general use)
