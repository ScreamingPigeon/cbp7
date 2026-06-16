import numpy as np
import plotly.graph_objects as go

# ============================================================
# CBP VFS PARAMETERS
# ============================================================

ALPHA = 1.625
BETA = 4 * ALPHA / (ALPHA - 1)**2
GAMMA = 2 / (ALPHA - 1)

MU = 0.05

# Reference predictor
IPC_STAR = 0.6
CPI_STAR = 0.05

WPI_STAR = IPC_STAR * CPI_STAR

LAMBDA = 1 / (1 + WPI_STAR / 2) - MU

# ============================================================
# VFS MODEL
# ============================================================

def compute_vfs(x, y, z):
    """
    x = IPC / IPC*
    y = CPI / CPI*
    z = EPI / EPI*
    """

    WPI = WPI_STAR * x * y

    speedup = (
        x
        * (1 + WPI_STAR)
        / (1 + WPI)
    )

    normalized_epi = (
        (
            MU * z
            + LAMBDA * speedup**GAMMA
        )
        * (1 + WPI / 2)
    )

    vfs = (
        speedup
        * ALPHA
        * (
            1
            - 2
            / (
                1
                + np.sqrt(
                    1
                    + BETA /
                    (
                        speedup *
                        normalized_epi
                    )
                )
            )
        )
    )

    return vfs


# ============================================================
# DESIGN SPACE
# ============================================================

x_vals = np.linspace(0.7, 1.3, 20)
y_vals = np.linspace(0.6, 1.4, 20)
z_vals = np.linspace(0.8, 1.5, 15)

X, Y, Z = np.meshgrid(
    x_vals,
    y_vals,
    z_vals,
    indexing="ij"
)

V = compute_vfs(X, Y, Z)

print("VFS range")
print(V.min(), V.max())

# ============================================================
# FLATTEN
# ============================================================

x_flat = X.flatten()
y_flat = Y.flatten()
z_flat = Z.flatten()

vfs_flat = V.flatten()

# ============================================================
# REFERENCE POINT
# ============================================================

ref_vfs = compute_vfs(1.0, 1.0, 1.0)

# ============================================================
# LOCAL SENSITIVITY
# ============================================================

eps = 0.01

dIPC = (
    compute_vfs(1 + eps, 1, 1)
    - ref_vfs
) / eps

dCPI = (
    compute_vfs(1, 1 + eps, 1)
    - ref_vfs
) / eps

dEPI = (
    compute_vfs(1, 1, 1 + eps)
    - ref_vfs
) / eps

print()
print("Local sensitivity around reference:")
print(f"dVFS/dIPC ≈ {dIPC:.4f}")
print(f"dVFS/dCPI ≈ {dCPI:.4f}")
print(f"dVFS/dEPI ≈ {dEPI:.4f}")

# ============================================================
# PLOT
# ============================================================

fig = go.Figure()

# VFS cloud

fig.add_trace(
    go.Scatter3d(
        x=x_flat,
        y=y_flat,
        z=z_flat,
        mode="markers",
        marker=dict(
            size=3,
            color=vfs_flat,
            colorscale="Viridis",
            opacity=0.55,
            colorbar=dict(
                title="VFS"
            ),
        ),
        name="VFS space"
    )
)

# Reference predictor

fig.add_trace(
    go.Scatter3d(
        x=[1],
        y=[1],
        z=[1],
        mode="markers+text",
        text=[f"Reference\nVFS={ref_vfs:.3f}"],
        textposition="top center",
        marker=dict(
            size=10,
            color="red"
        ),
        name="Reference"
    )
)

# Sensitivity vector

vec = np.array([
    dIPC,
    dCPI,
    dEPI
])

vec = vec / np.linalg.norm(vec)

scale = 0.35

fig.add_trace(
    go.Scatter3d(
        x=[1, 1 + scale * vec[0]],
        y=[1, 1 + scale * vec[1]],
        z=[1, 1 + scale * vec[2]],
        mode="lines+markers",
        line=dict(width=8),
        marker=dict(size=4),
        name="Steepest VFS increase"
    )
)

fig.update_layout(
    title="VFS Design Space",
    scene=dict(
        xaxis_title="IPC / IPC*",
        yaxis_title="CPI / CPI*",
        zaxis_title="EPI / EPI*",
        aspectmode="cube"
    ),
    margin=dict(
        l=0,
        r=0,
        b=0,
        t=50
    )
)

fig.show()
