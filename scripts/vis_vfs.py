import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters (same as before, but now z is variable)
alpha = 1.625
beta = 16.64
gamma = 3.2

# Reference values (replace with real ones!)
IPC_cbp_star = 0.6
CPI_cbp_star = 0.05
WPI_star = IPC_cbp_star * CPI_cbp_star
mu = 0.05
lambda_val = 1/(1 + WPI_star/2) - mu

# Define grid ranges
x_vals = np.linspace(0.7, 1.3, 15)
y_vals = np.linspace(0.6, 1.4, 15)
z_vals = np.linspace(0.8, 1.5, 10)   # EPI_cbp ratio

# Create meshgrid
X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')

# Compute VFS for each point
VFS = np.zeros_like(X)

for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        for k in range(X.shape[2]):
            x = X[i, j, k]
            y = Y[i, j, k]
            z = Z[i, j, k]
            
            # IPC ratio
            IPC_ratio = (1 + WPI_star) / (1/x + WPI_star * y)
            # EPI ratio (core)
            EPI_ratio = (lambda_val * (IPC_ratio ** gamma) + mu * z) * (1 + (WPI_star/2) * x * y)
            # VFS
            sqrt_term = np.sqrt(1 + beta / EPI_ratio)   # because EPI_star/EPI = 1/EPI_ratio
            VFS[i, j, k] = IPC_ratio * alpha * (1 - 2 / (1 + sqrt_term / IPC_ratio))

# Flatten for scatter plot
points = np.array([X.flatten(), Y.flatten(), Z.flatten()]).T
vfs_flat = VFS.flatten()

# Plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(points[:,0], points[:,1], points[:,2], c=vfs_flat, cmap='viridis', s=10, alpha=0.7)
ax.set_xlabel('x = IPC_cbp / IPC_cbp*')
ax.set_ylabel('y = CPI_cbp / CPI_cbp*')
ax.set_zlabel('z = EPI_cbp / EPI_cbp*')
fig.colorbar(sc, label='VFS')
plt.title('VFS as color in (x, y, z) space')
plt.show()
