import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Initial parameters (per instruction, N_cp = 1)
init_L1 = 1.0
init_L2 = 1.0
init_Lpipe = 10
init_short_rate = 0.02   # N_short / N_cp
init_long_rate  = 0.01   # N_long / N_cp
init_block_rate = 0.2    # N_block / N_cp
init_coincide_frac = 0.8 # N_coincide / N_short

fig, ax = plt.subplots(figsize=(7, 5))
plt.subplots_adjust(left=0.15, bottom=0.45)
ax.axis('off')
text_display = ax.text(0.05, 0.95, '', transform=ax.transAxes, va='top', fontfamily='monospace')

# Sliders
sliders = {}
params = [
    ('L1 (cycles)', init_L1, 0.5, 5.0),
    ('L2 (cycles)', init_L2, 0.5, 5.0),
    ('Lpipe (cycles)', init_Lpipe, 5, 20),
    ('short_rate', init_short_rate, 0.0, 0.1),
    ('long_rate', init_long_rate, 0.0, 0.05),
    ('block_rate', init_block_rate, 0.05, 0.5),
    ('coincide_frac', init_coincide_frac, 0.0, 1.0)
]

for i, (label, init_val, vmin, vmax) in enumerate(params):
    ax_slider = plt.axes([0.15, 0.35 - i*0.035, 0.65, 0.02])
    sliders[label] = Slider(ax_slider, label, vmin, vmax, valinit=init_val, valfmt='%.3f')

def update(val):
    L1 = sliders['L1 (cycles)'].val
    L2 = sliders['L2 (cycles)'].val
    Lpipe = sliders['Lpipe (cycles)'].val
    short_rate = sliders['short_rate'].val
    long_rate = sliders['long_rate'].val
    block_rate = sliders['block_rate'].val
    coincide_frac = sliders['coincide_frac'].val

    # Normalized per instruction (N_cp = 1)
    N_short = short_rate
    N_long = long_rate
    N_block = block_rate
    N_coincide = coincide_frac * N_short

    # Correct‑path cycles (Eq. from paper, assuming L1 < L2)
    T_cp = (N_block - N_coincide) * max(1, L1) + N_short * L2

    # Wrong‑path cycles
    T_wp = N_long * (L2 + Lpipe)

    # Derived metrics
    IPC_cbp = 1.0 / T_cp
    CPI_cbp = T_wp   # because N_cp = 1
    total_cycles_per_inst = T_cp + T_wp
    IPC_overall = 1.0 / total_cycles_per_inst

    # Display
    text_str = (f"L1 = {L1:.2f}   L2 = {L2:.2f}   Lpipe = {Lpipe:.1f}\n"
                f"N_short = {N_short:.4f}   N_long = {N_long:.4f}   N_block = {N_block:.3f}\n"
                f"T_cp = {T_cp:.4f} cycles/instruction\n"
                f"T_wp = {T_wp:.4f} cycles/instruction\n"
                f"IPC_cbp = {IPC_cbp:.4f}\n"
                f"CPI_cbp = {CPI_cbp:.4f}\n"
                f"Overall IPC = {IPC_overall:.4f}")
    text_display.set_text(text_str)
    fig.canvas.draw_idle()

for slider in sliders.values():
    slider.on_changed(update)

update(None)
plt.show()
