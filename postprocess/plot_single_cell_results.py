"""
Create plots of action potential and ion channel
outputs from single cell simulation.

Authors: Xixian Jiang, Barton Wei, Jorik Stoop
"""
import sys
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")  # Save figures without opening windows.
import matplotlib.pyplot as plt
import re

###########  USER OPTIONS ################
PLOT_CURRENTS = False
PLOT_STATES = False
PLOT_ALL_STATES = False
##########################################

def get_cl(filename):
    """Extract cycle length from a filename containing _CL_750.00."""
    match = re.search(
        r"(?:^|_)CL_(\d+(?:\.\d+)?)(?=_|$)",
        Path(filename).stem,
        flags=re.IGNORECASE,
    )

    if not match:
        raise ValueError(f"Missing CL value in filename: {filename}")

    return float(match.group(1))

# Load data. Supply a filename on the command line, or enter it at the prompt.
filename = sys.argv[1] if len(sys.argv) > 1 else input('Data file: ').strip()
data = np.loadtxt(filename, ndmin=2)
if data.shape[1] != 9:
    raise ValueError('Expected columns: t, v, IKr, ICaL, C1, C2, C3, I, O')
CL = get_cl(filename)  # cycle length (ms)

# Extract values (Python indices start at zero).
t, v, IKr, ICaL, closed1, closed2, closed3, inactivated, open_state = data.T
closed = closed1 + closed2 + closed3

if PLOT_CURRENTS:
    fig, axes = plt.subplots(1, 3, num=1, figsize=(9, 4))
    for ax, values, label in zip(axes, [v, IKr, ICaL],
                               ['Membrane Potential (mV)', 'IKr (pA/pF)', 'ICaL (pA/pF)']):
        ax.plot(t, values)
        ax.set_ylabel(label)
        ax.set_xlabel('Time (ms)')
    fig.tight_layout()

if PLOT_STATES:
    fig, ax = plt.subplots(num=2)
    ax.set_title('IKr States')
    ax.plot(t, open_state, 'r-', label='Open')
    ax.set_ylabel('Open State Occupancy')
    right = ax.twinx()
    right.plot(t, closed, 'g-', label='Closed')
    right.plot(t, inactivated, 'b-', label='Inactivated')
    right.set_ylabel('Closed / Inactivated State Occupancy')
    ax.set_xlabel('Time (ms)')
    fig.legend()

if PLOT_ALL_STATES:
    mask = t < 600  # first 600 ms; use actual timestamps
    fig, axes = plt.subplots(7, 1, figsize=(4, 10), sharex=True)
    for ax, values, label in zip(axes, [v, IKr, closed1, closed2, closed3, open_state, inactivated],
                               ['AP', 'IKr', 'C1', 'C2', 'C3', 'O', 'I']):
        ax.plot(t[mask], values[mask])
        ax.set_ylabel(label)
        ax.set_ylim((-100, 50) if label == 'AP' else (-0.1, 1))
    axes[-1].set_xlabel('Time (ms)')
    fig.tight_layout()

# Default: plot the full voltage trace, matching the active MATLAB code.
mask = t < CL
plt.figure(4)
plt.plot(t, v)
# plt.plot(t[mask], v[mask], linewidth=2)  # first AP only
plt.ylabel('Membrane Potential (mV)')
plt.xlabel('Time (ms)')

source = Path(filename)
plot_dir = source.parent / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)

for number in plt.get_fignums():
    fig = plt.figure(number)
    output = plot_dir / f"{source.stem}_figure_{number}.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print(f"Saved {output}")

plt.close("all")
