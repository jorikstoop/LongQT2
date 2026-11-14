"""
Create plots comparing single cell action potentials
for condition without vs with drug treatment.

Authors: Xixian Jiang, Barton Wei, Jorik Stoop
"""
import sys
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import re

###########  USER OPTIONS ################
PLOT_CURRENTS = False
PLOT_STATES = False
PLOT_ALL_STATES = False
##########################################

from pathlib import Path
import re


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


def get_label(filename):
    """Extract model and treatment(s) for plot legends."""
    tokens = Path(filename).stem.upper().split("_")

    model = next(
        (token for token in tokens if token in {"WT", "T474I", "R56Q"}),
        "",
    )

    treatment_names = {
        "ISO": "Isoproterenol",
        "KSUPP": "K+ supplementation",
        "VERA": "Verapamil",
        "NODRUG": "No drug",
    }

    treatments = [
        treatment_names[token]
        for token in tokens
        if token in treatment_names
    ]

    if not treatments:
        raise ValueError(f"Missing treatment tag in filename: {filename}")

    treatment = " + ".join(treatments)
    return f"{model} — {treatment}" if model else treatment


# File 1 -- no drug; file 2 -- drug application.
file1 = sys.argv[1] if len(sys.argv) > 1 else input('No-drug data file: ').strip()
file2 = sys.argv[2] if len(sys.argv) > 2 else input('Drug data file: ').strip()
data1 = np.loadtxt(file1, ndmin=2)
data2 = np.loadtxt(file2, ndmin=2)
if data1.shape[1] != 9 or data2.shape[1] != 9:
    raise ValueError('Each file must contain t, v, IKr, ICaL, C1, C2, C3, I, O')
t, v1, IKr1, ICaL1, closed11, closed21, closed31, inactivated1, open1 = data1.T
t2, v2, IKr2, ICaL2, closed12, closed22, closed32, inactivated2, open2 = data2.T
closed1 = closed11 + closed21 + closed31
closed2 = closed12 + closed22 + closed32

labels = [get_label(file1), get_label(file2)]


if PLOT_CURRENTS:
    fig, axes = plt.subplots(1, 2, num=1, figsize=(9, 4))
    for ax, y1, y2, label in zip(axes, [v1, IKr1], [v2, IKr2],
                                ['Membrane Potential (mV)', 'IKr (pA/pF)']):
        ax.plot(t, y1, 'b', label=labels[0])
        ax.plot(t2, y2, 'r', label=labels[1])
        ax.set_ylabel(label)
        ax.set_xlabel('Time (ms)')
        ax.legend()
    fig.tight_layout()

if PLOT_STATES:
    fig, ax = plt.subplots(num=2)
    ax.set_title('IKr States')
    ax.plot(t, open1, 'r-', label='Control Open')
    ax.plot(t2, open2, 'r--', label='Drug Open')
    ax.set_ylabel('Open State Occupancy')
    right = ax.twinx()
    right.plot(t, closed1, 'g-', label='Control Closed')
    right.plot(t2, closed2, 'g--', label='Drug Closed')
    right.plot(t, inactivated1, 'b-', label='Control Inactivated')
    right.plot(t2, inactivated2, 'b--', label='Drug Inactivated')
    right.set_ylabel('Closed / Inactivated State Occupancy')
    ax.set_xlabel('Time (ms)')
    fig.legend()

# AP: first cycle. Each file uses its own timestamps.
CL1 = get_cl(file1)
CL2 = get_cl(file2)
mask = t < CL1
mask2 = t2 < CL2
plt.figure(3)
plt.plot(t[mask], v1[mask], linewidth=2, label=labels[0])
plt.plot(t2[mask2], v2[mask2], linewidth=2, label=labels[1])
plt.ylabel('Membrane Potential (mV)')
plt.xlabel('Time (ms)')
plt.legend()

if PLOT_ALL_STATES:
    mask = t < 600
    mask2 = t2 < 600
    fig, axes = plt.subplots(7, 1, figsize=(4, 10), sharex=True)
    for ax, y1, y2, label in zip(axes,
            [v1, IKr1, closed11, closed21, closed31, open1, inactivated1],
            [v2, IKr2, closed12, closed22, closed32, open2, inactivated2],
            ['AP', 'IKr', 'C1', 'C2', 'C3', 'O', 'I']):
        ax.plot(t[mask], y1[mask], label=labels[0])
        ax.plot(t2[mask2], y2[mask2], label=labels[1])
        ax.set_ylabel(label)
        ax.set_ylim((-100, 50) if label == 'AP' else (-0.1, 1))
        ax.legend()
    axes[-1].set_xlabel('Time (ms)')
    fig.tight_layout()

# save figure
plot_dir = Path(file1).parent / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)
prefix = f"{Path(file1).stem}_vs_{Path(file2).stem}"

for number in plt.get_fignums():
    fig = plt.figure(number)
    output = plot_dir / f"{prefix}_figure_{number}.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print(f"Saved {output}")

plt.close("all")
