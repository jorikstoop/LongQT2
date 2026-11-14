"""
Create plots comparing single cell action potentials,
ion channel currents, and Markov states for two conditions.

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
PLOT_AP = True
PLOT_IKR = True
PLOT_MARKOV_STATES = True
PLOT_ALL = True
##########################################

# helper functions to extract info from filename
def get_cl(filename):
    """Extract cycle length from filename"""
    match = re.search(r"(?:^|_)CL_(\d+(?:\.\d+)?)(?=_|$)",
                      Path(filename).stem,
                     flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"Missing CL value in filename: {filename}")
    return float(match.group(1))


def get_label(filename):
    """Extract model and treatment(s) from filename"""
    tokens = Path(filename).stem.upper().split("_")
    model = next((token for token in tokens
                  if token in {"WT", "T474I", "R56Q"}),
                  "")
    treatment_names = {
        "ISO": "Isoproterenol",
        "KSUPP": "K+_supplement",
        "VERA": "Verapamil",
        "NODRUG": "No_drug",
    }
    treatments = [treatment_names[token]
                  for token in tokens
                  if token in treatment_names]
    if not treatments:
        raise ValueError(f"Missing treatment tag in filename: {filename}")
    treatment = "+".join(treatments)
    return f"{model}—{treatment}" if model else treatment

# read files
file1 = sys.argv[1] if len(sys.argv) > 1 else input('Enter file1 path: ').strip()
file2 = sys.argv[2] if len(sys.argv) > 2 else input('Enter file2 path: ').strip()
data1 = np.loadtxt(file1, ndmin=2)
data2 = np.loadtxt(file2, ndmin=2)
if data1.shape[1] != 9 or data2.shape[1] != 9:
    raise ValueError('Each file must contain t, v, IKr, ICaL, C1, C2, C3, I, O')
# extract data
t, v1, IKr1, ICaL1, closed11, closed21, closed31, inactivated1, open1 = data1.T
t2, v2, IKr2, ICaL2, closed12, closed22, closed32, inactivated2, open2 = data2.T
# sum closed states
closed1 = closed11 + closed21 + closed31
closed2 = closed12 + closed22 + closed32

labels = [get_label(file1), get_label(file2)]

if PLOT_AP:
    # plot first cycle
    CL1 = get_cl(file1)
    CL2 = get_cl(file2)
    mask = t < CL1
    mask2 = t2 < CL2
    plt.figure(num=1)
    plt.plot(t[mask], v1[mask], linewidth=2, label=labels[0])
    plt.plot(t2[mask2], v2[mask2], linewidth=2, label=labels[1])
    plt.ylabel('Membrane Potential (mV)')
    plt.xlabel('Time (ms)')
    plt.legend()

if PLOT_IKR:
    # plot IKr channel current
    fig, axes = plt.subplots(1, 2, num=2, figsize=(9, 4))
    for ax, y1, y2, ylabel in zip(
                axes, [v1, IKr1], [v2, IKr2],
                ['Membrane Potential (mV)', 'IKr (pA/pF)']):
        ax.plot(t, y1, color='C0', label=labels[0])
        ax.plot(t2, y2, color='C1', label=labels[1])
        ax.set_ylabel(ylabel)
        ax.set_xlabel('Time (ms)')
    # make legend
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, legend_labels,
        loc='lower center', ncol=2, frameon=False,
    )
    fig.tight_layout(rect=[0, 0.10, 1, 1])

if PLOT_MARKOV_STATES:
    # plot markov model open, inactive and closed states
    # make one panel for each condition
    fig, axes = plt.subplots(2, 1, num=3, figsize=(10, 10))
    state_handles = None

    for ax, time, opened, closed, inactivated, condition in zip(
            axes, [t, t2], [open1, open2], [closed1, closed2],
            [inactivated1, inactivated2], labels):
        ax.set_title(condition)
        open_line, = ax.plot(time, opened, color='r', label='Open')
        ax.set_ylabel('Open State Occupancy')
        ax.set_xlabel('Time (ms)')
        right = ax.twinx()
        closed_line, = right.plot(time, closed,
                                  color='g', label='Closed')
        inactive_line, = right.plot(time, inactivated,
                                    color='b', label='Inactivated')
        right.set_ylabel('Closed / Inactivated State Occupancy')
        if state_handles is None:
            state_handles = [open_line, closed_line, inactive_line]
    # add legend
    fig.legend(handles=state_handles,
               loc='lower center', ncol=3, frameon=False)


if PLOT_ALL:
    mask = t < CL1
    mask2 = t2 < CL2
    fig, axes = plt.subplots(7, 1, figsize=(4, 10), sharex=True)
    for ax, y1, y2, label in zip(axes,
            [v1, IKr1, closed11, closed21, closed31, open1, inactivated1],
            [v2, IKr2, closed12, closed22, closed32, open2, inactivated2],
            ['AP', 'IKr', 'C1', 'C2', 'C3', 'O', 'I']):
        ax.plot(t[mask], y1[mask], color='C0', label=labels[0])
        ax.plot(t2[mask2], y2[mask2], color='C1', label=labels[1])
        ax.set_ylabel(label)
        ax.set_ylim((-100, 50) if label == 'AP' else (-0.1, 1))
    axes[-1].set_xlabel('Time (ms)')
    # create legend
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, legend_labels,
               loc='lower center', ncol=1, frameon=False)
    fig.tight_layout(rect=[0, 0.06, 1, 1])

# save figure
plot_dir = Path(file1).parent / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)
prefix = f"{labels[0]}_vs_{labels[1]}"

for number in plt.get_fignums():
    fig = plt.figure(number)
    output = plot_dir / f"{prefix}_figure_{number}.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print(f"Saved {output}")

plt.close("all")
