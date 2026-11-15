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
from helpers import get_cl, get_label

###########  USER OPTIONS ################
PLOT_AP = True
PLOT_ION_CURRENTS = True
PLOT_MARKOV_STATES = True
PLOT_ALL = True
##########################################

# read files
# either supply a filename on the command line, or enter when prompted
if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = input('Enter data filepath: ').strip()
data = np.loadtxt(filename, ndmin=2)
if data.shape[1] != 9:
    raise ValueError('Expected columns: t, v, IKr, ICaL, C1, C2, C3, I, O')
# extract data
CL = get_cl(filename)  # cycle length (ms)
t, v, IKr, ICaL, closed1, closed2, closed3, inactivated, open_state = data.T
# sum closed states
closed = closed1 + closed2 + closed3

if PLOT_AP:
    # Default: plot the full voltage trace
    mask = t < CL
    plt.figure(num=1)
    plt.plot(t, v) 
    # plt.plot(t[mask], v[mask], linewidth=2)  # uncomment to plot first AP only
    plt.ylabel('Membrane Potential (mV)')
    plt.xlabel('Time (ms)')

if PLOT_ION_CURRENTS:
    # plot AP and IKr, ICaL channel currents
    fig, axes = plt.subplots(1, 3, num=2, figsize=(9, 4))
    for ax, values, label in zip(
                axes, [v, IKr, ICaL],
                ['Membrane Potential (mV)', 'IKr (pA/pF)', 'ICaL (pA/pF)']):
        ax.plot(t, values)
        ax.set_ylabel(label)
        ax.set_xlabel('Time (ms)')
    fig.tight_layout()

if PLOT_MARKOV_STATES:
    # plot markov model open, inactive and closed states
    # make one panel for each condition
    fig, ax = plt.subplots(num=3)
    ax.set_title('IKr States')
    ax.plot(t, open_state, 'r-', label='Open')
    ax.set_ylabel('Open State Occupancy')
    right = ax.twinx()
    right.plot(t, closed, 'g-', label='Closed')
    right.plot(t, inactivated, 'b-', label='Inactivated')
    right.set_ylabel('Closed / Inactivated State Occupancy')
    ax.set_xlabel('Time (ms)')
    fig.legend(loc='lower center', ncol=3, frameon=False)
    fig.tight_layout(rect=[0, 0.06, 1, 1])

if PLOT_ALL:
    mask = t < CL
    fig, axes = plt.subplots(7, 1, figsize=(4, 10), sharex=True)
    for ax, values, label in zip(axes, [v, IKr, closed1, closed2, closed3, open_state, inactivated],
                               ['AP', 'IKr', 'C1', 'C2', 'C3', 'O', 'I']):
        ax.plot(t[mask], values[mask])
        ax.set_ylabel(label)
        ax.set_ylim((-100, 50) if label == 'AP' else (-0.1, 1))
    axes[-1].set_xlabel('Time (ms)')
    fig.tight_layout()

# save figure
plot_dir = Path(filename).parent / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)

for number in plt.get_fignums():
    fig = plt.figure(number)
    output = plot_dir / f"SingleCell_{get_label(filename)}_figure_{number}.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print(f"Saved {output}")

plt.close("all")
