"""
Create pseudo ECG figure and/or movie of action potential
from oneD transmural simulation.

Authors: Xixian Jiang, Barton Wei, Jorik Stoop
"""
import sys
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")  # Save figures without opening windows.
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
import re

###########  USER OPTIONS ################
PLOT_HEATMAP = False
PLOT_SURFACE = False
PLOT_VOLTAGE = False
PLOT_PECG = False
SAVE_GIF = False
SAVE_VIDEO = False
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

# Load data: time in the first column, then one voltage per cell.
filename = sys.argv[1] if len(sys.argv) > 1 else input('1D data file: ').strip()
data = np.loadtxt(filename, ndmin=2)
if data.shape[1] < 2 or len(data) < 2:
    raise ValueError('Expected at least two rows: time followed by cell voltages')
CL = get_cl(filename)  # cycle length (ms)

# Optional trimming:
# data = data[(data[:, 0] >= 1000) & (data[:, 0] <= 3000)]
t = data[:, 0]
voltages = data[:, 1:]
numCells = voltages.shape[1]
vmin, vmax = np.min(voltages), np.max(voltages)

if PLOT_HEATMAP:
    fig, ax = plt.subplots(num=1)
    image = ax.imshow(voltages.T, extent=(t[0], t[-1], numCells + 0.5, 0.5),
                    aspect='auto', origin='upper', cmap='jet', vmin=vmin, vmax=vmax)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Position')
    ax.set_title('Voltage over Time in Transmural Strand')
    fig.colorbar(image, ax=ax)

if PLOT_SURFACE:
    fig = plt.figure(2)
    ax = fig.add_subplot(111, projection='3d')
    times, positions = np.meshgrid(t, np.arange(1, numCells + 1))
    surface = ax.plot_surface(times, positions, voltages.T, cmap='jet',
                              vmin=vmin, vmax=vmax, linewidth=0)
    ax.invert_yaxis()
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Position')
    ax.set_title('Voltage over Time in Transmural Strand')
    fig.colorbar(surface, ax=ax)

if PLOT_VOLTAGE:
    cells = [10]
    if any(cell < 1 or cell > numCells for cell in cells):
        raise ValueError(f'Cell numbers must be between 1 and {numCells}')
    vm_plot = np.zeros((len(t), len(cells)))
    for i, cell in enumerate(cells):
        vm_plot[:, i] = voltages[:, cell - 1]
    numPlots = len(cells)
    mask = (t > 2000) & (t < 2500)
    # mask = t > (t[-1] - CL)  # alternatively: final cycle
    fig, axes = plt.subplots(numPlots, 1, num=3, squeeze=False)
    for i in range(numPlots):
        ax = axes[i, 0]
        ax.plot(t[mask], vm_plot[mask, i], linewidth=2)
        ax.set_title(f'Cell {cells[i]}')
        ax.set_ylabel('Vm (mV)')
        if i == 0:
            ax.set_xlabel('Time (ms)')
    fig.tight_layout()


if PLOT_PECG:
    Ci = 0.16
    Ce = 0.21
    a = 0.001
    dx = 0.015
    pECG = np.zeros(len(t))
    first_j, last_j = 25, 176
    if data.shape[1] < last_j:
        raise ValueError('Pseudo ECG range requires at least 176 data columns')
    for i in range(len(t)):
        for j in range(first_j, last_j + 1):
            r = 2 + (100 - j) * dx
            Vm = data[i, j - 1]
            Vm_1 = data[i, j - 2]
            pECG[i] += ((a**2 * Ci) / (4 * Ce)) * (Vm_1 - Vm) * (1 / r) * dx
    mask = t > (t[-1] - CL)
    plt.figure(4)
    plt.plot(t[mask], pECG[mask], linewidth=2)
    plt.xlabel('Time (ms)')
    plt.ylabel('Pseudo ECG / Extracellular unipolar potentials (mV)')

# Pillow writes GIF; MP4 also requires FFmpeg.
if SAVE_GIF or SAVE_VIDEO:
    def make_movie(output, compression_level, fps, writer):
        fig, ax = plt.subplots(figsize=(10, 2))
        strip = np.repeat(voltages[0:1], 3, axis=0)
        image = ax.imshow(strip, cmap='jet', vmin=vmin, vmax=vmax, aspect='auto')
        fig.colorbar(image, ax=ax, orientation='horizontal')
        ax.get_yaxis().set_visible(False)
        title = ax.set_title('')

        def update(i):
            image.set_data(np.repeat(voltages[i:i+1], 3, axis=0))
            title.set_text(f'Voltage in Transmural Strand (mV), t={t[i] / 1000:.3f}s')
            return image, title

        movie = FuncAnimation(fig, update, frames=range(0, len(t), compression_level),
                              interval=1000 / fps, blit=False)
        movie.save(output, writer=writer(fps=fps))
        plt.close(fig)

    if SAVE_GIF:
        make_movie('Animated.gif', 3, 10, PillowWriter)
    if SAVE_VIDEO:
        if not FFMpegWriter.isAvailable():
            raise RuntimeError('MP4 export requires FFmpeg installed and available on PATH')
        make_movie('movie.mp4', 5, 100, FFMpegWriter)

source = Path(filename)
plot_dir = source.parent / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)

for number in plt.get_fignums():
    fig = plt.figure(number)
    output = plot_dir / f"{source.stem}_figure_{number}.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print(f"Saved {output}")

plt.close("all")
