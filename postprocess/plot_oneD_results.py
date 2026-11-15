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
from matplotlib.animation import FuncAnimation, FFMpegWriter
from helpers import get_cl, get_label

###########  USER OPTIONS ################
PLOT_VOLTAGE_HEATMAP = True
PLOT_VOLTAGE_SELECT_CELLS = True
SELECT_CELLS = [20, 100, 180]
PLOT_PECG = True
SAVE_VIDEO = True
##########################################

# read files
if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = input('Enter 1D data filepath: ').strip()
data = np.loadtxt(filename, ndmin=2)
if data.shape[1] < 2 or len(data) < 2:
    raise ValueError('Expected at least two rows: time followed by cell voltages')
# extract data
CL = get_cl(filename)  # cycle length (ms)
# Optional trimming:
# data = data[(data[:, 0] >= 1000) & (data[:, 0] <= 3000)]
t = data[:, 0]
voltages = data[:, 1:]
numCells = voltages.shape[1]
vmin, vmax = np.min(voltages), np.max(voltages)

if PLOT_VOLTAGE_HEATMAP:
    fig, ax = plt.subplots(num=1)
    image = ax.imshow(voltages.T, extent=(t[0], t[-1], numCells + 0.5, 0.5),
                    aspect='auto', origin='upper', cmap='jet', vmin=vmin, vmax=vmax)
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Position')
    ax.set_title('Voltage over Time in Transmural Strand')
    fig.colorbar(image, ax=ax)

if PLOT_VOLTAGE_SELECT_CELLS:
    # plot voltage for selected cells in 1D chain
    cells = SELECT_CELLS
    if any(cell < 1 or cell > numCells for cell in cells):
        raise ValueError(f'Cell numbers must be between 1 and {numCells}')
    vm_plot = np.zeros((len(t), len(cells)))
    for i, cell in enumerate(cells):
        vm_plot[:, i] = voltages[:, cell - 1]
    numPlots = len(cells)
    # mask = (t > 2000) & (t < 2500)
    mask = t > (t[-1] - CL)  # final cycle
    fig, axes = plt.subplots(numPlots, 1, num=2, squeeze=False)
    for i in range(numPlots):
        ax = axes[i, 0]
        ax.plot(t[mask], vm_plot[mask, i], linewidth=2)
        ax.set_title(f'Cell {cells[i]}')
        ax.set_ylabel('Vm (mV)')
        if i == 0:
            ax.set_xlabel('Time (ms)')
    fig.tight_layout()

if PLOT_PECG:
    # plot pseudo ECG
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
            Vm = voltages[i, j - 1]
            Vm_1 = voltages[i, j - 2]
            pECG[i] += ((a**2 * Ci) / (4 * Ce)) * (Vm_1 - Vm) * (1 / r) * dx
    mask = t > (t[-1] - CL)
    plt.figure(4)
    plt.plot(t[mask], pECG[mask], linewidth=2)
    plt.xlabel('Time (ms)')
    plt.ylabel('Pseudo ECG / Extracellular unipolar potentials (mV)')

# save figures
plot_dir = Path(filename).parent / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)

for number in plt.get_fignums():
    fig = plt.figure(number)
    output = plot_dir / f"OneD_{get_label(filename)}_figure_{number}.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print(f"Saved {output}")

if SAVE_VIDEO:
    def make_movie(output, frame_step, fps, writer):
        output = Path(output)
        frames = range(0, len(t), frame_step)
        total_frames = len(frames)
        report_every = max(1, total_frames // 100)
        print(f"Starting {output.name}: {total_frames} frames", flush=True)

        fig, ax = plt.subplots(figsize=(10, 2))
        strip = np.repeat(voltages[0:1], 3, axis=0)
        image = ax.imshow(strip, cmap='jet', vmin=vmin,
                          vmax=vmax, aspect='auto')
        fig.colorbar(image, ax=ax, orientation='horizontal')
        ax.get_yaxis().set_visible(False)
        title = ax.set_title('')

        def update(i):
            image.set_data(np.repeat(voltages[i:i+1], 3, axis=0))
            title.set_text(f'Voltage in Transmural Strand (mV), t={t[i] / 1000:.3f}s')
            return image, title

        def report_progress(current_frame, _total):
            completed = current_frame + 1
            if (completed == 1
                or completed % report_every == 0
                or completed == total_frames
            ):
                percent = 100 * completed / total_frames
                print(
                    f"\r{output.name}: "
                    f"{completed}/{total_frames} frames "
                    f"({percent:5.1f}%)",
                    end="",
                    flush=True,
                )
            if completed == total_frames:
                print("\nFinalizing file...", flush=True)

        movie = FuncAnimation(fig, update, frames=frames,
                              interval=1000 / fps, blit=False,
                              cache_frame_data=False)
        movie.save(str(output), writer=writer(fps=fps),
                   progress_callback=report_progress)
        print(f"Saved {output}", flush=True)
        plt.close(fig)

    if not FFMpegWriter.isAvailable():
        raise RuntimeError('MP4 export requires FFmpeg installed and available on PATH')
    make_movie(plot_dir / f"{get_label(filename)}_movie.mp4",
                frame_step=10,
                fps=200,
                writer=FFMpegWriter)

plt.close("all")
