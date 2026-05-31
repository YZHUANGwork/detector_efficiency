import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from scipy.interpolate import interp1d

# ── configuration ─────────────────────────────────────────────────────────────
EFFICIENCY_folder = './detector_efficiency'
FIGURES_folder    = './figures'
os.makedirs(FIGURES_folder, exist_ok=True)

detectors    = ['LZ', 'G3', 'XENONnT', 'Xe1t', 'Xe100t-5', 'Argon-Darkside']
interactions = ['beta', 'gammaRay', 'NR']
colors       = ['red', 'blue', 'green', 'purple']   # 4th color for Darkside measured
signal_type  = 'S1S2'

E_plot = np.logspace(-3, np.log10(1000), 1000)  # keV


def _read(filepath):
    """Read a two-column (E, efficiency) .txt file, skipping non-numeric header lines."""
    with open(filepath) as f:
        lines = f.readlines()
    skiprows = 0
    for line in lines:
        try:
            float(line.split()[0])
            break
        except (ValueError, IndexError):
            skiprows += 1
    data = np.loadtxt(filepath, skiprows=skiprows)
    return data[:, 0], data[:, 1]


# ── plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, len(detectors),
                         figsize=(4 * len(detectors), 3),
                         sharey=True, sharex=True)

for d, detector in enumerate(detectors):
    ax = axes[d]

    if detector == 'Argon-Darkside':
        # measured efficiency from data, not simulated
        file = os.path.join(EFFICIENCY_folder, f'det_eff_{detector}_NR.txt')
        if os.path.exists(file):
            E, eff = _read(file)
            f = interp1d(E, eff, fill_value=(0, eff[-1]), bounds_error=False)
            ax.semilogx(E_plot, f(E_plot), lw=2, label='NR (measured)', color=colors[3])
    else:
        for i, interaction in enumerate(interactions):
            file = os.path.join(EFFICIENCY_folder, f'det_eff_{detector}_{interaction}.csv')
            if not os.path.exists(file):
                continue

            df  = pd.read_csv(file, header=0, index_col=0)
            E   = df['E_center [keV]']
            eff = df[f'eff {signal_type}']

            f = interp1d(E, eff, fill_value=(0, eff.iloc[-1]), bounds_error=False)
            ax.semilogx(E_plot, f(E_plot), lw=2, label=interaction, color=colors[i])

    # formatting
    ax.set_xlim(0.1, 100)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel(r'$E_{r}$ [keV]', fontsize=13)
    ax.set_xticks([0.1, 1, 10, 100])
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%1g'))
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%1g'))
    ax.tick_params(labelsize=13)
    ax.grid(True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.text(0.1, 1.0, detector, transform=ax.transAxes,
            fontsize=15, va='top', ha='left')

axes[0].set_ylabel(r'Acceptance Rate $\epsilon(E_{r})$', fontsize=13)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, bbox_to_anchor=(0.4, 0.97, 0.5, .102),
           loc='lower left', ncol=3, mode='expand', fontsize=13)
fig.suptitle('S1S2 signals', fontsize=15, x=0.1, y=1.05, ha='left')
fig.tight_layout()

save_path = os.path.join(FIGURES_folder, 'detector_efficiency.png')
fig.savefig(save_path, bbox_inches='tight')
print(f'saved -> {save_path}')