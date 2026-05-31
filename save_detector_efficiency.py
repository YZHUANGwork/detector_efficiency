# Detector Efficiency from NEST Simulation
#
# For each simulated event, a recoil energy is drawn uniformly and passed to
# nestpy (Python bindings for NEST) with a chosen interaction type
# (NR, beta, or gammaRay). nestpy returns S1 and S2 signal values.
# A positive S1 or S2 means the signal was detected; a negative value means
# it was not. Detection efficiency is then computed by binning events in
# true recoil energy and counting the fraction with S1 > 0, S2 > 0, or both.

import numpy as np
import os
import pandas as pd

# ── configuration ────────────────────────────────────────────────────────────
binnum        = 100
RAWDATA_folder        = '../nest_sim_raw'       # folder with raw simulation .txt files
EFFICIENCY_folder = './detector_efficiency'  # output folder for .csv files

os.makedirs(EFFICIENCY_folder, exist_ok=True)

# ── helper functions ──────────────────────────────────────────────────────────
def read_file_dataline(file):
    """Read tab-separated simulation output, skipping header comment lines."""
    data = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('s') or '\t' not in line:
                continue
            vals = line.strip().split('\t')
            vals = [v for v in vals if v != '']
            try:
                data.append([float(v) for v in vals])
            except ValueError:
                continue
    data = np.array(data).T
    return data


def get_efficiency(eff_S1, eff_S2):
    """Combined S1+S2 efficiency (either signal counts)."""
    return eff_S1 + (1 - eff_S1) * eff_S2


def generate_eff(file_name, binnum, folder=RAWDATA_folder,
                 logbin=True, s1phe_boundary=[0, 0], s2phe_boundary=[0, 0]):
    """
    Read a raw simulation .txt file and compute detection efficiency per energy bin.

    Returns
    -------
    bins_E_centers : array  bin center energies [keV]
    bins_E         : array  bin edges [keV]
    eff            : array  combined S1+S2 efficiency
    eff_S1         : array  S1-only efficiency
    eff_S2         : array  S2-only efficiency
    eff_S1S2       : array  S1 AND S2 efficiency
    """
    file = os.path.join(folder, file_name + '.txt')
    data = read_file_dataline(file)

    dataf = pd.DataFrame({
        'E[keV]':    data[0],
        'cS1[phe]':  data[1],
        'cS2[phe]':  data[2],
        'Etrue[keV]': data[3]
    })

    counts, bins = np.histogram(dataf['Etrue[keV]'], binnum, density=True)

    if logbin:
        bins_E = np.logspace(
            np.log10(min(bins) - min(bins) / 10),
            np.log10(max(bins) + max(bins) / 100),
            binnum
        )
    else:
        bins_E = bins

    bins_E_centers = (bins_E[1:] + bins_E[:-1]) / 2
    eff_S1   = np.zeros(bins_E_centers.shape)
    eff_S2   = np.zeros(bins_E_centers.shape)
    eff_S1S2 = np.zeros(bins_E_centers.shape)

    for b in range(len(bins_E) - 1):
        good = (dataf['Etrue[keV]'] >= bins_E[b]) & (dataf['Etrue[keV]'] < bins_E[b + 1])
        n = good.sum()

        if n == 0:
            continue

        if sum(s1phe_boundary) > 0 and sum(s2phe_boundary) > 0:
            s1_lo, s1_hi = min(s1phe_boundary), max(s1phe_boundary)
            s2_lo, s2_hi = min(s2phe_boundary), max(s2phe_boundary)
            s1_vals = dataf.loc[good, 'cS1[phe]']
            s2_vals = dataf.loc[good, 'cS2[phe]']
            good_S1   = ((s1_vals > s1_lo) & (s1_vals < s1_hi)).sum()
            good_S2   = ((s2_vals > s2_lo) & (s2_vals < s2_hi)).sum()
            good_S1S2 = (((s1_vals > s1_lo) & (s1_vals < s1_hi)) &
                         ((s2_vals > s2_lo) & (s2_vals < s2_hi))).sum()
        else:
            s1_vals = dataf.loc[good, 'cS1[phe]']
            s2_vals = dataf.loc[good, 'cS2[phe]']
            good_S1   = (s1_vals > 0).sum()
            good_S2   = (s2_vals > 0).sum()
            good_S1S2 = ((s1_vals > 0) & (s2_vals > 0)).sum()

        eff_S1[b]   = good_S1   / n
        eff_S2[b]   = good_S2   / n
        eff_S1S2[b] = good_S1S2 / n

    eff_S1[np.isnan(eff_S1)]     = 0
    eff_S2[np.isnan(eff_S2)]     = 0
    eff_S1S2[np.isnan(eff_S1S2)] = 0

    return bins_E_centers, bins_E, get_efficiency(eff_S1, eff_S2), eff_S1, eff_S2, eff_S1S2


# ── main: loop over detectors and interactions ────────────────────────────────
detectors    = ['LZ', 'G3', 'XENONnT', 'Xe1t', 'Xe100t-5']#
interactions = ['beta', 'gammaRay', 'NR']

for detector in detectors:
    for interaction in interactions:
        FILE_NAME = f'{detector}_S1S2_{interaction}_Xenon_1e+06_0-100keV'
        print(f'Processing {FILE_NAME} ...')

        bins_E_centers, bins_E, eff, eff_S1, eff_S2, eff_S1S2 = generate_eff(
            FILE_NAME, binnum, folder=RAWDATA_folder, logbin=True,
            s1phe_boundary=[0, 0], s2phe_boundary=[0, 0]
        )

        df_eff = pd.DataFrame({
            'E_center [keV]': bins_E_centers,
            'E [keV]':        [f'{bins_E[i]}-{bins_E[i+1]}' for i in range(len(bins_E) - 1)],
            'trigger eff':    eff,
            'eff S1':         eff_S1,
            'eff S2':         eff_S2,
            'eff S1S2':       eff_S1S2,
        })

        out_file = os.path.join(EFFICIENCY_folder, f'det_eff_{detector}_{interaction}.csv')
        df_eff.to_csv(out_file)
        print(f'  saved -> {out_file}')
