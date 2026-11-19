#!/usr/bin/python
"""
Plots for the particle-population check on candidate (a) of the jet grid
(single_Igor_impingingJet_fixedT_particles_caseA.py): early particles alone,
with the late (bi-modal Fe2O3) population, and with an empty late population
(control), plus the no-particle run for the gas Fe reference.

    python plot_caseA_populations.py

Writes run/single_Igor/caseA_early.png (early particles) and
run/single_Igor/caseA_late.png (late particles) and prints the metrics.
"""
import glob
import h5py
import numpy as np
import cantera as ct
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RUN = 'run/single_Igor/caseA_'
gas = ct.Solution('opt-compact-mech-OF2_OFkinetics_FeTransp.yaml')
W = gas.molecular_weights
WFe = gas.atomic_weight('Fe')
feFrac = np.array([gas.n_atoms(k, 'Fe') for k in range(gas.n_species)]) * WFe / W
kFe2O3 = gas.species_index('FE2O3(s)')
NA = ct.avogadro

# reference palette (categorical slots 1-3 in fixed order, neutral for the reference)
INK, INK2, GRID = '#1f1f1e', '#5f5e58', '#e1e0d9'
STYLE = {
    'early': dict(color='#2a78d6', ls='-', lw=2.0, label='early only (case a)'),
    'late': dict(color='#eb6834', ls='-', lw=2.0, label='early + late particles'),
    'lateoff': dict(color='#1baf7a', ls=(0, (4, 2)), lw=2.0,
                    label='early + empty late population (control)'),
    'none': dict(color='#8a8980', ls=(0, (1, 1.5)), lw=2.0, label='no particles'),
}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
                     'axes.spines.top': False, 'axes.spines.right': False})


def load(mode):
    files = glob.glob(RUN + mode + '/prof*.h5')
    f = max((h5py.File(fn, 'r') for fn in files), key=lambda f: float(np.array(f['t'])))
    d = {k: np.array(f[k]) for k in f}
    d['feGas'] = (feFrac[:, None] * d['Y']).sum(axis=0)  # Fe in the gas [kg/kg]
    if 'moments_early' in d:
        N, mM, mO, H = d['moments_early']
        d['massEarly'] = mM + mO
    return d


runs = {m: load(m) for m in ('early', 'late', 'lateoff', 'none')}
ref = runs['none']


def extinction(d):
    m = d['massEarly']; j = np.argmax(m)
    return d['x'][j + np.nonzero(m[j:] < 0.01 * m[j])[0][0]]


def metrics(d):
    x = d['x']; xEnd = extinction(d)
    sel = (x >= 2e-3) & (x <= 8e-3)
    F = d['feGas'] / np.interp(x, ref['x'], ref['feGas'])
    up = x <= xEnd
    # F mean as a plain average over the grid points (as in the doc) and
    # weighted by distance (independent of the adaptive grid)
    return (np.median((d['particleTemperature_early'] - d['T'])[sel]), 1e3 * xEnd,
            F[up].min(), F[up].mean(),
            np.trapezoid(F[up], x[up]) / (x[up][-1] - x[up][0]))


rows = [('early', 'early only'), ('lateoff', 'empty late (control)'), ('late', 'early + late')]
print('%-22s %7s  %13s  %10s  %11s  %s' % ('', 't [s]', 'Tp-T 2-8 mm', 'gone at',
                                           'F min/mean', 'F mean x-weighted'))
print('%-22s %7s  %13s  %10s  %11s' % ('doc, case (a)', '', '607 K', '8.6 mm', '0.28 / 0.62'))
table = []
for m, name in rows:
    table.append((name, float(runs[m]['t'])) + metrics(runs[m]))
    print('%-22s %7.4f  %11.1f K  %7.2f mm  %.3f / %.3f  %.3f' % table[-1])


def style(ax, xlabel=True):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel('distance from burner x [mm]')


# ---------------------------------------------------------------- early particles
fig, axs = plt.subplots(3, 3, figsize=(12, 8.0), constrained_layout=True,
                        height_ratios=[1, 1, 0.32])
gsTable = axs[2, 0].get_gridspec()
for ax in axs[2]:
    ax.remove()
axTable = fig.add_subplot(gsTable[2, :])
axTable.axis('off')
axs = axs[:2]
xmax = 11.0
for m in ('early', 'lateoff', 'late'):
    d = runs[m]; x = 1e3 * d['x']; s = STYLE[m]
    live = d['massEarly'] > 1e-3 * d['massEarly'].max()  # T_p, d_p only where particles exist
    axs[0, 0].plot(x[live], (d['particleTemperature_early'] - d['T'])[live], **s)
    axs[0, 1].plot(x, 1e3 * d['massEarly'], **s)
    axs[0, 2].plot(x, d['numberDensity_early'], **s)
    axs[1, 0].plot(x[live], 1e9 * d['particleDiameter_early'][live], **s)
    axs[1, 1].plot(x, d['feGas'] / np.interp(d['x'], ref['x'], ref['feGas']), **s)
    if m != 'early':
        e = runs['early']
        diff = (np.interp(e['x'], d['x'], d['massEarly']) - e['massEarly']) / e['massEarly'].max()
        axs[1, 2].plot(1e3 * e['x'], 100 * diff, **s)
axs[0, 0].set_ylabel('particle overheating $T_p - T$ [K]')
axs[0, 1].set_ylabel('early particle mass $m_{Fe}+m_O$ [g / kg gas]')
axs[0, 2].set_ylabel('early number density [1/m³]')
axs[0, 2].set_yscale('log'); axs[0, 2].set_ylim(1e14, None)
axs[1, 0].set_ylabel('early particle diameter $d_p$ [nm]')
axs[1, 1].set_ylabel('gas-phase Fe / no-particle run')
axs[1, 1].axhline(1.0, color=INK2, lw=0.8)
axs[1, 2].set_ylabel('early mass − early-only run [% of peak]')
axs[1, 2].axhline(0.0, color=INK2, lw=0.8)
for ax, letter in zip(axs.flat, 'abcdef'):
    style(ax); ax.set_xlim(0, xmax)
    ax.set_title('(%s)' % letter, loc='left', fontsize=9, color=INK2)
for ax in (axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]):
    xe = 1e3 * extinction(runs['early'])
    ax.axvline(xe, color=STYLE['early']['color'], lw=0.8, ls=':')
cells = [['607 K', '8.6 mm', '0.28', '0.62', '–', '–']]
cells += [['%.1f K' % dT, '%.2f mm' % xe, '%.3f' % a, '%.3f' % b, '%.3f' % c, '%.4f s' % t]
          for n, t, dT, xe, a, b, c in table]
tab = axTable.table(cellText=cells,
                    rowLabels=['doc, case (a)'] + [n for n, *_ in table],
                    colLabels=['$T_p - T$, median 2–8 mm', 'particles gone at',
                               'gas Fe ratio F, min', 'F mean (grid points)',
                               'F mean (distance-weighted)', 'end time'],
                    loc='center', cellLoc='center', edges='horizontal')
tab.auto_set_font_size(False)
tab.set_fontsize(8.5)
tab.scale(1, 1.35)
for cell in tab.get_celld().values():
    cell.set_edgecolor(GRID)
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='outside lower center', ncol=3, frameon=False)
fig.suptitle('Early particles, candidate (a): unchanged by the late population except near '
             'the evaporation front', fontsize=11, color=INK, x=0.01, ha='left')
fig.savefig(RUN + 'early.png', dpi=150)

# ---------------------------------------------------------------- late particles
d, e = runs['late'], runs['early']
x = 1e3 * d['x']
fig, axs = plt.subplots(2, 3, figsize=(12, 6.6), constrained_layout=True)
LATE = '#eb6834'
seq = ['#86b6ef', '#2a78d6', '#184f95']  # one-hue ramp for the three diameters
axs[0, 0].plot(x, d['numberDensity1_late'], color=seq[0], lw=2, label='mode 1 (monomers, 0.46 nm)')
axs[0, 0].plot(x, d['numberDensity2_late'], color=seq[2], lw=2, label='mode 2 (aggregates)')
axs[0, 0].set_yscale('log'); axs[0, 0].set_ylim(1e12, None)
axs[0, 0].set_ylabel('late number density [1/m³]')
axs[0, 0].legend(frameon=False)
# mode-2 geometry only where there are aggregates (> 1e13 per m^3): below that
# it is integrator noise (N2, M2 near zero at the burner)
mode2 = d['numberDensity2_late'] > 1e13
for key, c, lab in (('massMeanDiameter_late', seq[0], 'mass-mean (both modes)'),
                    ('volumeDiameter2_late', seq[1], 'mode-2 volume-equivalent'),
                    ('collisionDiameter_late', seq[2], 'mode-2 collision (fractal)')):
    v = d[key] if key == 'massMeanDiameter_late' else np.where(mode2, d[key], np.nan)
    axs[0, 1].plot(x, 1e9 * v, color=c, lw=2, label=lab)
axs[0, 1].set_yscale('log'); axs[0, 1].set_ylim(0.3, None)
axs[0, 1].set_ylabel('late particle diameter [nm]')
axs[0, 1].legend(frameon=False)
axs[0, 2].plot(x, np.where(mode2, d['primaryCount_late'], np.nan), color=LATE, lw=2)
axs[0, 2].set_yscale('log'); axs[0, 2].set_ylim(1, None)
axs[0, 2].set_ylabel('primaries per aggregate $n_p$ (no sintering)')
for m in ('early', 'late', 'none'):
    r = runs[m]
    axs[1, 0].plot(1e3 * r['x'], r['Y'][kFe2O3], **STYLE[m])
axs[1, 0].set_yscale('log'); axs[1, 0].set_ylim(1e-10, None)
axs[1, 0].set_ylabel('gas-phase Fe$_2$O$_3$ mass fraction')
axs[1, 0].legend(frameon=False, loc='lower right')
N1, N2, M2, A2 = d['moments_late']
m1 = W[kFe2O3] / NA
feLate = (N1 * NA * m1 + M2) * 2 * WFe / W[kFe2O3]
axs[1, 1].stackplot(x, 1e3 * d['feGas'], 1e3 * d['moments_early'][1], 1e3 * feLate,
                    colors=['#8a8980', '#2a78d6', LATE], edgecolor='white', linewidth=0.5,
                    labels=['gas', 'early particles', 'late particles'])
axs[1, 1].axhline(1e3 * ref['feGas'][0], color=INK, lw=1, ls='--',
                  label='gas Fe at x = 0, no particles')
axs[1, 1].set_ylabel('Fe [g / kg gas]')
axs[1, 1].legend(frameon=False, loc='center left')
axs[1, 2].plot(x, d['T'], color=INK2, lw=2)
axs[1, 2].set_ylabel('gas temperature T [K] (prescribed)')
for ax, letter in zip(axs.flat, 'abcdef'):
    style(ax); ax.set_xlim(0, x[-1])
    ax.set_title('(%s)' % letter, loc='left', fontsize=9, color=INK2)
axs[1, 1].annotate('late mass piles up at the wall\n(no deposition, $D_p$ ~ 1e-9 m²/s)',
                   xy=(0.97, 0.97), xycoords='axes fraction', ha='right', va='top', color=INK2)
fig.suptitle('Late particles (bi-modal, Fe$_2$O$_3$ monomer conversion, Df = 1.8, no '
             'sintering), stopped at t = %.2f s' % float(d['t']),
             fontsize=11, color=INK, x=0.01, ha='left')
fig.savefig(RUN + 'late.png', dpi=150)
print('wrote', RUN + 'early.png', RUN + 'late.png')
