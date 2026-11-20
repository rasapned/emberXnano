#!/usr/bin/python
"""
Comparison of Sellmann2019_DFF_flameA.py with Sellmann et al., Proc.
Combust. Inst. 37 (2019) 1241, flame A, down-firing flame (digitised data
in Sellmann2019/): gas temperature (input), axial velocity and residence
time (Fig. 4), mass-averaged particle diameter (Fig. 7: unperturbed
simulation and PMS) against distance and against residence time, and the
late particle number densities.

    python plot_Sellmann2019_DFF_flameA.py

Writes run/Sellmann2019/DFF_flameA_fixedU.png and prints the metrics.
"""
import glob
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RUN = 'run/Sellmann2019/DFF_flameA_fixedU'
REF = 'Sellmann2019/'

INK, INK2, GRID = '#1f1f1e', '#5f5e58', '#e1e0d9'
EMBER, PAPER, EXP = '#2a78d6', '#1f1f1e', '#eb6834'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
                     'axes.spines.top': False, 'axes.spines.right': False})


def csv(name):
    return np.genfromtxt(REF + name, delimiter=',', comments='#').T


files = glob.glob(RUN + '/prof*.h5')
f = max((h5py.File(fn, 'r') for fn in files), key=lambda f: float(np.array(f['t'])))
d = {k: np.array(f[k]) for k in f}
x = d['x']
u = d['V'] / d['rho']

xT, TT = csv('T_of_x-DFF_flameA.csv')
xU, UU = csv('Ux_DFF_flameA.csv')
xD, DD = csv('dav_sim_DFF_flameA.csv')
xP, DP, DPlo, DPhi = csv('dav_PMS_DFF_flameA.csv')
# profiles fed to ember (written by the case script)
xIn, _, _, _, uIn = np.genfromtxt(RUN + '/inputProfiles.csv', delimiter=',', comments='#').T


def residenceTime(x, u, xMax):
    """t(x) = int_0^x dx/u, up to xMax (u -> 0 at the stagnation point)."""
    sel = (x <= xMax) & (u > 0)
    xs, us = x[sel], u[sel]
    return xs, np.r_[0.0, np.cumsum(0.5 * (1 / us[1:] + 1 / us[:-1]) * np.diff(xs))]


xMax = 44e-3  # paper U_x = 0 at 44.6 mm, ember wall at 46 mm
xtE, tE = residenceTime(x, u, xMax)
xtP, tP = residenceTime(xU, UU, xMax)
dE = d['massMeanDiameter_late']

# metrics
print('t = %.4f s, %d points, a = %.1f 1/s, inlet mass flux %.5f kg/m2/s, inlet u %.3f m/s'
      % (float(d['t']), len(x), float(d['a']), d['V'][0], u[0]))
print('%10s %10s %10s %10s %10s' % ('x [mm]', 't_res ember', 't_res paper', 'd_av ember', 'd_av paper'))
for xx in (20, 25, 30, 35, 40, 42, 43.4):
    print('%10.1f %9.1f ms %9.1f ms %8.2f nm %8.2f nm' % (
        xx, 1e3 * np.interp(xx * 1e-3, xtE, tE), 1e3 * np.interp(xx * 1e-3, xtP, tP),
        1e9 * np.interp(xx * 1e-3, x, dE), np.interp(xx * 1e-3, xD, DD)))
up = x <= 44e-3
j = np.argmax(dE)
print('d_av peak: ember %.2f nm at %.2f mm (up to 44 mm: %.2f nm); paper %.2f nm at %.1f mm'
      % (1e9 * dE[j], 1e3 * x[j], 1e9 * dE[up].max(), DD.max(), 1e3 * xD[np.argmax(DD)]))


def style(ax, xlabel='distance from burner x [mm]'):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel(xlabel)


fig, axd = plt.subplot_mosaic([['a', 'b', 'c'], ['d', 'e', 'f']],
                               figsize=(12, 7), constrained_layout=True)
ax = axd['a']
ax.plot(1e3 * xT, TT, color=PAPER, lw=2.5, ls=(0, (1, 1.5)), label='input (paper Figs. 4, 9)')
ax.plot(1e3 * x, d['T'], color=EMBER, lw=1.5, label='ember')
ax.set_ylabel('gas temperature T [K]')
ax.legend(frameon=False, loc='lower center')
ax = axd['b']
ax.plot(1e3 * xU, UU, color=PAPER, marker='v', ms=4, lw=1, ls='none', mfc='none', label='paper (Fig. 4)')
ax.plot(1e3 * xIn, uIn, color=PAPER, lw=2.5, ls=(0, (1, 1.5)), label='from input U, $\\rho$')
ax.plot(1e3 * x, u, color=EMBER, lw=1.5, label='ember (wall at 46 mm)')
ax.axhline(0, color=INK2, lw=0.8)
ax.set_ylabel('axial velocity $U_x$ [m/s]')
ax.legend(frameon=False)
ax = axd['c']
ax.plot(1e3 * xtP, 1e3 * tP, color=PAPER, lw=2, ls=(0, (4, 2)), label='paper (Fig. 4)')
ax.plot(1e3 * xtE, 1e3 * tE, color=EMBER, lw=2, label='ember')
ax.set_ylabel('residence time [ms]')
ax.set_ylim(0, 100)
ax.legend(frameon=False)
ax = axd['d']
ax.plot(1e3 * xD, DD, color=PAPER, lw=2, ls=(0, (4, 2)), label='paper sim. (Fig. 7)')
ax.errorbar(1e3 * xP, DP, yerr=[DP - DPlo, DPhi - DP], color=EXP, fmt='o', mfc='none', capsize=3,
            label='PMS (Fig. 7)')
ax.plot(1e3 * x, 1e9 * dE, color=EMBER, lw=2, label='ember, mass-mean')
ax.set_ylabel('mass-averaged particle diameter $d_{av}$ [nm]')
ax.set_xlim(0, 62)
ax.legend(frameon=False, loc='upper left')
ax = axd['e']
ax.plot(1e3 * tP, np.interp(xtP, xD, DD), color=PAPER, lw=2, ls=(0, (4, 2)), label='paper, unperturbed sim.')
ax.plot(1e3 * tE, 1e9 * np.interp(xtE, x, dE), color=EMBER, lw=2, label='ember')
ax.set_ylabel('$d_{av}$ [nm]')
ax.set_xlim(0, 150)
ax.legend(frameon=False, loc='upper left')
style(ax, 'residence time [ms] (up to x = 44 mm)')
ax = axd['f']
ax.plot(1e3 * x, d['numberDensity1_late'], color='#86b6ef', lw=2, label='mode 1 (monomers)')
ax.plot(1e3 * x, d['numberDensity2_late'], color='#184f95', lw=2, label='mode 2 (aggregates)')
ax.set_yscale('log')
ax.set_ylim(1e12, None)
ax.set_ylabel('number density [1/m³]')
ax2 = ax.twinx()
ax2.plot(1e3 * x, np.where(d['numberDensity2_late'] > 1e13, d['primaryCount_late'], np.nan),
         color=EXP, lw=1.5, ls=(0, (4, 2)), label='primaries per aggregate')
ax2.set_ylabel('primaries per aggregate $n_p$ (sintering)', color=EXP)
ax2.spines['right'].set_visible(True)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, frameon=False, loc='lower center')
for letter, a in axd.items():
    if letter != 'e':
        style(a)
    a.set_title('(%s)' % letter, loc='left', fontsize=9, color=INK2)
for letter in 'abcf':
    axd[letter].set_xlim(0, 1e3 * x[-1])
fig.suptitle('Sellmann et al. (2019), flame A, down-firing: Fe$_2$O$_3$ particles (bi-modal, sintering '
             '$A$ = 8e16 s/(K m$^4$), $T_a$ = 30000 K), ember at t = %.3f s' % float(d['t']),
             fontsize=11, color=INK, x=0.01, ha='left')
fig.savefig(RUN + '.png', dpi=150)
print('wrote', RUN + '.png')
