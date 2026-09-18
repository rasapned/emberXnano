#!/usr/bin/python
"""
Diagnostic plots for burnerFla_test.py: load a saved state (profNow.h5,
preAdapt.h5, postAdapt.h5, a numbered profXXXXXX.h5, or one of the
err_*.h5 dumps written when a split-operator stage throws) and plot the
fields needed to see *where* something has gone wrong and *which*
split-operator term (convection/diffusion/production/cross) is responsible.

Usage:
    python plot_burner_diagnostics.py [output_dir] [filename]

Defaults to run/burnerFla_test/profNow.h5 -- pass a different filename
(e.g. err_convectionIntegration.h5, preAdapt.h5, prof000003.h5) to inspect
a different snapshot. All plots are saved as PNGs into output_dir.
"""
import sys
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct

from ember import utils

# Species order must match whatever mechanism generated the data (only used
# for labeling/selecting a handful of species to plot -- doesn't need to be
# exact if it's slightly off, just point mechanism at the right .yaml).
mechanism = 'opt-compact-mech-OF2.yaml'

# Species singled out because they showed up going unphysical in earlier
# runs. Edit this list once you know what to look for in your own case.
species_watch = ['FEC5O5', 'FE', 'H2', 'O2', 'H2O', 'AR', 'OH', 'H']


def plot_diagnostics(output_dir, filename='profNow.h5'):
    path = os.path.join(output_dir, filename)
    struct = utils.load(path)
    x = struct.x
    dx = np.diff(x)

    gas = ct.Solution(mechanism)
    names = gas.species_names
    idx = {name: names.index(name) for name in species_watch if name in names}

    print(f"Loaded {path}")
    print(f"  {len(x)} grid points, x in [{x[0]:.5g}, {x[-1]:.5g}] m")
    print(f"  grid spacing dx in [{dx.min():.3g}, {dx.max():.3g}] m "
          f"(ratio {dx.max()/dx.min():.1f})")
    jTmin = int(np.argmin(struct.T))
    print(f"  T in [{struct.T.min():.1f}, {struct.T.max():.1f}] K "
          f"(min at j={jTmin}, x={x[jTmin]:.5g})")
    for name, k in idx.items():
        Yk = struct.Y[k, :]
        jmin, jmax = int(np.argmin(Yk)), int(np.argmax(Yk))
        print(f"  Y[{name}] in [{Yk.min():.4g}, {Yk.max():.4g}] "
              f"(min at j={jmin} x={x[jmin]:.5g}, "
              f"max at j={jmax} x={x[jmax]:.5g})")

    base = os.path.splitext(filename)[0]

    # --- Temperature and grid spacing --------------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
    ax1.plot(x, struct.T, '.-', ms=3)
    ax1.set_ylabel('T [K]')
    ax1.axhline(struct.T[0], color='gray', ls=':', lw=0.8,
                label=f'T[0] = {struct.T[0]:.1f} K (inlet)')
    ax1.legend(loc='best', fontsize=8)
    ax2.semilogy(0.5*(x[1:]+x[:-1]), dx, '.-', ms=3, color='C1')
    ax2.set_ylabel('grid spacing dx [m]')
    ax2.set_xlabel('Position [m]')
    fig.suptitle(f'{filename}: temperature & grid spacing')
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f'{base}_T_grid.png'))
    plt.close(fig)

    # --- Watched species mass fractions -------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    for name, k in idx.items():
        ax.plot(x, struct.Y[k, :], '.-', ms=2, label=name)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlabel('Position [m]')
    ax.set_ylabel('Mass fraction')
    ax.legend(loc='best', fontsize=8, ncol=2)
    fig.suptitle(f'{filename}: watched species')
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f'{base}_species.png'))
    plt.close(fig)

    # --- Split-operator breakdown for T --------------------------------
    if 'dTdtConv' in struct:
        fig, ax = plt.subplots(figsize=(8, 5))
        for key, label in [('dTdtConv', 'convection'), ('dTdtDiff', 'diffusion'),
                            ('dTdtProd', 'production'), ('dTdtCross', 'cross')]:
            if key in struct:
                ax.plot(x, struct[key], '.-', ms=2, label=label)
        ax.axhline(0, color='k', lw=0.8)
        ax.set_xlabel('Position [m]')
        ax.set_ylabel('dT/dt contribution [K/s]')
        ax.legend(loc='best', fontsize=8)
        fig.suptitle(f'{filename}: split-operator dT/dt breakdown')
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, f'{base}_dTdt.png'))
        plt.close(fig)

    # --- Split-operator breakdown for each watched species -------------
    if 'dYdtConv' in struct:
        for name, k in idx.items():
            fig, ax = plt.subplots(figsize=(8, 5))
            for key, label in [('dYdtConv', 'convection'), ('dYdtDiff', 'diffusion'),
                                ('dYdtProd', 'production'), ('dYdtCross', 'cross')]:
                if key in struct:
                    ax.plot(x, struct[key][k, :], '.-', ms=2, label=label)
            ax.axhline(0, color='k', lw=0.8)
            ax.set_xlabel('Position [m]')
            ax.set_ylabel(f'dY[{name}]/dt contribution [1/s]')
            ax.legend(loc='best', fontsize=8)
            fig.suptitle(f'{filename}: split-operator dY[{name}]/dt breakdown')
            fig.tight_layout()
            fig.savefig(os.path.join(output_dir, f'{base}_dYdt_{name}.png'))
            plt.close(fig)

    print(f"Plots written to {output_dir}/{base}_*.png")


if __name__ == '__main__':
    output_dir = sys.argv[1] if len(sys.argv) > 1 else 'run/burnerFla_test'
    filename = sys.argv[2] if len(sys.argv) > 2 else 'profNow.h5'
    plot_diagnostics(output_dir, filename)
