#!/usr/bin/python
"""
Premixed jet impinging on a cold wall, the counterpart of Cantera's
ImpingingJet: plug-flow inlet at x = 0 with a prescribed mass flux, inert
no-slip wall at x = L. The same case is solved with Cantera, whose solution is
used as the reference in the comparison plot.

Cantera solves for the pressure curvature Lambda that gives the inlet mass
flux. In ember Lambda = -rhou*a^2/beta^2 is set by the strain rate a, which
the ImpingingJet(massFlux=...) controller adjusts until the inlet mass flux
is reached (the converged a is written to out.h5 together with mdot).
"""

from ember import *
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct

output = 'run/ex_impingingJet'
L = 0.01  # nozzle-to-wall distance [m]
Twall = 400.0  # [K]
uIn = 3.0  # inlet velocity [m/s]

# Cantera reference
gas = ct.Solution('h2o2.yaml', transport_model='mixture-averaged')
gas.set_equivalence_ratio(0.7, 'H2:1.0', 'O2:1.0, N2:3.76')
gas.TP = 300.0, ct.one_atm
rhou = gas.density
mdot = rhou * uIn

conf = Config(
    Paths(outputDir=output),
    Chemistry(mechanismFile='h2o2.yaml', transportModel='Mix'),
    General(flameGeometry='disc',
            nThreads=4),
    InitialCondition(fuel='H2:1.0',
                     oxidizer='O2:1.0, N2:3.76',
                     equivalenceRatio=0.7,
                     Tu=300.0,
                     xRight=L,
                     centerWidth=0.003,
                     slopeWidth=0.001),
    ImpingingJet(wallTemperature=Twall,
                 massFlux=mdot),
    # Starting guess for a; the controller takes it from here
    StrainParameters(initial=600.0,
                     final=600.0),
    Times(globalTimestep=1e-5),
    TerminationCondition(tEnd=0.05,
                         measurement=None),
)

if __name__ == '__main__':
    conf.run()

    jet = ct.ImpingingJet(gas=gas, width=L)
    jet.inlet.mdot = mdot
    jet.surface.T = Twall
    jet.set_initial_guess(products='equil')
    jet.set_refine_criteria(ratio=3, slope=0.05, curve=0.1)
    jet.solve(loglevel=0, auto=True)
    aCantera = 2.0 * np.sqrt(-jet.L[0] / rhou)  # beta = 2 (disc)

    struct = utils.load(output + '/profNow.h5')
    print('Inlet mass flux: Cantera {:.4f}, ember {:.4f} kg/m^2/s'.format(
        mdot, struct.V[0]))
    print('Equivalent strain rate: Cantera {:.1f}, ember {:.1f} 1/s'.format(
        aCantera, struct.a))

    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    ax[0].plot(jet.grid, jet.T, 'k--', label='Cantera')
    ax[0].plot(struct.x, struct.T, label='ember')
    ax[1].plot(jet.grid, jet.velocity, 'k--')
    ax[1].plot(struct.x, struct.V / struct.rho)
    ax[2].plot(jet.grid, jet.spread_rate, 'k--')
    ax[2].plot(struct.x, struct.U)
    ax[0].set_ylabel('Temperature [K]')
    ax[1].set_ylabel('Axial velocity [m/s]')
    ax[2].set_ylabel('Radial velocity / r [1/s]')
    for a in ax:
        a.set_xlabel('Position [m]')
    ax[0].legend()
    plt.tight_layout()
    plt.savefig(output + '/comparison.png')
    plt.close()
