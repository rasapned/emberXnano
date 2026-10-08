#!/usr/bin/python
"""
Igor's burner (lean, 500 ppm FEC5O5, 3000 Pa) as a premixed jet impinging on
the stagnation plate, instead of single_Igor.py's flame opposing a cold inert
stream. Burner data as in drobnostki/mechCompare/mechCompare.py (case L500):
36 mm nozzle, flow rates in sccm, nozzle-to-plate distance and the burner /
plate temperatures from the measured profile T_of_x-L500.csv.

The inlet mass flux is prescribed; ember's ImpingingJet controller adjusts the
strain rate (pressure curvature) to reach it. The same case is solved with
Cantera's ImpingingJet for comparison.
"""

from ember import *
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct

output = 'run/single_Igor_impingingJet'
mech = 'Iron_elte_Syngas-newTransp.yaml'

# Lean case, 500 ppm FEC5O5, volume flow rates in sccm
X_FEC5O5 = 0.0005
VS_total = (400 + 400 + 600) / (1 - X_FEC5O5)
X_H2 = 400 / VS_total
X_O2 = 400 / VS_total
X_AR = 1.0 - (X_H2 + X_O2 + X_FEC5O5)
react = f"H2:{X_H2:.8e}, AR:{X_AR:.8e}, FEC5O5:{X_FEC5O5:.8e}, O2:{X_O2:.8e}"

pressure = 3000.0  # [Pa]
L = 0.0444  # nozzle-to-plate distance [m]
Tburner = 426.8  # [K]
Tplate = 1146.2  # [K]

# Inlet mass flux: standard (1 bar, 273.15 K) volume flow through the nozzle
gas = ct.Solution(mech, transport_model='mixture-averaged')
gas.TPX = 273.15, 1e5, react
mdot = gas.density * VS_total / 60e6 / (0.25 * np.pi * 0.036**2)
gas.TPX = Tburner, pressure, react
rhou = gas.density

conf = Config(
    Paths(outputDir=output),
    Chemistry(mechanismFile=mech,
              transportModel='Mix'),
    General(flameGeometry='disc',
            nThreads=16,
            chemistryIntegrator='cvode',
            splittingMethod='balanced'),
    InitialCondition(reactants=react,
                     Tu=Tburner,
                     pressure=pressure,
                     xRight=L,
                     # flame front at ~6 mm (from the Cantera solution)
                     centerWidth=0.032,
                     slopeWidth=0.003),
    ImpingingJet(wallTemperature=Tplate,
                 massFlux=mdot),
    # Starting guess for a; the controller takes it from here
    StrainParameters(initial=100.0,
                     final=100.0),
    CvodeTolerances(relativeTolerance=1e-7,
                    speciesAbsTol=1e-12,
                    energyAbsTol=1e-8,
                    minimumTimestep=1e-16),
    # Stop once the heat release rate is steady; a period of ~10 ms covers
    # the slow final approach of the mass flux controller
    TerminationCondition(tEnd=0.3,
                         tMin=0.02,
                         steadyPeriod=0.01),
)

if __name__ == '__main__':
    conf.run()

    jet = ct.ImpingingJet(gas=gas, width=L)
    jet.inlet.mdot = mdot
    jet.surface.T = Tplate
    jet.set_initial_guess(products='equil')
    jet.set_refine_criteria(ratio=3, slope=0.05, curve=0.1, prune=0.02)
    jet.solve(loglevel=0, auto=True)
    aCantera = 2.0 * np.sqrt(-jet.L[0] / rhou)  # beta = 2 (disc)

    struct = utils.load(output + '/profNow.h5')
    print('Inlet mass flux: Cantera {:.5f}, ember {:.5f} kg/m^2/s'.format(
        mdot, struct.V[0]))
    print('Equivalent strain rate: Cantera {:.2f}, ember {:.2f} 1/s'.format(
        aCantera, struct.a))
    print('Maximum temperature: Cantera {:.1f}, ember {:.1f} K'.format(
        jet.T.max(), struct.T.max()))

    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    ax = ax.ravel()
    ax[0].plot(jet.grid, jet.T, 'k--', label='Cantera')
    ax[0].plot(struct.x, struct.T, label='ember')
    ax[0].set_ylabel('Temperature [K]')
    ax[1].plot(jet.grid, jet.velocity, 'k--')
    ax[1].plot(struct.x, struct.V / struct.rho)
    ax[1].set_ylabel('Axial velocity [m/s]')
    ax[2].plot(jet.grid, jet.spread_rate, 'k--')
    ax[2].plot(struct.x, struct.U)
    ax[2].set_ylabel('Radial velocity / r [1/s]')
    for i, species in enumerate(['H2O', 'FEC5O5', 'FE']):
        k = gas.species_index(species)
        ax[3+i].plot(jet.grid, jet.Y[k], 'k--')
        ax[3+i].plot(struct.x, struct.Y[k])
        ax[3+i].set_ylabel('Y ' + species)
    for a in ax:
        a.set_xlabel('Position [m]')
    ax[0].legend()
    plt.suptitle(f'Impinging jet, p={pressure:.0f} Pa, mdot={mdot:.4f} kg/m²s')
    plt.tight_layout()
    plt.savefig(output + '/comparison.png')
    plt.close()
