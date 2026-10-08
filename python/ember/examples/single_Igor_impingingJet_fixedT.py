#!/usr/bin/python
"""
Igor's burner (lean, 500 ppm FEC5O5, 3000 Pa) as an impinging jet, solved
first with Cantera's ImpingingJet and then with ember started from it.

Cantera's converged pressure curvature Lambda gives the equivalent ember
strain rate, a = beta*sqrt(-Lambda/rhou), and its profiles are the initial
condition for ember, so ember starts close to its own steady state. The mass
flux controller then only has to correct the small remaining difference.

With fixedTemperature = True the energy equation is not solved in either code:
T is the measured profile T_of_x-L500.csv (burner and plate temperatures are
its end values), as in drobnostki/mechCompare/mechCompare.py.
"""

from ember import *
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct

fixedTemperature = True

output = 'run/single_Igor_impingingJet_' + ('fixedT' if fixedTemperature else 'energy')
mech = 'Iron_elte_Syngas-newTransp.yaml'

# Lean case, 500 ppm FEC5O5, volume flow rates in sccm
X_FEC5O5 = 0.0005
VS_total = (400 + 400 + 600) / (1 - X_FEC5O5)
X_H2 = 400 / VS_total
X_O2 = 400 / VS_total
X_AR = 1.0 - (X_H2 + X_O2 + X_FEC5O5)
react = f"H2:{X_H2:.8e}, AR:{X_AR:.8e}, FEC5O5:{X_FEC5O5:.8e}, O2:{X_O2:.8e}"

pressure = 3000.0  # [Pa]

# Measured temperature profile: nozzle at x = 0, plate at x = L
xProfile, TProfile = np.genfromtxt('T_of_x-L500.csv', delimiter=',', comments='#').T
L = xProfile[-1]
Tburner = TProfile[0]
Tplate = TProfile[-1]

# Inlet mass flux: standard (1 bar, 273.15 K) volume flow through the nozzle
gas = ct.Solution(mech, transport_model='mixture-averaged')
gas.TPX = 273.15, 1e5, react
mdot = gas.density * VS_total / 60e6 / (0.25 * np.pi * 0.036**2)
gas.TPX = Tburner, pressure, react
rhou = gas.density

# *** Cantera ***
jet = ct.ImpingingJet(gas=gas, width=L)
jet.inlet.mdot = mdot
jet.surface.T = Tplate
jet.set_refine_criteria(ratio=3, slope=0.05, curve=0.1, prune=0.02)
if fixedTemperature:
    jet.flame.set_fixed_temp_profile(xProfile / L, TProfile)
    jet.energy_enabled = False
    jet.solve(loglevel=0, refine_grid=True)
else:
    jet.set_initial_guess(products='equil')
    jet.solve(loglevel=0, auto=True)
aCantera = 2.0 * np.sqrt(-jet.L[0] / rhou)  # beta = 2 (disc)
print('Cantera: {} points, equivalent strain rate a = {:.3f} 1/s'.format(
    len(jet.grid), aCantera))

# *** ember, started from the Cantera solution ***
# ember takes the inlet stream composition from the first point of the given
# profile, so that point must hold the supplied mixture, not Cantera's inlet
# node (depleted in H2, which diffuses upstream against the inlet flux)
Y0 = jet.Y.copy()
gas.TPX = Tburner, pressure, react
Y0[:, 0] = gas.Y

conf = Config(
    Paths(outputDir=output),
    Chemistry(mechanismFile=mech,
              transportModel='Mix'),
    General(flameGeometry='disc',
            nThreads=16,
            chemistryIntegrator='cvode'),
    InitialCondition(reactants=react,
                     Tu=Tburner,
                     pressure=pressure,
                     xRight=L,
                     haveProfiles=True,
                     x=jet.grid,
                     T=jet.T,
                     U=jet.spread_rate,
                     V=jet.density * jet.velocity,
                     Y=Y0),
    ImpingingJet(wallTemperature=Tplate,
                 massFlux=mdot,
                 temperatureProfile=((xProfile, TProfile) if fixedTemperature
                                     else None)),
    StrainParameters(initial=aCantera,
                     final=aCantera),
    CvodeTolerances(relativeTolerance=1e-7,
                    speciesAbsTol=1e-12,
                    energyAbsTol=1e-8,
                    minimumTimestep=1e-16),
    # Stop once the heat release rate is steady over ~10 ms
    TerminationCondition(tEnd=0.3,
                         tMin=0.01,
                         steadyPeriod=0.01),
)

if __name__ == '__main__':
    conf.run()

    struct = utils.load(output + '/profNow.h5')
    print('Stopped at t = {:.4f} s'.format(struct.t))
    print('Inlet mass flux: Cantera {:.5f}, ember {:.5f} kg/m^2/s'.format(
        mdot, struct.V[0]))
    print('Equivalent strain rate: Cantera {:.3f}, ember {:.3f} 1/s'.format(
        aCantera, struct.a))
    print('Maximum temperature: Cantera {:.1f}, ember {:.1f} K'.format(
        jet.T.max(), struct.T.max()))

    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    ax = ax.ravel()
    ax[0].plot(jet.grid, jet.T, 'k--', label='Cantera')
    ax[0].plot(struct.x, struct.T, label='ember')
    if fixedTemperature:
        ax[0].plot(xProfile, TProfile, 'r:', label='T_of_x-L500.csv')
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
    plt.suptitle('Impinging jet, p={:.0f} Pa, mdot={:.4f} kg/m²s, {}'.format(
        pressure, mdot, 'fixed T profile' if fixedTemperature else 'energy equation'))
    plt.tight_layout()
    plt.savefig(output + '/comparison.png')
    plt.close()
