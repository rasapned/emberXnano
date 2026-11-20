#!/usr/bin/python
"""
Test against Sellmann et al., Proc. Combust. Inst. 37 (2019) 1241-1248,
flame A (H2/O2/Ar 400/400/600 sccm, 200 ppm Fe(CO)5, 3000 Pa, T_in = 300 K)
in the down-firing configuration (DFF): only the iron oxide particles of
the recombination zone, with their bi-modal model (Jeong & Choi 2003):
every gas-phase Fe2O3 molecule becomes a mode-1 particle (monomer
conversion), coagulation within and between the modes (Fuchs kernel,
fractal collision diameter, Df = 1.8), sintering of the mode-2 aggregates,
size-dependent diffusion. No early (Fe cluster) particles, as in the paper.

The DFF's detached stagnation point (buoyancy) is modelled as the wall of
an impinging jet at L = 46 mm (paper, Sec. 4.2); the gas temperature is the
paper's simulated centre-line profile, digitised from Figs. 4 and 9
(Sellmann2019/T_of_x-DFF_flameA.csv). The flow is prescribed too: the radial
strain U(x) follows from the paper's centre-line axial mass flux rho*u
(u from Fig. 4, Sellmann2019/Ux_DFF_flameA.csv) by continuity, so ember
has the paper's velocity and residence time. Cantera's ImpingingJet (gas
only, same fixed T) gives the initial composition and the density.

Differences to the paper: the monomer diameter follows from the Fe2O3
mass and density (0.46 nm; paper 0.4 nm); the paper's u crosses zero at
44.5 mm, so the last 1.5 mm to the wall are stagnant (U = V = 0);
particles do not cross the wall (zero flux), so late mass piles up there;
the run stops at tEnd.

    python Sellmann2019_DFF_flameA.py
    python plot_Sellmann2019_DFF_flameA.py
"""
from ember import *
import os
import numpy as np
import cantera as ct
from scipy.ndimage import gaussian_filter1d

output = 'run/Sellmann2019/DFF_flameA_fixedU'
mech = 'opt-compact-mech-OF2_OFkinetics_FeTransp.yaml'

# Flame A (Table 1): 200 ppm FEC5O5, volume flow rates in sccm
X_FEC5O5 = 0.0002
VS_total = (400 + 400 + 600) / (1 - X_FEC5O5)
X_H2 = 400 / VS_total
X_O2 = 400 / VS_total
X_AR = 1.0 - (X_H2 + X_O2 + X_FEC5O5)
react = f"H2:{X_H2:.8e}, AR:{X_AR:.8e}, FEC5O5:{X_FEC5O5:.8e}, O2:{X_O2:.8e}"

pressure = 3000.0  # [Pa]

# Simulated temperature profile of the paper: burner at x = 0, stagnation
# point (ember wall) at x = L = 46 mm
xProfile, TProfile = np.genfromtxt('Sellmann2019/T_of_x-DFF_flameA.csv',
                                   delimiter=',', comments='#').T
L = xProfile[-1]
Tburner = TProfile[0]
Tplate = TProfile[-1]

# Inlet mass flux: standard (1 bar, 273.15 K) volume flow through the
# 36 mm burner matrix
gas = ct.Solution(mech, transport_model='mixture-averaged')
gas.TPX = 273.15, 1e5, react
mdot = gas.density * VS_total / 60e6 / (0.25 * np.pi * 0.036**2)
gas.TPX = Tburner, pressure, react
rhou = gas.density

# *** Cantera (gas phase only), fixed temperature profile ***
jet = ct.ImpingingJet(gas=gas, width=L)
jet.inlet.mdot = mdot
jet.surface.T = Tplate
jet.set_refine_criteria(ratio=3, slope=0.05, curve=0.1, prune=0.02)
jet.flame.set_fixed_temp_profile(xProfile / L, TProfile)
jet.energy_enabled = False
jet.solve(loglevel=0, refine_grid=True)
aCantera = 2.0 * np.sqrt(-jet.L[0] / rhou)  # beta = 2 (disc)
print('Cantera: {} points, equivalent strain rate a = {:.3f} 1/s, inlet velocity {:.3f} m/s'.format(
    len(jet.grid), aCantera, jet.velocity[0]))

# Prescribed flow: axial mass flux m = rho*u of the paper (u from Fig. 4,
# rho of the gas-only solution at the same T), m = 0 beyond the paper's
# stagnation point; radial strain from continuity dm/dx = -beta*rho*U
# (beta = 2, disc). The inlet mass flux (m at x = 0) replaces mdot.
# m, not u, is interpolated between the markers (2 mm apart): u jumps with
# the steep T rise near the burner, m does not. U is then smoothed
# (Gaussian, 0.5 mm): its steps at the markers and its corner at the
# stagnation point made the grid adaptation refine to the minimum spacing
# in hot gas, which cut the diffusion time step ~4x (u at the markers
# stays within 2 %). m follows from continuity, integrated from the wall.
xU, uU = np.genfromtxt('Sellmann2019/Ux_DFF_flameA.csv', delimiter=',', comments='#').T
xFlow = np.linspace(0, L, int(round(L / 1e-4)) + 1)
rhoProfile = np.interp(xFlow, jet.grid, jet.density)
mU = np.interp(xU, jet.grid, jet.density) * np.maximum(uU, 0)
UProfile = gaussian_filter1d(-np.gradient(np.interp(xFlow, xU, mU), xFlow) / (2 * rhoProfile),
                             5.0, mode='nearest')
mProfile = np.r_[np.cumsum((2 * rhoProfile * UProfile * np.gradient(xFlow))[:0:-1])[::-1], 0.0]
print('prescribed flow: inlet mass flux {:.5f} kg/m2/s (mdot {:.5f}), U from {:.1f} to {:.1f} 1/s'.format(
    mProfile[0], mdot, UProfile.min(), UProfile.max()))

# ember takes the inlet stream composition from the first point of the given
# profile, so that point must hold the supplied mixture
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
                     U=np.interp(jet.grid, xFlow, UProfile),
                     V=np.interp(jet.grid, xFlow, mProfile),
                     Y=Y0),
    ImpingingJet(wallTemperature=Tplate,
                 temperatureProfile=(xProfile, TProfile),
                 strainProfile=(xFlow, UProfile)),
    StrainParameters(initial=aCantera,     # not used with strainProfile
                     final=aCantera),
    CvodeTolerances(relativeTolerance=1e-7,
                    speciesAbsTol=1e-12,
                    energyAbsTol=1e-8,
                    minimumTimestep=1e-16),
    # Fixed end time: the late mass keeps piling up at the wall, so the
    # 'moments' criterion would not settle (and its validation does not
    # yet accept a bimodal-only case)
    TerminationCondition(measurement=None,
                         tEnd=0.2),
    Particles(
        name='late',
        model='bimodal',
        metal='Fe',
        phases=[(1.5, 5240.0)],               # Fe2O3 [kg/m3]
        nucleation=NucleationChannel(
            conversionSpecies=['FE2O3(s)'],   # gas-phase Fe2O3 of the mechanism
            conversionRate=1e5,               # [1/s]: every molecule a monomer
        ),
        coagulation=True,
        fractalDimension=1.8,                 # Poliak et al. (2015)
        # tau = A T d_pp^4 exp(Ta/T): A [s/(K m^4)], Ta [K]
        sintering=(8e16, 30000.0),
    ),
)



def writeInputProfiles():
    """Profiles fed to ember, for plot_Sellmann2019_DFF_flameA.py"""
    os.makedirs(output, exist_ok=True)
    np.savetxt(output + '/inputProfiles.csv',
               np.c_[xFlow, np.interp(xFlow, xProfile, TProfile), UProfile, mProfile, mProfile / rhoProfile],
               delimiter=',', header='x [m], T [K], U [1/s], rho*u [kg/m2/s], u [m/s]')


if __name__ == '__main__':
    writeInputProfiles()
    conf.run()
