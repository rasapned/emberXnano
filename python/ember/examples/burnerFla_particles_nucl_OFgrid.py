#!/usr/bin/python
"""
Nucleation-only burner flame (as run/burnerFla_particles_nucl: the
burnerFla_particles_cond.py setup with coagulation=False and no
condensation) on the fixed OpenFOAM mesh, to check whether the remaining
ember-OpenFOAM differences in the nucleation case come from the grids.

Grid: the face positions of the OpenFOAM blockMeshDict (e.g.
OFcomparison/obsolete/compRes2/constant/polyMesh/blockMeshDict) -- 627
cells over 160 mm, simpleGrading 9 (last/first cell) -- so 628 nodes, the
first one on the inlet face as OpenFOAM's inlet patch. Grid adaptation and
boundary regridding are disabled (regrid intervals never reached), so the
grid stays exactly this one for the whole run.
"""
from ember import *
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct

output = 'run/burnerFla_particles_nucl_OFgrid'

# Same reactant composition as single_Igor_nucl.py
X_FEC5O5 = 0.0003
VS_total = 400 + 400 + 600
VS_total /= (1 - X_FEC5O5)
X_H2 = 400 / VS_total
X_O2 = 400 / VS_total
X_AR = 1.0 - (X_H2 + X_O2 + X_FEC5O5)
react = f"H2:{X_H2:.8e}, AR:{X_AR:.8e}, FEC5O5:{X_FEC5O5:.8e}, O2:{X_O2:.8e}"
print("Reactants:", react)

mechanism = 'opt-compact-mech-OF2_OFkinetics.yaml'
pressure = 3000.0
Tu = 300.0
u_in = 1.5

# --- OpenFOAM mesh: 627 cells over 160 mm, geometric grading 9 -------------
nCells, xRight, grading = 627, 0.16, 9.0
r = grading**(1 / (nCells - 1))
x = np.concatenate(([0.0], np.cumsum(r**np.arange(nCells))))
x *= xRight / x[-1]
xLeft = x[0]
print(f"OpenFOAM grid: {len(x)} points, dx = {np.diff(x)[0]*1e6:.1f} .. "
      f"{np.diff(x)[-1]*1e6:.1f} um")

# --- Reference solution: Cantera BurnerFlame, same mechanism/composition ---
gas = ct.Solution(mechanism)
gas.TPX = Tu, pressure, react
rho_u = gas.density
mdot = rho_u * u_in
print(f"rho_u = {rho_u:.6f} kg/m3 -> fixed V (mass flux) = {mdot:.6e} kg/m2/s")

# Wider than ember's domain so the initial profile never needs extrapolating
flame = ct.BurnerFlame(gas, width=xRight+0.02)
flame.burner.mdot = mdot
flame.transport_model = 'mixture-averaged'
flame.radiation_enabled = False
flame.set_refine_criteria(ratio=3.0, slope=0.06, curve=0.12, prune=0.01)
flame.solve(loglevel=1, refine_grid=True, auto=True)

print(f"Cantera BurnerFlame: unburned velocity = {flame.velocity[0]:.4f} m/s "
      f"(u_in = {u_in} m/s), Tb = {flame.T[-1]:.1f} K")

# --- Initial guess: Cantera solution on the OpenFOAM grid ------------------
T0 = np.interp(x, flame.grid, flame.T)
Y0 = np.array([np.interp(x, flame.grid, flame.Y[k])
               for k in range(gas.n_species)]).T
Y0 = np.clip(Y0, 0.0, None)
Y0 /= Y0.sum(axis=1, keepdims=True)          # shape (nPoints, nSpecies)
# Ember takes the inlet stream from the first point of the initial profile;
# put the actual reactants there, not Cantera's burner-face composition
# (see burnerFla_particles_cond.py)
gas.TPX = Tu, pressure, react
T0[0] = Tu
Y0[0] = gas.Y
V0 = np.full(len(x), mdot)                   # constant mass flux (a=0 planar continuity)
U0 = np.zeros(len(x))                        # no strain -> no radial velocity gradient

# Collisional nucleation: clusters up to FE7 in the mechanism; collisions
# reaching FE8 or larger form a particle nucleus.
maxClusterSize = 7
clusterSpecies = ['FE'] + [f'FE{n}' for n in range(2, maxClusterSize + 1)]
collisionSpeciesA, collisionSpeciesB = [], []
for i in range(1, maxClusterSize + 1):
    for j in range(i, maxClusterSize + 1):
        if i + j >= 8:
            collisionSpeciesA.append(clusterSpecies[i-1])
            collisionSpeciesB.append(clusterSpecies[j-1])

conf = Config(
    Paths(outputDir=output),
    Chemistry(mechanismFile=mechanism,
              transportModel='Mix'),
    General(twinFlame=False,
            flameGeometry='planar',
            nThreads=16,
            chemistryIntegrator='cvode',
            splittingMethod='balanced',
            continuityBC='fixedLeft',
            unburnedLeft=True,
            fixedBurnedVal=False,
            fixedLeftLocation=True,
    ),
    InitialCondition(reactants=react,
                      pressure=pressure,
                      Tu=Tu,
                      xLeft=xLeft,
                      xRight=xRight,
                      haveProfiles=True,
                      x=x, T=T0, U=U0, V=V0, Y=Y0,
                      ),
    StrainParameters(initial=0.0,
                      final=0.0),
    # Fixed grid: never regrid/adapt (FlameSolver::step only adapts once
    # regridTimeInterval or regridStepInterval is reached)
    Times(regridTimeInterval=1e9,
          regridStepInterval=10**9),
    CvodeTolerances(
        relativeTolerance=1e-7,
        speciesAbsTol=1e-12,
        energyAbsTol=1e-8,
        minimumTimestep=1e-16,
    ),
    TerminationCondition(
        measurement='dTdt',
        dTdtTol=0.03,
        tMin=0.02,
        tEnd=0.2,
    ),
    Particles(
        nMoments=3,                 # N, Fe mass, O mass
        momentBCLeft=0.0,
        metal='Fe',
        # (O/metal atomic ratio, density [kg/m3]): Fe, FeO, Fe3O4, Fe2O3
        phases=[(0.0, 7874.0), (1.0, 5745.0), (4/3, 5170.0), (1.5, 5240.0)],
        nucleation=NucleationChannel(
            collisionSpeciesA=collisionSpeciesA,
            collisionSpeciesB=collisionSpeciesB,
        ),
        coagulation=False,
        minParticleDiameter=0.4e-9,
    ),
)

if __name__ == '__main__':
    conf.run()

    struct = utils.load(output + '/profNow.h5')
    print(f"Final grid: {len(struct.x)} points (initial {len(x)}), "
          f"max |x - x_OF| = {np.abs(struct.x - x).max() if len(struct.x) == len(x) else np.nan:.2e} m")

    # Overview: T, particle fields and Fe-O-H / radical mole fractions
    Xember = struct.Y / gas.molecular_weights[:, None]
    Xember /= Xember.sum(axis=0)
    fields = [
        ('T', 'T [K]', struct.T, flame.T),
        ('N', 'N [1/m³]', struct.numberDensity, None),
        ('Y_p', 'Particle mass fraction [-]', struct.moments[1:3].sum(axis=0), None),
        ('d_p', 'd_p [nm]', struct.particleDiameter * 1e9, None),
    ]
    for name in ('FE', 'FEO', 'FEO2', 'FEOH', 'FEO2H2', 'O', 'H', 'OH'):
        k = gas.species_index(name)
        fields.append((f'x_{name}', 'Mole fraction', Xember[k], flame.X[k]))

    fig, axes = plt.subplots(3, 4, figsize=(18, 10), sharex=True)
    for ax, (title, ylabel, yE, yC) in zip(axes.flat, fields):
        ax.plot(struct.x, yE, lw=2, label='Ember (OpenFOAM grid)')
        if yC is not None:
            ax.plot(flame.grid, yC, '--', lw=2, label='Cantera')
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.legend(loc='best', fontsize=8)
    for ax in axes[-1, :]:
        ax.set_xlabel('Position [m]')
    ax.set_xlim(0, 0.1)
    fig.suptitle('Ember (nucleation, OpenFOAM grid) vs Cantera BurnerFlame')
    fig.tight_layout()
    fig.savefig(output + '/Particles_overview.png')
    plt.close(fig)
