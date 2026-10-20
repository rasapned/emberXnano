#!/usr/bin/python
"""
No-particle burner flame (as burnerFla_test.py) with ember AND Cantera both
on the fixed OpenFOAM mesh, to compare all three codes at the same resolution.

Grid: the face positions of the OpenFOAM blockMeshDict (627 cells over
160 mm, simpleGrading 9 -> 628 nodes, the first on the inlet face).
- ember: grid adaptation and boundary regridding disabled (regrid intervals
  never reached), Cantera start, forced to run t >= tMin.
- Cantera: BurnerFlame on that grid, refinement disabled; initial guess from
  a refined Cantera solution. Saved to <output>/cantera_OFgrid.csv.
"""
from ember import *
import numpy as np
import cantera as ct

output = 'run/burnerFla_OFmech_OFgrid'

# Same reactant composition as single_Igor_nucl.py
X_FEC5O5 = 0.0003
VS_total = 400 + 400 + 600
VS_total /= (1 - X_FEC5O5)
X_H2 = 400 / VS_total
X_O2 = 400 / VS_total
X_AR = 1.0 - (X_H2 + X_O2 + X_FEC5O5)
react = f"H2:{X_H2:.8e}, AR:{X_AR:.8e}, FEC5O5:{X_FEC5O5:.8e}, O2:{X_O2:.8e}"

mechanism = 'opt-compact-mech-OF2_OFkinetics_FeTransp.yaml'
pressure = 3000.0
Tu = 300.0
u_in = 1.5

# --- OpenFOAM mesh: 627 cells over 160 mm, geometric grading 9 -------------
nCells, xRight, grading = 627, 0.16, 9.0
r = grading**(1 / (nCells - 1))
x = np.concatenate(([0.0], np.cumsum(r**np.arange(nCells))))
x *= xRight / x[-1]
xLeft = x[0]

gas = ct.Solution(mechanism)
gas.TPX = Tu, pressure, react
mdot = gas.density * u_in

# --- Cantera: refined solution, then the same flame on the OpenFOAM grid ----
ref = ct.BurnerFlame(gas, width=xRight)
ref.burner.mdot = mdot
ref.transport_model = 'mixture-averaged'
ref.radiation_enabled = False
ref.set_refine_criteria(ratio=3.0, slope=0.06, curve=0.12, prune=0.0)
ref.solve(loglevel=0, refine_grid=True, auto=True)

gas.TPX = Tu, pressure, react
flame = ct.BurnerFlame(gas, grid=x)
flame.burner.mdot = mdot
flame.transport_model = 'mixture-averaged'
flame.radiation_enabled = False
zRef = ref.grid / ref.grid[-1]
flame.set_profile('velocity', zRef, ref.velocity)
flame.set_profile('T', zRef, ref.T)
for k, name in enumerate(gas.species_names):
    flame.set_profile(name, zRef, ref.Y[k])
flame.solve(loglevel=0, refine_grid=False)
assert len(flame.grid) == len(x)
print(f"Cantera on OpenFOAM grid: {len(flame.grid)} points, "
      f"T(100 mm) = {np.interp(0.1, flame.grid, flame.T):.1f} K")

# --- ember initial guess: the Cantera solution on the same grid --------------
T0 = flame.T.copy()
Y0 = np.clip(flame.Y.T.copy(), 0.0, None)
Y0 /= Y0.sum(axis=1, keepdims=True)
gas.TPX = Tu, pressure, react
T0[0] = Tu                      # inlet stream = reactants (see burnerFla_test.py)
Y0[0] = gas.Y
V0 = np.full(len(x), mdot)
U0 = np.zeros(len(x))

conf = Config(
    Paths(outputDir=output),
    Chemistry(mechanismFile=mechanism,
              transportModel='Mix'),
    General(twinFlame=False,
            flameGeometry='planar',
            nThreads=8,
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
    # Fixed grid: never regrid/adapt
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
        tMin=0.08,          # several flow-throughs, so slow minor species settle
        tEnd=0.2,
    ),
)

if __name__ == '__main__':
    conf.run()
    # Cantera solution on the OpenFOAM grid, for the comparison plots
    cols = ['x', 'T', 'u'] + [f'Y_{s}' for s in gas.species_names] + \
           [f'X_{s}' for s in gas.species_names]
    data = np.column_stack([flame.grid, flame.T, flame.velocity, flame.Y.T, flame.X.T])
    np.savetxt(output + '/cantera_OFgrid.csv', data, delimiter=',',
               header=','.join(cols), comments='')
    struct = utils.load(output + '/profNow.h5')
    print(f"ember final grid: {len(struct.x)} points (initial {len(x)}), "
          f"T(100 mm) = {np.interp(0.1, struct.x, struct.T):.1f} K")
