#!/usr/bin/python
"""
Burner-stabilized flame in ember: single premixed inlet at a prescribed mass
flux (no opposing stream), unstrained (a=0), free/adiabatic outlet.

Ember's initial guess is a crude, smooth tanh ramp from unburned to
equilibrium -- deliberately NOT the fully converged Cantera BurnerFlame
solution (solved separately below, purely as a reference to validate
against). Ember integrates the transient PDE via Strang-split sub-stepping
(convection/diffusion/production evaluated as separate, sequentially-coupled
sub-steps, unlike Cantera's fully-implicit steady Newton solve), so handing
it an already-converged profile means its sharply-peaked, near-zero-balanced
trace radicals (O, H, HO2, ...) get evaluated by a single operator alone for
a sub-step, before any splitConst cross-term estimate even exists -- easily
enough to drive them slightly negative. A crude ramp has nothing precariously
balanced yet, so it can't undershoot through zero this way; ember is designed
to develop the real balance over physical time.

Because ember's automatic profile generator ties the fixed-left mass flux to
the strain rate (u = a*z potential-flow relation), it has no formula for V
when a=0 (see the premixed branch of ConcreteConfig.generateInitialProfiles
in input.py -- U ends up identically 0, which makes its finite-difference
build of V identically 0 too). So we build the initial profile by hand with
haveProfiles=True.
"""
from ember import *
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct

output = 'run/burnerFla_test'

# Same reactant composition as single_Igor_nucl.py
X_FEC5O5 = 0.0005
VS_total = 400 + 400 + 600
VS_total /= (1 - X_FEC5O5)
X_H2 = 400 / VS_total
X_O2 = 400 / VS_total
X_AR = 1.0 - (X_H2 + X_O2 + X_FEC5O5)
react = f"H2:{X_H2:.8e}, AR:{X_AR:.8e}, FEC5O5:{X_FEC5O5:.8e}, O2:{X_O2:.8e}"
print("Reactants:", react)

mechanism = 'Iron_elte_Syngas-newTransp.yaml'
pressure = 3000.0
Tu = 300.0
u_in = 1.1
xLeft = 0.0
xRight = 0.065
nPoints = 200

# --- Reference solution: Cantera BurnerFlame, same mechanism/composition ---
gas = ct.Solution(mechanism)
gas.TPX = Tu, pressure, react
rho_u = gas.density
mdot = rho_u * u_in
print(f"rho_u = {rho_u:.6f} kg/m3 -> fixed V (mass flux) = {mdot:.6e} kg/m2/s")

flame = ct.BurnerFlame(gas, width=xRight)
flame.burner.mdot = mdot
flame.transport_model = 'mixture-averaged'
flame.radiation_enabled = False
flame.set_refine_criteria(
    ratio=3.0,     # default 4.0 (smaller -> more aggressive)
    slope=0.06,    # default 0.04
    curve=0.12,    # default 0.08
    prune=0.01     # default 0.0 (prune disabled)
)
flame.solve(loglevel=1, refine_grid=True, auto=True)

print(f"Cantera BurnerFlame: unburned velocity = {flame.velocity[0]:.4f} m/s "
      f"(u_in = {u_in} m/s), Tb = {flame.T[-1]:.1f} K")

# --- Build ember's initial guess: a crude, smooth tanh ramp (see module
# docstring for why this -- not the converged Cantera profile -- is what
# gets handed to Ember) --------------------------------------------------
x_flame_guess = xLeft + 0.15 * (xRight - xLeft)
flameThickness_guess = 0.02 * (xRight - xLeft)

gas.TPX = Tu, pressure, react
Yu = gas.Y.copy()
gas.equilibrate('HP')
Tb = gas.T
Yb = gas.Y.copy()
print(f"Adiabatic flame temperature: {Tb:.1f} K")

x = np.linspace(xLeft, xRight, nPoints)
s = 0.5 * (1 + np.tanh((x - x_flame_guess) / flameThickness_guess))
T0 = Tu + (Tb - Tu) * s
Y0 = np.outer(1 - s, Yu) + np.outer(s, Yb)   # shape (nPoints, nSpecies)
V0 = np.full(nPoints, mdot)                  # constant mass flux (a=0 planar continuity)
U0 = np.zeros(nPoints)                       # no strain -> no radial velocity gradient

conf = Config(
    Paths(outputDir=output),
    Chemistry(mechanismFile=mechanism,
              transportModel='Mix'),
    General(twinFlame=False,
            flameGeometry='planar',
            nThreads=1,
            chemistryIntegrator='cvode',
            splittingMethod='strang',
            continuityBC='fixedLeft',
            unburnedLeft=True,
            fixedBurnedVal=False,
            # x=0 is the burner face here, a hard wall the domain must never
            # extend past -- same requirement as a twin/cylindrical flame's
            # symmetry plane, even though this isn't one. Without this,
            # OneDimGrid::addLeft() treats any steep near-inlet gradient
            # (expected right after a Dirichlet-pinned burner surface) as a
            # reason to keep extending the domain to negative x.
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
    # No boundaryTol/boundaryTolRm overrides here: for a=0, FlameSolver
    # switches to OneDimGrid::regridUnstrained (grid.cpp), which sizes the
    # domain from the qdot (heat release) profile via
    # Grid.unstrainedDownstreamWidth rather than boundaryTol -- loosening
    # boundaryTol doesn't touch the code path that's actually active here.
    Grid(
        # Tightened from 0.15/0.25: the FEC5O5 -> FE decomposition front is
        # much narrower than the main flame front and was going under-
        # resolved for long enough that OneDimGrid::addPoint's cubic-spline
        # interpolation (fit through all current points, no bound on
        # overshoot) produced a wildly unphysical FE spike when inserting a
        # new point into the gap. Refining sooner leaves less room for that
        # overshoot. 0.08/0.12 matches the "high accuracy" values documented
        # in input.py's Grid options.
        vtol=0.08,
        dvtol=0.12,
        gridMax=3e-4,
        addPointCount=2,
    ),
    CvodeTolerances(
        relativeTolerance=1e-7,
        speciesAbsTol=1e-12,
        energyAbsTol=1e-8,
        minimumTimestep=1e-16,
    ),
    TerminationCondition(
        tolerance=3e-3,     # relative RMS heat-release tolerance
        steadyPeriod=0.01,  # average over a longer window
        tMin=0.01,          # must exceed steadyPeriod and allow the flame to relax
        tEnd=0.6,           # several residence times (xRight/u_in ~ 0.06 s here)
    ),
)

if __name__ == '__main__':
    conf.run()

    struct = utils.load(output + '/profNow.h5')

    plt.figure()
    plt.plot(struct.x, struct.T, label='Ember')
    plt.plot(flame.grid, flame.T, '--', label='Cantera BurnerFlame')
    plt.xlabel('Position [m]')
    plt.ylabel('Temperature [K]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output + '/FinalTemperature.png')
    plt.close()
