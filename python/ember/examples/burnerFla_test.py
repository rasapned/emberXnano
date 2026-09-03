#!/usr/bin/python
"""
Burner-stabilized flame in ember, built from EXAMPLE_SINGLE / your
stagnation-point case: single premixed inlet at a prescribed velocity
(no opposing stream), unstrained (a=0), free/adiabatic outlet.

Because ember's automatic profile generator ties the fixed-left mass flux
to the strain rate (u = a*z potential-flow relation), it has no formula for
V when a=0. So we build the initial profile by hand with haveProfiles=True
and set V = rho_u * u_in directly at every point (consistent with a=0
planar continuity: rho*u = const through the whole domain).
"""
from ember import *
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct

output = 'run/test'

# --- Same reactant composition as your stagnation-point case ---
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
u_in = 1.1                 # m/s, prescribed inlet (burner) velocity
xLeft = 0.0
xRight = 0.065
nPoints = 200

# Initial guess for where the flame front sits and how thick it is.
# Not critical to get exact -- the transient solve will relax it to the
# physically correct standoff position -- but a wildly wrong guess can
# hurt convergence robustness.
x_flame_guess = xLeft + 0.15 * (xRight - xLeft)
flameThickness_guess = 0.02 * (xRight - xLeft)

# --- Reference states via Cantera, used only to build the initial guess ---
gas = ct.Solution(mechanism)
gas.TPX = Tu, pressure, react
rho_u = gas.density
Yu = gas.Y.copy()

mdot = rho_u * u_in
print(f"rho_u = {rho_u:.6f} kg/m3 -> fixed V (mass flux) = {mdot:.6e} kg/m2/s")

gas.equilibrate('HP')
Tb = gas.T
Yb = gas.Y.copy()
print(f"Adiabatic flame temperature: {Tb:.1f} K")

# u_in must be below the laminar flame speed or the flame blows off
# downstream; well above it and it flashes back into the burner.
gas.TPX = Tu, pressure, react
flame = ct.FreeFlame(gas, width=0.02)
#flame.solve(loglevel=0, auto=True)
S_L = flame.velocity[0]
print(f"Laminar flame speed S_L = {S_L:.4f} m/s (u_in = {u_in} m/s)")

# --- Hand-built initial profile (tanh transition, unburned -> burned) ---
x = np.linspace(xLeft, xRight, nPoints)
s = 0.5 * (1 + np.tanh((x - x_flame_guess) / flameThickness_guess))

T0 = Tu + (Tb - Tu) * s
Y0 = np.outer(1 - s, Yu) + np.outer(s, Yb)   # shape (nPoints, nSpecies)
V0 = np.full(nPoints, mdot)                  # constant mass flux (a=0 -> conserved)
U0 = np.zeros(nPoints)                       # no strain -> no radial velocity gradient

conf = Config(
    Paths(outputDir=output),
    Chemistry(mechanismFile=mechanism,
              transportModel='Mix'),
    General(twinFlame=False,
            flameGeometry='planar',
            nThreads=1,
     #       chemistryIntegrator='cvode',
     #       splittingMethod='strang',
            continuityBC='fixedLeft',   # fixes V at the inlet -- this is the default,
                                         # but explicit here since your original used
                                         # 'stagnationPoint' for the opposed-flow case
            unburnedLeft=True,
            fixedBurnedVal=False,       # right boundary is zero-gradient / free outflow
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
    Grid(
    vtol=0.15,
    dvtol=0.25,
    gridMax=3e-4,
    boundaryTol=2e-4,     # was 5e-5 default — this is the one driving unbounded growth
    boundaryTolRm=5e-5,
    addPointCount=2,
    ),
    CvodeTolerances(
        relativeTolerance=1e-7,
        speciesAbsTol=1e-12,
        energyAbsTol=1e-8,
        minimumTimestep=1e-16,
    ),
    TerminationCondition(
        tolerance=3e-3,     # relative RMS heat-release tolerance (was 5e-4, too tight)
        steadyPeriod=0.02,  # average over a longer window (default 0.002 too short here)
        tMin=0.1,           # must exceed steadyPeriod and allow the flame to relax
        tEnd=0.6,           # several residence times (xRight/u_in ~ 0.06 s here)
    ),
)

if __name__ == '__main__':
    #conf.run()

    struct = utils.load(output + '/profNow.h5')

    plt.figure()
    plt.plot(struct.x, struct.T)
    plt.xlabel('Position [m]')
    plt.ylabel('Temperature [K]')
    plt.savefig(output + '/FinalTemperature.png')
    plt.close()