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

output = 'run/burnerFla_particles'


# Same reactant composition as single_Igor_nucl.py
X_FEC5O5 = 0.0005
VS_total = 400 + 400 + 600
VS_total /= (1 - X_FEC5O5)
X_H2 = 400 / VS_total
X_O2 = 400 / VS_total
X_AR = 1.0 - (X_H2 + X_O2 + X_FEC5O5)
react = f"H2:{X_H2:.8e}, AR:{X_AR:.8e}, FEC5O5:{X_FEC5O5:.8e}, O2:{X_O2:.8e}"
print("Reactants:", react)

#mechanism = 'Iron_elte_Syngas-newTransp.yaml'
mechanism = 'opt-compact-mech-OF2.yaml'
pressure = 3000.0
Tu = 300.0
u_in = 1.1
xLeft = 0.0
xRight = 0.1
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

# Collisional nucleation: the mechanism only tracks gas-phase Fe clusters up
# to FE7 (FE, FE2, ..., FE7); any collision whose combined size reaches FE8
# or larger is treated as leaving the gas phase to form a particle nucleus.
maxClusterSize = 7
clusterSpecies = ['FE'] + [f'FE{n}' for n in range(2, maxClusterSize + 1)]
collisionSpeciesA, collisionSpeciesB = [], []
for i in range(1, maxClusterSize + 1):
    for j in range(i, maxClusterSize + 1):
        if i + j >= 8:
            collisionSpeciesA.append(clusterSpecies[i-1])
            collisionSpeciesB.append(clusterSpecies[j-1])
print("Nucleation channels:",
      list(zip(collisionSpeciesA, collisionSpeciesB)))

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
        # Tightened from 0.15/0.25: didnt work as expected
        vtol=0.15,
        dvtol=0.25,
        gridMax=3e-4,
        gridMin=1e-5,
        addPointCount=2,
    ),
    CvodeTolerances(
        relativeTolerance=1e-7,
        speciesAbsTol=1e-12,
        energyAbsTol=1e-8,
        minimumTimestep=1e-16,
    ),
    TerminationCondition(
        # Previous settings (kept for reference): the 'Q' (bulk heat-release)
        # measurement declared burnerFla_test.py "steady" while xFlame was
        # still visibly drifting at ~0.84 m/s -- 'Q' isn't sensitive to a
        # moving-but-otherwise-unchanging front. Matters even more here:
        # nucleation/coagulation rates depend on local residence time at
        # each temperature, so basing particle moments on a flame that
        # hasn't actually reached steady state would corrupt exactly the
        # numbers being compared against OpenFOAM.
        #tolerance=5e-3,     # relative RMS heat-release tolerance
        #steadyPeriod=0.01,  # average over a longer window
        #tMin=0.01,          # must exceed steadyPeriod and allow the flame to relax
        #tEnd=0.015,         # several residence times (xRight/u_in ~ 0.06 s here)
        # 'dTdt' checks whether the temperature FIELD itself is still
        # changing anywhere in space -- validated in burnerFla_test.py.
        measurement='dTdt',
        dTdtTol=0.6,        # tightened from default 10.0 [1/s]
        tMin=0.02,          # let the initial-guess transient clear first
        tEnd=0.2,           # generous headroom
    ),
    Particles(
        nMoments=2,
        momentBCLeft=0.0,
        particleDensity=7874.0,
        nucleation=NucleationChannel(
            collisionSpeciesA=collisionSpeciesA,
            collisionSpeciesB=collisionSpeciesB,
        ),
        coagulation=True,
    ),
)

if __name__ == '__main__':
    conf.run()

    struct = utils.load(output + '/profNow.h5')

    plt.figure()
    plt.plot(struct.x, struct.T, lw=2, label='Ember')
    plt.plot(flame.grid, flame.T, '--', lw=2, label='Cantera BurnerFlame')
    plt.xlabel('Position [m]')
    plt.ylabel('Temperature [K]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output + '/FinalTemperature.png')
    plt.close()

    # Plot the number particles and volume of the particles
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(struct.x, struct.numberDensity, 'b-', label='Number of particles')
    ax2.plot(struct.x, struct.particleDiameter, 'r-', label='Particle diameter')
    ax1.set_xlabel('Position [m]')
    ax1.set_ylabel('Particle number density / m⁻³', color='b')
    ax2.set_ylabel('Particle diameter / m', color='r')
    ax1.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    #plt.title(f'Ember planar, p=3000 Pa, a={a} s⁻¹')
    plt.tight_layout()
    plt.savefig(output + '/Particles_nucl.png')
    plt.close()

    # --- Species comparison: ember vs the Cantera BurnerFlame reference ---
    # One panel per species since the scales differ by orders of magnitude
    # (FE is a trace species from the 500 ppm FEC5O5; H2O/O2/H2 are majors).
    # Both solutions use the same mechanism, so species indices match.
    speciesPlot = ['H2O', 'FE', 'OH', 'H2', 'O2', 'H']
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True)
    for ax, name in zip(axes.flat, speciesPlot):
        k = gas.species_index(name)
        ax.plot(struct.x, struct.Y[k, :], lw=2, label='Ember')
        ax.plot(flame.grid, flame.Y[k, :], '--', lw=2, label='Cantera')
        ax.set_title(name)
        ax.set_ylabel('Mass fraction')
        ax.axhline(0, color='k', lw=0.5)
        ax.legend(loc='best', fontsize=8)
    for ax in axes.flat[len(speciesPlot):]:
        ax.axis('off')
    for ax in axes[-1, :]:
        ax.set_xlabel('Position [m]')
    fig.suptitle('Ember vs Cantera BurnerFlame: species mass fractions')
    fig.tight_layout()
    fig.savefig(output + '/Species_comparison.png')
    plt.close(fig)

    # Same thing on a log scale, so the trace species (FE, OH) are readable
    # next to the majors and any negative excursions show up as gaps.
    fig, ax = plt.subplots(figsize=(9, 6))
    for i, name in enumerate(speciesPlot):
        k = gas.species_index(name)
        color = f'C{i}'
        ax.semilogy(struct.x, struct.Y[k, :], color=color, lw=2, label=f'{name} (Ember)')
        ax.semilogy(flame.grid, flame.Y[k, :], '--', color=color, lw=2,
                    label=f'{name} (Cantera)')
    ax.set_xlabel('Position [m]')
    ax.set_ylabel('Mass fraction')
    ax.legend(loc='best', fontsize=8, ncol=2)
    fig.suptitle('Ember vs Cantera BurnerFlame: species (log scale)')
    fig.tight_layout()
    fig.savefig(output + '/Species_comparison_log.png')
    plt.close(fig)
