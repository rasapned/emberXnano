#!/usr/bin/python
"""
Burner-stabilized flame in ember: single premixed inlet at a prescribed mass
flux (no opposing stream), unstrained (a=0), free/adiabatic outlet.

Gas-phase setup as in burnerFla_test.py, plus iron particle moments
(nucleation from Fe clusters, coagulation) -- as burnerFla_particles.py, but
with particle-phase oxygen (nMoments=3) and heterogeneous condensation of Fe,
the Fe clusters and FeO onto the particles at the kinetic collision rate. No
gas-particle surface reactions (compare burnerFla_particles_oxid.py, which
adds them): particle oxygen comes from FeO condensation only.

Ember's initial guess is the converged Cantera BurnerFlame solution (solved
below, also used as the reference to compare against), which converges
fastest. burnerFla_test.py checked that shifted and tanh starts reach the
same steady state. The particle moments start from zero everywhere and build
up from nucleation; the Cantera gas profile has no particle sink for the Fe
clusters, so the gas phase relaxes a little as the particles form.

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

output = 'run/burnerFla_particles_cond'


# Same reactant composition as single_Igor_nucl.py
X_FEC5O5 = 0.0003
VS_total = 400 + 400 + 600
VS_total /= (1 - X_FEC5O5)
X_H2 = 400 / VS_total
X_O2 = 400 / VS_total
X_AR = 1.0 - (X_H2 + X_O2 + X_FEC5O5)
react = f"H2:{X_H2:.8e}, AR:{X_AR:.8e}, FEC5O5:{X_FEC5O5:.8e}, O2:{X_O2:.8e}"
print("Reactants:", react)

#mechanism = 'Iron_elte_Syngas-newTransp.yaml'
mechanism = 'opt-compact-mech-OF2_OFkinetics.yaml'
#mechanism = 'opt-compact-mech-OF2.yaml'
pressure = 3000.0
Tu = 300.0
u_in = 1.1
xLeft = 0.0
xRight = 0.1
nPoints = 200
gridMin = 1e-5

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
flame.set_refine_criteria(
    ratio=3.0,     # default 4.0 (smaller -> more aggressive)
    slope=0.06,    # default 0.04
    curve=0.12,    # default 0.08
    prune=0.01     # default 0.0 (prune disabled)
)
flame.solve(loglevel=1, refine_grid=True, auto=True)

print(f"Cantera BurnerFlame: unburned velocity = {flame.velocity[0]:.4f} m/s "
      f"(u_in = {u_in} m/s), Tb = {flame.T[-1]:.1f} K")

# --- Build ember's initial guess from the Cantera solution -----------------
# Merge Cantera's grid into the uniform one rather than resampling onto the
# uniform grid alone: Cantera's cells near the burner are a few um wide, and
# a uniform 200-point grid would smear out that gradient (and with it the
# burner heat loss). The uniform points keep dx <= gridMax downstream, where
# Cantera's grid is coarse; points closer than gridMin to their predecessor
# are dropped.
xAll = np.union1d(flame.grid[flame.grid < xRight],
                  np.linspace(xLeft, xRight, nPoints))
x = [xAll[0]]
for xi in xAll[1:]:
    if xi - x[-1] >= gridMin:
        x.append(xi)
x[-1] = xRight
x = np.array(x)
T0 = np.interp(x, flame.grid, flame.T)
Y0 = np.array([np.interp(x, flame.grid, flame.Y[k])
               for k in range(gas.n_species)]).T
Y0 = np.clip(Y0, 0.0, None)
Y0 /= Y0.sum(axis=1, keepdims=True)          # shape (nPoints, nSpecies)
# Ember takes the inlet stream (Tleft, Yleft) from the first point of the
# initial profile (FlameSolver::loadProfile). Cantera's x=0 value is the
# burner-face composition, already partly burned by back-diffusion (~25% of
# the H is in H2O there), so feeding it in as the inlet stream starves the
# flame of hydrogen. Put the actual reactants at x=0; the face relaxes to its
# flux balance within microseconds.
gas.TPX = Tu, pressure, react
T0[0] = Tu
Y0[0] = gas.Y
V0 = np.full(len(x), mdot)                   # constant mass flux (a=0 planar continuity)
U0 = np.zeros(len(x))                        # no strain -> no radial velocity gradient

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
            nThreads=16,
            chemistryIntegrator='cvode',
            # Strang leaves the cross terms out of drhodt, which leaks mass
            # at the burner (V ~ +7.6%); balanced splitting doesn't
            splittingMethod='balanced',
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
        gridMin=gridMin,
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
        dTdtTol=1.0,        # tightened from default 10.0 [1/s], as burnerFla_test.py
        tMin=0.02,         # let the start-up transient clear first, as burnerFla_test.py
        tEnd=0.2,           # generous headroom
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
        coagulation=True,
        condensationSpecies=clusterSpecies + ['FEO'],
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

    # --- Fe clusters: the only intended difference from Cantera is the Fe that
    # nucleation moves from the gas phase into the particle moments. The last
    # panel is the total Fe (element) mass fraction left in the gas phase;
    # Cantera has no particle sink, so the gap between the curves is the Fe in
    # particles.
    def gasFe(Y):
        kFe = gas.element_index('Fe')
        WFe = gas.atomic_weights[kFe]
        nFe = np.array([gas.n_atoms(k, kFe) for k in range(gas.n_species)])
        return (nFe * WFe / gas.molecular_weights) @ Y

    fig, axes = plt.subplots(2, 4, figsize=(17, 7), sharex=True)
    for ax, name in zip(axes.flat, clusterSpecies):
        k = gas.species_index(name)
        ax.plot(struct.x, struct.Y[k, :], lw=2, label='Ember')
        ax.plot(flame.grid, flame.Y[k, :], '--', lw=2, label='Cantera')
        ax.set_title(name)
        ax.set_ylabel('Mass fraction')
        ax.axhline(0, color='k', lw=0.5)
        ax.legend(loc='best', fontsize=8)
    ax = axes.flat[len(clusterSpecies)]
    ax.plot(struct.x, gasFe(struct.Y), lw=2, label='Ember')
    ax.plot(flame.grid, gasFe(flame.Y), '--', lw=2, label='Cantera')
    ax.set_title('Total Fe in gas phase (element)')
    ax.set_ylabel('Mass fraction')
    ax.legend(loc='best', fontsize=8)
    for ax in axes[-1, :]:
        ax.set_xlabel('Position [m]')
    fig.suptitle('Ember vs Cantera BurnerFlame: Fe clusters')
    fig.tight_layout()
    fig.savefig(output + '/FeClusters_comparison.png')
    plt.close(fig)

    # --- Overview: T, particle fields and mole fractions of the Fe-O-H and
    # radical species. Cantera has no particle phase, so N, Y_p and d_p are
    # ember only. moments[1:] are particle-phase Fe (and O) mass per mass of
    # gas [kg/kg], so their sum is the particle mass fraction.
    Xember = struct.Y / gas.molecular_weights[:, None]
    Xember /= Xember.sum(axis=0)
    fields = [
        ('T', 'T [K]', struct.T, flame.T),
        ('N', 'N [1/m³]', struct.numberDensity, None),
        ('Y_p', 'Particle mass fraction [-]', struct.moments[1:].sum(axis=0), None),
        ('d_p', 'd_p [nm]', struct.particleDiameter * 1e9, None),
    ]
    for name in ('FE', 'FEO', 'FEO2', 'FEOH', 'FEO2H2', 'O', 'H', 'OH'):
        k = gas.species_index(name)
        fields.append((f'x_{name}', 'Mole fraction', Xember[k], flame.X[k]))

    fig, axes = plt.subplots(3, 4, figsize=(18, 10), sharex=True)
    for ax, (title, ylabel, yE, yC) in zip(axes.flat, fields):
        ax.plot(struct.x, yE, lw=2, label='Ember')
        if yC is not None:
            ax.plot(flame.grid, yC, '--', lw=2, label='Cantera')
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.legend(loc='best', fontsize=8)
    for ax in axes[-1, :]:
        ax.set_xlabel('Position [m]')
    fig.suptitle('Ember (with particles) vs Cantera BurnerFlame')
    fig.tight_layout()
    fig.savefig(output + '/Particles_overview.png')
    plt.close(fig)

    # --- Particle composition: Fe and O mass in the particle phase and the
    # O/Fe atomic ratio (capped at 1.5 = Fe2O3, the last entry of phases).
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharex=True)
    axes[0].plot(struct.x, struct.moments[1], lw=2)
    axes[0].set_title('Particle-phase Fe')
    axes[0].set_ylabel('Mass per mass of gas [-]')
    axes[1].plot(struct.x, struct.moments[2], lw=2)
    axes[1].set_title('Particle-phase O')
    axes[1].set_ylabel('Mass per mass of gas [-]')
    axes[2].plot(struct.x, struct.particleOxygenRatio, lw=2)
    axes[2].axhline(1.5, color='k', ls='--', lw=0.8, label='Fe2O3')
    axes[2].axhline(4/3, color='0.5', ls=':', lw=0.8, label='Fe3O4')
    axes[2].axhline(1.0, color='0.7', ls=':', lw=0.8, label='FeO')
    axes[2].set_title('Particle O/Fe atomic ratio')
    axes[2].legend(loc='best', fontsize=8)
    for ax in axes:
        ax.set_xlabel('Position [m]')
    fig.tight_layout()
    fig.savefig(output + '/Particles_composition.png')
    plt.close(fig)
