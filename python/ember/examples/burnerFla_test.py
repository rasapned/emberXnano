#!/usr/bin/python
"""
Burner-stabilized flame in ember: single premixed inlet at a prescribed mass
flux (no opposing stream), unstrained (a=0), free/adiabatic outlet.

Ember's initial guess is chosen by initMode: by default the converged
Cantera BurnerFlame solution (solved below, also used as the reference to
validate against), which converges fastest. 'shifted' and 'tanh' start
elsewhere on purpose, to check that the steady state doesn't depend on the
initial condition; their plots also show the 'cantera' ember run
(referenceOutput) for direct comparison. Before the burner-face fixes, the
tanh start flashed back to the burner, was quenched and washed out (seen in
run/burnerFla_diffBC_corr).

Caveat with a converged start: ember integrates the transient PDE via
Strang-split sub-stepping, so on the first step the sharply-peaked,
near-zero-balanced trace radicals (O, H, HO2, ...) are advanced by each
operator alone before any splitConst estimate exists, which can push them
slightly negative. The Cantera profile is clipped to Y >= 0 and
renormalized; if a first-step undershoot still causes trouble, that is
where to look.

Because ember's automatic profile generator ties the fixed-left mass flux to
the strain rate (u = a*z potential-flow relation), it has no formula for V
when a=0 (see the premixed branch of ConcreteConfig.generateInitialProfiles
in input.py -- U ends up identically 0, which makes its finite-difference
build of V identically 0 too). So we build the initial profile by hand with
haveProfiles=True.
"""
from ember import *
import os
import glob
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct

output = 'run/burnerFla'

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
xRight = 0.08
nPoints = 200
gridMin = 1e-5

# Continue the run in the output folder (chosen by initMode below) from its
# last saved state, profNow.h5, in the same folder: the clock resumes at the
# saved time (Times.tStart) and the profNNNNNN.h5 numbering continues after
# the last file. tMin/tEnd in TerminationCondition are absolute times, so
# raise tEnd past the saved time to actually run longer. The profile is
# loaded here in Python rather than via InitialCondition(restartFile=...) so
# that its first point can be reset to the reactants: that point is the
# burner-face composition, and ember would otherwise take it as the inlet
# stream (see the note at the Cantera start below). out.h5 is rewritten from
# scratch by every run, so the earlier part is kept as out_until_<t>s.h5.
# Particle moments are not carried over.
restart = False

# --- Reference solution: Cantera BurnerFlame, same mechanism/composition ---
gas = ct.Solution(mechanism)
gas.TPX = Tu, pressure, react
rho_u = gas.density
mdot = rho_u * u_in
print(f"rho_u = {rho_u:.6f} kg/m3 -> fixed V (mass flux) = {mdot:.6e} kg/m2/s")

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

# --- Build ember's initial guess (see module docstring) ------------------
# 'cantera': converged Cantera BurnerFlame (fastest convergence)
# 'shifted': the Cantera profile moved initShift downstream, fresh reactants
#            filling the gap -- the flame must travel back to its standoff,
#            so the result shows the steady state doesn't depend on the start
# 'tanh':    crude ramp from reactants to adiabatic equilibrium (the hardest
#            start: wrong position, thickness and temperature)
initMode = 'tanh'
initShift = 4e-3  # [m], used by 'shifted' only
if initMode == 'shifted':
    output += f'_shifted{initShift*1e3:g}mm'
elif initMode == 'tanh':
    output += '_tanh'
# Converged run to compare the other starts against in the plots below
referenceOutput = 'run/burnerFla_canteraInit'

if initMode in ('cantera', 'shifted'):
    shift = initShift if initMode == 'shifted' else 0.0
    gridC = flame.grid + shift
    # Merge Cantera's grid into the uniform one rather than resampling onto
    # the uniform grid alone: Cantera's cells near the burner are a few um
    # wide, and a uniform 200-point grid would smear out that gradient (and
    # with it the burner heat loss). The uniform points keep dx <= gridMax
    # downstream, where Cantera's grid is up to ~7.5 mm coarse; points closer
    # than gridMin to their predecessor are dropped.
    xAll = np.union1d(gridC[gridC < xRight],
                      np.linspace(xLeft, xRight, nPoints))
    x = [xAll[0]]
    for xi in xAll[1:]:
        if xi - x[-1] >= gridMin:
            x.append(xi)
    x[-1] = xRight
    x = np.array(x)
    # Upstream of a shifted profile: fresh reactants at Tu
    gas.TPX = Tu, pressure, react
    T0 = np.interp(x, gridC, flame.T, left=Tu)
    Y0 = np.array([np.interp(x, gridC, flame.Y[k], left=gas.Y[k])
                   for k in range(gas.n_species)]).T
    Y0 = np.clip(Y0, 0.0, None)
    Y0 /= Y0.sum(axis=1, keepdims=True)      # shape (nPoints, nSpecies)
    # Ember takes the inlet stream (Tleft, Yleft) from the first point of the
    # initial profile (FlameSolver::loadProfile). Cantera's x=0 value is the
    # burner-face composition, already partly burned by back-diffusion
    # (~25% of the H is in H2O there), so feeding it in as the inlet stream
    # starves the flame of hydrogen. Put the actual reactants at x=0; the
    # face relaxes to its flux balance within microseconds.
    gas.TPX = Tu, pressure, react
    T0[0] = Tu
    Y0[0] = gas.Y
else:
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

V0 = np.full(len(x), mdot)                   # constant mass flux (a=0 planar continuity)
U0 = np.zeros(len(x))                        # no strain -> no radial velocity gradient

tStart = 0.0
firstFileNumber = 0
if restart:
    restartFile = output + '/profNow.h5'
    prof = utils.load(restartFile)
    tStart = float(prof.t)
    # Next number after the files actually present (profNow's own fileNumber
    # field can lag behind the last numbered file)
    firstFileNumber = 1 + max(int(os.path.basename(f)[4:10])
                              for f in glob.glob(output + '/prof[0-9]*.h5'))
    print(f"Restarting from {restartFile}: t = {tStart:.4f} s, "
          f"next file prof{firstFileNumber:06d}")
    x = prof.x
    T0 = prof.T.copy()
    Y0 = prof.Y.T.copy()                     # file stores (nSpecies, nPoints)
    V0 = prof.V
    U0 = prof.U
    gas.TPX = Tu, pressure, react
    T0[0] = Tu
    Y0[0] = gas.Y

conf = Config(
    Paths(outputDir=output),
    Chemistry(mechanismFile=mechanism,
              transportModel='Mix'),
    General(twinFlame=False,
            flameGeometry='planar',
            nThreads=1,
            chemistryIntegrator='cvode',
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
                      #restartFile=restartFile,
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
        # measurement declared this "steady" at t=0.012s while xFlame was
        # still visibly drifting at ~0.84 m/s -- 'Q' isn't sensitive to a
        # moving-but-otherwise-unchanging front.
        #tolerance=5e-3,     # relative RMS heat-release tolerance
        #steadyPeriod=0.01,  # average over a longer window
        #tMin=0.01,          # must exceed steadyPeriod and allow the flame to relax
        #tEnd=0.015,         # several residence times (xRight/u_in ~ 0.06 s here)
        # 'dTdt' checks whether the temperature FIELD itself is still
        # changing anywhere in space, which is what we actually care about.
        measurement='dTdt',
        dTdtTol=1.0,        # tightened from default 10.0 [1/s]
        tMin=0.025,          # let the initial-guess transient clear first; 0.02
                            # stopped the tanh-start run mid-transient
        tEnd=0.2,           # generous headroom; previous run needed >0.012s
                            # and was still drifting steadily
    ),
    # Dump a state file after every split sub-stage (conv1/diff1/diff2/prod/
    # diff3/conv2/diff4, see splitSolver.cpp) within this time window, so we
    # can see exactly which operator first introduces the bad values, rather
    # than only catching the coarser periodic profNow/prof000NNN snapshots
    # (which land between crashes and miss the actual divergence). Narrow
    # window bracketing the ~0.0044s crash time seen in the last run --
    # adjust once you know a new crash time.
    Times(tStart=tStart),
    OutputFiles(firstFileNumber=firstFileNumber,
                debugIntegratorStages=False),
    Debug(startTime=0.0043, stopTime=0.0046),
)

if __name__ == '__main__':
    if restart and os.path.exists(output + '/out.h5'):
        # out.h5 only holds the current run's time series; keep the earlier part
        os.replace(output + '/out.h5', output + f'/out_until_{tStart:.4f}s.h5')
    conf.run()

    struct = utils.load(output + '/profNow.h5')

    # For the non-default starts, compare against the converged 'cantera'
    # ember run: independent starts must reach the same steady state
    ref = None
    if output != referenceOutput and os.path.exists(referenceOutput + '/profNow.h5'):
        ref = utils.load(referenceOutput + '/profNow.h5')
        xc = ref.x[ref.x <= min(ref.x[-1], struct.x[-1])]
        dT = np.interp(xc, struct.x, struct.T) - np.interp(xc, ref.x, ref.T)
        print(f"\nvs reference ({referenceOutput}): max |dT| = {abs(dT).max():.2f} K "
              f"at x = {xc[abs(dT).argmax()]*1e3:.2f} mm")
        for name in ('H2', 'H', 'OH', 'H2O'):
            k = gas.species_index(name)
            print(f"  {name:4} at burner face: {struct.Y[k, 0]:.4e} vs {ref.Y[k, 0]:.4e} "
                  f"(rel. diff {struct.Y[k, 0]/ref.Y[k, 0] - 1:+.2e})")

    plt.figure()
    plt.plot(struct.x, struct.T, lw=2, label='Ember')
    if ref is not None:
        plt.plot(ref.x, ref.T, ':', lw=2, label='Ember (Cantera start)')
    plt.plot(flame.grid, flame.T, '--', lw=2, label='Cantera BurnerFlame')
    plt.xlabel('Position [m]')
    plt.ylabel('Temperature [K]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output + '/FinalTemperature.png')
    plt.close()
    
    plt.figure()
    plt.plot(struct.x, struct.V/struct.rho, lw=2, label='Ember')
    plt.plot(flame.grid, flame.velocity, '--', lw=2, label='Cantera BurnerFlame')
    plt.xlabel('Position [m]')
    plt.ylabel('Velocity [m/s]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output + '/FinalVelocity.png')
    plt.close()

    # --- Species comparison: ember vs the Cantera BurnerFlame reference ---
    # One panel per species since the scales differ by orders of magnitude
    # (FE is a trace species from the 500 ppm FEC5O5; H2O/O2/H2 are majors).
    # Both solutions use the same mechanism, so species indices match.
    speciesPlot = ['H2', 'O2', 'AR', 'H', 'O', 'H2O', 'FE', 'OH']
    fig, axes = plt.subplots(2, 4, figsize=(17, 7), sharex=True)
    for ax, name in zip(axes.flat, speciesPlot):
        k = gas.species_index(name)
        ax.plot(struct.x, struct.Y[k, :], lw=2, label='Ember')
        if ref is not None:
            ax.plot(ref.x, ref.Y[k, :], ':', lw=2, label='Ember (Cantera start)')
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
