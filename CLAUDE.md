# CLAUDE.md

Behavioral guidelines to reduce common LLM coding mistakes. Merge with project-specific instructions as needed.

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.

# Ember (emberXnano fork)

Quasi-1D unsteady strained-flame solver (C++ core + Cython/Python interface), extended
in this fork with a particle "moments" module for gas-to-particle nucleation
(two-way coupled with the gas phase). Upstream: https://github.com/speth/ember

## Build system: SCons, not pip/CMake

Everything is driven by `SConstruct` at the repo root. Do not invoke `pip install .`,
`setup.py`, or `cmake` directly — `python/setup.py`/`setup.cfg` are generated
artifacts consumed *by* the SCons build, not standalone entry points.

```
scons build      # compiles core C++ lib + cythonizes/compiles python/ember/_ember.pyx
scons test        # C++ and python test suites
scons install     # builds wheel + pip-installs it
```

Run these from the repo root (not from `python/`).

- `scons build` recompiles `src/*.cpp` into the core lib, and separately
  cythonizes `python/ember/_ember.pyx` → `_ember.cpp` → compiles the extension
  module, writing it **in place** to `python/ember/_ember*.so`.
- If the active Python environment has `ember` installed in **editable/develop
  mode** (pointing at this source tree), `scons build` alone is enough — the
  new `.so` is picked up immediately.
- If it's a non-editable install (a copied `.so` in `site-packages`), you also
  need `scons install` to rebuild the wheel and `pip install` it, or you'll
  silently keep running the old binary.
- After editing `python/ember/_ember.pyx` specifically, the cythonize step
  re-triggers automatically via SCons's dependency tracking; you don't need to
  run Cython by hand.

## Layout

- `src/` — C++ solver core.
  - `flameSolver.cpp/h` — top-level driver; owns the `state` matrix
    (`nVars × nPoints`), Strang-split orchestration, thermo/transport property
    updates.
  - `splitSolver.cpp/h` — generic operator-split integration (diffusion →
    convection → diffusion → production → diffusion → convection → diffusion,
    each stage measured as a `delta` and turned into a `splitConst` for the
    other operators).
  - `convectionSystem.cpp/h`, `diffusionSystem.cpp/h`, `sourceSystem.cpp/h` —
    the three split operators. Diffusion solves
    `ẏ = B·d/dx(D·dy/dx) + C` per row of `state` (`B` = prefactor,
    typically `1/ρ` since state variables are mass-specific; `D` = diffusion
    coefficient). Convection: one solver per scalar (species, moments), T/U
    coupled together. Production: one coupled ODE per grid point (chemistry +
    nucleation).
  - `readConfig.h` — `ConfigOptions` (mirrors `python/ember/input.py`), plus
    `kN`/`kV` moment-index constants.
- `python/ember/` — Python package.
  - `_ember.pyx` — Cython interface layer wrapping the C++ `FlameSolver`;
    read-only properties (`T`, `Y`, `moments`, `rho`, ...) and derived
    quantities that are pure post-processing (no C++ change needed) belong
    here.
  - `input.py` — `ConfigOptions` and sub-option classes (`Particles`,
    `NucleationChannel`, etc.) with docstrings — the canonical description of
    every option's meaning/units.
  - `output.py` — `StateWriter`/`TimeSeriesWriter`, what gets written to
    `.h5`/`.npz` output files.
  - `examples/` — example case scripts (e.g. `single_Igor_nucl.py`).
- `test/` — C++ (`test_*.cpp`) and Python (`test/python`) test suites, run via
  `scons test`.

## Particle moments module (nucleation / surface reactions / transport)

Extra scalar fields transported alongside `T`, `U`, `Y` — rows
`kMoments = nSpec+2` .. `kMoments+nMoments-1` of `state`. Up to three
moments (`src/readConfig.h`): `kN` (number, index 0), `kM` (particle-phase
metal mass, index 1) and `kO` (particle-phase O mass, index 2; only with
`nMoments = 3`, otherwise particles are pure metal). The metal is set per
case (`Particles.metal`, an element name of the mechanism, e.g. `'Fe'`);
nothing in the code is specific to iron — only the defaults are.

**Unit convention — easy to get wrong:** moments are transported exactly like
mass fractions (`B = 1/ρ` in diffusion, same convection-solver treatment as a
species), so they are **mass-specific** quantities, not raw densities:
- `kN` = kmol of particles **per kg of gas** (a molar count, consistent with
  the kmol convention used everywhere else for species — see
  `sourceSystem.cpp::computeNucleationRates()`). To get a true number density
  [particles/m³]: `n = kN · ρ · Nₐ`, where `Nₐ` must be Cantera's
  **kmol-based** Avogadro constant (`Cantera::Avogadro` in C++,
  `cantera.avogadro` in Python, ≈6.022e26/kmol) — the everyday 6.022e23/mol
  value is 1000× off. This conversion is done at the Python interface layer
  (`FlameSolver.numberDensity` property in `_ember.pyx`), not inside the
  transport equations.
- `kM`, `kO` = kg of particle-phase metal / O **per kg of gas**. All source
  terms (nucleation, surface reactions) are built from species molecular
  weights so they conserve mass exactly between gas and particles
  (`sum(Y) + mM + mO` is invariant under the source terms). Never multiply
  these by Avogadro's number.
- Particle **size** is derived, not transported: mass per particle
  `(mM+mO)/(N·Nₐ)` divided by a density interpolated linearly in the
  O/metal ratio between the case-defined `Particles.phases`
  (`(ratio, density)` pairs from pure metal up to the most oxidized phase,
  whose ratio is also the oxidation cap) — `src/particleUtils.cpp`.
- Surface reactions (`SurfaceReaction` in `input.py`,
  `computeSurfaceReactionRates()` in `sourceSystem.cpp`): Arrhenius rate ×
  particle surface area × reactant concentration × available-site fraction;
  the particle change (Δmetal, ΔO) is derived from reactant − product
  composition. Below `minParticleDiameter` mass removal also removes
  particles from N (disintegration).
- Heterogeneous condensation (`Particles.condensationSpecies`,
  `computeCondensationRates()` in `sourceSystem.cpp`): whole gas molecules
  (metal, clusters, oxides) added to particles at the hard-sphere
  gas-particle collision rate, collision diameters from the mechanism's
  transport data, as for collision nucleation.
- Cross terms (Soret/thermal-diffusion coupling) are explicitly zeroed for
  moments (`flameSolver.cpp`, `updateCrossTerms`) — only plain Fickian
  diffusion applies.
- Particle diffusivity (`updateParticleDiffusivity()` in `flameSolver.cpp`) is
  computed per grid point from the local mean particle size
  (`d_p` from the size closure above) via Stokes-Einstein with Cunningham
  slip correction, falling back to the constant `options.momentDiffusivity`
  wherever `N` or `mM` ≤ 0 (no particles present yet). Because
  `D_p ∝ 1/d_p²` in the free-molecular regime, this ratio is numerically
  sensitive where `N, mM` are both near zero (e.g. domain edges) — treat
  large diffusivities there with suspicion; likely a mass/`N` noise artifact
  rather than real physics.

## Impinging jet (`ImpingingJet` option)

Counterpart of Cantera's `ImpingingJet` (Inlet1D + AxisymmetricFlow +
Surface1D), for `flameGeometry` `'disc'` (= Cantera) or `'planar'`, premixed
only. Example/validation: `python/ember/examples/example_impingingJet.py`.
- Fixed domain `[0, xRight]` (xLeft forced to 0): boundary regridding is
  skipped in `flameSolver.cpp` (same branch as quasi2d); interior adaptation
  still runs.
- Left = `BoundaryCondition::InletFlux` (burner-style flux balance for Y and
  moments, T fixed) plus U = 0 (plug flow). Right =
  `BoundaryCondition::Wall`: T = `wallTemperature`, U = 0, V = 0, zero-flux
  half control volume for Y/moments in `diffusionSystem.cpp` (T and U rows are
  overridden to FixedValue in `prepareIntegrators`), conservative wall term for
  the jCorr/Soret cross fluxes, and no chemistry at the wall node (as at the
  inlet). Particles: zero flux at the wall (no deposition/thermophoresis).
- Continuity is integrated from the wall (`ContinuityBoundaryCondition::Wall`,
  `rV[jj] = 0`), so the inlet mass flux `V[0]` is an output.
- Pressure curvature: Cantera solves for Λ; ember prescribes it through the
  strain rate, `Λ = -rhou·a²/β²` (so `a = β·sqrt(-Λ/rhou)` maps a Cantera
  solution to ember). With `ImpingingJet(massFlux=...)` a PI controller on
  log(a) (`FlameSolver::updateMassFluxControl`, `ControlledFunction` in
  `scalarFunction.h`) adjusts `a` until `V[0]` matches — steady state only,
  transients are not physical. `a` and `mdot` are in `out.h5`.
- The thin plug-flow layer at the inlet limits the split timestep: the default
  `globalTimestep = 2e-5` gave NaNs at a ≈ 900 1/s, 1 atm; 1e-5 was stable
  (2e-5 is fine at 3000 Pa, `single_Igor_impingingJet.py`).
- Use `splittingMethod='balanced'` (default) for steady comparisons: with
  `'strang'` the Igor case converged to an `a` 7.6% above Cantera's; balanced
  matched it to 0.02%.
- Prescribed temperature (`ImpingingJet(temperatureProfile=(x, T))`, Cantera's
  `set_fixed_temp_profile` + `energy_enabled = False`): `options.fixedTemperature`
  resets T to the interpolated profile in `setupStep` and zeroes dT/dt in all
  operators (convection UTW, source systems, energy diffusion B = 0, cross
  term, energy ddt rows). Example: `single_Igor_impingingJet_fixedT.py`
  (measured profile `T_of_x-L500.csv`, ember started from Cantera's solution).
  With T prescribed, `TerminationCondition(measurement='dTdt')` is meaningless
  (stops at tMin); use `measurement='moments'` (`momentsTol` [1/s]: RMS rate
  of every particle moment relative to its peak, `checkTerminationCondition`)
  for particle cases, e.g. `single_Igor_impingingJet_fixedT_particles.py`.
- Starting from given profiles (`haveProfiles=True`) the inlet stream
  composition `Yleft` is taken from `Y[:, 0]`; overwrite that point with the
  supplied mixture, since a converged burner/jet inlet node is depleted in
  fast diffusers (H2).
- `InletFlux` with the flame on the burner (strong back-diffusion through the
  inlet face) needed two conservation fixes, which also apply to burner
  flames: the UTW system's Wmx[0] now gets the inflow term (otherwise rho[0],
  hence the inlet velocity used by the species, drifted within each step and
  ~2.5% of the H atoms were lost), and the jCorr/Soret cross flux across the
  first face is credited to node 0. Check with an element balance: inflow
  `mdot*Z_in` vs radial outflow `∫ β·rho·U·Z dx`.

## Working conventions

- Commit/push only when explicitly asked; this session's git user is
  `sppicwie`, default branch `main`.
- Prefer small, targeted edits matching existing style (this codebase mixes
  Eigen `dvec`/`dmatrix`, raw loops with TBB `parallel_for`, and Cantera's
  kmol-based unit convention throughout — match whichever a given file
  already uses).
