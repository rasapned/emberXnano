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

## Particle moments module (nucleation / passive-scalar transport)

Extra scalar fields transported alongside `T`, `U`, `Y` — rows
`kMoments = nSpec+2` .. `kMoments+nMoments-1` of `state`. Two moments are
currently used: `kN` (number, index 0) and `kV` (volume, index 1),
`src/readConfig.h`.

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
- `kV` = m³ of particle volume **per kg of gas**, deliberately built so that
  `particleDensity · V` equals the particle-phase mass fraction, preserving
  the invariant `sum(Y) + particleDensity·V == 1`. Never multiply this by
  Avogadro's number.
- Cross terms (Soret/thermal-diffusion coupling) are explicitly zeroed for
  moments (`flameSolver.cpp`, `updateCrossTerms`) — only plain Fickian
  diffusion applies.
- Particle diffusivity (`updateParticleDiffusivity()` in `flameSolver.cpp`) is
  computed per grid point from the local mean particle size
  (`d_p` from `V/(N·Nₐ)`) via Stokes-Einstein with Cunningham slip correction,
  falling back to the constant `options.momentDiffusivity` wherever `N` or
  `V` ≤ 0 (no particles present yet). Because `D_p ∝ 1/d_p²` in the
  free-molecular regime, this ratio is numerically sensitive where `N, V` are
  both near zero (e.g. domain edges) — treat large diffusivities there with
  suspicion; likely a `V/N` noise artifact rather than real physics.

## Working conventions

- Commit/push only when explicitly asked; this session's git user is
  `sppicwie`, default branch `main`.
- Prefer small, targeted edits matching existing style (this codebase mixes
  Eigen `dvec`/`dmatrix`, raw loops with TBB `parallel_for`, and Cantera's
  kmol-based unit convention throughout — match whichever a given file
  already uses).
