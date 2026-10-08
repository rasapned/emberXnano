# Impinging jet in ember (`ImpingingJet` option)

Ember can now simulate a premixed jet impinging on a wall, the counterpart of Cantera's
`ImpingingJet`. Before, this case had to be approximated as a flame opposing an inert stream,
which has no wall: no no-slip, a stagnation plane that drifts, heat lost into the inert stream
instead of a wall at fixed temperature, and inert gas diffusing into the flame.

This note covers what was added, how it maps onto Cantera, two conservation bugs in the
burner inlet condition that were found and fixed on the way, and the validation runs.

## Usage

```python
ImpingingJet(wallTemperature=400.0,          # [K]
             massFlux=0.0272,                # target inlet mass flux [kg/m^2/s], or None
             temperatureProfile=(x, T))      # optional: prescribed T(x), no energy equation
```

- Premixed, `unburnedLeft`, `flameGeometry='disc'` (axisymmetric, as Cantera) or `'planar'`.
- The domain is fixed: `[0, InitialCondition.xRight]`, with `xLeft` forced to 0. Boundary
  regridding is skipped; interior grid adaptation still runs.
- Without `massFlux`, the strain rate from `StrainParameters` is used as given and the inlet
  mass flux is an output.

## How it maps onto Cantera

Cantera's `ImpingingJet` is `Inlet1D` + `AxisymmetricFlow` + `Surface1D` on a fixed interval:

| | Inlet (x = 0) | Wall (x = L) |
|---|---|---|
| Mass flux | = ṁ | = 0 |
| Radial velocity (ember `U`) | 0 (plug flow) | 0 (no slip) |
| T | inlet T | wall T |
| Species, particle moments | flux balance `ρuY + j = ṁ·Y_in` (`InletFlux`) | zero flux |
| Chemistry at the node | none | none |

In ember:

- Left boundary: the existing `BoundaryCondition::InletFlux`, plus U held at 0.
- Right boundary: new `BoundaryCondition::Wall`. T and U are fixed; species and moments use a
  zero-flux half control volume (`diffusionSystem.cpp`). The Soret and correction fluxes
  across the last face are credited to the wall node, so they are conserved too.
- Continuity is integrated from the wall (`ContinuityBoundaryCondition::Wall`, `rV[jj] = 0`),
  so the inlet mass flux `V[0]` follows from the solution.
- Particles cannot cross the wall (no deposition, no thermophoresis).

### The pressure term and the mass flux controller

Cantera solves for the pressure curvature Λ, a single constant that lets both mass-flux
conditions hold. Ember prescribes it through the strain rate:

```
Λ = −ρu·a²/β²      ⇔      a = β·sqrt(−Λ/ρu)     (β = 2 disc, 1 planar)
```

So a Cantera solution maps onto ember exactly through `a`. With `ImpingingJet(massFlux=...)`,
a PI controller on log(a) (`FlameSolver::updateMassFluxControl`) adjusts `a` until `V[0]`
reaches the target, starting from `StrainParameters.initial`. Only the converged state is
physical, not the path to it. `a` and `mdot` are written to `out.h5`.

The cold-flow estimate `a ≈ ṁ/(ρu·L)` is far too low once there is a flame (3–4× here). Fastest
is to take `a` from a Cantera run, as the Igor case scripts do.

### Prescribed temperature profile

`ImpingingJet(temperatureProfile=(x, T))` matches Cantera's `set_fixed_temp_profile` with
`energy_enabled = False`. T is reset to the interpolated profile every step, and dT/dt is zero
in every operator (convection, diffusion, source terms, cross terms).

### New termination measure: `'moments'`

With a prescribed T, `TerminationCondition(measurement='dTdt')` stops at `tMin` (dT/dt is zero
everywhere), and the heat release (`'Q'`) barely depends on the particles. The new
`measurement='moments'` stops once every particle moment changes by less than `momentsTol`
[1/s], as an RMS over the grid relative to the moment's peak value.

## Bugs found and fixed in the burner inlet (`InletFlux`)

The fixed-temperature Igor case puts the flame on the burner, so H2 diffuses back through the
inlet face strongly. Ember first lost about 2.5% of the hydrogen atoms there, which showed up
as a 0.6% error in the strain rate. An element balance (inflow `ṁ·Z_in` against outflow
`∫ β·ρ·U·Z dx`) located the loss at the inlet node. Two causes:

1. **Mean molecular weight at the inlet node (main cause).** In the convection system,
   `Wmx[0]` had no inflow term for `InletFlux` (the `ControlVolume` and `WallFlux` conditions
   have one). With strong diffusion through the face, `Wmx[0]`, and so `ρ[0]`, drifted within
   each step. The species then saw an inlet velocity about 14% too low, and part of the inflow
   was the depleted node-0 mixture instead of the supplied stream. Fixed in
   `ConvectionSystemUTW::f`.
2. **Cross fluxes at the first face.** The correction and Soret fluxes leaving node 1 toward
   the inlet were never credited to node 0. Fixed in `FlameSolver::updateCrossTerms`.

Both fixes also apply to burner-stabilized flames. Rerunning `burnerFla_test.py` with the
old and the new code (same Cantera start) shows a small effect:

| | Cantera BurnerFlame | ember, old | ember, new |
|---|---|---|---|
| T_max [K] | 2178.5 | 2190.0 | 2194.4 |
| Flame position (T = 1000 K) | 4.069 mm | 4.116 mm | 4.093 mm |
| H at outlet / inlet | 0.9967 | 0.9937 | 0.9985 |
| H element fraction at 10 mm | 0.9519 | 0.9426 | 0.9486 |

The new code is closer to Cantera in hydrogen; earlier burner results shift by under 1%. The
new reference run is in `run/burnerFla_canteraInit`.

## Pitfalls

- **Use `splittingMethod='balanced'` (the default).** With `'strang'`, copied from
  `single_Igor.py`, the Igor jet converged to a strain rate 7.6% above Cantera's. Balanced
  splitting matched it to 0.02%. Existing `single_Igor` runs use `'strang'`.
- **Starting from given profiles (`haveProfiles=True`).** Ember takes the inlet stream
  composition from the first point, `Y[:, 0]`. A converged Cantera or ember inlet node is
  depleted in H2, so overwrite that point with the supplied mixture (all new scripts do this;
  `burnerFla_particles*.py` already did).
- **Timestep.** The thin plug-flow layer at the inlet limits the split timestep. The default
  2e-5 s gave NaNs at a ≈ 900 1/s and 1 atm (1e-5 was stable); 2e-5 is fine at 3000 Pa.

## Validation

**Lean H2/air at 1 atm**, φ = 0.7, L = 1 cm, 400 K wall, 3 m/s inlet, disc
(`example_impingingJet.py`):

| | Cantera | ember, prescribed `a` | ember, ṁ controller | ember, finer grid |
|---|---|---|---|---|
| Inlet ṁ [kg/m²s] | 2.7727 | 2.7668 | 2.7727 (target) | 2.7646 |
| `a` [1/s] | 922.6 | 922.6 (input) | 924.0 | 922.6 (input) |
| T_max [K] | 2016 | 1997 | 1997 | 2016 |
| Flame position | 5.016 mm | 5.075 mm | 5.081 mm | 5.035 mm |

The remaining gap closes with grid refinement. The planar geometry runs and its controller
converges (no Cantera counterpart).

**Igor's burner, lean 500 ppm FEC5O5, 3000 Pa** (36 mm nozzle, L = 44.4 mm, burner 426.8 K,
plate 1146.2 K, from `T_of_x-L500.csv`):

| | Cantera | ember |
|---|---|---|
| Energy solved: `a` [1/s] | 113.62 | 113.65 |
| Energy solved: T_max [K] | 1947.7 | 1955.5 |
| Measured T profile: `a` [1/s] | 105.79 | 105.53 |
| Measured T profile: H outflow / inflow | 0.9984 | 0.9985 |

All velocity, temperature and species profiles overlay Cantera's.

## Igor case with particles

`single_Igor_impingingJet_fixedT_particles.py` adds the full particle model of
`burnerFla_particles_oxid.py`:

- 3 moments (N, Fe mass, O mass);
- Fe-cluster collision nucleation and coagulation;
- OpenFOAM surface reactions (oxidation, reduction, etching);
- condensation of Fe, FE2–FE7 and FeO.

It uses the measured temperature profile and terminates on `'moments'` (converged by about
10 ms, stopped at 20 ms).

- Particles form 1–5 mm from the burner (peak N ≈ 8×10¹⁸ m⁻³) and are completely gone by
  about 8 mm.
- They stay at 0.55–0.68 nm, close to `minParticleDiameter` (0.5 nm), so they are etched
  and disintegrated almost as fast as they grow. O/Fe is 0.8–0.95.
- Compared with Cantera's gas-only solution, O, H and OH are about 25%, 28% and 10% lower
  across the whole domain. This is probably radical consumption by the surface reactions
  on the large particle surface near the burner; it is worth checking the surface-reaction
  rates in this regime.

## Files

| File | Content |
|---|---|
| `python/ember/examples/example_impingingJet.py` | H2/air validation against Cantera |
| `python/ember/examples/single_Igor_impingingJet.py` | Igor burner, energy equation solved |
| `python/ember/examples/single_Igor_impingingJet_fixedT.py` | Igor burner, Cantera first, measured T profile (switch `fixedTemperature`) |
| `python/ember/examples/single_Igor_impingingJet_fixedT_particles.py` | as above, with the full particle model |
| `python/ember/examples/T_of_x-L500.csv` | measured temperature profile, lean 500 ppm case |

## Open points

- The `ControlVolume` left boundary (twin and cylindrical flames) has the same kind of
  cross-flux gap at its first face as the one fixed for `InletFlux`; it was not changed.
- The C++ test `DiffusionSystemTest.CylindricalCoordinates` fails; it already failed before
  this work. All other tests pass.
- No particle deposition or thermophoresis at the wall.
