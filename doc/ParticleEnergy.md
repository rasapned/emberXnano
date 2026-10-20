# Particle temperature (particle enthalpy moment)

Before this work, the particles in ember's moments model were at the gas temperature by
construction. The heat released when gas molecules condense onto them, or react on their
surface, was not accounted for anywhere. The gas energy equation does not see particle
processes either.

The particle enthalpy is now transported as an extra moment, and the particle temperature T_p
follows from it. Particles can therefore be hotter or colder than the gas, and the
surface-reaction rates respond to T_p.

## Usage

```python
Particles(nMoments=3,              # N, metal mass, O mass
          energy=True,             # + particle enthalpy -> particle temperature
          phases=[(0.0, 7874.0, ['Fe(a)', 'Fe(c)', 'Fe(d)', 'Fe(L)']),
                  (1.0, 5745.0, ['FeO(s)', 'FeO(L)']),
                  (4/3, 5170.0, ['Fe3O4(s)']),
                  (1.5, 5240.0, ['Fe2O3(s)'])],
          thermoFile='nasa_condensed.yaml',  # condensed-phase thermo (default)
          thermalAccommodation=0.3,          # gas-particle conduction
          radiationAbsorption=0.3,           # E(m); 0 switches radiation off
          radiationTsurr=300.0,              # surroundings [K]
          surfaceReactions=[SurfaceReaction(reactant='O', product=None,
                                            A=..., Ta=..., deltaH=-578000.0), ...],
          ...)
```

- **`energy=False` (default):** there is no enthalpy moment, T_p equals T_gas, and the code
  paths are unchanged.
- **`energy=True`:** requires `nMoments=3`. Internally the enthalpy is one more moment row
  (index `kH = 3`).
- **Output:** `particleTemperature` is written to the profile files. Where there are no
  particles it equals the gas temperature.
- **Integrator:** use `chemistryIntegrator='cvode'`. The particle temperature relaxes on
  microsecond time scales, which the explicit QSS integrator handles poorly.

### Surface reaction enthalpies

`SurfaceReaction.deltaH` is the reaction enthalpy in **J/mol** of events for
`reactant(g) + particle → product(g) + particle'`.

- **Sign:** negative = exothermic; the heat −deltaH goes into the particle.
- **Composition change:** the particle's own change (for example Fe → FeO) is booked at T_p
  from the condensed-phase thermo. Only −deltaH changes the particle temperature.
- **OpenFOAM values:** these are heat released into the particle, so enter them as
  `deltaH = −value`. The Igor case script prints the 298 K standard reaction enthalpies next to
  the values in use, for comparison.

## Model

### State variable

H is the **absolute** particle enthalpy [J per kg of gas], on the same reference as the
gas-phase NASA thermo (formation enthalpies included).

- **Why enthalpy and not temperature:** enthalpy is additive, so convection, diffusion,
  mixing and coagulation need no extra terms. It is transported exactly like the other
  moments, with the particle diffusivity.
- **Boundaries:** H = 0 at the inlet (no particles); zero flux at an impinging-jet wall.

### Particle thermo and temperature

- Each phase lists its condensed species; for iron these are Cantera's `nasa_condensed.yaml`
  entries for Fe(s, l), FeO(s, l), Fe₃O₄ and Fe₂O₃.
- **Tables:** at each temperature (200–5000 K, 2 K steps) the species with the lowest Gibbs
  energy among those whose data range covers it is used. Polymorphs and melting therefore
  appear as enthalpy jumps. Outside all data ranges, h is extrapolated linearly with c_p.
- **Composition:** the enthalpy per kmol of metal is interpolated in the O/metal ratio between
  the bracketing phases (lever rule, as for the density).
- **Inversion:** T_p is obtained by bisection on the monotone table:
  `particleTemperature()` in `src/particleUtils.cpp`. The tables are built in Python by
  `_phaseEnthalpyTable()` in `python/ember/_ember.pyx`.
- **Consistency check:** the mechanisms' gas-phase iron thermo agrees with the condensed data.
  Sublimation heats are 415.5 kJ/mol for Fe (literature ≈ 416) and 523 kJ/mol for FeO.

### Sources of particle enthalpy (`src/sourceSystem.cpp`)

| Process | Contribution to dH/dt |
|---|---|
| Nucleation (collision channels, CNT) | + enthalpy of the consumed gas molecules at T_gas, `Σ ν_k·h_k(T_gas)·J / ρ` |
| Condensation | + `R·h_k(T_gas) / ρ` for each condensing species |
| Surface reactions | + `R·(ΔH_comp − deltaH) / ρ`, with `ΔH_comp = (h − x·h′)·ΔM + h′·ΔO` (composition change at T_p, h and h′ = ∂h/∂x per metal atom); the Arrhenius factor uses **T_p** |
| Conduction to the gas | − `N·N_A·q`, with q from the free-molecular limit `α·πd²·(p/8)·c̄·(γ+1)/(γ−1)·(T_p/T_g − 1)` combined harmonically with the continuum limit `2πd·λ_g·(T_p − T_g)` |
| Radiation (Rayleigh limit) | − `N·N_A·P`, with `P = 199.1·π³·d³·E(m)·k_B⁵·(T_p⁵ − T_surr⁵) / (h⁴c³)` |
| Coagulation, disintegration | none (H is additive; removed mass is already in the reactions) |

The latent heat therefore follows directly from the mechanism. When gas molecules condense,
the particle gains their full gas-phase enthalpy, and its own (much lower) condensed-phase
enthalpy at that composition shows the excess as a higher T_p.

Collision rates for nucleation, condensation and coagulation stay at T_gas.

## Implementation

| Where | What |
|---|---|
| `src/readConfig.h` | `kH`, phase enthalpy tables, accommodation, radiation inputs, `surfDeltaH` |
| `src/particleUtils.cpp/.h` | `particleEnthalpyPerMetal()`, `particleTemperature()` |
| `src/sourceSystem.cpp/.h` | T_p per right-hand-side call, enthalpy sources, `computeParticleHeatTransfer()` |
| `src/flameSolver.cpp/.h` | `particleTemperature` field for output |
| `src/convectionSystem.cpp`, `src/sourceSystem.cpp` | absolute tolerances of the H row (scaled ×10⁷); QSS lower bound disabled for H, which can be negative |
| `python/ember/input.py` | `Particles.energy`, `thermoFile`, `thermalAccommodation`, `radiationAbsorption`, `radiationTsurr`, phase species, `SurfaceReaction.deltaH` |
| `python/ember/_ember.pyx/.pxd` | enthalpy tables, option transfer, `particleTemperature` property |
| `python/ember/output.py` | writes `particleTemperature` |

## Verification

Case: Igor's burner as an impinging jet with the measured gas temperature profile
(`single_Igor_impingingJet_fixedT_particlesEnergy.py`).

| Test | Expected | Result |
|---|---|---|
| `energy=False` (regression) | identical to the particle case without energy | moments and strain rate unchanged (rel. diff. ≤ 4·10⁻⁷) |
| Limit test: α = 1000, E(m) = 0 | T_p → T_gas | T_p − T_gas within −0.1 … +0.9 K; N and particle Fe match the energy-off run to 5 digits |
| deltaH = 0 (only latent heat) | particles slightly hotter than the gas | T_p − T_gas = +2 … +67 K, largest where condensation is strongest |
| 298 K standard deltaH (placeholders) | net endothermic in the H-rich zone | T_p − T_gas = −72 … −139 K: reduction by H and etching by OH dominate |

**Why particles run hotter than the gas without reaction heat.** Relaxation is fast
(τ = m·c_p/(hA) ≈ 1 µs for a 1 nm particle at 3000 Pa), but the particle settles at the
temperature where the heat input balances the loss.

- Loss: each gas collision removes about α·k·ΔT·(γ+1)/(2(γ−1)) ≈ 0.6·k·ΔT.
- Gain: a fraction x of the collisions are Fe-bearing molecules that stick, each delivering
  L ≈ 4–5 eV.
- Balance: ΔT ≈ x·L / (0.6·k), so x = 10⁻⁴ … 5·10⁻⁴ gives about 9 … 45 K, in line with the
  runs.

Radiation is negligible at these sizes (≈ 2 K equivalent at 1 nm and 2000 K), and conduction
carries the heat back to the gas (up to ≈ 0.3 MW/m³).

The results with the OpenFOAM reaction enthalpies (all reactions heating the particle) are
not part of this note; that run is set up in the case script.

## Pitfalls

- **Plots:** the particle mass is `moments[1:3]`, not `moments[1:]`, which would add H. Size
  and O/metal ratio are 0/0 noise where there are no particles; the case scripts mask them
  where the particle mass is below 10⁻³ of its peak.
- **Termination:** with a prescribed gas temperature use
  `TerminationCondition(measurement='moments')`. The `dTdt` check stops at `tMin`.
- **Collision cross-sections:** condensation uses π(d_p + d_k)²/4, conduction only πd_p². For
  sub-nm particles this overstates the heat input relative to the loss by about 2–3×. Using
  (d_p + d_gas)² in the conduction term would be more consistent.
- **Surface reaction rates:** the OpenFOAM pre-factors exceed the gas-particle collision rate
  (sticking probability above 1), which inflates every reaction heat.

## Not included yet

- **Evaporation/sublimation** of Fe and FeO from hot particles (Kelvin-enhanced vapour
  pressure). It is likely the main limit on overheating, and it is also missing as a mass
  process.
- **Coalescence heat:** surface-energy release when particles merge, which is significant for
  particles of a few nm.
- **Two-way coupling** to the gas energy equation. It has no effect when the gas temperature
  is prescribed.
- **A heat-partition factor** between particle and desorbing product.
- **Particle-size effects on the thermo:** bulk data is used even for sub-nm clusters.
