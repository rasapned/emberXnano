# Particle temperature in the impinging-jet case: what was checked (2026-10-22)

Case: `single_Igor_impingingJet_fixedT_particlesEnergy.py`, run with the measured gas
temperature profile (T_max 1701 K at 21 mm). Mechanism `opt-compact-mech-OF2_OFkinetics_FeTransp.yaml`,
particle energy on, branch `feat/particle-thermo-evaporation`.

Variant runs `run/jetSens_*` (deleted 2026-11-17, superseded by the jet grid below). Details and numbers:
`python/ember/examples/caseLog/CASE_LOG.md`.

## What was checked

1. **OpenFOAM reaction heats as prescribed.** They gave T_p ≈ 5000 K. Several have the wrong
   sign compared with the gas and condensed-phase thermochemistry. Reduction by H2 (+21 kJ/mol)
   and H (+83), and OH etching (+61), are endothermic; OpenFOAM sets them exothermic.
   Oxidation–reduction cycles then accumulate heat in the particle without limit.
2. **Heat from the thermochemistry** (committed in 97eaf4f). T_p is bounded, but the particles
   are 20–140 K *colder* than the gas: the endothermic reduction and etching reactions plus
   radiation outweigh the heating. This held in every sensitivity variant: evaporation, OH
   etching ×0 / ×0.1, coagulation ×3, nucleation efficiency 0.1 / 0.3, extra Fe carriers.
3. **Each OpenFOAM heat on its own** (added as an extra on top of the thermochemistry, evaporation on):
   - H→OH alone causes the runaway to the 5000 K table limit and crashes the run.
   - O2→O, H2O→H2 and H2→H2O make the particles hotter than the gas: peaks of about +205,
     +90 and +125 K respectively (O2 and H2 near 2 mm, H2O over 5–10 mm).
   - O→particle and OH→FeOH have no visible effect.
4. **Extra heat on H2O→H2 tuned for the measured 400–800 K overheating:**
   - With full OH etching: particles evaporate by ~4 mm (−300 kJ/mol extra), or the run runs
     away (−500).
   - With OH etching ×0.1 and −175 kJ/mol extra: T_p − T = 375–470 K from 3 to 10 mm,
     particles gone by 11.4 mm. That extra is about 9× the thermochemical H2O→H2 heat
     (−21 kJ/mol).
   - The switch is sharp: −150 gives +320–380 K and particles survive to 16 mm; −200 gives
     +540–560 K and particles are gone by 7 mm.
5. **OpenFOAM heats kept, OH etching ×0.1, H→OH heat scanned.**
   - With all six OpenFOAM heats, the H→OH heat is not the lever: even at its thermochemical
     value (+83 kJ/mol) the particles run ~830 K above the gas at 1 mm and evaporate by 6 mm,
     and the runs do not converge.
   - The cause is the H2O→H2 / H2→H2O pair. In OpenFOAM both are exothermic (−87 and
     −52 kJ/mol), so every oxidation–reduction cycle leaves 139 kJ/mol in the particle, where
     it should leave none.
   - With the pair made heat-neutral (H2→H2O = +87 kJ/mol) and the other OpenFOAM heats kept,
     H→OH controls the overheating smoothly:

     | H→OH total heat | T_p − T at 3/5/8/10/13/15 mm [K] | Particles | d_p, 3 → 21 mm | Gas Fe / no particles at 3/8/13/21 mm | Run |
     |---|---|---|---|---|---|
     | +83 (thermochemistry) | +167/+192/+206/+191/+151/+124 | whole domain | 1.0 → 1.7 nm | 0.27/0.17/0.30/0.40 | converged |
     | **+60** | **+264/+296/+302/+276/+212/+172** | **whole domain** | 1.0 → ~1.6 nm | 0.27/0.21/0.41/0.54 | converged |
     | +40 | +346/+382/+376/+336/+245/+191 | partly evaporate from ~10 mm | ~1–1.4 nm | 0.28/0.29/0.61/0.80 | not converged (oscillates to 0.5 s) |
     | 0 | +503/+536 (3/5 mm) | gone by 7.7 mm | ≤ 1.2 nm | 0.30/0.88/0.99/1.00 | converged |

6. **Gas-phase Fe species with and without particles** (mole-fraction ratio to the
   no-particle run, at 3/5/10/15/20/30/40 mm):

   | Case | Fe | FeO | FeOH | FeO2 |
   |---|---|---|---|---|
   | thermo heats, OH ×1 | 0.94/0.92/0.78/0.74/0.67/0.68/0.78 | 0.65/0.99/0.98/0.92/0.82/0.77/0.81 | 0.77/**1.40**/**1.31**/**1.18**/0.98/0.85/0.86 | 0.52/0.87/0.99/1.07/1.02/0.93/0.91 |
   | thermo heats, OH ×0.1 | 0.68/0.30/**0.18**/0.27/0.45/**1.18**/**1.83** | 0.43/0.29/**0.18**/0.24/0.33/0.66/0.97 | 0.44/0.34/0.20/0.21/0.25/0.39/0.41 | 0.37/0.27/0.18/0.22/0.25/0.31/0.31 |
   | OF heats\*, H→OH +83 | 0.67/0.34/0.23/0.37/0.57/**1.08**/**1.49** | 0.37/0.31/0.25/0.36/0.48/0.73/1.00 | 0.38/0.39/0.29/0.35/0.41/0.51/0.56 | 0.30/0.28/0.25/0.34/0.40/0.42/0.41 |
   | OF heats\*, H→OH +60 | 0.70/0.36/0.29/0.48/0.64/0.94/**1.21** | 0.36/0.33/0.33/0.49/0.60/0.76/1.00 | 0.37/0.43/0.40/0.51/0.57/0.62/0.70 | 0.29/0.30/0.33/0.49/0.56/0.55/0.55 |

   \* OpenFOAM heats with the H2O/H2 pair made heat-neutral, OH etching ×0.1.
   - With full OH etching, Fe and FeO barely drop, and FeOH even rises (etching produces it).
   - With OH ×0.1, Fe and FeO drop to 0.18–0.33 around 10 mm, close to the measured 1/6–1/5,
     but recover towards the flame. Beyond ~28 mm atomic Fe is *above* the no-particle level
     (up to 1.8×), i.e. the particles release Fe near the plate. The hotter the particles, the
     smaller the drop.

## Reaction enthalpies: thermochemistry vs. OpenFOAM vs. tested [kJ/mol]

Negative = heat into the particle. Thermochemistry: mechanism gas thermo + `nasa_condensed.yaml`
for the particle, Fe → FeO stage. In the code the thermochemical heat is always included at
the actual temperatures; the tested values were imposed as an extra
deltaH = target − thermochemistry(1000 K), so they hold exactly at 1000 K only.

| Surface reaction | Thermo 300 K | Thermo 1000 K | OpenFOAM | Item 4 (extra on H2O→H2) | Item 5 (OF heats\*) |
|---|---|---|---|---|---|
| O2 → O + particle-O | −23 | −16 | −80 | −16 (thermo) | −80 |
| O → particle-O | −521 | −522 | −578 | −522 (thermo) | −578 |
| H2O → H2 + particle-O | −30 | −21 | −87 | −171 / −196 / −221 / −271 / −321 | −87 |
| H2 + particle-O → H2O | +30 | +21 | −52 | +21 (thermo) | **+87** (pair heat-neutral) |
| H + particle-O → OH | +91 | +83 | −53 | +83 (thermo) | **0 / +40 / +60 / +83** |
| OH + particle-Fe → FeOH | +68 | +61 | −157 | +61 (thermo) | −157 (rate ×0.1) |

Item 4 totals correspond to extras of −150 / −175 / −200 / −250 / −300 kJ/mol on top of the
thermochemical −21.

Figures (script `python/ember/examples/caseLog/plot_particleTemperature.py`):
- [figures/particleTemperature_OFheats.png](figures/particleTemperature_OFheats.png): heat from
  the thermochemistry only vs. OpenFOAM heats one at a time (item 3).
- [figures/particleTemperature_extraHeat.png](figures/particleTemperature_extraHeat.png):
  OH etching ×0.1 and the extra H2O→H2 heat scan (item 4). Grey bands are the measured ranges.
  The small steps and the plateau in T_p near 1667 K are the Fe γ→δ transition in the
  condensed-phase tables.
- [figures/particleTemperature_HOHscan.png](figures/particleTemperature_HOHscan.png): OpenFOAM
  heats with the H2O/H2 pair heat-neutral, H→OH scanned (item 5), whole domain. The steep d_p
  rise beyond ~40 mm and the jump at 43 mm in the thermo-heat case are wall artefacts.
- [figures/gasFeSpecies.png](figures/gasFeSpecies.png): Fe, FeO, FeOH, FeO2 mole fractions and
  their ratio to the no-particle run (item 6).

## Main conclusions

- **Missing heat source.** The measured overheating can only be reached with heats that do not
  conserve energy. If the measurement is right, the model lacks a heat source. One candidate
  is H- or O-atom recombination on the particle surface (H + H → H2 releases 436 kJ/mol);
  this has not been checked against the literature.
- **OpenFOAM heats.** The 5000 K comes from two things: the H→OH heat (exothermic in OpenFOAM,
  endothermic by thermochemistry) and the exothermic H2O→H2 / H2→H2O pair, which heats the
  particle on every oxidation–reduction cycle. With the pair made heat-neutral, the remaining
  OpenFOAM heats with H→OH at +60 kJ/mol give a converged solution with particles ~300 K
  above the gas over 3–10 mm, surviving the whole domain and growing to ~1.6 nm.
- **Evaporation ceiling.** Once ~1 nm Fe/FeO particles are 400 K or more above the gas, they
  pass about 1500–1650 K and evaporate within a few mm (Kelvin effect). This held for both
  ways of adding heat (items 4 and 5); near the ceiling the runs oscillate instead of
  converging.
- **Conflict with the other measurements.** In the model, hot particles evaporate early, then
  gas Fe returns to the no-particle level and the particles stay small (1.1–1.3 nm). This
  contradicts the measured 5–6× lower gas Fe over the whole profile and the 2–3 nm sizes.
  Matching all three would need weaker evaporation:
  - an evaporation coefficient below 1,
  - a lower vapour pressure over the higher oxides (Fe3O4/Fe2O3),
  - a smaller Kelvin effect at nm sizes.

  These are the open points in `doc/EvaporationModel.md`.
- **Gas-phase Fe** responds mainly to OH etching, not to the particle temperature. With
  OH ×0.1, Fe and FeO drop to ~0.2 of the no-particle level around 10 mm, close to the
  measurement, but not over the whole profile: they recover towards the flame, and near the
  plate atomic Fe even exceeds the no-particle level. Hotter particles reduce the drop.

---

# Sensitivity grid with prescribed (OpenFOAM) surface heats (2026-11-12)

Code: commit `7daa663` (`Particles.surfaceReactionHeat = 'prescribed'`: the only reaction heat
is −ΔH, as OpenFOAM; sticking molecules bring (1 − α) of their sensible heat, so nothing is
counted twice with conduction) plus `TerminationCondition.momentsFloor = 1e-3` (commit `4fa93b4`).
Same jet case, particle energy and evaporation on. 133 runs `run/jetGrid_2026-11/jetGrid_<variant>`;
driver, variant lists, logs, the full metrics table and the plotting scripts in
`python/ember/examples/run/jetGrid_2026-11/`.

Targets: T_p − T ≈ 500–800 K, particles gone at 8–13 mm, gas Fe below the no-particle run.
Metrics: median T_p − T over 2–8 mm; x where the particle mass falls below 1 % of its peak;
F = gas Fe / no-particle gas Fe (min and mean up to that point).

## Results

1. **Reference (all OpenFOAM heats and rates):** 614 K, particles gone at 6.3 mm, F = 1.04
   (gas Fe not depleted), d_p 0.54 nm. T_p is pinned at ~1530–1560 K, the evaporation limit of
   sub-nm Fe/FeO particles. O2→O, H2O→H2, H2→H2O and H→OH each supply more heat than
   conduction removes; the surplus evaporates the particles.
2. **One change at a time (24 runs: each heat → 0 or its thermochemical value, each rate ×0.1 or
   ×3): nothing moves the result** (520–660 K, 5.6–6.9 mm, F ≈ 1.04). This reproduces October:
   with the H2O/H2 pair at the OpenFOAM heats, the H→OH heat has no effect.
3. **The H2O→H2 / H2→H2O pair is the gatekeeper.** Both are exothermic in OpenFOAM (−87, −52),
   so each oxidation–reduction cycle leaves 139 kJ/mol in the particle. Made heat-neutral
   (H2→H2O = −ΔH(H2O→H2)), the particles grow and deplete gas Fe. Without it no other change
   helps; on its own it is not enough (565 K, 4.6 mm, F 0.97).
4. **OH etching ×0.1** matters only once the particles survive; then it controls the Fe drop.
5. **Trade-off:** more heat → more overheating but earlier extinction. With the pair neutral and
   OH ×0.1, scaling all heats by f: f = 0.5 → 440 K, 14.2 mm; 0.55 → 477 K, 10.8 mm;
   0.6 → 509 K, 8.3 mm; 0.65 → 541 K, 6.9 mm; 0.7 → 574 K, 6.0 mm.
6. **Oxidation state shifts the trade-off.** Faster oxidation (O2→O, O, H2O→H2 rates) or slower
   reduction by H raises O/Fe above 1. In the lever-rule evaporation model the pure-Fe phase then
   has zero activity, so only FeO evaporates and the particles survive to T_p ≈ 1780–1880 K
   (instead of ~1550 K): more overheating *and* later extinction. The particles stay ≤ 1.3 nm.

## Candidates meeting all three targets

All keep non-conserving prescribed heats.

| | Changes from OpenFOAM | T_p − T | Gone at | F min / mean |
|---|---|---|---|---|
| (a) | pair heat-neutral, OH rate ×0.1, H→OH rate ×0.55 | 607 K | 8.7 mm | 0.28 / 0.62 |
| (b) | all heats ×0.6, pair heat-neutral, OH rate ×0.1 | 509 K | 8.3 mm | 0.28 / 0.64 |
| (c) | heats ×0.8, pair heat-neutral, OH ×0.1, O2→O rate ×2.5, H2O→H2 rate ×3 | 684 K | 8.6 mm | 0.25 / 0.62 |
| (c), heats ×0.78 | as (c), heats ×0.78 | 670 K | 9.4 mm | 0.24 / 0.62 |

- (a) changes the fewest parameters but sits on a knife-edge: H→OH ×0.5 gives 587 K / 13.2 mm,
  ×0.6 gives 618 K / 6.7 mm.
- (c) reaches the highest overheating, but O2→O ×2.5 is near the collision limit. The unmodified
  OpenFOAM rates are already close to it: reactions per collision at full site coverage are
  1.1–1.7 for O, 0.6–0.8 for O2, 0.9 for H, at 800–2000 K.
- Above ~700 K overheating no setting keeps the particles beyond ~7 mm.
- (c) and most other runs with raised O2/H2O/O rates exceed the collision limit; see the next
  section for the grid redone within it.

Figures: `python/ember/examples/run/jetGrid_2026-11/jetGrid_summary.png` (all runs: overheating vs extinction,
coloured by F), `jetGrid_best.png` (profiles of the candidates), `jetGrid_phase3_best.png` (same folder).

## Within the collision limit (2026-11-16)

The uncapped grid used surface rates above the gas–particle collision rate: per collision on an
available site (hard-sphere kernel, as condensation), O→p reaches 1.3 in (a)–(c), O2→O and
H2O→H2 1.6 in (c). Commit `3cb1c1c` adds `Particles.surfaceCollisionLimit`, which caps each rate
coefficient at that collision rate. 25 more runs: ref and (a)–(c) capped, a re-tune of all
three routes (phase 9), and the Kelvin nucleation limit on ref and (a)–(c) (no effect:
nucleation happens at 0–3 mm, where Fe is far supersaturated).

| | Changes from OpenFOAM, collision cap on | T_p − T | Gone at | F min / mean |
|---|---|---|---|---|
| ref | — | 614 K | 6.3 mm | 0.91 / 1.04 |
| (a) | pair heat-neutral, OH rate ×0.1, H→OH rate ×0.55 | 607 K | 8.6 mm | 0.28 / 0.62 |
| (b) | all heats ×0.6, pair heat-neutral, OH rate ×0.1 | 510 K | 8.3 mm | 0.28 / 0.64 |
| (c), capped | as (c) | 663 K | 6.7 mm | 0.31 / 0.65 |
| (c′) | heats ×0.7, pair heat-neutral, OH ×0.1, O2→O and H2O→H2 rates ×1.5 | 604 K | 8.6 mm | 0.26 / 0.61 |

1. **Below the cap nothing changes.** Heats, H→OH and OH rates give the same results capped and
   uncapped (≤ 3 K, ≤ 0.1 mm). O→p is capped over ~58 % of the particle zone in (a) and (b)
   without effect: it is not rate-limiting.
2. **Oxidation boosts saturate.** At OpenFOAM rates O2→O and H2O→H2 already use ~half the
   collision rate, so at most ~×2 is effective. Extinction vs O2/H2O factor ×1 / ×1.5 /
   ×2.5–3 / ×5: 6.0 / 8.6 / 9.6 / 9.6 mm at heats ×0.7; – / 7.0 / 7.9 / 7.9 mm at ×0.75;
   – / 5.7 / 6.7 / 6.7 mm at ×0.8. Route (c) at heats ≥ ×0.75 is closed; (c′) replaces it,
   with O2 and H2O below the cap and only O→p capped, as in (a) and (b).
3. **9 of the 15 uncapped all-target runs relied on boosts above the cap**; most of them also
   regrew particles near the wall (4–23 % of the peak mass). (a), (b) and its variants stand.
4. **The window is narrow.** Heats: 500 K needs f ≥ 0.6, 8 mm needs f ≤ 0.6. H→OH rate:
   ×0.65 / 0.6 / 0.55 / 0.5 → 5.5 / 6.7 / 8.6 / 11.8 mm (×0.45 never extinguishes), i.e. 10 % in
   the rate moves extinction by 2–3 mm. Both (a) and (c′) end near 605 K / 8.6 mm.
5. **Only 8–9 mm extinction is clean.** Runs surviving to ~9.5–12 mm regrow particles near the
   wall (H→OH ×0.5: 4.9 %, heats ×0.5: 5.7 %, heats ×0.7 with ×2.5/×3 boosts: 4.3 %).
6. F (mean 0.60–0.72) is set by the OH etching rate and is not the limiting target.

Figure: `python/ember/examples/run/jetGrid_2026-11/jetGrid_candidates_capped.png`; all 158 runs
with a regrowth column in `grid_table.json` (same folder).

## Caveats

- The zero Fe activity above O/Fe = 1 (lever rule between FeO and Fe3O4) is a model choice that
  the candidates rely on; real Fe activity in an oxidised particle is small but not zero.
- Downstream of the evaporation front N and the O moment carry residual negative noise (N up to
  several times the small particle peak, mass ~1e-3 of its peak); gas Fe is not affected.
- The gas F stays depleted only up to the extinction; beyond it gas Fe returns to the
  no-particle level, as in October.
