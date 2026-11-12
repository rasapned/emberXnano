# Physics & numerics changelog

One entry per implementation commit, newest first, written by the `physics-changelog` skill
(`.claude/skills/physics-changelog/SKILL.md`). The full current model is in
[`PhysicsAndNumerics.md`](PhysicsAndNumerics.md). Commits up to `5837ad2` (2026-10-13) are
summarised in its §7.

## 2026-11-12 — feat: surface-reaction heat modes; sticking molecules' sensible heat split with conduction (staged)

**Area:** particles: energy (surface reactions, condensation, conduction)

**What changed (physics):** (1) `Particles.surfaceReactionHeat` selects one of two exclusive
reaction-heat models: `'thermo'` (default, from the thermochemistry; `deltaH` must be 0) or
`'prescribed'` (only $-\Delta H_r$ = `SurfaceReaction.deltaH`, as OpenFOAM). (2) Conduction counts
every gas collision with the accommodation $\alpha_T$; a molecule that sticks (reacts or condenses)
fully accommodates, so the source terms now add only the remaining $(1-\alpha_T)$ of its sensible
heat. Before, `'thermo'` counted it fully *and* with $\alpha_T$ in conduction. (3) The $d_{min}$
floor of the surface sources is also applied in the heat transfer.

**Equations** (per event, [J/kmol]; $S_k = h_k(T) - h_k(T_p)$):

$$
\Delta H_p^{thermo} = h_R(T_p) - h_P(T_p) + (1-\alpha_T)S_R
\qquad\text{(was } h_R(T) - h_P(T_p)\text{)}
$$

$$
\Delta H_p^{prescribed} = (h_p - x_O h_p')\,\Delta n_M + h_p'\,\Delta n_O - \Delta H_r + (1-\alpha_T)S_R,
\qquad h_p' = \partial h_p/\partial x_O
$$

$$
\Delta H_p^{cond} = h_k(T_p) + (1-\alpha_T)S_k \qquad\text{(was } h_k(T)\text{)}
$$

$h_p(x_O,T_p)$ [J/kmol metal] the lever-rule particle enthalpy; in `'prescribed'` the composition
term carries the gained/lost atoms at their enthalpy in the particle, so $T_p$ changes only through
$-\Delta H_r$ and $S_R$. $\dot H \mathrel{+}= R\,\Delta H_p/\rho$. Nucleation ($h_k(T)$) and
evaporation ($-h_k(T_p)$) unchanged. Heat transfer below $d_{min}$: $n \to n(d_p/d_{min})^3$,
$d_p \to d_{min}$ (radiation $\propto n d_p^3$ unchanged, conduction $\propto n d_p^2$ reduced).

**Discretisation / numerics:** production operator, per point; no change to the integrator.

**Conservation & units:** `'thermo'` with $\alpha_T = 1$ is the previous energy-conserving form
minus the double count; for $\alpha_T<1$ the sticking molecule's sensible heat is split between
conduction ($\alpha_T$, mixture-averaged) and the source term ($1-\alpha_T$, species $h$).
`'prescribed'` is not energy conserving. Steady solution affected only with `energy=True`.

**Code:** `sourceSystem.cpp::computeSurfaceReactionRates`, `computeCondensationRates`,
`computeParticleHeatTransfer`; `_ember.pyx` (`surfPrescribedHeat`, error for `deltaH` ≠ 0 in
`'thermo'`); option `Particles.surfaceReactionHeat`.

**Verification:** 15 Python tests pass. Impinging jet without particle energy bit-identical.
Energy jet, no evaporation: `'thermo'` $T_p - T$ −129/−97/−107/−88 K at 1/3/5/10 mm (was
−123/−88/−97/−81); `'prescribed'`, $\Delta H = 0$: +6/+8/+8/+14 K (was +8/+11/+12/+19); gas Fe and
$d_p$ within 0.01. `'thermo'` before this change reproduced bit-identically by the mode switch.

**Limits / open points:** conduction is mixture-averaged ($\bar c$, $\gamma$ at $T$), the
source-term part is species-specific, so the split is exact only for a single-species gas; the
correction uses the free-molecular $\alpha_T$ also when the continuum limit dominates (not the case
at nm sizes). $\alpha_T = 0.3$ has no cited source. OpenFOAM's H + particle-O rate exceeds the
H collision rate above ~1400 K at full O coverage.

## 2026-11-12 — feat: Kelvin-limited collision nucleation (staged)

**Area:** particles: nucleation

**What changed (physics):** optional gate on collision nucleation (`NucleationChannel(kelvinLimit=True)`,
default off). A channel whose nucleus is below the critical size at the local monomer pressure
forms no particles; above it the rate is reduced by the fraction of the monomer flux not
balanced by evaporation from the nucleus. The nucleus is treated as a bulk pure-metal sphere
(capillarity approximation, same $\sigma(T)$ and thermo as evaporation), at the gas temperature.

**Equations:** for channel $c$ with $n_c$ metal atoms per nucleus,

$$
J_c \leftarrow J_c\,\max\!\left(0,\;1-\frac{p_{eq}(n_c,T)}{p_1}\right),\qquad
p_1 = \frac{Y_1\rho}{W_1}R_uT,\qquad
p_{eq} = p_{sat,1}(T)\,\exp\!\left(\frac{4\sigma(T)v_1}{d_n k_BT}\right),
$$

$$
p_{sat,1} = p^\circ\exp\!\left(\frac{g_{0}(T)-g_1^\circ(T)}{R_uT}\right),\quad
v_1 = \frac{W_1}{\rho_{p,0}N_A},\quad d_n = \left(\frac{6n_cv_1}{\pi}\right)^{1/3},\quad
\sigma = \sigma_A - \sigma_BT
$$

$p_1$ [Pa] monomer (species of composition {metal: 1}) partial pressure; $g_0$ [J/kmol] molar
Gibbs energy of the pure-metal phase (`phases[0]` condensed species, lowest $g$ per $T$),
$g_1^\circ$ that of the gas monomer, $p^\circ$ the reference pressure; $\rho_{p,0}$ the phase-0
density; $n_c = \mathrm{round}(\nu_A n_{M,A} + \nu_B n_{M,B})$.

**Discretisation / numerics:** production operator (`computeNucleationRates`), evaluated per
point; $\ln p_{sat,1}$ tabulated on the 2 K `particleThermoT` grid, linear interpolation (shared
with the condensation $p_{sat}$ table). The factor is not differentiated separately; it enters
the CVODE RHS like the rest of $J$.

**Conservation & units:** only $J$ is scaled, so the gas/particle mass and enthalpy transfer of
each event is unchanged and $\sum Y_k + M_M + M_O$ is still conserved. $R_u$, $N_A$ kmol-based.
Affects the steady solution.

**Code:** `sourceSystem.cpp::computeNucleationRates`, `particleUtils.cpp::nucleusEquilibriumPressure`,
`_ember.pyx` (tables, $n_c$, monomer); option `NucleationChannel.kelvinLimit` (needs condensed
species in `Particles.phases`).

**Verification:** $p_{eq}(8, 1300\,\mathrm K) = 1.49$ Pa, same as an independent Python
evaluation. Burner nucl-only and cond cases: gate = 1 below ~1100 K, 0.70 at 1190 K, 0 above
~1250 K; $N$ −0.5 to −1 % downstream, $d_p$ and $T$ unchanged; with evaporation the residual
negative moments halve (min $N$/peak −1.8 % → −0.86 %). Impinging jet (no particle energy,
evaporation): no visible change.

**Limits / open points:** bulk $\sigma$ and density for 8-atom nuclei; gate uses $T_{gas}$, not
$T_p$; only the monomer pressure enters (clusters Fe$_2$–Fe$_7$ are not counted as supersaturation);
the mechanism's own Fe$_n$ clusters are less stable than this model, so the gate is a consistency
limit, not a calibrated rate.

## 2026-11-12 — fix: surface sources of sub-d_min particles evaluated at d_min (staged)

**Area:** particles: surface reactions, condensation / evaporation, disintegration

**What changed (numerics/physics):** below $d_{min}$ (`minParticleDiameter`) particles
disintegrate instead of shrinking, so their surface sources are now evaluated for the same
mass in particles of diameter $d_{min}$. Before, $d_p \to 0$ with $N$ finite kept the
evaporation sink finite (Kelvin factor floored, $n_p$ unchanged) while the particle mass
vanished: mass and $N$ crossed zero in finite time, the $N, M_M \le 0$ guard then switched all
sinks off, and the mean particle mass became $\ll$ one atom ($\dot N/N \sim 10^{17}$ 1/s),
which made CVODE fail at the evaporation front.

**Equations** ($m_p = (M_M+M_O)/(NN_A)$, $n_p = \rho N N_A$): for $d_p < d_{min}$

$$
d_p \to d_{min},\qquad n_p \to n_p\left(\frac{d_p}{d_{min}}\right)^3,\qquad
m_p \to m_p\left(\frac{d_{min}}{d_p}\right)^3
$$

(same $\rho_p(x_O)$, so $n_p m_p$ is kept), used in $\pi d_p^2 n_p$ of the surface reactions
and in $\beta_{kp}(d_p, m_p)\,n_p$ of condensation and evaporation; the Kelvin diameter is this
$d_p$ (replaces $\max(d_p, d_{min})$). Mass sinks are then $\propto M_M+M_O$, and the
disintegration term $\dot N = -N\dot m_{loss}/(M_M+M_O)$ (ramp = 1) becomes
$\dot N = -kN$ with the same $k$ as $\dot m = -k\,(M_M+M_O)$: $N$ and mass decay together.
Continuous at $d_p = d_{min}$; no change for $d_p \ge d_{min}$.

**Conservation & units:** unchanged; mass-specific moments, kmol-based $N_A$.

**Code:** `particleUtils.cpp::floorParticleSize`, called in
`sourceSystem.cpp::computeSurfaceReactionRates` and `computeCondensationRates`. Coagulation
and the particle diffusivity (already floored in `flameSolver.cpp`) are unchanged.

**Verification:** `burnerFla_particles_cond` (no evaporation) unchanged (max $\Delta N/N_{max}$
1.4e-5, $\Delta T$ 4e-4 K). `burnerFla_particles_cond_evap` ($d_{min}$ = 0.4 nm): t = 0.08 s
with 0 CVODE failures (before: abort at 0.031 s after 500); particles evaporate at
15–22 mm (1400–1570 K). Python tests pass; `DiffusionSystemTest.CylindricalCoordinates`
fails as before (diffusion code unchanged).

**Limits / open points:** affects the steady solution only where $d_p < d_{min}$. Residual
negative moments downstream of the evaporation front (N to −3 % of its peak, mass to −0.5 %):
sub-tolerance $N$ (~1e-14 kmol/kg vs abstol 1e-12) flips sign in the source stage, and the
balanced-splitting carry-over of the production rate acts in the transport stages.

## 2026-10-22 — test: example cases on OpenFOAM kinetics mechanism, local OpenFOAM thermo (staged)

No physics/numerics change: burner particle cases and the impinging-jet cases use
`opt-compact-mech-OF2_OFkinetics_FeTransp.yaml`, OH etching A −2.568e26, and read Sutherland
coefficients from the local copy `examples/OF_thermo.compressibleGas`.

## 2026-10-22 — feat: reversible condensation with Kelvin evaporation; coagulation and nucleation multipliers (staged)

**Area:** particles: condensation / evaporation, coagulation, nucleation

**What changed (physics):** With `Particles(evaporation=True)` heterogeneous condensation is
reversible: each condensing species whose O/metal ratio matches a particle phase also
evaporates, at the same collision rate times its saturation concentration over the bulk phase,
the surface activity of that phase (lever-rule fraction) and the Kelvin factor. Species
without a matching phase (e.g. FeO$_2$) only condense. Evaporated mass below $d_{min}$ removes
particles (same disintegration ramp as surface reactions). Two sensitivity multipliers:
`Particles.coagulationEnhancement` ($\beta_c$) and `NucleationChannel.efficiency` ($J$).

**Equations** ($\beta_{kp}$ hard-sphere gas–particle collision kernel at $T$ [m³/s], $n_p$
[1/m³], $c_k$ [kmol/m³], $a_i$ lever-rule weight of phase $i$, $v_{m,k} = W_k/(\rho_i N_A)$
[m³], $\sigma = a_\sigma - b_\sigma T_p$ [J/m²], $p^\circ$ gas reference pressure, $g$ standard
molar Gibbs energies [J/kmol], $\Delta_M$ metal atoms in $k$, $g_i$ per kmol metal):

$$
R_k = \beta_{kp}\,n_p\left[c_k f_{cap} - a_i\,c_{sat,k}(T_p)\,
\exp\!\frac{4\sigma v_{m,k}}{\max(d_p,d_{min})\,k_BT_p}\right],\qquad
c_{sat,k} = \frac{p^\circ}{R_uT}\exp\frac{\Delta_M g_i(T) - g_k(T)}{R_uT}
$$

Evaporation part: $\dot M_{M,O} \mathrel{-}= \Delta_{M,O}W_{M,O}R^{ev}_k/\rho$,
$\dot Y_k \mathrel{+}= W_kR^{ev}_k/\rho$, $\dot H \mathrel{-}= R^{ev}_k h_k(T_p)/\rho$;
$\dot N = -\text{ramp}\,N\,\dot m_{ev}/(M_M+M_O)$ below $1.1d_{min}$.
Multipliers: $\beta_c \to f_{coag}\beta_c$, $J \to \eta J$.

**Discretisation / numerics:** production operator; $\ln p_{sat,k}$ tabulated in Python on the
`particleThermoT` grid (200–5000 K, 2 K) from the gas thermo and the lowest-Gibbs condensed
species of `Particles.phases` (`_phaseEnthalpyTable` now also returns $g_i$), interpolated
linearly in $T_p$. The Kelvin diameter is floored at $d_{min}$: the factor grows faster than
$d_p^2$ shrinks, and without the floor the first step diverged.

**Conservation & units:** mass exactly ($W_k = \Delta_MW_M+\Delta_OW_O$); energy with $h_k(T_p)$
leaving. Fe element flux (inflow vs radial outflow, gas + particles) closes to 0.995–0.998 in
the impinging-jet runs (gas-only 0.998).

**Code:** `sourceSystem.cpp::computeCondensationRates`, `addDisintegration` (shared with surface
reactions), `computeCoagulationRates`, `computeNucleationRates`;
`particleUtils.cpp::phaseActivity`, `condensationSaturationPressure`; options
`Particles.evaporation`, `surfaceTension`, `coagulationEnhancement`,
`NucleationChannel.efficiency`.

**Verification:** $p_{sat}$(Fe) = 3.5 Pa at 1811 K (literature ≈ 3.6 Pa). `evaporation=False`
reproduces the previous jet result bit for bit. Jet case: below 1400 K unchanged, particles
evaporate from ~13 mm and are gone by 15–18 mm (T ≈ 1600–1700 K). Python tests pass; C++
`DiffusionSystemTest.CylindricalCoordinates` fails, unrelated (diffusion code unchanged from
main).

**Limits / open points:** activity = lever-rule fraction (surface-coverage picture, not
two-phase equilibrium); bulk $\sigma$ and thermo at 1 nm; clusters (Fe$_n$) evaporate as whole
molecules; $N$ slightly negative at the wall node in some runs.

## 2026-10-22 — feat: thermo-consistent surface-reaction heat for the particle energy (staged)

**Area:** particles: energy (surface reactions)

**What changed (physics):** The heat a surface reaction deposits in the particle is no longer
prescribed (`SurfaceReaction.deltaH`, previously the OpenFOAM values) but follows from the
thermochemistry: the particle receives the absolute enthalpy of the impinging gas molecule and
loses that of the released one. Since $H$ is the absolute particle enthalpy and $T_p$ comes from
the condensed-phase tables, the reaction heat (including the enthalpy of the atom attached to or
removed from the Fe$_x$O$_y$ particle) is implied by gas + condensed thermo, and energy is
conserved exactly. `deltaH` remains as an optional extra heat (default 0, not conserving).

**Equations** ($R_r$ [kmol/m³/s], $h_k(T)$ absolute molar enthalpy of gas species $k$
[J/kmol] from the mechanism, $\rho$ gas density):

$$
\dot H\big|_{surf} = \sum_r \frac{R_r\,\big[h_{react}(T) - h_{prod}(T_p) - \Delta H_r^{extra}\big]}{\rho},
\qquad h_{prod} = 0 \text{ if no gas product}
$$

Replaced: $\dot H|_{surf} = \sum_r R_r\big[(h - x_O h')\Delta_M + h'\Delta_O - \Delta H_r\big]/\rho$
(composition change at $T_p$ plus user heat $-\Delta H_r$).

**Discretisation / numerics:** production operator (CVODE); $h_{prod}(T_p)$ from the species'
NASA polynomials (`CanteraGas::getSpeciesEnthalpy`).

**Conservation & units:** gas enthalpy flux of the reactant/product + $\dot H$ balance exactly;
J/kg gas for $H$.

**Code:** `sourceSystem.cpp::computeSurfaceReactionRates`, `chemistry0d.cpp::getSpeciesEnthalpy`;
`SurfaceReaction.deltaH` (meaning changed).

**Verification:** impinging-jet energy case (`single_Igor_impingingJet_fixedT_particlesEnergy.py`):
with the OpenFOAM deltaH table $T_p$ reached ~5000 K; now $T_p$ is bounded,
$T_p - T$ = −16 to −123 K. 300 K reaction enthalpies (thermo, kJ/mol): O$_2$→O −23, O→p −521,
H$_2$O→H$_2$ −30, H$_2$→H$_2$O +30, H→OH +91, OH→FeOH +68 (OpenFOAM: all exothermic).

**Limits / open points:** bulk condensed thermo for nm particles (no surface energy);
particles end up colder than the gas because reduction and etching are endothermic, while the
experiment reports hotter particles; heat budget not yet analysed.

## 2026-10-21 — fix: coagulation kernel guard constant (staged)

**Area:** particles: coagulation

**What changed (physics):** No model change. The Brownian coagulation kernel now evaluates the
equations it was meant to: the division guard `ROOTVSMALL` was 1e-18 (OpenFOAM's
single-precision value) and dominated the particle mass and squared diameter it was added to,
which made coagulation ~300× too slow, i.e. effectively off.

**Equations** (unchanged intent, §4.4c of PhysicsAndNumerics.md):

$$
\dot N = -\tfrac12\beta_c\,\rho N_A N^2,\quad
\beta_c = \left(\beta_{fm}^{-1}+\beta_{c\ell}^{-1}\right)^{-1},\quad
\beta_{fm} = 4d_p^2\sqrt{\frac{\pi k_BT}{m_p}},\quad
\beta_{c\ell} = \frac{4k_BTd_p}{3\mu}\left(\frac{2}{d_p} + \frac{5.028\,\lambda_g}{d_p^2}\right)
$$

Previously implemented: $m_p \to m_p + 10^{-18}$ kg and $d_p^2 \to d_p^2 + 10^{-18}$ m² inside
the kernel. For iron nuclei $m_p \sim 10^{-24}$ kg and $d_p^2 \sim 3\times10^{-19}$ m², so
$\beta_{fm}$ was reduced by $\approx\sqrt{10^{-18}/m_p}$. Now the guard is $10^{-150}$
(OpenFOAM's double-precision `ROOTVSMALL`, as used by the reference `monoCalcBeta`).

**Discretisation / numerics:** none; production operator (per grid point, CVODE).

**Conservation & units:** unchanged; coagulation changes only $N$ [kmol/kg gas], particle
mass moments conserved.

**Code:** `particleUtils.cpp::monodisperseCoagulationKernel`; option `Particles.coagulation`.

**Verification:** at 50 mm in `burnerFla_particles_cond` ($d_p$ = 0.565 nm,
$m_p = 7.4\times10^{-25}$ kg, $T$ = 1738 K) the logged kernel was $1.35\times10^{-18}$ m³/s,
reproduced exactly with the old guard; correct value $4.06\times10^{-16}$ (×301). Before the
fix, cond with and without coagulation differed by 1.7 % in $N$ at 100 mm. After: vs the
OpenFOAM reference (same kernel), $N$ within 0.2–3 %, $d_p$ 1.639 vs 1.633 nm at 100 mm.

**Limits / open points:** all earlier runs with `coagulation=True` (cond, oxid cases) are
affected and must be rerun.

## 2026-10-20 — feat: OpenFOAM Sutherland viscosity/conductivity option (staged)

**Area:** gas transport

**What changed (physics):** Optional replacement of Cantera's mixture viscosity and thermal
conductivity by OpenFOAM 5's `sutherland` transport with `reactingMixture` mixing: species
Sutherland coefficients are averaged by mass fraction and the conductivity follows a modified
Eucken relation. Species diffusion is unchanged (still from `transportModel`). Intended for
code-to-code comparison with OpenFOAM, not as a better transport model.

**Equations** ($Y_k$ mass fractions as normalised by Cantera; $A_{s,k}$ [Pa·s·K$^{-1/2}$],
$T_{s,k}$ [K]; $c_v$ mixture mass-specific heat [J/kg/K]; $R = R_u/W$ [J/kg/K]):

$$
A_s = \sum_k Y_k A_{s,k},\quad T_s = \sum_k Y_k T_{s,k},\quad
\mu = \frac{A_s\sqrt{T}}{1 + T_s/T},\quad
\lambda = \mu\, c_v\left(1.32 + 1.77\,\frac{R}{c_v}\right)
$$

Previously (and still with `'cantera'`): $\mu$, $\lambda$ from the Cantera transport object
(mixture-averaged: Wilke / Mathur, mole-fraction based). Species without given coefficients:
least-squares fit $\sqrt{T}/\mu_k = 1/A_{s,k} + (T_{s,k}/A_{s,k})/T$ to Cantera's pure-species
viscosity at 28 temperatures in 300–3000 K.

**Discretisation / numerics:** none new; $\lambda$ enters the energy diffusion operator
($B = 1/(\rho c_p)$, $D = \lambda$), $\mu$ the momentum diffusion, the grid damping length and,
with particles, the mean free path and Stokes–Einstein particle diffusivity. Evaluated per grid
point in `updateChemicalProperties`. Affects the steady solution.

**Conservation & units:** no change to conservation; coefficients in SI as in OpenFOAM's
`thermo.compressibleGas`.

**Code:** `chemistry0d.cpp::CanteraGas::getViscosity/getThermalConductivity`,
`_ember.pyx::_sutherlandCoefficients`, `utils.readFoamSutherland`; options
`Chemistry.conductivityModel` (`'cantera'` | `'sutherland'`), `Chemistry.sutherlandCoefficients`.

**Verification:** ember $\mu$, $\lambda$ equal an independent evaluation of the formulas to
$6\times10^{-16}$; $\lambda/\lambda_{Cantera}$ = 0.58 (0.5 mm) … 0.99 (100 mm) in the
H$_2$/O$_2$/Ar burner flame. `burnerFla_test_OFgrid.py` vs OpenFOAM without radiation (same
grid): $T$ within 4–9 K (was up to 41 K near the burner), mean RMS deviation of $T$ + 11 species
2.89 % → 1.57 %.

**Limits / open points:** mass-fraction mixing strongly under-weights H$_2$ (13–42 % lower
$\lambda$ in the preheat zone than Wilke/Mathur); diffusion coefficients do not use the
Sutherland data; remaining ember–OpenFOAM differences come from the correction-flux and
energy-equation formulations.

