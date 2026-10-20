# Physics & numerics changelog

One entry per implementation commit, newest first, written by the `physics-changelog` skill
(`.claude/skills/physics-changelog/SKILL.md`). The full current model is in
[`PhysicsAndNumerics.md`](PhysicsAndNumerics.md). Commits up to `5837ad2` (2026-10-13) are
summarised in its §7.

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

