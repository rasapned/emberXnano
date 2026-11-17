# Evaporation of iron/iron-oxide nanoparticles in emberXnano: the implemented model

*Short note for checking against the literature. Code: `src/sourceSystem.cpp::computeCondensationRates`,
`src/particleUtils.cpp`, `python/ember/_ember.pyx` (saturation tables). Option:
`Particles(evaporation=True)`, branch `feat/particle-thermo-evaporation`.*

## Context

emberXnano is a 1D flame solver with a monodisperse particle model. It transports per kg of
gas the particle number $N$ [kmol/kg], the particle metal mass $M_M$ and oxygen mass $M_O$
[kg/kg], and optionally the particle enthalpy $H$ [J/kg], which gives the particle temperature
$T_p$. The mean diameter $d_p$ follows from mass/number and a bulk density interpolated in
the O/Fe ratio $x_O$ between the phases Fe, FeO, Fe$_3$O$_4$ and Fe$_2$O$_3$.

Gas molecules $k$ that contain Fe (Fe, Fe$_2$…Fe$_7$, FeO, FeO$_2$) condense at the hard-sphere
collision rate. Evaporation makes this process reversible.

## Model

Net rate for gas species $k$ [kmol m⁻³ s⁻¹]:

$$
R_k = \beta_{kp}\,n_p\Big[\,c_k \;-\; a_i\,c_{sat,k}(T_p)\,\exp\!\Big(\frac{4\sigma(T_p)\,v_{m,k}}{\tilde d_p\,k_B T_p}\Big)\Big]
$$

| Symbol | Meaning |
|---|---|
| $\beta_{kp} = \frac{\pi}{4}(d_p+d_k)^2\sqrt{8k_BT/(\pi\mu_{kp})}$ | free-molecular gas–particle collision kernel at the **gas** temperature $T$; $\mu_{kp}$ reduced mass, $d_k$ from the mechanism's transport data; sticking/evaporation coefficient = 1 |
| $n_p$ | particle number density [m⁻³] |
| $c_k$ | gas concentration [kmol/m³] |
| $i$ | particle phase with the same O/Fe ratio as $k$: Fe$_n$ → Fe, FeO → FeO. Species without a matching phase (FeO$_2$) only condense |
| $a_i$ | **activity of phase $i$ at the surface = its lever-rule fraction** in the particle; e.g. at $x_O = 0.3$, $a_{Fe} = 0.7$, $a_{FeO} = 0.3$ |
| $c_{sat,k} = p_{sat,k}/(R_uT_p)$ | saturation concentration over the bulk (flat) phase |
| $\sigma = 2.4152 - 2.85\times10^{-4}\,T_p$ [J/m²] | surface tension (iron; same form as the OpenFOAM `surfCoeffs`) |
| $v_{m,k} = W_k/(\rho_i N_A)$ | volume of molecule $k$ in phase $i$ (so Fe$_7$ has 7× the Fe volume) |
| $\tilde d_p = \max(d_p, d_{min})$ | Kelvin diameter, floored at the disintegration diameter $d_{min}$ (0.5 nm) |

**Saturation pressure** of the reaction $k(g) \rightleftharpoons \Delta_M\,(\text{phase } i)$,
where $\Delta_M$ is the number of Fe atoms in $k$:

$$
p_{sat,k}(T) = p^\circ \exp\!\Big(\frac{\Delta_M\,g_i(T) - g_k^\circ(T)}{R_uT}\Big)
$$

- $g_k^\circ$ is the standard molar Gibbs energy of the gas species, from the mechanism's NASA
  polynomials.
- $g_i$ is the Gibbs energy of the condensed phase per kmol of Fe. It comes from Cantera's
  `nasa_condensed.yaml`, using at each $T$ the lowest-Gibbs species of the phase (Fe α/γ/δ/liquid;
  FeO solid/liquid).
- $p_{sat,k}$ is tabulated on 200–5000 K and interpolated linearly at $T_p$.

**Source terms.** The condensation part adds to $M_M$, $M_O$ and $H$, carrying $h_k(T)$. The
evaporation part:

$$
\dot M_M \mathrel{-}= \frac{\Delta_M W_{Fe} R^{ev}_k}{\rho},\quad
\dot M_O \mathrel{-}= \frac{\Delta_O W_O R^{ev}_k}{\rho},\quad
\dot Y_k \mathrel{+}= \frac{W_k R^{ev}_k}{\rho},\quad
\dot H \mathrel{-}= \frac{R^{ev}_k\,h_k(T_p)}{\rho}
$$

- **Mass:** conserved exactly.
- **Energy:** the evaporated molecules leave with their gas enthalpy at $T_p$, so the particle
  loses the latent heat $h_k(T_p) - h_{cond}(T_p)$ per molecule. There is no coupling back to the
  gas energy equation.

**Number.** Evaporation leaves $N$ unchanged, except below $1.1\,d_{min}$. There, removed mass
also removes particles, $\dot N = -\text{ramp}\cdot N\,\dot m_{ev}/(M_M+M_O)$, so the
particles disintegrate instead of shrinking to zero.

## What was checked

- **Saturation pressure:** $p_{sat}$(Fe) = 3.5 Pa at the melting point (1811 K); the literature
  gives about 3.6 Pa. At 1600 K: Fe 0.11 Pa, FeO 3×10⁻⁴ Pa.
- **Kelvin factor:** about $e^{8}$ for Fe at 0.5 nm and 1700 K.
- **Default off:** with `evaporation=False` the results are unchanged bit for bit.
- **Fe conservation:** the Fe element balance closes to 0.5 %.
- **Impinging-jet case** (3000 Pa, $T_{max}$ 1700 K): below ~1400 K nothing changes. Particles
  (d ≈ 0.5–1 nm) evaporate from ~13 mm (≈1550 K) and are gone by 15–18 mm.

## Points to check against the literature

1. **Surface activity as lever-rule fraction.** In bulk two-phase equilibrium (e.g. Fe + FeO)
   both phases have activity 1. Here the surface is treated as a patchwork, with the evaporating
   area proportional to the phase fraction. Is this reasonable for nm particles? Are there
   core–shell or oxide-layer models for this?
2. **Kinetic prefactor.** Evaporation uses the condensation kernel at gas $T$ (detailed balance
   with the condensation term). The Hertz–Knudsen flux would use $T_p$, i.e.
   $\propto\sqrt{T_p/m_k}$. Also: should evaporation/condensation coefficients < 1 be used, and
   which values?
3. **Kelvin equation at 0.5–2 nm.** Bulk $\sigma$ and bulk thermo are used, with no Tolman or
   size-dependent correction. Which $\sigma(T)$ for Fe and FeO (liquid/solid)? The same $\sigma$
   is currently used for all phases.
4. **Cluster evaporation.** Fe$_2$…Fe$_7$ evaporate as whole molecules, each with its own
   $p_{sat}$ (thermodynamically consistent, but tiny). Physically evaporation is
   monomer-dominated. Is the cluster channel negligible, or should it be removed?
5. **Floor at $d_{min}$.** Without it the Kelvin factor makes the rate diverge as $d_p \to 0$.
   Below $d_{min}$, disintegration takes over. Compare with how sectional or moment models treat
   evaporation of the smallest clusters.
6. **Monodisperse closure.** Real size distributions evaporate the small particles first.
   With one mean size, all particles evaporate together.
7. **Missing channels.** There is no evaporation for FeO$_2$, FeOH or Fe(OH)$_2$ (no matching
   phase), and no decomposition/sublimation of Fe$_3$O$_4$/Fe$_2$O$_3$ to FeO + O$_2$.

Suggested starting points:
- Friedlander, *Smoke, Dust and Haze* (Kelvin effect, condensation/evaporation kernels).
- Seinfeld & Pandis (transition-regime corrections, accommodation coefficients).
- Literature on iron/iron-oxide nanoparticle synthesis in Fe(CO)$_5$-doped flames, for
  measured particle evaporation temperatures and the Fe vapour pressure over nanoparticles.
