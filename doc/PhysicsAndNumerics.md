# Ember (emberXnano): physics, equations and numerics

A technical reference for the solver as it stands, written from the code (state: commit
`5837ad2`, 2026-10-13). Each item has a short description, the equation the code solves, and
the file that implements it. Changes made after this snapshot are logged in
[`PhysicsChangelog.md`](PhysicsChangelog.md). Background and validation for individual features
are in [`InletFluxBC.md`](InletFluxBC.md), [`ImpingingJet.md`](ImpingingJet.md) and
[`ParticleEnergy.md`](ParticleEnergy.md).

---

## 1. What kind of solver this is

- **Transient.** Ember integrates the time-dependent equations, $\partial\phi/\partial t = \dots$.
  Steady flames (burner, impinging jet, steady strained flame) are reached by time-marching
  until a termination criterion holds (§6.4). There is no Newton solver for the steady problem,
  as Cantera has.
- **PDEs reduced to ODEs (method of lines).** The model is a system of 1D parabolic PDEs in
  $(x,t)$, the quasi-1D similarity reduction of the 2D/axisymmetric reacting Navier–Stokes
  equations. Space is discretised by finite differences / finite volumes on an adaptive,
  non-uniform grid. Time is then integrated with operator splitting. Each split operator is a
  system of ODEs, integrated by its own stiff or non-stiff integrator:
  - **diffusion:** one tridiagonal linear ODE system per variable;
  - **convection:** one ODE system per scalar, with $T$, $U$ and $W$ coupled;
  - **production (chemistry and particles):** one independent ODE system per grid point.
- **Low Mach number.** Pressure $p$ is constant in space and time. Only the pressure curvature
  of the similarity solution enters, through the strain rate $a(t)$.

## 2. Coordinates, state vector, geometry

The independent variables are the normal coordinate $x$ and time $t$, with $r = x^{\alpha}$ the
"radius" weighting.

| Geometry (`flameGeometry`) | $\alpha$ | $\beta$ | Configuration |
|---|---|---|---|
| `planar` | 0 | 1 | planar strained flame / planar jet |
| `disc` | 0 | 2 | axisymmetric opposed / impinging jet (as Cantera) |
| `cylindrical` | 1 | 1 | curved (tubular) flame |

Each grid point $j$ carries the column of the state matrix `state` ($n_{vars}\times n_{points}$,
`flameSolver.cpp`):

$$
\boldsymbol\phi_j = \big(U,\; T,\; Y_1,\dots,Y_K,\; \underbrace{N,\; M_M,\; M_O,\; H}_{\text{particle moments}}\big)_j
$$

- $U = v_r/r$ (or $v_y/y$ in the planar case) is the tangential-velocity gradient. In the
  potential-flow region $U \to a/\beta$.
- $V = \rho v$ is the normal mass flux, and $rV$ is the "radius-weighted" flux. $V$ is not a
  state variable; continuity gives it.
- $W$ (mean molecular weight) is carried as an extra variable inside the convection operator
  only.

## 3. Gas-phase governing equations

Mixture-averaged transport, ideal gas, $\rho = pW/(R_u T)$.

**Viscosity and conductivity** (`CanteraGas::getViscosity/getThermalConductivity`).
`Chemistry.conductivityModel = 'cantera'` (default): $\mu$, $\lambda$ from the Cantera transport
model. `'sutherland'`: OpenFOAM's `sutherland` + `reactingMixture` model, with per-species
coefficients $A_{s,k}$, $T_{s,k}$ (`Chemistry.sutherlandCoefficients`, e.g.
`utils.readFoamSutherland`; missing species fitted to Cantera's pure-species viscosity):

$$
\mu = \frac{A_s\sqrt{T}}{1 + T_s/T},\quad A_s = \sum_k Y_k A_{s,k},\ T_s = \sum_k Y_k T_{s,k},\qquad
\lambda = \mu\, c_v\left(1.32 + 1.77\,\frac{R}{c_v}\right),\ R = \frac{R_u}{W}
$$

Species diffusion coefficients always come from `transportModel`.

**Continuity** (`ConvectionSystemUTW::f`). $\partial\rho/\partial t$ is assembled from the other
operators' rates, $\dot\rho = -\rho\,(\dot T/T + W\sum_k \dot Y_k/W_k)$ (`prepareIntegrators`),
and $rV$ is obtained by quadrature:

$$
\frac{\partial\rho}{\partial t} + \beta\rho U + \frac{1}{r}\frac{\partial (rV)}{\partial x} = 0
\quad\Longrightarrow\quad
(rV)_{j+1} = (rV)_j - h_j\, r_{j+\frac12}\big(\dot\rho_j + \beta\rho_j U_j\big)
$$

The integration constant is fixed by `continuityBC`: `Left` (given $rV_0$), `Zero` (at the
stagnation point), `Temp`, `Qdot`, and `Wall` (new, §5.3).

**Tangential momentum**, which carries the imposed pressure curvature through the strain rate
$a(t)$:

$$
\rho\frac{\partial U}{\partial t} + V\frac{\partial U}{\partial x}
= \frac{\partial}{\partial x}\!\left(\mu\frac{\partial U}{\partial x}\right)
- \rho U^2 + \rho_u\!\left(\frac{a^2}{\beta^2} + \frac{1}{\beta}\frac{da}{dt}\right)
$$

so that the pressure-curvature term is $\Lambda = -\rho_u a^2/\beta^2$ for constant $a$.

**Energy**, with $j_k = -\rho D_k\,\partial Y_k/\partial x - (D_k^T/T)\,\partial T/\partial x + Y_k j_{corr}$:

$$
\rho c_p\frac{\partial T}{\partial t} + \rho c_p v\frac{\partial T}{\partial x}
= \frac{1}{r}\frac{\partial}{\partial x}\!\left(r\lambda\frac{\partial T}{\partial x}\right)
- \sum_k \dot\omega_k h_k
- \Big(\sum_k \frac{c_{p,k}}{W_k} j_k\Big)\frac{\partial T}{\partial x} - \dot q_{loss}
$$

**Species**:

$$
\rho\frac{\partial Y_k}{\partial t} + \rho v\frac{\partial Y_k}{\partial x}
= -\frac{1}{r}\frac{\partial (r j_k)}{\partial x} + \dot\omega_k W_k + \rho\,S^{p}_{Y_k}
$$

- $S^{p}_{Y_k}$ is the gas-phase sink/source from the particles (§4), which is the two-way
  coupling of mass.
- The Fickian part $-\rho D_k\nabla Y_k$ is treated implicitly in the diffusion operator.
- The Soret flux and the mass-conserving correction flux
  $j_{corr} = -\sum_k(j_k^{Fick}+j_k^{Soret})$ are explicit "cross terms" (`updateCrossTerms`),
  held constant over a global step. $j_{corr}$ is smoothed with one implicit Le = 1 diffusion
  step for stability.

## 4. Particle phase: moment equations

### 4.1 Transported moments and units

The moments are mass-specific, like $Y_k$. All of them obey the same transport equation (§4.2).

| Row | Symbol | Meaning | Unit |
|---|---|---|---|
| `kN`=0 | $N$ | particle number | kmol particles / kg gas |
| `kM`=1 | $M_M$ | particle-phase metal mass (metal = `Particles.metal`) | kg / kg gas |
| `kO`=2 | $M_O$ | particle-phase O mass (only with `nMoments=3`) | kg / kg gas |
| `kH`=3 | $H$ | absolute particle enthalpy (only with `energy=True`) | J / kg gas |

- Number density: $n_p = N\rho N_A$, with Cantera's $N_A = 6.022\times10^{26}\,\mathrm{kmol^{-1}}$.
- Exact invariant of all source terms: $\sum_k Y_k + M_M + M_O = 1$.

### 4.2 Transport (`flameSolver.cpp::prepareIntegrators`, `updateParticleDiffusivity`)

$$
\rho\frac{\partial \phi}{\partial t} + \rho v\frac{\partial \phi}{\partial x}
= \frac{1}{r}\frac{\partial}{\partial x}\!\left(r\rho D_p\frac{\partial\phi}{\partial x}\right)
+ \rho\,\dot S_\phi ,\qquad \phi\in\{N, M_M, M_O, H\}
$$

- Plain Fickian diffusion only: no Soret, no correction flux, no thermophoresis.
- Stokes–Einstein–Cunningham diffusivity, evaluated per grid point from the local mean size:

$$
D_p = \frac{k_B T\, C_c(Kn)}{3\pi\mu d_p},\quad
C_c = 1 + Kn\big(1.257 + 0.4\,e^{-1.1/Kn}\big),\quad
Kn = \frac{2\lambda_g}{d_p},\quad
\lambda_g = \frac{\mu}{p}\sqrt{\frac{\pi R_u T}{2W}}
$$

  with $d_p \ge d_{min}$, and $D_p =$ `momentDiffusivity` where no particles exist.
- *Fix 9206616:* the diffusion coefficient is $\rho D_p$, with prefactor $B = 1/\rho$. Before
  the fix it was $D_p/\rho$, which is 30–200× too fast at 3000 Pa.

### 4.3 Size closure (`particleUtils.cpp`)

Particle size is derived from the moments, not transported. Monodisperse spheres:

$$
x_O = \frac{M_O/W_O}{M_M/W_M},\qquad
\rho_p(x_O)=\text{piecewise-linear interpolation of } (x_i,\rho_i)\in\texttt{Particles.phases}
$$

$$
v_p = \frac{M_M + M_O}{N N_A\,\rho_p(x_O)},\qquad d_p = \left(\frac{6v_p}{\pi}\right)^{1/3},\qquad
m_p = \frac{M_M+M_O}{N N_A}
$$

The ratio of the last phase, $x_{max}$, is the oxidation cap,
$f_{cap} = \max(1 - x_O/x_{max},\,0)$.

### 4.4 Source terms (`sourceSystem.cpp`, evaluated per grid point)

The particle source terms are added to the chemistry ODE right-hand side:

$$
\frac{d\phi}{dt} = Q_\phi - D_\phi,\qquad
\frac{dY_k}{dt} \mathrel{+}= Q_{Y_k} - D_{Y_k}
$$

(production $Q$ and destruction $D$ are kept separate for the QSS integrator).

**(a) Collision nucleation** (`computeNucleationRates`). A hard-sphere A + B collision forms one
particle:

$$
J = \tfrac12^{\delta_{AB}}\;\pi d_{AB}^2\sqrt{\frac{8k_BT}{\pi\mu_{AB}}}\;N_A\,c_A c_B
\quad[\mathrm{kmol\,m^{-3}s^{-1}}],\qquad d_{AB}=\tfrac12(d_A+d_B)
$$

$$
\dot N = \frac{J}{\rho},\quad
\dot M_M = \frac{J(\nu_A W_A+\nu_B W_B)}{\rho},\quad
\dot Y_{A} = -\frac{\nu_A J W_{A}}{\rho},\quad
\dot Y_{B} = -\frac{\nu_B J W_{B}}{\rho},\quad
\dot H = \frac{J(\nu_A h_A+\nu_B h_B)}{\rho}
$$

Here $c_A, c_B$ are concentrations [kmol/m³], and $d$, $\mu_{AB}$ come from the mechanism's
transport data.

**(b) Classical nucleation theory** (optional, `precursorSpecies`):

$$
S = \frac{p_{mon}}{p_{sat}},\quad p_{sat}=10^{A-B/(C+T-273.15)}\cdot10^5,\quad
\Theta = \frac{\pi d_1^2\sigma}{k_BT},\quad g^* = \left(\frac{2\Theta}{3\ln S}\right)^3
$$

$$
J = n_1^2 S v_1\sqrt{\frac{2\sigma}{\pi m_1}}\exp\!\left(\Theta - \frac{4\Theta^3}{27\ln^2S}\right),
\quad \dot N = J/\rho,\quad \dot M_M = J g^* W_1/\rho
$$

> ⚠ **Unverified branch, likely unit errors** (no example case uses CNT):
> 1. $p_{mon}$ is computed as $c\cdot1000\cdot R_u T$ with Cantera's $R_u$ in J/(kmol K), which
>    is 1000× too large.
> 2. $J$ here is in particles/(m³ s), but it is added to $N$ [kmol/kg] and to $M_M$ via $W$
>    [kg/kmol] without dividing by $N_A$.
>
> Check both before using CNT.

**(c) Brownian coagulation** (`computeCoagulationRates`). Monodisperse, conserves mass, so only
$N$ changes. The kernel interpolates harmonically (Fuchs) between the free-molecular and
continuum-slip limits:

$$
\dot N = -\tfrac12\beta_c\,\rho N_A N^2,\quad
\beta_c = \left(\beta_{fm}^{-1}+\beta_{c\ell}^{-1}\right)^{-1},\quad
\beta_{fm} = 4d_p^2\sqrt{\frac{\pi k_BT}{m_p}},\quad
\beta_{c\ell} = \frac{4k_BTd_p}{3\mu}\left(\frac{2}{d_p} + \frac{5.028\,\lambda_g}{d_p^2}\right)
$$

**(d) Surface reactions** (`computeSurfaceReactionRates`, commit 7b00e2a). One reaction is
$\text{reactant}(g) + \text{particle} \to \text{product}(g) + \text{particle}'$:

$$
R_r = A_r\,e^{-T_{a,r}/T_p}\;\theta_r\;\big(\pi d_p^2 n_p\big)\;c_{k_r}\quad[\mathrm{kmol\,m^{-3}s^{-1}}]
$$

- **Site fraction $\theta_r$:** $X_M = 1/(1+x_O)$ for oxidation and etching; $X_O = 1 - X_M$
  for reduction; 1 for metal-bearing deposition. It is multiplied by $f_{cap}$ if the reaction
  raises O/metal.
- **Particle change** $(\Delta_M, \Delta_O)$ per event, from the reactant − product
  composition: $\dot M_M = \Delta_M W_M R/\rho$, $\dot M_O = \Delta_O W_O R/\rho$.
- **Gas species:** $\dot Y_{react} = -W_{react}R/\rho$, $\dot Y_{prod} = +W_{prod}R/\rho$
  (mass conserved by construction).
- **Disintegration** below $d_{min}$ (ramped linearly over $[d_{min}, 1.1d_{min}]$): mass loss
  also removes particles, $\dot N = -\text{ramp}\cdot N\,\dot m_{loss}/(M_M+M_O)$, so the mean
  particle mass stays fixed.

**(e) Heterogeneous condensation** (`computeCondensationRates`, commit df9cb16). Whole gas
molecules $k$ (metal, clusters, oxides) stick on collision, with sticking probability 1 and
$N$ unchanged:

$$
R_k = \frac{\pi}{4}(d_p+d_k)^2\sqrt{\frac{8k_BT}{\pi\mu_{kp}}}\;n_p\,c_k\,(\times f_{cap}),\qquad
\dot M_M = \Delta_M W_M R_k/\rho,\;\; \dot M_O = \Delta_O W_O R_k/\rho,\;\;
\dot Y_k = -W_kR_k/\rho,\;\; \dot H = R_k h_k(T)/\rho
$$

### 4.5 Particle energy and temperature (commit b98ce11, `energy=True`)

$H$ is the absolute particle enthalpy, on the same reference as the gas NASA thermo.

**Particle temperature.** $T_p$ is obtained by inverting tabulated phase enthalpies
$h_i(T)$ [J/kmol metal]:

- The tables are built in `_ember.pyx::_phaseEnthalpyTable` from `nasa_condensed.yaml`, using
  at each $T$ the species with the lowest Gibbs energy.
- The enthalpy is interpolated in O/metal by the lever rule:

$$
\frac{H}{M_M/W_M} = (1-w)\,h_i(T_p) + w\,h_{i+1}(T_p),\qquad w = \frac{x_O-x_i}{x_{i+1}-x_i}
$$

- $T_p$ is found by bisection on the table (`particleTemperature`).

**Enthalpy balance**:

$$
\dot H = \underbrace{\sum_{\text{nucl.}}\frac{J\sum\nu h_k(T)}{\rho} + \sum_k\frac{R_k h_k(T)}{\rho}}_{\text{gas molecules incl. latent heat}}
+ \sum_r \frac{R_r\big[(h - x_O h')\Delta_M + h'\Delta_O - \Delta H_r\big]}{\rho}
- N N_A\,(q + P)
$$

with $h' = \partial h/\partial x$.

- **Conduction:** the free-molecular and continuum limits are combined harmonically,
  $q = q_{fm}q_c/(q_{fm}+q_c)$:

$$
q_{fm} = \alpha_T\,\pi d_p^2\,\frac{p}{8}\,\bar c\,\frac{\gamma+1}{\gamma-1}\left(\frac{T_p}{T}-1\right),\qquad
q_c = 2\pi d_p\lambda (T_p - T)
$$

- **Radiation, Rayleigh limit:**

$$
P = \frac{199.1\,\pi^3 d_p^3 E(m)\,k_B^5}{h^4c^3}\,(T_p^5 - T_{surr}^5)
$$

- **Coupling:** the heat goes to the particle only. It is **not** fed back to the gas energy
  equation. The surface-reaction Arrhenius rates use $T_p$; collision rates use $T$.

## 5. Boundary conditions

### 5.1 Overview

| BC (`BoundaryCondition::`) | Where | $Y_k$, moments | $T$ | $U$ |
|---|---|---|---|---|
| `FixedValue` | left/right | Dirichlet | Dirichlet | from momentum eq. |
| `ZeroGradient` | left/right | $\phi_0=\phi_1$ substituted in stencil | same | same |
| `ControlVolume` | left (twin/cylindrical) | upwind inflow + diffusion into CV $[0,x_1]$ | same | — |
| `WallFlux` | left (cylindrical) | — | Robin: $\lambda\,\partial_x T = K_{wall}(T-T_\infty)$ | — |
| `InletFlux` (new) | left, burner/jet inlet | **flux balance** | Dirichlet $T_{left}$ | free (burner) / 0 (jet) |
| `Wall` (new) | right, impinging jet | **zero flux**, half CV | Dirichlet $T_{wall}$ | 0 (no slip) |

### 5.2 `InletFlux` — burner face (commit 9bf4b55, conservation fixes in 37240fd)

This is Cantera's `Inlet1D` condition. The **total** flux of species/moment $\phi$ is imposed,
not its face value:

$$
\dot m\,\phi_{in} = \dot m\,\phi_0 + j_\phi(0)
$$

Discretely, on the cell $[0,x_1]$ with $V_c = x_1^{\alpha+1}/(\alpha+1)$:

$$
\left.\frac{d\phi_0}{dt}\right|_{conv} = -\frac{\max(rV_0,0)}{\rho_0 V_c}(\phi_0-\phi_{in}),\qquad
\left.\frac{d\phi_0}{dt}\right|_{diff} = B_0\,(\alpha+1)\,\frac{\overline{\rho D}_{1/2}\,(\phi_1-\phi_0)}{h_0^2}
$$

- **Flux conservation at the first face:** the conductance $c_{2,0}$ is set so that node 1
  receives exactly what node 0 loses.
- **Mean molecular weight at the inlet:** $W_0$ gets the matching inflow term,
  $\dot W_0 = \frac{rV_0}{\rho_0V_c}W_0\frac{W_{in}-W_0}{W_{in}}$.
- **Cross fluxes:** the Soret and correction fluxes across face $\frac12$ are credited to node 0.
- **Temperature:** $T_0 = T_{left}$, and its rates are zeroed, so they do not enter $\dot\rho$.
- **No chemistry at node 0.**

### 5.3 `Wall` — impinging-jet wall (commit 37240fd)

The node $N-1$ is a half control volume of width $h/2$ whose only flux crosses the inner face:

$$
\frac{d\phi_{N-1}}{dt} = \frac{B_{N-1}}{r_{N-1}\,h_{N-2}/2}\;r_{N-\frac32}\,\overline{\rho D}\,\frac{\phi_{N-2}-\phi_{N-1}}{h_{N-2}}
$$

- **Wall values:** $T = T_{wall}$, $U = 0$ (no slip), and $U(0) = 0$ (plug-flow inlet).
- **Continuity** is integrated backward from $(rV)_{N-1}=0$, so the inlet mass flux $V_0$ is
  an output.
- No chemistry, no particle deposition, and no thermophoresis at the wall.

### 5.4 Strain-rate / mass-flux control (commit 37240fd)

With `ImpingingJet(massFlux=ṁ*)`, a PI controller on $\ln a$ is applied each global step:

$$
a = a_0\exp\!\Big(K_P\,e + K_I\!\int e\,dt\Big),\qquad e = \frac{\dot m^* - V_0}{\dot m^*}
$$

The mapping to Cantera is $a = \beta\sqrt{-\Lambda/\rho_u}$. Only the converged state is
physical.

### 5.5 Prescribed temperature (`ImpingingJet(temperatureProfile=…)`)

$T(x)$ is reset to the interpolated profile each step, and $\partial T/\partial t \equiv 0$ in
every operator (convection, diffusion with $B=0$, source, cross terms). This matches Cantera's
`energy_enabled = False`.

## 6. Numerics

### 6.1 Spatial discretisation (`diffusionSystem.cpp`, `convectionSystem.cpp`, `grid.cpp`)

The grid is non-uniform, with $h_j = x_{j+1}-x_j$ and $\Delta_j = (h_{j-1}+h_j)/2$.

**Diffusion** is second-order central, in conservative flux form:

$$
\frac{B_j}{r_j}\frac{\partial}{\partial x}\!\left(rD\frac{\partial\phi}{\partial x}\right)_j \approx
\frac{B_j}{r_j\Delta_j}\left[r_{j+\frac12}\bar D_{j+\frac12}\frac{\phi_{j+1}-\phi_j}{h_j} - r_{j-\frac12}\bar D_{j-\frac12}\frac{\phi_j-\phi_{j-1}}{h_{j-1}}\right]
$$

with $\bar D_{j+\frac12} = \tfrac12(D_j+D_{j+1})$. For each variable, $(B, D)$ is:

| Variable | $B$ | $D$ |
|---|---|---|
| $U$ | $1/\rho$ | $\mu$ |
| $T$ | $1/(\rho c_p)$ | $\lambda$ |
| $Y_k$ | $1/\rho$ | $\rho D_k$ |
| moments | $1/\rho$ | $\rho D_p$ |

**Convection** is first-order upwind:

$$
v_j\frac{\partial\phi}{\partial x}\Big|_j \approx v_j\frac{\phi_{j+1}-\phi_j}{h_j}\;(v_j<0),\qquad
v_j\frac{\phi_j-\phi_{j-1}}{h_{j-1}}\;(v_j\ge0)
$$

The outflow boundary uses a one-sided difference. Continuity uses first-order left-point
quadrature (§3).

**Adaptive grid**, applied between global steps:

- Points are inserted where $|\phi_{j+1}-\phi_j| > \texttt{vtol}\cdot\text{range}(\phi)$ or
  $|\phi'_{j+1}-\phi'_j| > \texttt{dvtol}\cdot\text{range}(\phi')$, and also to meet
  `gridMax` and `uniformityTol` (adjacent-spacing ratio).
- Points are removed with hysteresis (`rmTol`). Variables with range < `absvtol` are ignored.
- New values come from interpolation clamped to the bracketing points.
- Boundary regridding (extending or trimming the domain) is skipped for the fixed impinging-jet
  domain.

### 6.2 Operator splitting (`splitSolver.cpp`)

One global step $\Delta t$ (`globalTimestep`, fixed, default $2\times10^{-5}$ s) is symmetric
(Strang-type) in the order

$$
\mathcal D_{\Delta t/4}\;\mathcal C_{\Delta t/2}\;\mathcal D_{\Delta t/4}\;\mathcal P_{\Delta t}\;\mathcal D_{\Delta t/4}\;\mathcal C_{\Delta t/2}\;\mathcal D_{\Delta t/4}
$$

where $\mathcal D$ is diffusion, $\mathcal C$ convection and $\mathcal P$ production.

Each operator $i$ integrates $\dot\phi = f_i(\phi) + k_i$, where the split constant $k_i$
carries the other operators' rates from the previous step,
$\dot\phi_i = \Delta\phi_i/\Delta t - k_i$.

- **`balanced`** (default, Speth et al. 2013):

$$
k_D = \tfrac14(\dot\phi_P+\dot\phi_C+\dot\phi_X-3\dot\phi_D),\quad
k_C = \tfrac14(\dot\phi_P+\dot\phi_D+\dot\phi_X-3\dot\phi_C),\quad
k_P = \tfrac12(\dot\phi_C+\dot\phi_D+\dot\phi_X-\dot\phi_P)
$$

  At steady state the full residual $\sum_i\dot\phi_i = 0$ is satisfied exactly; there is no
  $O(\Delta t)$ splitting error in the steady solution.
- **`strang`:** $k_D = \dot\phi_X$, $k_C = k_P = 0$. The steady state depends on $\Delta t$.
  Measured: +7.6 % error in $a$ for the impinging jet, and a mass leak through $\dot\rho$ at
  burner inlets.

### 6.3 Time integrators per operator

| Operator | System | Integrator |
|---|---|---|
| Diffusion | linear tridiagonal, one per variable, $\dot\phi = A\phi + k$ | implicit BDF1 (8 substeps in the first step), then BDF2; Thomas algorithm; $\Delta t_{diff} = \min(\texttt{mult}\cdot\Delta_j^2/(BD),\,\Delta t)$ |
| Convection | $T$–$U$–$W$ coupled (3$N$ ODEs); each $Y_k$ and moment separately ($N$ ODEs) | SUNDIALS CVODE, Adams (non-stiff), adaptive |
| Production | per grid point: $U, T, Y_k$, moments ($K+2+n_m$ ODEs) | CVODE BDF with dense analytic/FD Jacobian (`chemistryIntegrator='cvode'`), or explicit QSS predictor–corrector (`'qss'`) |
| Cross terms | Soret, $j_{corr}$, energy enthalpy flux | explicit, frozen over $\Delta t$ |

The particle enthalpy relaxes on ~µs scales at nm sizes, so `energy=True` needs `'cvode'`.
Production is parallel over grid points and diffusion over variables (TBB).

### 6.4 Termination (`TerminationCondition.measurement`)

Checks start at `tMin`.

- `'Q'`: the integrated heat release is steady to `tolerance` (RMS) over `steadyPeriod`.
- `'dTdt'`: $\|T^{-1}\partial T/\partial t\|_2/\sqrt{N_{pts}} <$ `dTdtTol`.
- `'moments'` (new, 37240fd), for a prescribed $T$:
  $\max_m \|\dot\phi_m\|_2/\sqrt{N_{pts}}/\max_j|\phi_m| < \texttt{momentsTol}$ [1/s].

## 7. Recent implementations (2026-09 → 2026-10)

| Commit | Kind | Physics / math change |
|---|---|---|
| 9bf4b55 | BC | `InletFlux`: flux-balance inlet $\dot m\phi_{in} = \dot m\phi_0 + j_\phi$; conservative first diffusion face ($c_{2,0}$); fixed $T_0$ with rates zeroed in $\dot\rho$; no chemistry at face |
| 9206616 | fix | particle diffusion coefficient $\rho D_p$ (was $D_p/\rho$) |
| 7b00e2a | particles | volume moment → mass moments $M_M, M_O$; density closure $\rho_p(x_O)$; Arrhenius surface reactions with site fractions and oxidation cap; disintegration below $d_{min}$ |
| df9cb16 | particles | heterogeneous condensation at hard-sphere gas–particle collision rate |
| 37240fd | BC/config | impinging jet: `Wall` BC (zero-flux half CV, no slip, fixed $T$), continuity from the wall, PI control of $\ln a$ on $\dot m$, prescribed $T(x)$, `'moments'` termination; inlet $W_0$ inflow term and first-face cross-flux credit |
| b98ce11 | particles | particle enthalpy moment $H$, $T_p$ from condensed-phase tables, latent heat, reaction heat, conduction and Rayleigh radiation |
| 2026-10-20 | transport | `conductivityModel='sutherland'`: OpenFOAM's mass-fraction-mixed Sutherland $\mu$ and modified-Eucken $\lambda$ (code-to-code comparison) |

## 8. Known limitations / open points

- No particle → gas energy coupling; no evaporation/sublimation; no coalescence heat.
- Monodisperse size closure (no size distribution); bulk thermo for sub-nm clusters.
- No thermophoresis, no wall deposition.
- CNT branch unit issues (§4.4b).
- Convection is first-order upwind, so grid adaptation carries the accuracy; refinement closed
  the remaining 1 % gap to Cantera in the impinging-jet validation.
- `ControlVolume` left BC still lacks the first-face cross-flux credit given to `InletFlux`.
- `conductivityModel='sutherland'` reproduces OpenFOAM, not better physics: mass-fraction
  mixing under-weights H$_2$ ($\lambda$ 13–42 % below Wilke/Mathur in H$_2$-rich preheat zones).
  Remaining differences to OpenFOAM: correction flux ($\propto Y_k$ in ember, inert species in
  OpenFOAM) and OpenFOAM's energy equation omitting the inert species' enthalpy terms.
