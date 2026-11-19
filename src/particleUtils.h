#pragma once

#include <cstddef>

struct ParticleOptions;

//! Small standalone helpers for the particle moments (N, mM, mO) module.
//!
//! These are shared between FlameSolver (particle diffusivity, used in the
//! diffusion split operator) and SourceSystem (nucleation, coagulation and
//! surface reactions, used in the production split operator) so the
//! underlying particle-size and gas-transport formulas live in exactly one
//! place. The particle material (metal, phases) comes from `opts`.

//! Particle-phase O/metal atomic ratio implied by the metal and O mass
//! moments [kg / kg gas]. Returns 0 where no metal is present (mM <= 0).
double particleOxygenRatio(double mM, double mO, const ParticleOptions& opts);

//! Particle bulk density [kg/m^3], interpolated linearly in the O/metal
//! ratio x between the phases (opts.phaseRatio, opts.phaseDensity). x is
//! clamped to the range of opts.phaseRatio.
double particleBulkDensity(double mM, double mO, const ParticleOptions& opts);

//! Mean single-particle volume [m^3] implied by the local N [kmol particles
//! / kg gas], mM and mO [kg / kg gas] moments. All three carry the same
//! "per kg of gas" normalization, so rho cancels out of the ratio. Returns
//! 0 where no particles are present (N or mM <= 0).
double particleVolumeFromMoments(double N, double mM, double mO,
                                 const ParticleOptions& opts);

//! Particle diameter [m] for a given single-particle volume [m^3] (sphere-
//! equivalent). Returns 0 for vParticle <= 0.
double particleDiameterFromVolume(double vParticle);

//! Mean single-particle diameter [m] implied by the local N, mM, mO
//! moments. Returns 0 where no particles are present (N or mM <= 0).
double computeParticleDiameter(double N, double mM, double mO,
                               const ParticleOptions& opts);

//! Particles below dmin (opts.minParticleDiameter) disintegrate instead of
//! shrinking, so their surface sources (surface reactions, condensation,
//! evaporation) are evaluated for the same mass in particles of diameter
//! dmin: raises dp [m] to dmin and scales the number density nParticles and
//! the single-particle mass mp (if given) so that the total mass is kept.
//! Mass sinks then stay proportional to the particle mass, and N and mass
//! go to zero together instead of crossing it. No change for dp >= dmin.
void floorParticleSize(double& dp, double& nParticles, double* mp,
                       const ParticleOptions& opts);

//! Particle enthalpy [J per kmol of metal] at temperature T [K] for the
//! O/metal atomic ratio x: lever rule between the bracketing phases (as for
//! the density), linear in T on the opts.particleThermoT grid (T clamped to
//! its range). If dhdx is given, it receives the slope dh/dx at fixed T.
double particleEnthalpyPerMetal(double T, double x, const ParticleOptions& opts,
                                double* dhdx=nullptr);

//! Particle temperature [K] from the enthalpy moment H [J / kg gas] and the
//! mass moments mM, mO [kg / kg gas], by inverting
//! particleEnthalpyPerMetal() (clamped to the table range). Returns Tgas
//! where no particles are present (mM <= 0).
double particleTemperature(double H, double mM, double mO, double Tgas,
                           const ParticleOptions& opts);

//! Lever-rule weight of phase `phase` (index into opts.phaseRatio) in a
//! particle of O/metal ratio x: 1 - w and w for the two phases bracketing x,
//! 0 for the others (1 if there is only one phase). Used as the activity of
//! the phase at the particle surface for evaporation.
double phaseActivity(double x, int phase, const ParticleOptions& opts);

//! Saturation pressure [Pa] of condensation species c (opts.condSpecies[c])
//! over its bulk condensed phase at temperature T [K]: exp of
//! opts.condLnPsat, interpolated linearly on the opts.particleThermoT grid
//! (T clamped to its range).
double condensationSaturationPressure(size_t c, double T, const ParticleOptions& opts);

//! Equilibrium pressure [Pa] of the metal monomer over a nucleus of nAtoms
//! metal atoms at temperature T [K]: p_sat over the pure-metal phase
//! (opts.nucLnPsatMonomer) times the Kelvin factor
//! exp(4 sigma v_1 / (d_n k_B T)), d_n = (6 nAtoms v_1 / pi)^(1/3).
double nucleusEquilibriumPressure(int nAtoms, double T, const ParticleOptions& opts);

//! Gas mean free path [m] from kinetic theory (Chapman-Enskog), evaluated
//! from the local dynamic viscosity [Pa*s], temperature [K], mixture
//! molecular weight [kg/kmol] and pressure [Pa].
double gasMeanFreePath(double mu, double T, double Wmx, double pressure);

//! Cunningham slip correction factor (Allen & Raabe, 1985 coefficients),
//! given the particle Knudsen number Kn = 2*meanFreePath/dp.
double cunninghamSlipCorrection(double Kn);

//! Monodisperse Brownian coagulation kernel [m^3/s]: the free-molecular and
//! continuum (slip-corrected) limits combined by harmonic mean (Fuchs
//! interpolation). dc = particle diameter [m], mp = single-particle mass
//! [kg], T = temperature [K], mu = dynamic viscosity [Pa*s], meanFreePath =
//! gas mean free path [m]. Returns 0 for dc <= 0 (no particles).
double monodisperseCoagulationKernel(double dc, double mp, double T,
                                      double mu, double meanFreePath);

//! Stokes-Einstein-Cunningham particle diffusivity [m^2/s] for a particle
//! of (mobility) diameter dp [m] at temperature T [K], gas viscosity mu
//! [Pa*s] and mean free path [m].
double particleDiffusivityFromDiameter(double dp, double T, double mu,
                                       double meanFreePath);

//! Fuchs interpolation kernel [m^3/s] for Brownian coagulation between
//! particles i and j (Seinfeld & Pandis, Table 13.1), valid from the
//! free-molecular to the continuum regime: collision diameters d [m],
//! masses m [kg] and diffusivities D [m^2/s] of both particles, gas
//! temperature T [K]. Tends to (pi/4) (d_i+d_j)^2 sqrt(c_i^2+c_j^2) in the
//! free-molecular and 2 pi (D_i+D_j)(d_i+d_j) in the continuum limit.
double fuchsCoagulationKernel(double di, double mi, double Di,
                              double dj, double mj, double Dj, double T);

//! Geometry of an aggregate of volume v [m^3] and surface area a [m^2]
//! made of equal spherical primary particles: primary diameter dpp = 6 v/a,
//! number of primaries np = a^3 / (36 pi v^2) and collision diameter
//! dc = dpp np^(1/Df) (Kruis et al. 1993), Df = fractal dimension. An area
//! below that of the volume-equivalent sphere is taken as the sphere
//! (np = 1, dc = dpp = volume-equivalent diameter).
void aggregateGeometry(double v, double a, double Df,
                       double& dpp, double& np, double& dc);

//! Bimodal coagulation (Jeong & Choi 2003; Sellmann et al. 2019, Eqs. 1-4)
//! per unit volume and time [1/m^3/s], for mode-1 (monomer, volume v1) and
//! mode-2 (volume v2 = r v1) number densities n1, n2 [1/m^3] and kernels
//! b11, b12, b22 [m^3/s]. A collision of two monomers forms 1/(r-1) mode-2
//! particles and puts back (r-2)/(r-1) monomers, which conserves number
//! and volume; r is raised to at least 2.
struct BimodalCoagulationRates
{
    double n1Loss; //!< mode-1 particles lost
    double n2Gain; //!< mode-2 particles formed (monomer-monomer collisions)
    double n2Loss; //!< mode-2 particles lost (mode 2 - mode 2 collisions)
    double monomersToMode2; //!< monomer units moved into mode 2 (== n1Loss)
};
BimodalCoagulationRates bimodalCoagulation(double n1, double n2, double r,
                                           double b11, double b12, double b22);

//! Characteristic sintering time [s] of the aggregate primaries (Kruis et
//! al. 1993 form): tau = A T dpp^4 exp(Ta / T), A [s/(K m^4)], Ta [K].
double sinteringTime(double T, double dpp, double A, double Ta);

//! Collision diameter [m] of single-mode particles of volume-equivalent
//! diameter dv [m]: with a prescribed primary diameter dpp
//! (opts.primaryDiameter > 0) and dv > dpp, the fractal aggregate value
//! dpp (dv/dpp)^(3/Df), Df = opts.fractalDimension; dv otherwise.
double fractalCollisionDiameter(double dv, const ParticleOptions& opts);

//! Mode-2 aggregates of a bimodal population (opts.bimodal) from its
//! moments N2 [kmol / kg gas], M2 [kg / kg gas], A2 [m^2 / kg gas] at gas
//! density rho [kg/m^3]. Returns false if there are none: less than one
//! aggregate per m^3 (vanishing N2 with finite M2 in the integrator's trial
//! states would mean absurdly large aggregates) or no mass. Otherwise sets
//! the aggregate mass mp [kg] (at least two monomers), volume v [m^3] and,
//! via aggregateGeometry(), dpp, np and dc; the area per aggregate is
//! bounded by those of the volume-equivalent sphere and of its unsintered
//! monomers, so that dpp >= opts.monomerDiameter.
bool bimodalAggregate(double N2, double M2, double A2, double rho,
                      const ParticleOptions& opts, double& mp, double& v,
                      double& dpp, double& np, double& dc);
