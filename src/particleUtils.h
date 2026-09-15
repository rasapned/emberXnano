#pragma once

//! Small standalone helpers for the particle moments (N, V) module.
//!
//! These are shared between FlameSolver (particle diffusivity, used in the
//! diffusion split operator) and SourceSystem (nucleation and coagulation,
//! used in the production split operator) so the underlying particle-size
//! and gas-transport formulas live in exactly one place.

//! Mean single-particle volume [m^3] implied by the local N, V moments.
//! N [kmol particles / kg gas] and V [m^3 particle volume / kg gas] carry
//! the same "per kg of gas" normalization, so rho cancels out of the ratio.
//! Returns 0 where no particles are present (N or V <= 0).
double particleVolumeFromMoments(double N, double V);

//! Particle diameter [m] for a given single-particle volume [m^3] (sphere-
//! equivalent). Returns 0 for vParticle <= 0.
double particleDiameterFromVolume(double vParticle);

//! Mean single-particle diameter [m] implied by the local N, V moments.
//! Returns 0 where no particles are present (N or V <= 0).
double computeParticleDiameter(double N, double V);

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
