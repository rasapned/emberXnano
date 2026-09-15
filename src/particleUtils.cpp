#include "particleUtils.h"
#include "chemistry0d.h"

#include <cmath>

double particleVolumeFromMoments(double N, double V)
{
    if (!(N > 0) || !(V > 0)) {
        return 0.0;
    }
    // Avogadro's number here is Cantera's kmol-based constant, consistent
    // with computeNucleationRates() in sourceSystem.cpp.
    return V / (N * Cantera::Avogadro); // [m^3]
}

double particleDiameterFromVolume(double vParticle)
{
    if (!(vParticle > 0)) {
        return 0.0;
    }
    return std::cbrt(6.0 * vParticle / M_PI);
}

double computeParticleDiameter(double N, double V)
{
    return particleDiameterFromVolume(particleVolumeFromMoments(N, V));
}

double gasMeanFreePath(double mu, double T, double Wmx, double pressure)
{
    const double Ru_ = Cantera::GasConstant;
    return mu / pressure * std::sqrt(M_PI * Ru_ * T / (2.0 * Wmx));
}

double cunninghamSlipCorrection(double Kn)
{
    return 1.0 + Kn * (1.257 + 0.4 * std::exp(-1.1 / Kn));
}

double monodisperseCoagulationKernel(double dc, double mp, double T,
                                      double mu, double meanFreePath)
{
    if (!(dc > 0)) {
        return 0.0;
    }

    const double kB_ = Cantera::Boltzmann;
    const double ROOTVSMALL = 1e-18; // guards against division by zero

    // Free-molecular limit (kinetic theory of gases applied to particle-
    // particle collisions).
    double fm = std::sqrt(M_PI * kB_ * T / (mp + ROOTVSMALL)) * 4.0 * dc * dc;

    // Continuum (slip-flow) limit: Brownian-diffusion-driven collision rate,
    // with a leading-order Cunningham slip correction.
    double sf = 4.0 * kB_ * T * dc / (3.0 * mu + ROOTVSMALL) *
        (2.0 / (dc + ROOTVSMALL) +
         1.257 * (2.0 * meanFreePath) * 2.0 / (dc * dc + ROOTVSMALL));

    // Harmonic mean (Fuchs interpolation) between the two limits.
    return 1.0 / (1.0 / (ROOTVSMALL + fm) + 1.0 / (ROOTVSMALL + sf));
}
