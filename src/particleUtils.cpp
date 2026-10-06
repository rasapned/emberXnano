#include "particleUtils.h"
#include "chemistry0d.h"
#include "readConfig.h"

#include <algorithm>
#include <cmath>

double particleOxygenRatio(double mM, double mO, const ConfigOptions& opts)
{
    if (!(mM > 0) || !(mO > 0)) {
        return 0.0;
    }
    return (mO / opts.oxygenWeight) / (mM / opts.metalWeight);
}

double particleBulkDensity(double mM, double mO, const ConfigOptions& opts)
{
    const std::vector<double>& xs = opts.phaseRatio;
    const std::vector<double>& rhos = opts.phaseDensity;
    size_t n = xs.size();

    double x = std::min(std::max(particleOxygenRatio(mM, mO, opts), xs[0]), xs[n-1]);
    for (size_t i = 0; i + 1 < n; i++) {
        if (x <= xs[i+1]) {
            double w = (x - xs[i]) / (xs[i+1] - xs[i]);
            return rhos[i] + w * (rhos[i+1] - rhos[i]);
        }
    }
    return rhos[n-1];
}

double particleVolumeFromMoments(double N, double mM, double mO,
                                 const ConfigOptions& opts)
{
    if (!(N > 0) || !(mM > 0)) {
        return 0.0;
    }
    double m = mM + std::max(mO, 0.0);
    // Avogadro's number here is Cantera's kmol-based constant, consistent
    // with computeNucleationRates() in sourceSystem.cpp.
    return m / (N * Cantera::Avogadro) / particleBulkDensity(mM, mO, opts); // [m^3]
}

double particleDiameterFromVolume(double vParticle)
{
    if (!(vParticle > 0)) {
        return 0.0;
    }
    return std::cbrt(6.0 * vParticle / M_PI);
}

double computeParticleDiameter(double N, double mM, double mO,
                               const ConfigOptions& opts)
{
    return particleDiameterFromVolume(particleVolumeFromMoments(N, mM, mO, opts));
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
