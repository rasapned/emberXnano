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

namespace {
// Phase segment [i, i+1] bracketing the O/metal ratio x, and the lever-rule
// weight w of phase i+1 (w = 0 with i = last phase if there is only one).
void phaseSegment(double x, const ConfigOptions& opts, size_t& i, double& w)
{
    const std::vector<double>& xs = opts.phaseRatio;
    size_t n = xs.size();
    i = 0;
    w = 0.0;
    if (n < 2) {
        return;
    }
    x = std::min(std::max(x, xs[0]), xs[n-1]);
    while (i + 2 < n && x > xs[i+1]) {
        i++;
    }
    w = (x - xs[i]) / (xs[i+1] - xs[i]);
}

double phaseTableValue(size_t phase, size_t n, const ConfigOptions& opts)
{
    return opts.phaseEnthalpy[phase * opts.particleThermoT.size() + n];
}
}

double particleEnthalpyPerMetal(double T, double x, const ConfigOptions& opts,
                                double* dhdx)
{
    const std::vector<double>& Ts = opts.particleThermoT;
    size_t nT = Ts.size();
    T = std::min(std::max(T, Ts[0]), Ts[nT-1]);
    size_t n = std::min(static_cast<size_t>(
        std::upper_bound(Ts.begin(), Ts.end(), T) - Ts.begin()), nT - 1);
    n = std::max(n, static_cast<size_t>(1)) - 1; // Ts[n] <= T <= Ts[n+1]
    double s = (T - Ts[n]) / (Ts[n+1] - Ts[n]);

    size_t i;
    double w;
    phaseSegment(x, opts, i, w);
    double hi = (1 - s) * phaseTableValue(i, n, opts) + s * phaseTableValue(i, n+1, opts);
    if (opts.phaseRatio.size() < 2) {
        if (dhdx) {
            *dhdx = 0.0;
        }
        return hi;
    }
    double hj = (1 - s) * phaseTableValue(i+1, n, opts) + s * phaseTableValue(i+1, n+1, opts);
    if (dhdx) {
        *dhdx = (hj - hi) / (opts.phaseRatio[i+1] - opts.phaseRatio[i]);
    }
    return (1 - w) * hi + w * hj;
}

double particleTemperature(double H, double mM, double mO, double Tgas,
                           const ConfigOptions& opts)
{
    if (!(mM > 0)) {
        return Tgas;
    }
    const std::vector<double>& Ts = opts.particleThermoT;
    size_t nT = Ts.size();
    double hTarget = H / (mM / opts.metalWeight); // [J per kmol metal]

    size_t i;
    double w;
    phaseSegment(particleOxygenRatio(mM, mO, opts), opts, i, w);
    size_t i2 = (opts.phaseRatio.size() < 2) ? i : i + 1;
    auto h = [&](size_t n) {
        return (1 - w) * phaseTableValue(i, n, opts) + w * phaseTableValue(i2, n, opts);
    };

    // h(T) increases monotonically (the jumps at phase transitions are
    // resolved over one table interval), so bisect on the table index
    if (hTarget <= h(0)) {
        return Ts[0];
    } else if (hTarget >= h(nT-1)) {
        return Ts[nT-1];
    }
    size_t lo = 0, hi = nT - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (h(mid) <= hTarget) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return Ts[lo] + (hTarget - h(lo)) / (h(hi) - h(lo)) * (Ts[hi] - Ts[lo]);
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
