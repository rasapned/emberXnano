#include "particleUtils.h"
#include "chemistry0d.h"
#include "readConfig.h"

#include <algorithm>
#include <cmath>

double particleOxygenRatio(double mM, double mO, const ParticleOptions& opts)
{
    if (!(mM > 0) || !(mO > 0)) {
        return 0.0;
    }
    return (mO / opts.oxygenWeight) / (mM / opts.metalWeight);
}

double particleBulkDensity(double mM, double mO, const ParticleOptions& opts)
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
                                 const ParticleOptions& opts)
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
                               const ParticleOptions& opts)
{
    return particleDiameterFromVolume(particleVolumeFromMoments(N, mM, mO, opts));
}

void floorParticleSize(double& dp, double& nParticles, double* mp,
                       const ParticleOptions& opts)
{
    double dmin = opts.minParticleDiameter;
    if (!(dp < dmin)) {
        return;
    }
    double r3 = std::pow(dp / dmin, 3); // same density, so mass ~ dp^3
    nParticles *= r3;
    if (mp) {
        *mp /= r3;
    }
    dp = dmin;
}

namespace {
// Phase segment [i, i+1] bracketing the O/metal ratio x, and the lever-rule
// weight w of phase i+1 (w = 0 with i = last phase if there is only one).
void phaseSegment(double x, const ParticleOptions& opts, size_t& i, double& w)
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

double phaseTableValue(size_t phase, size_t n, const ParticleOptions& opts)
{
    return opts.phaseEnthalpy[phase * opts.particleThermoT.size() + n];
}
}

double particleEnthalpyPerMetal(double T, double x, const ParticleOptions& opts,
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

double phaseActivity(double x, int phase, const ParticleOptions& opts)
{
    if (opts.phaseRatio.size() < 2) {
        return (phase == 0) ? 1.0 : 0.0;
    }
    size_t i;
    double w;
    phaseSegment(x, opts, i, w);
    if (phase == static_cast<int>(i)) {
        return 1.0 - w;
    } else if (phase == static_cast<int>(i) + 1) {
        return w;
    }
    return 0.0;
}

namespace {
// exp of a ln(p) table on opts.particleThermoT, interpolated linearly at T
// (T clamped to the table range)
double interpolateLnTable(const double* lnP, double T, const ParticleOptions& opts)
{
    const std::vector<double>& Ts = opts.particleThermoT;
    size_t nT = Ts.size();
    T = std::min(std::max(T, Ts[0]), Ts[nT-1]);
    size_t n = std::min(static_cast<size_t>(
        std::upper_bound(Ts.begin(), Ts.end(), T) - Ts.begin()), nT - 1);
    n = std::max(n, static_cast<size_t>(1)) - 1; // Ts[n] <= T <= Ts[n+1]
    double s = (T - Ts[n]) / (Ts[n+1] - Ts[n]);
    return exp((1 - s) * lnP[n] + s * lnP[n+1]);
}
} // namespace

double condensationSaturationPressure(size_t c, double T, const ParticleOptions& opts)
{
    return interpolateLnTable(&opts.condLnPsat[c * opts.particleThermoT.size()], T, opts);
}

double nucleusEquilibriumPressure(int nAtoms, double T, const ParticleOptions& opts)
{
    double v1 = opts.nucMonomerVolume;
    double dn = std::cbrt(6.0 * nAtoms * v1 / M_PI); // nucleus diameter [m]
    double sigma = opts.surfaceTensionA - opts.surfaceTensionB * T;
    return interpolateLnTable(opts.nucLnPsatMonomer.data(), T, opts) *
           exp(4.0 * sigma * v1 / (dn * Cantera::Boltzmann * T));
}

double particleTemperature(double H, double mM, double mO, double Tgas,
                           const ParticleOptions& opts)
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
    // Guards against division by zero only (OpenFOAM's ROOTVSMALL for double).
    // It must stay far below the physical values: 1e-18 exceeded mp
    // (~1e-24 kg) and dc^2 (~1e-19 m^2) and cut the kernel ~300x.
    const double ROOTVSMALL = 1e-150;

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

double particleDiffusivityFromDiameter(double dp, double T, double mu,
                                       double meanFreePath)
{
    double Kn = 2.0 * meanFreePath / dp;
    double Cc = cunninghamSlipCorrection(Kn);
    return Cantera::Boltzmann * T * Cc / (3.0 * M_PI * mu * dp);
}

double fuchsCoagulationKernel(double di, double mi, double Di,
                              double dj, double mj, double Dj, double T)
{
    const double kB_ = Cantera::Boltzmann;
    // mean thermal speed, apparent mean free path and transition length
    auto transition = [&](double d, double m, double D, double& c) {
        c = std::sqrt(8.0 * kB_ * T / (M_PI * m));
        double l = 8.0 * D / (M_PI * c);
        return (std::pow(d + l, 3) - std::pow(d * d + l * l, 1.5)) /
               (3.0 * d * l) - d;
    };
    double ci, cj;
    double gi = transition(di, mi, Di, ci);
    double gj = transition(dj, mj, Dj, cj);
    double dSum = di + dj;
    double DSum = Di + Dj;
    double gij = std::sqrt(gi * gi + gj * gj);
    double cij = std::sqrt(ci * ci + cj * cj);
    return 2.0 * M_PI * DSum * dSum /
           (dSum / (dSum + 2.0 * gij) + 8.0 * DSum / (cij * dSum));
}

void aggregateGeometry(double v, double a, double Df,
                       double& dpp, double& np, double& dc)
{
    double aSphere = std::pow(M_PI, 1.0 / 3.0) * std::pow(6.0 * v, 2.0 / 3.0);
    a = std::max(a, aSphere);
    dpp = 6.0 * v / a;
    np = std::max(a * a * a / (36.0 * M_PI * v * v), 1.0);
    dc = dpp * std::pow(np, 1.0 / Df);
}

BimodalCoagulationRates bimodalCoagulation(double n1, double n2, double r,
                                           double b11, double b12, double b22)
{
    r = std::max(r, 2.0);
    double R11 = 0.5 * b11 * n1 * n1; // monomer-monomer collisions
    double R12 = b12 * n1 * n2; // monomer-aggregate collisions
    double R22 = 0.5 * b22 * n2 * n2; // aggregate-aggregate collisions
    BimodalCoagulationRates rates;
    rates.n1Loss = R11 * r / (r - 1.0) + R12;
    rates.n2Gain = R11 / (r - 1.0);
    rates.n2Loss = R22;
    rates.monomersToMode2 = rates.n1Loss;
    return rates;
}

double sinteringTime(double T, double dpp, double A, double Ta)
{
    return A * T * std::pow(dpp, 4) * std::exp(Ta / T);
}

double fractalCollisionDiameter(double dv, const ParticleOptions& opts)
{
    double dpp = opts.primaryDiameter;
    if (dpp > 0 && dv > dpp) {
        return dpp * std::pow(dv / dpp, 3.0 / opts.fractalDimension);
    }
    return dv;
}

bool bimodalAggregate(double N2, double M2, double A2, double rho,
                      const ParticleOptions& opts, double& mp, double& v,
                      double& dpp, double& np, double& dc)
{
    double nPerKg = N2 * Cantera::Avogadro; // [aggregates / kg gas]
    if (!(nPerKg * rho >= 1.0 && M2 > 0)) {
        mp = v = dpp = dc = 0.0;
        np = 1.0;
        return false;
    }
    const double m1 = opts.monomerMass, d1 = opts.monomerDiameter;
    mp = std::max(M2 / nPerKg, 2.0 * m1);
    v = mp / opts.phaseDensity[0];
    double aMax = mp / m1 * M_PI * d1 * d1; // unsintered monomers
    double a = std::min(std::max(A2, 0.0) / nPerKg, aMax);
    aggregateGeometry(v, a, opts.fractalDimension, dpp, np, dc);
    return true;
}
