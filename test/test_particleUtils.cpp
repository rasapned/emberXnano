#include "../src/particleUtils.h"
#include "../src/readConfig.h"

#include "gtest/gtest.h"
#include <cmath>

namespace {
const double kB = 1.380649e-23; // [J/K]
const double T = 1500.0; // [K]
const double mu = 6e-5; // gas viscosity [Pa*s]

double thermalSpeed(double m) { return std::sqrt(8 * kB * T / (M_PI * m)); }
}

// Free-molecular limit (tiny particles, long mean free path):
// beta -> pi/4 (d_i+d_j)^2 sqrt(c_i^2 + c_j^2)
TEST(ParticleUtils, FuchsKernelFreeMolecularLimit)
{
    double lambda = 1e-3; // [m], very low pressure
    double di = 1e-9, dj = 2e-9;
    double mi = 1e-24, mj = 8e-24;
    double Di = particleDiffusivityFromDiameter(di, T, mu, lambda);
    double Dj = particleDiffusivityFromDiameter(dj, T, mu, lambda);
    double beta = fuchsCoagulationKernel(di, mi, Di, dj, mj, Dj, T);
    double ci = thermalSpeed(mi), cj = thermalSpeed(mj);
    double betaFM = M_PI / 4 * std::pow(di + dj, 2) * std::sqrt(ci * ci + cj * cj);
    EXPECT_NEAR(beta / betaFM, 1.0, 1e-3);
}

// Continuum limit (large particles, short mean free path):
// beta -> 2 pi (D_i+D_j)(d_i+d_j)
TEST(ParticleUtils, FuchsKernelContinuumLimit)
{
    // the corrections decay only as (lambda/d)^(1/2): mm-sized particles
    double lambda = 1e-9; // [m]
    double di = 1e-3, dj = 3e-3;
    double rhoP = 5000.0;
    double mi = rhoP * M_PI / 6 * std::pow(di, 3);
    double mj = rhoP * M_PI / 6 * std::pow(dj, 3);
    double Di = particleDiffusivityFromDiameter(di, T, mu, lambda);
    double Dj = particleDiffusivityFromDiameter(dj, T, mu, lambda);
    double beta = fuchsCoagulationKernel(di, mi, Di, dj, mj, Dj, T);
    double betaC = 2 * M_PI * (Di + Dj) * (di + dj);
    EXPECT_NEAR(beta / betaC, 1.0, 3e-3);
}

TEST(ParticleUtils, FuchsKernelSymmetric)
{
    double lambda = 1e-6;
    double Di = particleDiffusivityFromDiameter(1e-9, T, mu, lambda);
    double Dj = particleDiffusivityFromDiameter(5e-9, T, mu, lambda);
    double a = fuchsCoagulationKernel(1e-9, 1e-24, Di, 5e-9, 1e-22, Dj, T);
    double b = fuchsCoagulationKernel(5e-9, 1e-22, Dj, 1e-9, 1e-24, Di, T);
    EXPECT_DOUBLE_EQ(a, b);
}

// A sphere has one primary particle; n spheres of diameter d in contact
// (area n pi d^2, volume n pi d^3 / 6) have n primaries of diameter d
TEST(ParticleUtils, AggregateGeometry)
{
    double d = 2e-9, Df = 1.8;
    double v1 = M_PI / 6 * std::pow(d, 3), a1 = M_PI * d * d;
    double dpp, np, dc;
    aggregateGeometry(v1, a1, Df, dpp, np, dc);
    EXPECT_NEAR(dpp / d, 1.0, 1e-12);
    EXPECT_NEAR(np, 1.0, 1e-12);
    EXPECT_NEAR(dc / d, 1.0, 1e-12);

    aggregateGeometry(10 * v1, 10 * a1, Df, dpp, np, dc);
    EXPECT_NEAR(dpp / d, 1.0, 1e-12);
    EXPECT_NEAR(np, 10.0, 1e-9);
    EXPECT_NEAR(dc / (d * std::pow(10.0, 1 / Df)), 1.0, 1e-12);

    // Area below the sphere's: taken as the volume-equivalent sphere
    aggregateGeometry(8 * v1, 0.0, Df, dpp, np, dc);
    EXPECT_NEAR(dpp / (2 * d), 1.0, 1e-12);
    EXPECT_NEAR(np, 1.0, 1e-12);
}

// Every coagulation event removes exactly one particle in total, and the
// monomer units leaving mode 1 all arrive in mode 2 (volume conservation)
TEST(ParticleUtils, BimodalCoagulationConservation)
{
    double n1 = 3e17, n2 = 2e15, b11 = 1e-15, b12 = 4e-15, b22 = 2e-14;
    for (double r : {1.5, 2.0, 7.0, 1e4}) {
        BimodalCoagulationRates c = bimodalCoagulation(n1, n2, r, b11, b12, b22);
        double events = 0.5 * b11 * n1 * n1 + b12 * n1 * n2 + 0.5 * b22 * n2 * n2;
        double dTotal = -c.n1Loss + c.n2Gain - c.n2Loss;
        EXPECT_NEAR(dTotal / events, -1.0, 1e-12) << "r = " << r;
        EXPECT_DOUBLE_EQ(c.monomersToMode2, c.n1Loss);
        // the monomer units in the new mode-2 particles: r per particle
        double rr = std::max(r, 2.0);
        double fromPairs = c.n2Gain * rr; // volume (in v1) of new particles
        EXPECT_NEAR(fromPairs + b12 * n1 * n2, c.monomersToMode2, 1e-9 * c.n1Loss);
    }
}

// Constant kernel: the total number follows dn/dt = -beta/2 n^2 for any
// split between the modes (Smoluchowski)
TEST(ParticleUtils, BimodalCoagulationConstantKernel)
{
    double beta = 3e-15, n1 = 1e17, n2 = 4e16;
    BimodalCoagulationRates c = bimodalCoagulation(n1, n2, 12.0, beta, beta, beta);
    double dTotal = -c.n1Loss + c.n2Gain - c.n2Loss;
    EXPECT_NEAR(dTotal / (-0.5 * beta * std::pow(n1 + n2, 2)), 1.0, 1e-12);
}

TEST(ParticleUtils, SinteringTime)
{
    // Kruis et al. (1993), TiO2
    double tau = sinteringTime(1500.0, 5e-9, 7.44e16, 31000.0);
    EXPECT_NEAR(tau / (7.44e16 * 1500.0 * std::pow(5e-9, 4) * std::exp(31000.0 / 1500.0)),
                1.0, 1e-12);
}

TEST(ParticleUtils, FractalCollisionDiameter)
{
    ParticleOptions p;
    p.fractalDimension = 1.8;
    EXPECT_DOUBLE_EQ(fractalCollisionDiameter(3e-9, p), 3e-9); // no primaries
    p.primaryDiameter = 1e-9;
    EXPECT_DOUBLE_EQ(fractalCollisionDiameter(0.5e-9, p), 0.5e-9); // below d_pp
    EXPECT_NEAR(fractalCollisionDiameter(3e-9, p) / (1e-9 * std::pow(27.0, 1 / 1.8)),
                1.0, 1e-12);
}
