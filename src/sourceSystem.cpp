#include "sourceSystem.h"
#include "readConfig.h"
#include "perfTimer.h"
#include "sundialsUtils.h"
#include "chemistry0d.h"
#include "scalarFunction.h"
#include "particleUtils.h"

#include <algorithm>
#include <boost/format.hpp>

SourceSystem::SourceSystem()
    : U(NaN)
    , T(NaN)
    , Tp(NaN)
    , debug(false)
    , options(NULL)
    , gas(NULL)
    , strainFunction(NULL)
    , rateMultiplierFunction(NULL)
    , heatLoss(NULL)
    , pop(NULL)
    , pn(0)
    , qLoss(0.0)
    , quasi2d(false)
{
}

void SourceSystem::updateThermo()
{
    thermoTimer->start();
    gas->getEnthalpies(hk);
    rho = gas->getDensity();
    cp = gas->getSpecificHeatCapacity();
    thermoTimer->stop();
}

double SourceSystem::getQdotIgniter(double t)
{
    ConfigOptions& opt = *options;
    if (t >= opt.ignition_tStart &&
        t < opt.ignition_tStart + opt.ignition_duration) {
        return opt.ignition_energy /
                (opt.ignition_stddev * sqrt(2 * M_PI) * opt.ignition_duration) *
                exp(-pow(x - opt.ignition_center, 2) /
                    (2 * pow(opt.ignition_stddev, 2)));
    } else {
        return 0.0;
    }
}

void SourceSystem::setOptions(ConfigOptions& opts)
{
    options = &opts;
}

void SourceSystem::initialize(size_t new_nSpec, size_t new_nMoments)
{
    nSpec = new_nSpec;
    nMoments = new_nMoments;

    Y.setConstant(nSpec, NaN);
    moments.setConstant(nMoments, NaN);
    cpSpec.resize(nSpec);
    splitConst.resize(nSpec + 2 + nMoments);
    hk.resize(nSpec);

    W.resize(gas->nSpec); // move this to initialize
    gas->getMolecularWeights(W);
}

void SourceSystem::computeNucleationRates
(dvec& momentsQ, dvec& momentsD, dvec& speciesQ, dvec& speciesD)
{
    // With the bimodal model the moments here are the single-mode layout
    // (N, mM, mO) of the nuclei (see computeBimodalRates), without H
    const bool energy = pn > kH && !pop->bimodal;

    // --- Collision-based nucleation channels ----------------------------
    // The hard-sphere kinetic collision prefactor (diameter/reduced-mass
    // terms) is precomputed once in Python; only sqrt(T) and the
    // concentrations need to be evaluated here.
    double pMonomer = 0.0; // partial pressure of the metal monomer [Pa]
    if (pop->nucKelvinLimit) {
        size_t k1 = pop->nucMonomerSpecies;
        pMonomer = Y[k1] * rho / W[k1] * Cantera::GasConstant * T;
    }
    for (size_t c = 0; c < pop->nucSpeciesA.size(); c++) {
        size_t kA = pop->nucSpeciesA[c];
        size_t kB = pop->nucSpeciesB[c];
        double concA = Y[kA] * rho / W[kA]; // [kmol/m^3]
        double concB = Y[kB] * rho / W[kB]; // [kmol/m^3]

        // Calculate the nucleation rate (prefactor has a unit of m^3/s/sqrt(K)/kmol (already * avogadro))
        // Unit of J: kmol/m^3/s
        double J = pop->nucEfficiency * pop->nucCollisionPrefactor[c] *
                   sqrt(T) * concA * concB;
        // Halve the nucleation for like molecules
        if (kA == kB) {
            // Each A-A collision is counted once, not twice.
            J *= 0.5;
        }
        // Only nuclei above the critical size grow: net monomer uptake of the
        // nucleus ~ (p_1 - p_eq), so no nucleation for p_1 <= p_eq
        if (pop->nucKelvinLimit) {
            double pEq = nucleusEquilibriumPressure(pop->nucAtomsPerEvent[c], T, *pop);
            J *= (pMonomer > pEq) ? 1.0 - pEq / pMonomer : 0.0;
        }
        J = std::max(J, 0.0);

        // Particle moments: pure production
        momentsQ[kN] += J / rho;
        double massO = pop->nucOxygenMassPerEvent[c];
        momentsQ[kM] += J / rho * (pop->nucMassPerEvent[c] - massO);
        if (massO != 0) {
            momentsQ[kO] += J / rho * massO;
        }

        // Gas species consumption: pure destruction
        speciesD[kA] += pop->nucStoichA[c] * J * W[kA] / rho;
        speciesD[kB] += pop->nucStoichB[c] * J * W[kB] / rho;

        // The new particle takes the enthalpy of the consumed molecules
        if (energy) {
            momentsQ[kH] += J * (pop->nucStoichA[c] * hk[kA] +
                                 pop->nucStoichB[c] * hk[kB]) / rho;
        }
    }

    // --- Monomer conversion ------------------------------------------------
    // Every molecule of the species becomes one particle (first order)
    for (size_t c = 0; c < pop->convSpecies.size(); c++) {
        size_t k = pop->convSpecies[c];
        double J = std::max(pop->convRate * Y[k] * rho / W[k], 0.0); // [kmol/m^3/s]
        double massO = pop->convOxygenMass[c];
        momentsQ[kN] += J / rho;
        momentsQ[kM] += J / rho * (W[k] - massO);
        if (massO != 0) {
            momentsQ[kO] += J / rho * massO;
        }
        speciesD[k] += J * W[k] / rho;
        if (energy) {
            momentsQ[kH] += J * hk[k] / rho;
        }
    }

    // --- Classical nucleation theory (CNT) -------------------------------
    // Monomer attachment only (no pairwise cluster collisions). precursorConc
    // still sums all cluster sizes toward the monomer-equivalent
    // supersaturation, but the rate and mass balance are referenced to the
    // single true monomer species (nucMonomerIndex, atom count == 1).
    size_t nPrecursors = pop->nucPrecursorSpecies.size();
    if (nPrecursors == 0 || pop->nucMonomerIndex < 0) {
        return;
    }

    const double kB_ = Cantera::Boltzmann;
    const double NA_ = Cantera::Avogadro;
    const double Ru_ = Cantera::GasConstant;

    // Saturation vapor pressure (Antoine equation, T in degC, P in bar).
    double Tc = T - 273.15;
    double pSat = pow(10.0, pop->nucAntoineA -
                      pop->nucAntoineB / (Tc + pop->nucAntoineC)) * 1e5; // [Pa]

    // Monomer-equivalent partial pressure: sum of precursor concentrations
    // weighted by their cluster size (atom count).
    double pMon = 0.0;
    for (size_t i = 0; i < nPrecursors; i++) {
        size_t k = pop->nucPrecursorSpecies[i];
        double conc = Y[k] * rho / W[k]; // [kmol/m^3]
        pMon += conc * 1000.0 * pop->nucPrecursorAtomCount[i] * Ru_ * T;
    }
    pMon = std::max(pMon, 0.0);
    double S = pMon / (pSat + 1e-300); // supersaturation ratio
    if (S <= 1.0) {
        return; // no nucleation below saturation
    }

    // Monomer mass/volume, derived from the mechanism (not user-specified),
    // using the pure-metal phase density.
    size_t kMon = pop->nucPrecursorSpecies[pop->nucMonomerIndex];
    double n_mon = pMon / (kB_ * T); // monomer-equivalent number density [1/m^3]
    double m1 = W[kMon] / NA_; // monomer mass [kg] -- unrelated to M_PI
    double v1 = W[kMon] / (NA_ * pop->phaseDensity[0]); // monomer volume [m^3]

    // Surface tension and dimensionless surface energy parameter Theta,
    // referenced to a spherical monomer of volume v1.
    double sigma = pop->nucSurfaceTensionA + pop->nucSurfaceTensionB * T;
    double d1 = pow(6.0 * v1 / M_PI, 1.0 / 3.0);
    double theta = M_PI * d1 * d1 * sigma / (kB_ * T);

    // Critical cluster size (classical nucleation theory)
    double lnS = log(S);
    double gStar = pow(2.0 * theta / (3.0 * lnS), 3.0);

    // Nucleation rate
    double J = n_mon * n_mon * S * v1 * pow(2.0 * sigma / M_PI / m1, 0.5) *
              exp(theta - (4.0*theta*theta*theta/27.0/(lnS*lnS)));
    J = std::max(J, 0.0);

    momentsQ[kN] += J / rho;
    momentsQ[kM] += J * gStar * W[kMon] / rho;
    speciesD[kMon] += J * gStar * W[kMon] / rho;
    if (energy) {
        momentsQ[kH] += J * gStar * hk[kMon] / rho;
    }
}

void SourceSystem::computeCoagulationRates(dvec& momentsD)
{
    if (pn < 2 || !pop->coagulation) {
        // Both N (kN) and mM (kM) are required to define a particle size.
        return;
    }

    double N = pm[kN]; // [kmol particles / kg gas]
    double mM = pm[kM]; // [kg particle-phase metal / kg gas]
    double mO = (pn > kO) ? pm[kO] : 0.0; // [kg particle-phase O / kg gas]

    double vParticle = particleVolumeFromMoments(N, mM, mO, *pop);
    if (!(vParticle > 0)) {
        return; // no particles present -- nothing to coagulate
    }
    double dc = particleDiameterFromVolume(vParticle); // [m]
    double mp = (mM + std::max(mO, 0.0)) / (N * Cantera::Avogadro); // single-particle mass [kg]
    dc = fractalCollisionDiameter(dc, *pop);

    double Wmx_ = gas->getMixtureMolecularWeight(); // [kg/kmol]
    double mu_ = gas->getViscosity(); // [Pa*s]
    double meanFreePath = gasMeanFreePath(mu_, T, Wmx_, gas->pressure);

    double beta = pop->coagulationEnhancement *
                  monodisperseCoagulationKernel(dc, mp, T, mu_, meanFreePath);

    // Monodisperse Brownian coagulation halves the physical number density
    // at rate dn/dt = -0.5*beta*n^2 (n = N*rho*Avogadro [particles/m^3]),
    // while leaving the total particle mass (mM, mO) unchanged --
    // particles merge, matter is conserved. Converting back to the
    // "kmol particles per kg gas" state-variable convention used for kN
    // (rho and one power of Avogadro cancel; see the unit note above):
    double dNdt_coag = 0.5 * N * N * rho * Cantera::Avogadro * beta;
    momentsD[kN] += dNdt_coag;

    if (debug) {
        logFile.write(format(
            "coagulation: j=%i x=%.4g | N=%.4e[kmol/kg] dc=%.4e[m] "
            "mp=%.4e[kg] meanFreePath=%.4e[m] beta=%.4e[m^3/s] | "
            "dN/dt=-%.4e[kmol/kg/s]") %
            j % x % N % dc % mp % meanFreePath % beta % dNdt_coag);
    }
}

void SourceSystem::computeSurfaceReactionRates
(dvec& momentsQ, dvec& momentsD, dvec& speciesQ, dvec& speciesD)
{
    size_t nReactions = pop->surfReactant.size();
    if (nReactions == 0 || pn < 2) {
        return;
    }

    double N = pm[kN]; // [kmol particles / kg gas]
    double mM = pm[kM]; // [kg particle-phase metal / kg gas]
    double mO = (pn > kO) ? pm[kO] : 0.0; // [kg particle-phase O / kg gas]

    double vParticle = particleVolumeFromMoments(N, mM, mO, *pop);
    if (!(vParticle > 0)) {
        return; // no particles present -- no surface to react on
    }
    double dp = particleDiameterFromVolume(vParticle); // [m]
    double mp = (mM + std::max(mO, 0.0)) / (N * Cantera::Avogadro); // single-particle mass [kg]
    double nParticles = N * rho * Cantera::Avogadro; // [particles / m^3]
    floorParticleSize(dp, nParticles, &mp, *pop);
    double area = M_PI * dp * dp * nParticles; // particle surface area [m^2 / m^3]

    // The surface composition is assumed equal to the particle composition:
    // atomic fractions of metal and O sites, plus a gradual cap on oxidation
    // that vanishes at the most oxidized phase (O/metal = xMax).
    double xO = particleOxygenRatio(mM, mO, *pop);
    double XM = 1.0 / (1.0 + xO);
    double XO = 1.0 - XM;
    double xMax = pop->phaseRatio.back();
    double capFactor = (xMax > 0) ? std::max(1.0 - xO / xMax, 0.0) : 0.0;

    const double WM = pop->metalWeight;
    const double WO = pop->oxygenWeight;

    // Prescribed heat: the atoms gained or lost carry the particle's own
    // enthalpy at Tp, from H = nM * h(x, Tp) with x = nO/nM:
    //     dH = (h - x*dh/dx) * dnM + dh/dx * dnO
    double hP = 0.0, dhdx = 0.0;
    if (pn > kH && pop->surfPrescribedHeat) {
        hP = particleEnthalpyPerMetal(Tp, xO, *pop, &dhdx);
    }

    double massLoss = 0.0; // particle mass removal rate [kg / kg gas / s]
    for (size_t r = 0; r < nReactions; r++) {
        int dM = pop->surfDeltaM[r];
        int dO = pop->surfDeltaO[r];

        // Fraction of the surface available to this reaction
        double theta;
        if (dM < 0) {
            theta = XM; // etching: needs metal sites
        } else if (dO < 0) {
            theta = XO; // reduction: needs O sites
        } else if (dM == 0) {
            theta = XM; // oxidation: binds to metal sites
        } else {
            theta = 1.0; // condensation of metal-bearing species
        }
        // Gradual cap at the most oxidized phase for any reaction that would
        // raise O/metal beyond xMax: adding material richer in O, or removing
        // material richer in metal
        bool raisesO = (dM > 0) ? (dO > xMax * dM)
                     : (dM < 0) ? (-dO < -xMax * dM)
                     : (dO > 0);
        if (raisesO) {
            theta *= capFactor;
        }

        size_t k = pop->surfReactant[r];
        double conc = Y[k] * rho / W[k]; // [kmol/m^3]
        // Arrhenius factor at the particle (surface) temperature
        double kSurf = pop->surfA[r] * exp(-pop->surfTa[r] / Tp); // [m/s]
        if (pop->surfCollisionLimit) {
            // At most one reaction per collision on an available site:
            // hard-sphere collision rate per unit particle area (as for
            // condensation), at the gas temperature
            double mk = W[k] / Cantera::Avogadro; // molecular mass [kg]
            double mu = mk * mp / (mk + mp); // reduced mass [kg]
            double ratio = 1.0 + pop->surfDiameter[r] / dp;
            double kColl = 0.25 * ratio * ratio *
                           sqrt(8.0 * Cantera::Boltzmann * T / (M_PI * mu)); // [m/s]
            kSurf = std::min(kSurf, kColl);
        }
        double R = kSurf * theta * area * conc; // [kmol/m^3/s]
        R = std::max(R, 0.0);

        // Particle composition change. Losses go to the destruction terms,
        // which vanish with theta as the corresponding element runs out.
        double dmM = dM * WM * R / rho; // [kg / kg gas / s]
        double dmO = dO * WO * R / rho;
        if (dmM >= 0) {
            momentsQ[kM] += dmM;
        } else {
            momentsD[kM] -= dmM;
            massLoss -= dmM;
        }
        if (dO != 0) {
            if (dmO >= 0) {
                momentsQ[kO] += dmO;
            } else {
                momentsD[kO] -= dmO;
                massLoss -= dmO;
            }
        }

        // Gas species: reactant consumed, product released. By construction
        // W[reactant] - W[product] == dM*WM + dO*WO, so mass is conserved.
        speciesD[k] += W[k] * R / rho;
        int p = pop->surfProduct[r];
        if (p >= 0) {
            speciesQ[p] += W[p] * R / rho;
        }

        if (pn > kH) {
            // The reacting molecule sticks, so it brings its whole sensible
            // heat h(T) - h(Tp); conduction (computeParticleHeatTransfer)
            // already transfers the fraction alpha of it for every collision,
            // so only the remaining (1 - alpha) is added here
            double hReactantTp = gas->getSpeciesEnthalpy(k, Tp);
            double sensible = (1.0 - pop->thermalAccommodation) *
                              (hk[k] - hReactantTp); // [J/kmol]
            if (pop->surfPrescribedHeat) {
                // Prescribed heat -surfDeltaH only (as OpenFOAM): the atoms
                // gained or lost carry their enthalpy in the particle at Tp
                double dHcomposition = (hP - xO * dhdx) * dM + dhdx * dO; // [J/kmol]
                momentsQ[kH] += R * (dHcomposition - pop->surfDeltaH[r] + sensible) / rho;
            } else {
                // Reaction heat from the thermochemistry at Tp: the particle
                // receives the absolute enthalpy of the reactant and loses that
                // of the product. H being the absolute particle enthalpy, the
                // condensed-phase tables that set Tp from H account for the
                // atoms attached to or removed from the particle.
                double hProduct = (p >= 0) ? gas->getSpeciesEnthalpy(p, Tp) : 0.0;
                momentsQ[kH] += R * (hReactantTp - hProduct + sensible) / rho;
            }
        }
    }

    double dNdt_dis = addDisintegration(massLoss, dp, momentsD);

    if (debug) {
        logFile.write(format(
            "surface reactions: j=%i x=%.4g | N=%.4e[kmol/kg] mM=%.4e mO=%.4e[kg/kg] "
            "dp=%.4e[m] O/M=%.4g | massLoss=%.4e[kg/kg/s] dN/dt=-%.4e[kmol/kg/s]") %
            j % x % N % mM % mO % dp % xO % massLoss % dNdt_dis);
    }
}

double SourceSystem::addDisintegration(double massLoss, double dp, dvec& momentsD)
{
    // Particle disintegration: once particles have shrunk to the threshold
    // size, removing mass also removes particles, so that the mean particle
    // mass (and size) stays fixed instead of shrinking further. The switch is
    // ramped linearly over [dmin, 1.1*dmin]: with a hard switch at dmin the
    // integrator chatters wherever growth and etching balance near dmin.
    double dmin = pop->minParticleDiameter;
    if (!(massLoss > 0) || dp >= 1.1 * dmin) {
        return 0.0;
    }
    double N = pm[kN];
    double mM = pm[kM];
    double mO = (pn > kO) ? pm[kO] : 0.0;
    double ramp = std::min((1.1 * dmin - dp) / (0.1 * dmin), 1.0);
    double dNdt_dis = ramp * N * massLoss / (mM + std::max(mO, 0.0));
    momentsD[kN] += dNdt_dis;
    return dNdt_dis;
}

void SourceSystem::computeCondensationRates(dvec& momentsQ, dvec& momentsD,
                                            dvec& speciesQ, dvec& speciesD)
{
    size_t nCond = pop->condSpecies.size();
    if (nCond == 0 || pn < 2) {
        return;
    }

    double N = pm[kN]; // [kmol particles / kg gas]
    double mM = pm[kM]; // [kg particle-phase metal / kg gas]
    double mO = (pn > kO) ? pm[kO] : 0.0; // [kg particle-phase O / kg gas]

    double vParticle = particleVolumeFromMoments(N, mM, mO, *pop);
    if (!(vParticle > 0)) {
        return; // no particles present -- nothing to condense on
    }
    const double NA_ = Cantera::Avogadro;
    double dp = particleDiameterFromVolume(vParticle); // [m]
    double mp = (mM + std::max(mO, 0.0)) / (N * NA_); // single-particle mass [kg]
    double nParticles = N * rho * NA_; // [particles / m^3]
    floorParticleSize(dp, nParticles, &mp, *pop);

    // Same gradual cap as for surface reactions: species richer in O than the
    // most oxidized phase stop condensing as the particle approaches it.
    double xO = particleOxygenRatio(mM, mO, *pop);
    double xMax = pop->phaseRatio.back();
    double capFactor = (xMax > 0) ? std::max(1.0 - xO / xMax, 0.0) : 0.0;

    double massLoss = 0.0; // particle mass evaporation rate [kg / kg gas / s]
    for (size_t c = 0; c < nCond; c++) {
        size_t k = pop->condSpecies[c];
        int dM = pop->condDeltaM[c];
        int dO = pop->condDeltaO[c];

        // Hard-sphere collision rate between molecules of k and particles
        double mk = W[k] / NA_; // molecular mass [kg]
        double mu = mk * mp / (mk + mp); // reduced mass [kg]
        double dkp = dp + pop->condDiameter[c];
        double beta = 0.25 * M_PI * dkp * dkp *
                      sqrt(8.0 * Cantera::Boltzmann * T / (M_PI * mu)); // [m^3/s]

        double conc = Y[k] * rho / W[k]; // [kmol/m^3]
        double R = std::max(beta * nParticles * conc, 0.0); // [kmol/m^3/s]
        if (dO > xMax * dM) {
            R *= capFactor;
        }

        // By construction W[k] == dM*WM + dO*WO, so mass is conserved.
        momentsQ[kM] += dM * pop->metalWeight * R / rho;
        if (dO != 0) {
            momentsQ[kO] += dO * pop->oxygenWeight * R / rho;
        }
        speciesD[k] += W[k] * R / rho;
        if (pn > kH) {
            // Enthalpy of the condensed molecules at Tp, plus the part
            // (1 - alpha) of their sensible heat not already transferred by
            // conduction (as for the surface reactions)
            double hkTp = gas->getSpeciesEnthalpy(k, Tp);
            momentsQ[kH] += R * (hkTp + (1.0 - pop->thermalAccommodation) *
                                 (hk[k] - hkTp)) / rho;
        }

        // Evaporation (detailed balance with the condensation above): the
        // collision rate at the saturation concentration over the bulk phase
        // at Tp, times the activity of that phase at the surface and the
        // Kelvin factor of the curved surface.
        int phase = pop->evaporation ? pop->condPhase[c] : -1;
        if (phase < 0) {
            continue;
        }
        double activity = phaseActivity(xO, phase, *pop);
        if (!(activity > 0)) {
            continue;
        }
        double sigma = pop->surfaceTensionA - pop->surfaceTensionB * Tp;
        // Kelvin factor of the curved surface; dp is at least
        // minParticleDiameter (floorParticleSize), otherwise it would grow
        // faster than d_p^2 shrinks and the rate would diverge as d_p -> 0.
        double kelvin = 4.0 * sigma * pop->condMolarVolume[c] /
                        (dp * Cantera::Boltzmann * Tp);
        double cSat = condensationSaturationPressure(c, Tp, *pop) /
                      (Cantera::GasConstant * Tp); // [kmol/m^3]
        double Rev = beta * nParticles * activity * cSat * exp(kelvin); // [kmol/m^3/s]

        double dmM = dM * pop->metalWeight * Rev / rho;
        double dmO = dO * pop->oxygenWeight * Rev / rho;
        momentsD[kM] += dmM;
        if (dO != 0) {
            momentsD[kO] += dmO;
        }
        massLoss += dmM + dmO;
        speciesQ[k] += W[k] * Rev / rho;
        if (pn > kH) {
            // The evaporated molecules leave with their enthalpy at Tp
            momentsQ[kH] -= Rev * gas->getSpeciesEnthalpy(k, Tp) / rho;
        }
    }

    addDisintegration(massLoss, dp, momentsD);
}

void SourceSystem::updateParticleTemperature()
{
    if (pn > kH && !pop->bimodal) {
        Tp = particleTemperature(pm[kH], pm[kM], pm[kO], T, *pop);
    } else {
        Tp = T;
    }
}

void SourceSystem::computeParticleHeatTransfer(dvec& momentsQ)
{
    if (pn <= kH) {
        return;
    }
    double N = pm[kN];
    double vParticle = particleVolumeFromMoments(N, pm[kM], pm[kO], *pop);
    if (!(vParticle > 0) || Tp == T) {
        return;
    }
    const double NA_ = Cantera::Avogadro;
    const double kB_ = Cantera::Boltzmann;
    double dp = particleDiameterFromVolume(vParticle); // [m]
    // Same mass in d_min particles below d_min, as for the surface sources
    double nPerKg = N * NA_; // [particles / kg gas]
    floorParticleSize(dp, nPerKg, nullptr, *pop);

    // Conduction: free-molecular limit (Filippov & Rosner 2000) with thermal
    // accommodation alpha, combined harmonically with the continuum limit
    double Wmx_ = gas->getMixtureMolecularWeight(); // [kg/kmol]
    double cBar = sqrt(8.0 * kB_ * T / (M_PI * Wmx_ / NA_)); // mean molecular speed [m/s]
    double gamma = cp / gas->thermo->cv_mass();
    double qFM = pop->thermalAccommodation * M_PI * dp * dp *
                 gas->pressure / 8.0 * cBar * (gamma + 1) / (gamma - 1) *
                 (Tp / T - 1.0); // [W]
    double qC = 2.0 * M_PI * dp * gas->getThermalConductivity() * (Tp - T); // [W]
    double q = qFM * qC / (qFM + qC);

    // Radiation, Rayleigh limit (d << wavelength): emission ~ d^3 T^5 with
    // the absorption function E(m), net exchange with surroundings at Tsurr
    const double c5 = 8.0 * 24.0 * 1.0369277551433699; // 4*2*Gamma(5)*zeta(5)
    double P = c5 * pow(M_PI, 3) * pow(dp, 3) * pop->radiationAbsorption *
               pow(kB_, 5) / (pow(Cantera::Planck, 4) * pow(Cantera::lightSpeed, 3)) *
               (pow(Tp, 5) - pow(pop->radiationTsurr, 5)); // [W]

    momentsQ[kH] -= nPerKg * (q + P);
}

void SourceSystem::computeBimodalRates(dvec& momentsQ, dvec& momentsD,
                                       dvec& speciesQ, dvec& speciesD)
{
    const double NA_ = Cantera::Avogadro;
    const double rhoP = pop->phaseDensity[0]; // particle material density
    const double m1 = pop->monomerMass; // mode-1 particle mass [kg]
    const double d1 = pop->monomerDiameter; // mode-1 particle diameter [m]

    // Nucleation into mode 1, evaluated in the single-mode layout (N, mM,
    // mO); the nucleated mass gives the number of nuclei of mass m1, so that
    // the gas mass consumed equals the mode-1 mass formed
    dvec q = dvec::Zero(3), d = dvec::Zero(3);
    computeNucleationRates(q, d, speciesQ, speciesD);
    momentsQ[kN1] += (q[kM] + q[kO]) / (m1 * NA_);

    // Mode-2 aggregates: mass and area per particle -> primary particle
    // diameter, number of primaries and collision diameter
    double N1 = std::max(pm[kN1], 0.0); // [kmol / kg gas]
    double N2 = pm[kN2];
    double n1 = N1 * rho * NA_; // [1/m^3]
    double mp2, v2, dpp, np, dc2;
    bool haveMode2 = bimodalAggregate(N2, pm[kM2], pm[kA2], rho, *pop,
                                      mp2, v2, dpp, np, dc2);
    double n2 = haveMode2 ? N2 * rho * NA_ : 0.0;

    double mu_ = gas->getViscosity(); // [Pa*s]
    double meanFreePath = gasMeanFreePath(mu_, T, gas->getMixtureMolecularWeight(),
                                          gas->pressure);

    // Brownian coagulation within and between the modes (Fuchs kernel with
    // the collision diameters)
    if (pop->coagulation && (n1 > 0 || n2 > 0)) {
        double D1 = particleDiffusivityFromDiameter(d1, T, mu_, meanFreePath);
        double f = pop->coagulationEnhancement;
        double b11 = f * fuchsCoagulationKernel(d1, m1, D1, d1, m1, D1, T);
        double b12 = 0.0, b22 = 0.0, r = 2.0;
        if (haveMode2) {
            double D2 = particleDiffusivityFromDiameter(dc2, T, mu_, meanFreePath);
            b12 = f * fuchsCoagulationKernel(d1, m1, D1, dc2, mp2, D2, T);
            b22 = f * fuchsCoagulationKernel(dc2, mp2, D2, dc2, mp2, D2, T);
            r = mp2 / m1;
        }
        BimodalCoagulationRates c = bimodalCoagulation(n1, n2, r, b11, b12, b22);
        double toMoments = 1.0 / (rho * NA_); // [1/m^3/s] -> [kmol/kg/s]
        momentsD[kN1] += c.n1Loss * toMoments;
        momentsQ[kN2] += c.n2Gain * toMoments;
        momentsD[kN2] += c.n2Loss * toMoments;
        // the monomers keep their mass and (unsintered) surface
        momentsQ[kM2] += c.monomersToMode2 * m1 / rho;
        momentsQ[kA2] += c.monomersToMode2 * M_PI * d1 * d1 / rho;
    }

    if (!haveMode2) {
        return;
    }

    // Sintering: the aggregate surface relaxes to that of the
    // volume-equivalent spheres. The relaxation acts from both sides: A2 <
    // As is unphysical, but the operator splitting can carry A2 below As
    // where tau is very short, and a one-sided sink would then leave it
    // there (and let it run away through the split constants).
    if (pop->sintering) {
        double aSphere = std::pow(M_PI, 1.0 / 3.0) * std::pow(6.0 * v2, 2.0 / 3.0);
        double As = N2 * NA_ * aSphere; // [m^2 / kg gas]
        double tau = sinteringTime(T, dpp, pop->sinterA, pop->sinterTa);
        double rate = (pm[kA2] - As) / tau; // [m^2 / kg gas / s]
        if (rate > 0) {
            momentsD[kA2] += rate;
        } else {
            momentsQ[kA2] -= rate;
        }
    }

    // Condensation onto the aggregates (hard-sphere collision rate with the
    // collision diameter); the surface grows as for spheres of the primary
    // diameter, dA = 4 dV / dpp
    for (size_t c = 0; c < pop->condSpecies.size(); c++) {
        size_t k = pop->condSpecies[c];
        double mk = W[k] / NA_; // molecular mass [kg]
        double mu = mk * mp2 / (mk + mp2); // reduced mass [kg]
        double dkp = dc2 + pop->condDiameter[c];
        double beta = 0.25 * M_PI * dkp * dkp *
                      sqrt(8.0 * Cantera::Boltzmann * T / (M_PI * mu)); // [m^3/s]
        double conc = Y[k] * rho / W[k]; // [kmol/m^3]
        double R = std::max(beta * n2 * conc, 0.0); // [kmol/m^3/s]
        double dM = W[k] * R / rho; // [kg / kg gas / s]
        momentsQ[kM2] += dM;
        momentsQ[kA2] += 4.0 / dpp * dM / rhoP;
        speciesD[k] += dM;
    }
}

void SourceSystem::computeParticleSources
(dvec& momentsQ, dvec& momentsD, dvec& speciesQ, dvec& speciesD)
{
    momentsQ.setZero(nMoments);
    momentsD.setZero(nMoments);
    speciesQ.setZero(nSpec);
    speciesD.setZero(nSpec);

    dvec Q, D;
    for (const ParticleOptions& p : options->particles) {
        if (p.nMoments == 0) {
            continue;
        }
        // The compute*Rates functions act on this population only
        pop = &p;
        pn = p.nMoments;
        pm = moments.segment(p.offset, pn);
        Q.setZero(pn);
        D.setZero(pn);

        updateParticleTemperature();
        if (p.bimodal) {
            computeBimodalRates(Q, D, speciesQ, speciesD);
        } else {
            computeNucleationRates(Q, D, speciesQ, speciesD);
            computeSurfaceReactionRates(Q, D, speciesQ, speciesD);
            computeCondensationRates(Q, D, speciesQ, speciesD);
            computeCoagulationRates(D);
            computeParticleHeatTransfer(Q);
        }

        momentsQ.segment(p.offset, pn) = Q;
        momentsD.segment(p.offset, pn) = D;
    }
}

void SourceSystem::setTimers
(PerfTimer* reactionRates, PerfTimer* thermo, PerfTimer* jacobian)
{
    reactionRatesTimer = reactionRates;
    thermoTimer = thermo;
    jacobianTimer = jacobian;
}

void SourceSystem::setPosition(size_t _j, double _x)
{
    j = static_cast<int>(_j);
    x = _x;
}

void SourceSystem::setupQuasi2d(std::shared_ptr<BilinearInterpolator> vz_int,
                                std::shared_ptr<BilinearInterpolator> T_int)
{
    quasi2d = true;
    vzInterp = vz_int;
    TInterp = T_int;
}

void SourceSystem::writeState(std::ostream& out, bool init)
{
    if (init) {
        out << "T = []" << std::endl;
        out << "Y = []" << std::endl;
        out << "t = []" << std::endl;
    }

    Eigen::IOFormat fmt(15, Eigen::DontAlignCols, ",", ",", "", "", "[", "]");

    out << "T.append(" << T << ")" << std::endl;
    out << "Y.append(" << Y.format(fmt) << ")" << std::endl;
    out << "t.append(" << time() << ")" << std::endl;
}

// ----------------------------------------------------------------------------

void SourceSystemCVODE::initialize(size_t new_nSpec, size_t new_nMoments)
{
    SourceSystem::initialize(new_nSpec, new_nMoments);
    dYdt.resize(nSpec);
    wDot.resize(nSpec);
    dMomentsdt.resize(nMoments);

    integrator.reset(new SundialsCvode(static_cast<int>(nSpec+2+nMoments)));
    integrator->setODE(this);
    integrator->linearMultistepMethod = CV_BDF;
    integrator->maxNumSteps = 1000000;
}

void SourceSystemCVODE::setOptions(ConfigOptions& opts)
{
    SourceSystem::setOptions(opts);
    integrator->abstol[kMomentum] = options->integratorMomentumAbsTol;
    integrator->abstol[kEnergy] = options->integratorEnergyAbsTol;
    for (size_t k=0; k<nSpec; k++) {
        integrator->abstol[kSpecies+k] = options->integratorSpeciesAbsTol;
    }
    for (size_t m=0; m<nMoments; m++) {
        integrator->abstol[kSpecies+nSpec+m] = options->integratorSpeciesAbsTol *
            options->momentAbsTolScale[m];
    }
    integrator->reltol = options->integratorRelTol;
    integrator->minStep = options->integratorMinTimestep;
    integrator->errorStopCount = options->errorStopCount;
}

int SourceSystemCVODE::f(const realtype t, const sdVector& y, sdVector& ydot)
{
    unroll_y(y, t);

    reactionRatesTimer->start();
    if (rateMultiplierFunction) {
        gas->setRateMultiplier(rateMultiplierFunction->a(t));
    }
    gas->thermo->setMassFractions_NoNorm(Y.data());
    gas->thermo->setState_TP(T, gas->pressure);
    gas->getReactionRates(wDot);
    reactionRatesTimer->stop();

    updateThermo();

    // *** Calculate the time derivatives
    double scale;
    if (!quasi2d) {
        scale = 1.0;
        dUdt = splitConst[kMomentum];
        qDot = - (wDot * hk).sum() + getQdotIgniter(t);
        if (heatLoss && options->alwaysUpdateHeatFlux) {
            qDot -= heatLoss->eval(x, t, U, T, Y);
        } else {
            qDot -= qLoss;
        }
        dTdt = qDot/(rho*cp) + splitConst[kEnergy];
    } else {
        scale = 1.0/vzInterp->get(x, t);
        dUdt = splitConst[kMomentum];
        dTdt = splitConst[kEnergy];
    }
    if (options->fixedTemperature) {
        dTdt = 0;
    }
    dYdt = scale * wDot * W / rho + splitConst.segment(kSpecies, nSpec);

    // Particle nucleation source terms (two-way coupled with the gas phase)
    dvec momentsQ, momentsD, speciesQ, speciesD;
    computeParticleSources(momentsQ, momentsD, speciesQ, speciesD);
    dMomentsdt = momentsQ - momentsD + splitConst.segment(kSpecies+nSpec, nMoments);
    dYdt += speciesQ - speciesD;

    roll_ydot(ydot);
    return 0;
}

int SourceSystemCVODE::denseJacobian(const realtype t, const sdVector& y,
                                     const sdVector& ydot, sdMatrix& J)
{
    // TODO: Verify that f has just been called so that we don't need to
    // unroll y and compute all the transport properties.

    fdJacobian(t, y, ydot, J);
    return 0;

    double a = strainFunction->a(t);
    double dadt = strainFunction->dadt(t);
    double A = a*a + dadt;

    // Additional properties not needed for the normal function evaluations:
    thermoTimer->start();
    gas->getSpecificHeatCapacities(cpSpec);
    Wmx = gas->getMixtureMolecularWeight();
    thermoTimer->stop();

    // The constant "800" here has been empirically determined to give
    // good performance for typical test cases. This value can have
    // a substantial impact on the convergence rate of the solver.
    double eps = sqrt(DBL_EPSILON)*800;

    // *** Derivatives with respect to temperature
    double TplusdT = T*(1+eps);

    double dwhdT = 0;
    dvec dwdT(nSpec);
    dvec wDot2(nSpec);

    reactionRatesTimer->start();
    gas->setStateMass(Y, TplusdT);
    gas->getReactionRates(wDot2);
    reactionRatesTimer->stop();

    for (size_t k=0; k<nSpec; k++) {
        dwdT[k] = (wDot2[k]-wDot[k])/(TplusdT-T);
        dwhdT += hk[k]*dwdT[k] + cpSpec[k]*wDot[k];
    }

    double drhodT = -rho/T;

    // *** Derivatives with respect to species concentration
    dmatrix dwdY(nSpec, nSpec);
    dvec hdwdY = dvec::Zero(nSpec);
    dvec YplusdY(nSpec);
    dvec drhodY(nSpec);

    double scale = (quasi2d) ? 1.0/vzInterp->get(x, t) : 1.0;

    for (size_t k=0; k<nSpec; k++) {
        YplusdY[k] = (abs(Y[k]) > eps/2) ? Y[k]*(1+eps) : eps;
        reactionRatesTimer->start();
        gas->thermo->setMassFractions_NoNorm(YplusdY.data());
        gas->thermo->setState_TP(T, gas->pressure);

        gas->getReactionRates(wDot2);
        reactionRatesTimer->stop();

        for (size_t i=0; i<nSpec; i++) {
            dwdY(i,k) = (wDot2[i]-wDot[i])/(YplusdY[k]-Y[k]);
            hdwdY[k] += hk[i]*dwdY(i,k);
        }

        drhodY[k] = rho*(W[k]-Wmx)/(W[k]*(1-Y[k]*(1-eps)));
    }

    for (size_t k=0; k<nSpec; k++) {
        for (size_t i=0; i<nSpec; i++) {
            // dSpecies/dY
            J(kSpecies+k, kSpecies+i) = scale *
                (dwdY(k,i)*W[k]/rho - wDot[k]*W[k]*drhodY[i]/(rho*rho));
        }
        if (!quasi2d) {
            // dSpecies/dT
            J(kSpecies+k, kEnergy) = dwdT[k]*W[k]/rho -
                wDot[k]*W[k]*drhodT/(rho*rho);

            // dEnergy/dY
            J(kEnergy, kSpecies+k) = -hdwdY[k]/(rho*cp) -
                qDot*drhodY[k]/(rho*rho*cp);

            // dMomentum/dY
            J(kMomentum, kSpecies+k) = -A*drhodY[k]/(rho*rho);
        }
    }

    if (!quasi2d) {
        // dEnergy/dT
        J(kEnergy, kEnergy) = -dwhdT/(rho*cp) - qDot*drhodT/(rho*rho*cp);

        // dMomentum/dU
        J(kMomentum, kMomentum) = 0;

        // dMomentum/dT
        J(kMomentum, kEnergy) = -A*drhodT/(rho*rho);
    }

    return 0;
}

int SourceSystemCVODE::fdJacobian(const realtype t, const sdVector& y,
                                  const sdVector& ydot, sdMatrix& J)
{
    jacobianTimer->start();
    sdVector yplusdy(y.length(), sunContext);
    sdVector ydot2(y.length(), sunContext);
    size_t nVars = nSpec+2+nMoments;
    double eps = sqrt(DBL_EPSILON);
    double atol = DBL_EPSILON;

    for (size_t i=0; i<nVars; i++) {
        for (size_t k=0; k<nVars; k++) {
            yplusdy[k] = y[k];
        }
        double dy = (abs(y[i]) > atol) ? abs(y[i])*(eps) : abs(y[i])*eps + atol;
        yplusdy[i] += dy;
        f(t, yplusdy, ydot2);
        for (size_t k=0; k<nVars; k++) {
            J(k,i) = (ydot2[k]-ydot[k])/dy;
        }
    }

    jacobianTimer->stop();

    return 0;
}

void SourceSystemCVODE::setState
(double tInitial, double uu, double tt, const dvec& yy, const dvec& mm)
{
    integrator->t0 = tInitial;
    integrator->y[kMomentum] = uu;
    integrator->y[kEnergy] = tt;
    Eigen::Map<dvec>(&integrator->y[kSpecies], nSpec) = yy;
    Eigen::Map<dvec>(&integrator->y[kSpecies+nSpec], nMoments) = mm;
    integrator->initialize();
    if (heatLoss && !options->alwaysUpdateHeatFlux) {
        qLoss = heatLoss->eval(x, tInitial, uu, tt, const_cast<dvec&>(yy));
    }
}

int SourceSystemCVODE::integrateToTime(double tf)
{
    try {
        return integrator->integrateToTime(integrator->t0 + tf);
    } catch (Cantera::CanteraError& err) {
        logFile.write(err.what());
        return -1;
    }
}

int SourceSystemCVODE::integrateOneStep(double tf)
{
    return integrator->integrateOneStep(integrator->t0 + tf);
}

double SourceSystemCVODE::time() const
{
    return integrator->tInt - integrator->t0;
}

void SourceSystemCVODE::unroll_y(const sdVector& y, double t)
{
    if (!quasi2d) {
        T = y[kEnergy];
        U = y[kMomentum];
    } else {
        T = TInterp->get(x, t);
        U = 0;
    }
    Y = Eigen::Map<dvec>(&y[kSpecies], nSpec);
    moments = Eigen::Map<dvec>(&y[kSpecies+nSpec], nMoments);
}

void SourceSystemCVODE::roll_y(sdVector& y) const
{
    y[kEnergy] = T;
    y[kMomentum] = U;
    Eigen::Map<dvec>(&y[kSpecies], nSpec) = Y;
    Eigen::Map<dvec>(&y[kSpecies+nSpec], nMoments) = moments;
}

void SourceSystemCVODE::roll_ydot(sdVector& ydot) const
{
    ydot[kEnergy] = dTdt;
    ydot[kMomentum] = dUdt;
    Eigen::Map<dvec>(&ydot[kSpecies], nSpec) = dYdt;
    Eigen::Map<dvec>(&ydot[kSpecies+nSpec], nMoments) = dMomentsdt;
}

std::string SourceSystemCVODE::getStats()
{
    return (format("%i") % integrator->getNumSteps()).str();
}

void SourceSystemCVODE::writeState(std::ostream& out, bool init)
{
    SourceSystem::writeState(out, init);
    if (init) {
        out << "dTdt = []" << std::endl;
        out << "dYdt = []" << std::endl;
        out << "splitConstT = []" << std::endl;
        out << "splitConstY = []" << std::endl;
        out << "moments = []" << std::endl;
        out << "dMomentsdt = []" << std::endl;
    }

    Eigen::IOFormat fmt(15, Eigen::DontAlignCols, ",", ",", "", "", "[", "]");

    out << "dTdt.append(" << dTdt << ")" << std::endl;
    out << "dYdt.append(" << dYdt.format(fmt) << ")" << std::endl;
    out << "splitConstT.append(" << splitConst[kEnergy] << ")" << std::endl;
    out << "splitConstY.append(" << splitConst.segment(kSpecies, nSpec).format(fmt) << ")" << std::endl;
    out << "moments.append(" << moments.format(fmt) << ")" << std::endl;
    out << "dMomentsdt.append(" << dMomentsdt.format(fmt) << ")" << std::endl;
}

void SourceSystemCVODE::writeJacobian(std::ostream& out)
{
    int N = integrator->y.length();
    double t = integrator->tInt;
    sdMatrix J(N,N, sunContext);
    sdVector ydot(N, sunContext);
    f(t, integrator->y, ydot);
    denseJacobian(t, integrator->y, ydot, J);

    out << "J = []" << std::endl;
    for (int i=0; i<N; i++) {
        out << "J.append([";
        for (int k=0; k<N; k++) {
            out << boost::format("%.5e, ") % J(i,k);
        }
        out << "])" << std::endl;
    }
}

// ----------------------------------------------------------------------------

SourceSystemQSS::SourceSystemQSS()
{
    integrator.setOde(this);
    dUdtQ = 0;
    dUdtD = 0;
    dTdtQ = 0;
    dTdtD = 0;
}

void SourceSystemQSS::initialize(size_t new_nSpec, size_t new_nMoments)
{
    SourceSystem::initialize(new_nSpec, new_nMoments);
    integrator.initialize(new_nSpec + 2 + new_nMoments);

    dYdtQ.setConstant(nSpec, 0);
    dYdtD.setConstant(nSpec, 0);
    wDotD.resize(nSpec);
    wDotQ.resize(nSpec);
    dMomentsdtQ.setConstant(nMoments, 0);
    dMomentsdtD.setConstant(nMoments, 0);

    integrator.enforce_ymin[kMomentum] = false;
}

void SourceSystemQSS::setOptions(ConfigOptions& opts)
{
    SourceSystem::setOptions(opts);
    integrator.epsmin = options->qss_epsmin;
    integrator.epsmax = options->qss_epsmax;
    integrator.dtmin = options->qss_dtmin;
    integrator.dtmax = options->qss_dtmax;
    integrator.itermax = options->qss_iterationCount;
    integrator.abstol = options->qss_abstol;
    integrator.stabilityCheck = options->qss_stabilityCheck;
    integrator.ymin.setConstant(nSpec + 2 + nMoments, options->qss_minval);
    integrator.ymin[kMomentum] = -1e4;
    for (size_t m=0; m<nMoments; m++) {
        // The particle enthalpy is signed (oxide formation enthalpies); its
        // whole net source is put in the "production" term
        integrator.enforce_ymin[kSpecies+nSpec+m] = options->momentEnforceMin[m];
    }
}

void SourceSystemQSS::setState
(double tStart, double uu, double tt, const dvec& yy, const dvec& mm)
{
    dvec yIn(nSpec + 2 + nMoments);
    yIn << uu, tt, yy, mm;
    integrator.setState(yIn, tStart);
    if (heatLoss && !options->alwaysUpdateHeatFlux) {
        qLoss = heatLoss->eval(x, tStart, uu, tt, const_cast<dvec&>(yy));
    }
}

void SourceSystemQSS::odefun(double t, const dvec& y, dvec& q, dvec& d,
                             bool corrector)
{
    tCall = t;
    unroll_y(y, corrector);

    // *** Update auxiliary data ***
    reactionRatesTimer->start();
    if (rateMultiplierFunction) {
        gas->setRateMultiplier(rateMultiplierFunction->a(t));
    }
    gas->setStateMass(Y, T);
    gas->getCreationRates(wDotQ);
    gas->getDestructionRates(wDotD);
    reactionRatesTimer->stop();

    if (!corrector) {
        updateThermo();
    }

    qDot = - ((wDotQ - wDotD) * hk).sum() + getQdotIgniter(t);

    // *** Calculate the time derivatives
    double scale;
    if (!quasi2d) {
        scale = 1.0;
        dUdtQ = splitConst[kMomentum];
        dUdtD = 0;
        dTdtQ = qDot/(rho*cp) + splitConst[kEnergy];
    } else {
        scale = 1.0/vzInterp->get(x, t);
        dUdtQ = 0;
        dUdtD = 0;
        dTdtQ = 0;
    }

    if (heatLoss && options->alwaysUpdateHeatFlux) {
        dTdtD = heatLoss->eval(x, t, U, T, const_cast<dvec&>(Y)) / (rho*cp);
    } else {
        dTdtD = qLoss / (rho*cp);
    }
    if (options->fixedTemperature) {
        dTdtQ = 0;
        dTdtD = 0;
    }

    dYdtQ = scale * wDotQ * W / rho + splitConst.segment(kSpecies, nSpec);
    dYdtD = scale * wDotD * W / rho;

    // Particle nucleation source terms (two-way coupled with the gas phase)
    dvec momentsQ, momentsD, speciesQ, speciesD;
    computeParticleSources(momentsQ, momentsD, speciesQ, speciesD);
    dMomentsdtQ = momentsQ + splitConst.segment(kSpecies+nSpec, nMoments);
    dMomentsdtD = momentsD;
    dYdtQ += speciesQ;
    dYdtD += speciesD;

    assert(rhou > 0);
    assert(rho > 0);
    assert(U > -1e100 && U < 1e100);
    assert(splitConst[kMomentum] > -1e100 && splitConst[kMomentum] < 1e100);

    assert(dUdtQ > -1e100 && dUdtQ < 1e100);
    assert(dUdtD > -1e100 && dUdtD < 1e100);

    roll_ydot(q, d);
}

void SourceSystemQSS::unroll_y(const dvec& y, bool corrector)
{
    if (!quasi2d) {
        if (!corrector) {
            T = y[kEnergy];
        }
        U = y[kMomentum];
    } else {
        if (!corrector) {
            T = TInterp->get(x, tCall);
        }
        U = 0;
    }
    Y = y.segment(kSpecies, nSpec);
    moments = y.segment(kSpecies+nSpec, nMoments);
}

void SourceSystemQSS::roll_y(dvec& y) const
{
    y << U, T, Y, moments;
}

void SourceSystemQSS::roll_ydot(dvec& q, dvec& d) const
{
    q << dUdtQ, dTdtQ, dYdtQ, dMomentsdtQ;
    d << dUdtD, dTdtD, dYdtD, dMomentsdtD;
}

std::string SourceSystemQSS::getStats()
{
    return (format("%i/%i") % integrator.gcount % integrator.rcount).str();
}
