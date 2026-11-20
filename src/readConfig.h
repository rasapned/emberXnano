#pragma once
#include "mathUtils.h"

#include <string>
#include <vector>

// indices of the respective equations / solution components
const size_t kMomentum = 0;
const size_t kEnergy = 1;
const size_t kSpecies = 2;
const size_t kWmx = 2; // never used in the same systems as kSpecies

// Particle moment indices, relative to the start of the moment block
// (i.e. row kSpecies+nSpec in the full state vector). Fixed convention:
// index 0 = N (kmol particles / kg gas), index 1 = particle-phase metal mass
// (kg / kg gas), index 2 = particle-phase O mass (kg / kg gas; only present
// when nMoments >= 3, otherwise particles are pure metal).
const size_t kN = 0;
const size_t kM = 1;
const size_t kO = 2;
// Index 3 = particle enthalpy H (J / kg gas, absolute, on the same reference
// as the gas-phase thermo; only present with particles.energy, i.e. when
// nMoments >= 4 here, where nMoments counts all moment rows). The particle
// temperature follows from H and the mass moments.
const size_t kH = 3;
// Bimodal model (particles.model = 'bimodal'), relative to the population
// block: mode-1 particles of fixed size (nuclei), N1 [kmol / kg gas], and
// mode-2 aggregates: number N2 [kmol / kg gas], mass M2 [kg / kg gas] and
// surface area A2 [m^2 / kg gas].
const size_t kN1 = 0;
const size_t kN2 = 1;
const size_t kM2 = 2;
const size_t kA2 = 3;

//! Scale of the H moment [J/kg] relative to the mass-like moments, used to
//! set its absolute integrator tolerances
const double particleEnthalpyTolScale = 1e7;

//! Possible boundary conditions for the continuity equations
namespace ContinuityBoundaryCondition {
    enum BC {
        Left, //!< V is computed from j = 0
        Right, //!< V is computed from j = jj
        Zero, //!< V = 0 at x = 0
        Qdot, //!< V is calculated from the location of maximum heat release rate
        Temp, //!< V is calculated from the temperature midpoint
        Wall //!< V = 0 at j = jj (impermeable wall; impinging jet)
    };
}

//! Options of one particle population [see input.py Particles]. Each
//! population owns the moment rows offset .. offset+nMoments-1 of the moment
//! block (state rows kSpecies+nSpec+offset ...); within them the indices
//! kN, kM, kO, kH above apply. Populations interact only through the gas.
struct ParticleOptions
{
    size_t offset = 0; //!< first row of this population in the moment block
    size_t nMoments = 0; //!< number of moment rows (incl. the enthalpy row H)

    //! [particles.model] 'bimodal': rows kN1, kN2, kM2, kA2 of a single
    //! material of density phaseDensity[0]; otherwise single mode (kN, kM,
    //! kO, kH)
    bool bimodal = false;
    //! Bimodal: mass [kg] and diameter [m] of a mode-1 particle (nucleus)
    double monomerMass = 0.0, monomerDiameter = 0.0;
    //! [particles.fractalDimension] of the aggregates (bimodal mode 2, or
    //! single mode with primaryDiameter)
    double fractalDimension = 3.0;
    //! [particles.sintering] Bimodal: A2 relaxes to the area of the
    //! volume-equivalent spheres with tau = sinterA T dpp^4 exp(sinterTa/T)
    bool sintering = false;
    double sinterA = 0.0, sinterTa = 0.0;
    //! [particles.primaryDiameter] Single mode: prescribed primary particle
    //! diameter [m] for the fractal collision diameter used in coagulation
    //! and diffusion; 0 = compact spheres
    double primaryDiameter = 0.0;
    double momentDiffusivity = 0.0; //!< [particles.momentDiffusivity] Diffusivity of moment scalars [m^2/s]
    double momentBCLeft = 0.0; //!< [particles.momentBCLeft] Fixed left boundary value
    double momentBCRight = 0.0; //!< [particles.momentBCRight] Fixed right boundary value
    // Particle nucleation source terms (two-way coupled with the gas phase)
    //! [particles.phases] Particle material phases as parallel arrays of
    //! O/metal atomic ratio (ascending, starting at 0 = pure metal) and bulk
    //! density [kg/m^3]. The particle density is interpolated linearly in the
    //! O/metal ratio between them (used only to convert particle mass into
    //! particle size); the last ratio is the maximum oxidation state.
    std::vector<double> phaseRatio, phaseDensity;

    //! Atomic weights [kg/kmol] of the particle metal ([particles.metal])
    //! and of oxygen.
    double metalWeight = 0.0, oxygenWeight = 0.0;

    //! [particles.coagulation] Enable the Brownian coagulation sink for the
    //! N moment (particle-particle collisions: N decreases, mass unchanged).
    bool coagulation = false;

    //! [particles.minParticleDiameter] Particle diameter [m] below which
    //! mass-removing surface reactions also remove particles from N, in
    //! proportion to the removed mass (particle disintegration).
    double minParticleDiameter = 0.0;

    //! Collision-based nucleation channels (parallel arrays, one entry per
    //! channel). Channel c consumes gas species nucSpeciesA[c] (stoichiometry
    //! nucStoichA[c]) and nucSpeciesB[c] (stoichiometry nucStoichB[c]; may
    //! equal speciesA for a homomolecular collision) to form
    //! nucMassPerEvent[c] [kg/kmol of events] of new particle-phase mass,
    //! of which nucOxygenMassPerEvent[c] is O (oxide nuclei) and the rest metal.
    //! nucMassPerEvent is ALWAYS derived (in Python) from the consumed
    //! species' molecular weights so that the particle metal mass production
    //! equals the gas mass consumption rate identically -- this is what
    //! makes sum(Y) + mM + mO == 1 an exact invariant of the coupled ODE.
    //! Only species of the metal and O may nucleate. nucCollisionPrefactor[c] is the
    //! hard-sphere kinetic collision prefactor (diameter/reduced-mass terms,
    //! independent of local state) precomputed once in Python from the
    //! mechanism's transport data; only sqrt(T) and the concentrations need
    //! to be evaluated at each grid point/timestep.
    std::vector<int> nucSpeciesA, nucSpeciesB;
    std::vector<int> nucStoichA, nucStoichB;
    std::vector<double> nucMassPerEvent, nucOxygenMassPerEvent;
    std::vector<double> nucCollisionPrefactor;
    //! Monomer conversion [nucleation.conversionSpecies, conversionRate]:
    //! each molecule of species convSpecies[c] becomes one particle, at the
    //! rate convRate * C_k [kmol/m^3/s]; convOxygenMass[c] [kg/kmol] is the
    //! O part of its mass W_k, the rest is metal.
    std::vector<int> convSpecies;
    std::vector<double> convOxygenMass;
    double convRate = 0.0;
    //! [nucleation.efficiency] Multiplier of the collision nucleation rate
    //! (sticking probability of the nucleating collisions).
    double nucEfficiency = 1.0;
    //! [nucleation.kelvinLimit] Multiply channel c by max(0, 1 - p_eq/p_1):
    //! p_1 partial pressure of the metal monomer nucMonomerSpecies, p_eq =
    //! p_sat exp(4 sigma v_1 / (d_n k_B T)) over a nucleus of
    //! nucAtomsPerEvent[c] metal atoms (diameter d_n from n v_1).
    //! nucLnPsatMonomer holds ln(p_sat [Pa]) of the monomer over the
    //! pure-metal phase on particleThermoT; nucMonomerVolume [m^3] is v_1 in
    //! that phase; sigma as for evaporation (surfaceTensionA/B).
    bool nucKelvinLimit = false;
    int nucMonomerSpecies = -1;
    std::vector<int> nucAtomsPerEvent;
    std::vector<double> nucLnPsatMonomer;
    double nucMonomerVolume = 0.0;

    //! Gas-particle surface reactions (parallel arrays, one entry per
    //! reaction): gas species surfReactant[r] impinges on a particle and
    //! leaves as surfProduct[r] (-1 = no gas product), changing the particle
    //! by surfDeltaM[r] metal atoms and surfDeltaO[r] O atoms per event (both
    //! derived in Python from the species compositions; positive = added to
    //! the particle). Rate [kmol/m^3/s] = surfA * exp(-surfTa/T) * theta *
    //! (particle surface area per m^3) * C_reactant, with surfA in [m/s].
    std::vector<int> surfReactant, surfProduct;
    std::vector<double> surfA, surfTa;
    std::vector<int> surfDeltaM, surfDeltaO;

    //! [particles.surfaceCollisionLimit] Cap surfA * exp(-surfTa/T_p) at the
    //! hard-sphere collision rate per unit particle area, with the
    //! reactant's collision diameter surfDiameter[r] [m].
    bool surfCollisionLimit = false;
    std::vector<double> surfDiameter;

    //! Heterogeneous condensation (parallel arrays, one entry per condensing
    //! gas species): species condSpecies[k] (collision diameter
    //! condDiameter[k] [m], containing condDeltaM[k] metal and condDeltaO[k]
    //! O atoms) is added to the particles at the kinetic (hard-sphere,
    //! free-molecular) gas-particle collision rate.
    std::vector<int> condSpecies;
    std::vector<double> condDiameter;
    std::vector<int> condDeltaM, condDeltaO;

    //! [particles.evaporation] Reversible condensation: species c with a
    //! matching particle phase condPhase[c] (index into phaseRatio, -1 =
    //! none: condensation only) also evaporates, at the collision rate times
    //! its saturation concentration over the bulk phase, times the activity
    //! of the phase (lever-rule weight) and the Kelvin factor
    //! exp(4 sigma v_m / (d_p k_B T_p)). condLnPsat (row-major, nCond x nT,
    //! on particleThermoT) holds ln(p_sat [Pa]); condMolarVolume [m^3 per
    //! molecule] is v_m in the phase; sigma = surfaceTensionA -
    //! surfaceTensionB * T_p [J/m^2] [particles.surfaceTension].
    bool evaporation = false;
    std::vector<int> condPhase;
    std::vector<double> condLnPsat, condMolarVolume;
    double surfaceTensionA = 0.0, surfaceTensionB = 0.0;

    //! [particles.coagulationEnhancement] Multiplier of the coagulation
    //! kernel (e.g. van der Waals enhancement, 2-3 at nm sizes).
    double coagulationEnhancement = 1.0;

    //! Particle energy (particles.energy, nMoments >= 4). particleThermoT [K] is a temperature
    //! grid; phaseEnthalpy (row-major, nPhases x nT) holds for each phase the
    //! enthalpy of its stable condensed species [J per kmol of metal], on the
    //! same reference as the gas-phase thermo [particles.phases,
    //! particles.thermoFile].
    std::vector<double> particleThermoT, phaseEnthalpy;
    double thermalAccommodation = 0.0; //!< [particles.thermalAccommodation]
    double radiationAbsorption = 0.0; //!< E(m) [particles.radiationAbsorption]
    double radiationTsurr = 0.0; //!< [particles.radiationTsurr] [K]

    //! Prescribed surface-reaction heat [J/kmol of events]; negative heats
    //! the particle [surfaceReaction.deltaH]. Used only if surfPrescribedHeat
    //! (particles.surfaceReactionHeat = 'prescribed'), otherwise all 0 and
    //! the heat follows from the thermochemistry.
    std::vector<double> surfDeltaH;
    bool surfPrescribedHeat = false;

    //! Classical nucleation theory (CNT), for a single monomer material
    //! whose cluster species (e.g. Fe, Fe2, Fe3, ...) are all listed in
    //! nucPrecursorSpecies, with nucPrecursorAtomCount giving the number of
    //! monomer units in each cluster (index nucMonomerIndex is the monomer
    //! itself, atom count == 1). The monomer mass/volume used by the CNT
    //! rate are derived at runtime from W[nucPrecursorSpecies[nucMonomerIndex]]
    //! and the pure-metal density phaseDensity[0] -- not a separate user input.
    std::vector<int> nucPrecursorSpecies;
    std::vector<int> nucPrecursorAtomCount;
    int nucMonomerIndex = -1;

    //! Antoine equation coefficients for the precursor's saturation vapor
    //! pressure: log10(Psat[bar]) = A - B/(C + T[degC]).
    double nucAntoineA = 0.0, nucAntoineB = 0.0, nucAntoineC = 0.0;

    //! Surface tension as a linear function of temperature [K]:
    //! sigma = nucSurfaceTensionA + nucSurfaceTensionB * T.
    double nucSurfaceTensionA = 0.0, nucSurfaceTensionB = 0.0;
};

//! Configuration options for the simulation.
/*!
 *  These are more fully documented in `input.py`. For variables that
 *  correspond to a definition in `input.py`, the name used for the variable
 *  in that file is given in square brackets.
 */
class ConfigOptions
{
public:
    ConfigOptions() {}

    //! Returns true if integrator stage data should be saved at the current
    //! time step.
    bool debugIntegratorStages(double t) const;

    //! Set the value of #continuityBC by name
    void setContinuityBC(const std::string& condition);

    std::string outputDir; //!< [paths.outputDir]

    //! [general.fixedBurnedVal] True if burned gas state is fixed at
    //! equilibrium conditions.
    bool fixedBurnedVal;

    //! [general.fixedLeftLocation] True if the position of the leftmost grid
    //! point is fixed.
    bool fixedLeftLoc;

    //! [general.twinFlame] Set to true for the twin flame configuration.
    //! If specified, the leftmost grid point will not extend beyond x = 0.
    bool twinFlame;

    //! [general.discFlame] Set to true for the axisymmetric jet flame configuration.
    bool discFlame;

    ContinuityBoundaryCondition::BC continuityBC;

    //! [general.cylindricalFlame] Set to true for the cylindrical flame configuration.
    //! If specified, the leftmost grid point will not extend beyond x = 0.
    bool cylindricalFlame;
    bool unburnedLeft; //!< [general.unburnedLeft]
    bool fuelLeft; //!< [general.fuelLeft]

    std::string flameType; //!< [initialCondition.flameType]

    int regridStepInterval; //!< [times.regridStepInterval]
    int outputStepInterval; //!< [times.outputStepInterval]
    int profileStepInterval; //!< [times.profileStepInterval]
    int currentStateStepInterval; //!< [times.currentStateStepInterval]
    int terminateStepInterval; //!< [times.terminateStepInterval]
    double regridTimeInterval; //!< [times.regridTimeInterval]
    double outputTimeInterval; //!< [times.outputTimeInterval]
    double profileTimeInterval; //!< [times.profileTimeInterval]
    double globalTimestep; //!< [times.globalTimestep]
    double diffusionTimestepMultiplier;  //!< [times.diffusionTimestepMultiplier]
    std::string splittingMethod; //!< [general.splittingMethod]

    std::string chemistryIntegrator; //!< [general.chemistryIntegrator]
    double integratorRelTol; //!< [cvodeTolerances.relativeTolerance]
    double integratorMomentumAbsTol; //!< [cvodeTolerances.momentumAbsTol]
    double integratorEnergyAbsTol; //!< [cvodeTolerances.energyAbsTol]
    double integratorSpeciesAbsTol; //!< [cvodeTolerances.speciesAbsTol]
    double integratorMinTimestep; //!< [cvodeTolerances.minimumTimestep]

    // QSS integrator parameters
    double qss_epsmax; //!< [qssTolerances.epsmax]
    double qss_epsmin; //!< [qssTolerances.epsmin]
    double qss_dtmin; //!< [qssTolerances.dtmin]
    double qss_dtmax; //!< [qssTolerances.dtmax]
    int qss_iterationCount; //!< [qssTolerances.iterationCount]
    double qss_abstol; //!< [qssTolerances.abstol]
    double qss_minval; //!< [qssTolerances.minval]
    bool qss_stabilityCheck; //!< [qssTolerances.stabilityCheck]

    std::string gasMechanismFile; //!< [chemistry.mechanismFile]
    std::string gasPhaseID; //!< [chemistry.phaseID]
    std::string kineticsModel; //!< [chemistry.kineticsModel]
    std::string transportModel; //!< [chemistry.transportModel]
    double transportThreshold; //!< [chemistry.threshold]
    std::string conductivityModel; //!< [chemistry.conductivityModel]
    //! Sutherland coefficients per mechanism species, As [Pa*s/K^0.5] and
    //! Ts [K] (conductivityModel = "sutherland")
    std::vector<double> sutherlandAs, sutherlandTs;
    std::string rateMultiplierFunctionType; //!< [see chemistry.rateMultiplierFunction]

    double pressure; //!< [initialConditoin.pressure]

    std::string strainFunctionType;
    double strainRateInitial; //!< [strainParameters.initial]
    double strainRateFinal; //!< [strainParameters.final]
    double strainRateDt; //!< [strainParameters.dt]
    double strainRateT0; //!< [strainParameters.tStart]

    // Initial profiles specified in configuration file
    dvec x_initial; //!< [initialCondition.x]
    dvec T_initial; //!< [initialCondition.T]
    dvec U_initial; //!< [initialCondition.U]
    dvec V_initial; //!< [initialCondition.V]
    dmatrix Y_initial; //!< [initialCondition.Y]

    //! If false, energy, continuity, and momentum equations are integrated
    //! normally. If true, then time is actually a surrogate for a spatial
    //! variable and the temperature and velocity fields are interpolated from
    //! a supplied data file.
    bool quasi2d;

    //! Data file containing 2D interpolation data [general.interpFile]
    std::string interpFile;

    // Wall Flux properties
    bool wallFlux; //<! True if wall flux is applied when x[0] == 0
    double Tinf; //<! Temperature used for computing wall heat flux [wallFlux.Tinf]
    double Kwall; //<! Thermal conductance of the wall at x = 0 [wallFlux.Kwall]

    // Ignition parameters
    double ignition_tStart; //!< [ignition.tStart]
    double ignition_duration; //!< [ignition.duration]
    double ignition_energy; //!< [ignition.energy]
    double ignition_center; //!< [ignition.center]
    double ignition_stddev; //!< [ignition.stddev]

    bool alwaysUpdateHeatFlux; //!< [externalHeatFlux.alwaysUpdate]

    // Tolerances for adaptation and regridding
    double vtol; //!< [grid.vtol]
    double dvtol; //!< [grid.dvtol]
    double rmTol; //!< [grid.rmTol]
    double dampConst; //!< [grid.dampConst]
    double gridMin; //!< [grid.gridMin]
    double gridMax; //!< [grid.gridMax]
    double uniformityTol; //!< [grid.uniformityTol]
    double absvtol; //!< [grid.absvtol]
    double boundaryTol; //!< [grid.boundaryTol]
    double boundaryTolRm; //!< [grid.boundaryTolRm]
    double unstrainedDownstreamWidth; //!< [grid.unstrainedDownstreamWidth]
    int addPointCount; //!< [grid.addPointCount]

    double tStart; //!< [times.tStart]
    double tEnd; //!< [terminationCondition.tEnd]
    double tEndMin; //!< [terminationCondition.tMin]
    bool haveTStart;

    // Controls which variables are included in the outXXXXXX.h5 files
    bool outputProfiles; //!< [outputFiles.saveProfiles]

    //! Grid point index marked for verbose output (cvodeSteps.py)
    //! [debug.sourcePoint]
    int debugSourcePoint;

    //! Integrator time at which to generate verbose output (then terminate)
    //! [debug.sourceTime]
    double debugSourceTime;

    //! Integrator time at which to begin saving output between integrator
    //! stages [debug.startTime]
    double debugStartTime;

    //! Integrator time at which to stop saving output between integrator
    //! stages [debug.stopTime]
    double debugStopTime;
    bool outputDebugIntegratorStages; //!< [outputFiles.debugIntegratorStages]

    //! If non-empty, code finishes when the specified property is constant to
    //! within some tolerance [terminationCondition.measurement]
    std::string terminationMeasurement;

    //! relative tolerance for termination [terminationCondition.tolerance]
    double terminationTolerance;

    //! absolute tolerance for termination using measurement == 'Q'
    //! [terminationCondition.abstol]
    double terminationAbsTol;

    //! tolerance for termination when measurement == `dTdt`
    //! [terminationCondition.dTdtTol]
    double termination_dTdtTol;

    //! tolerance [1/s] for termination when measurement == `moments`
    //! [terminationCondition.momentsTol]
    double termination_momentsTol;

    //! points where a moment is below this fraction of its peak are left
    //! out of the `moments` termination check
    //! [terminationCondition.momentsFloor]
    double termination_momentsFloor;

    //! period over which to require constant heat release rate
    //! [terminationCondition.steadyPeriod]
    double terminationPeriod;

    int errorStopCount; //!< [general.errorStopCount]
    bool stopIfError;
    int nThreads; //!< [general.nThreads]

    double centerGridMin; //!< [grid.centerGridMin]

    bool xFlameControl;
    double xFlameInitial; //!< [positionControl.xInitial]
    double xFlameFinal; //!< [positionControl.xFinal]
    double xFlameDt; //!< [positionControl.dt]
    double xFlameT0; //!< [positionControl.tStart]
    double xFlameIntegralGain; //!< [positionControl.integralGain]
    double xFlameProportionalGain; //!< [positionControl.proportionalGain]

    bool xStagControl;
    double xStag;

    // Impinging jet: premixed inlet at x = 0, impermeable no-slip wall at
    // the right end of a fixed domain [see input.py ImpingingJet]
    bool impingingJet; //!< True if an ImpingingJet section is given
    double wallTemperature; //!< [impingingJet.wallTemperature]

    //! True if the temperature is prescribed as T(x) (Tfixed_T at positions
    //! Tfixed_x) instead of being solved for [impingingJet.temperatureProfile]
    bool fixedTemperature;
    dvec Tfixed_x; //!< positions of the prescribed temperature profile [m]
    dvec Tfixed_T; //!< prescribed temperature profile [K]

    //! True if the radial strain U(x) is prescribed (Ufixed_U at positions
    //! Ufixed_x) instead of being solved for [impingingJet.strainProfile]
    bool fixedStrain;
    dvec Ufixed_x; //!< positions of the prescribed strain profile [m]
    dvec Ufixed_U; //!< prescribed radial strain profile U = v_r/r [1/s]

    //! True if the strain rate a(t) (i.e. the pressure curvature term) is set
    //! by a feedback controller to reach #massFluxTarget at the inlet
    bool massFluxControl;
    double massFluxTarget; //!< [impingingJet.massFlux] [kg/m^2*s]
    double massFluxProportionalGain; //!< [impingingJet.proportionalGain]
    double massFluxIntegralGain; //!< [impingingJet.integralGain] [1/s]

    //! Particle populations [see input.py Particles], in the order of their
    //! rows in the moment block
    std::vector<ParticleOptions> particles;
    int nMoments; //!< total number of moment rows over all populations
    //! per moment row: multiplier of integratorSpeciesAbsTol (CVODE)
    std::vector<double> momentAbsTolScale;
    //! per moment row: 1 if the QSS integrator keeps the value >= ymin, 0
    //! for signed rows (particle enthalpy)
    std::vector<int> momentEnforceMin;
};
