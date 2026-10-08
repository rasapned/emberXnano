#include "flameSolver.h"
#include "scalarFunction.h"

FlameSolver::FlameSolver()
    : U(0, 0, Stride1X(1 ,1))
    , T(0, 0, Stride1X(1 ,1))
    , Y(0, 0, 0, StrideXX(1 ,1))
    , moments(0, 0, 0, StrideXX(1 ,1))
    , jCorrSolver(jCorrSystem)
    , strainfunc(NULL)
    , rateMultiplierFunction(NULL)
    , stateWriter(NULL)
    , timeseriesWriter(NULL)
    , heatLossFunction(NULL)
    , vzInterp(new BilinearInterpolator)
    , vrInterp(new BilinearInterpolator)
    , TInterp(new BilinearInterpolator)
{
    #ifndef NDEBUG
        Cantera::printStackTraceOnSegfault();
    #endif
}

FlameSolver::~FlameSolver()
{
    delete strainfunc;
    delete rateMultiplierFunction;
}

void FlameSolver::setOptions(const ConfigOptions& _options)
{
    SplitSolver::setOptions(_options);
    tStart = options.tStart;
    tEnd = options.tEnd;

    gas.setOptions(_options);
    grid.setOptions(_options);
}

void FlameSolver::initialize(void)
{
    tbbTaskSched.reset(new tbb::global_control(
        tbb::global_control::max_allowed_parallelism, options.nThreads));
    delete strainfunc;
    strainfunc = newScalarFunction(options.strainFunctionType, options);

    delete rateMultiplierFunction;
    if (options.rateMultiplierFunctionType != "") {
        rateMultiplierFunction = newScalarFunction(options.rateMultiplierFunctionType,
                                                   options);
    } else {
        rateMultiplierFunction = NULL;
    }

    flamePosIntegralError = 0;
    terminationCondition = 1e10;

    if (options.massFluxControl) {
        // a(t) is set by updateMassFluxControl(), starting from the initial
        // strain rate
        delete strainfunc;
        strainfunc = new ControlledFunction(options.strainRateInitial);
    }
    massFluxIntegralError = 0;

    // Cantera initialization
    gas.initialize();
    nSpec = gas.nSpec;
    nMoments = options.nMoments;
    kMoments = kSpecies + nSpec;
    nVars = nSpec + 2 + nMoments;
    W.resize(nSpec);
    gas.getMolecularWeights(W);
    // Set boundary conditions for particles (vector of nMoments values, each equal to opitions.<...>)
    momentsLeft = dvec::Constant(nMoments, options.momentBCLeft);
    momentsRight = dvec::Constant(nMoments, options.momentBCRight);

    // Get Initial Conditions
    loadProfile();

    grid.setSize(x.size());
    convectionSystem.setGas(gas);
    convectionSystem.setLeftBC(Tleft, Yleft, momentsLeft);
    convectionSystem.setTolerances(options);

    for (size_t k=0; k<nVars; k++) {
        DiffusionSystem* term = new DiffusionSystem();
        TridiagonalIntegrator* integrator = new TridiagonalIntegrator(*term);
        integrator->resize(nPoints);
        diffusionTerms.push_back(term);
        diffusionSolvers.push_back(integrator);
    }
    if (options.wallFlux) {
        diffusionTerms[kEnergy].yInf = options.Tinf;
        diffusionTerms[kEnergy].wallConst = options.Kwall;
    }

    resizeAuxiliary();

    ddtConv.setZero();
    ddtDiff.setZero();
    ddtProd.setZero();

    updateChemicalProperties();
    calculateQdot();

    t = tStart;
    tOutput = t;
    tRegrid = t + options.regridTimeInterval;
    tProfile = t + options.profileTimeInterval;
    nTotal = 0;
    nRegrid = 0;
    nOutput = 0;
    nProfile = 0;
    nTerminate = 0;
    nCurrentState = 0;

    grid.updateValues();
    resizeAuxiliary();

    tFlamePrev = t;
    tMassFluxPrev = t;
    tNow = t;

    totalTimer.start();
}

void FlameSolver::setupStep()
{
    setupTimer.start();

    // Debug sanity check
    #ifndef NDEBUG
        bool error = false;
        for (size_t j=0; j<nPoints; j++) {
            if (T(j) < 295 || T(j) > 3000) {
                logFile.write(format(
                    "WARNING: Unexpected Temperature: T = %f at j = %i") % T(j) % j);
                error = true;
            }
        }
        if (error) {
            writeStateFile("err_setupStep", true, false);
        }
    #endif

    // Reset boundary conditions to prevent numerical drift
    if (grid.leftBC == BoundaryCondition::FixedValue) {
        T(0) = Tleft;
        Y.col(0) = Yleft;
        if (nMoments > 0) {
            moments.col(0) = momentsLeft;
        }
    } else if (grid.leftBC == BoundaryCondition::InletFlux) {
        // Only temperature is held at the inlet value: the burner face acts
        // as a heat sink. Species and moments are left alone so that the
        // convective + diffusive flux balance at the face determines them.
        T(0) = Tleft;
    }

    if (grid.rightBC == BoundaryCondition::FixedValue) {
        T(jj) = Tright;
        Y.col(jj) = Yright;
        if (nMoments > 0) {
            moments.col(jj) = momentsRight;
        }
    } else if (grid.rightBC == BoundaryCondition::Wall) {
        // Impinging jet: wall temperature and no slip are prescribed;
        // species and moments at the wall follow from the zero-flux balance.
        // The inlet velocity is plug flow (U = 0), as in Cantera.
        T(jj) = Tright;
        U(jj) = 0;
        U(0) = 0;
    }

    if (options.fixedTemperature) {
        // Prescribed temperature profile (energy equation not solved),
        // interpolated onto the current grid
        T = mathUtils::interp1(options.Tfixed_x, options.Tfixed_T, x, false);
    }

    updateChemicalProperties();

    updateBC();
    if (options.xFlameControl) {
        update_xStag(t, true); // calculate the value of rVzero
    }
    if (options.massFluxControl) {
        updateMassFluxControl(t);
    }
    convectionSystem.set_rVzero(rVzero);
    setupTimer.stop();

    // Set up solvers for split integration
    updateCrossTerms();
}

void FlameSolver::prepareIntegrators()
{
    splitTimer.resume();
    updateParticleDiffusivity();
    // Diffusion terms
    if (!options.quasi2d) {
        // Diffusion solvers: Energy and momentum
        diffusionTerms[kMomentum].B = rho.inverse();
        diffusionTerms[kEnergy].B = (rho * cp).inverse();

        diffusionTerms[kMomentum].D = mu;
        diffusionTerms[kEnergy].D = lambda;

        // Diffusion solvers: Species
        for (size_t k=0; k<nSpec; k++) {
            diffusionTerms[kSpecies+k].B =rho.inverse();
            diffusionTerms[kSpecies+k].D = rhoD.row(k);
        }

        // Diffusion solvers: Particle moments (passive scalars)
        // particleDiffusivity is a plain diffusivity [m^2/s]; like rhoD for
        // the species, the solver needs rho*D [kg/m*s] with B = 1/rho.
        for (size_t m=0; m<nMoments; m++) {
            DiffusionSystem& sys = diffusionTerms[kMoments+m];
            sys.B = rho.inverse();
            sys.D = rho * particleDiffusivity;
        }
    } else {
        // Diffusion solvers: Energy and momentum
        diffusionTerms[kMomentum].B.setZero(nPoints);
        diffusionTerms[kEnergy].B.setZero(nPoints);

        diffusionTerms[kMomentum].D.setZero(nPoints);
        diffusionTerms[kEnergy].D.setZero(nPoints);

        // Diffusion solvers: Species
        for (size_t k=0; k<nSpec; k++) {
            DiffusionSystem& sys = diffusionTerms[kSpecies+k];
            sys.D = rhoD.row(k);
            for (size_t j = 0; j <= jj; j++) {
                sys.B[j] = 1 / (rho[j] * vzInterp->get(x[j], tNow));
            }
        }

        // Diffusion solvers: Particle moments (passive scalars)
        for (size_t m=0; m<nMoments; m++) {
            DiffusionSystem& sys = diffusionTerms[kMoments+m];
            sys.D = rho * particleDiffusivity;
            for (size_t j = 0; j <= jj; j++) {
                sys.B[j] = 1 / (rho[j] * vzInterp->get(x[j], tNow));
            }
        }
    }

    // Burner-stabilized inlet: T(0) is held at Tleft (see setupStep), so the
    // energy equation must not diffuse there. Species and moments keep the
    // InletFlux stencil so their face values follow from the flux balance.
    if (grid.leftBC == BoundaryCondition::InletFlux) {
        diffusionTerms[kEnergy].grid.leftBC = BoundaryCondition::FixedValue;
    }

    // Prescribed temperature profile: no heat conduction
    if (options.fixedTemperature) {
        diffusionTerms[kEnergy].B.setZero(nPoints);
    }

    // Impinging jet: T and U are prescribed at the wall and U at the inlet;
    // species and moments keep the zero-flux Wall stencil.
    if (grid.rightBC == BoundaryCondition::Wall) {
        diffusionTerms[kEnergy].grid.rightBC = BoundaryCondition::FixedValue;
        diffusionTerms[kMomentum].grid.rightBC = BoundaryCondition::FixedValue;
        diffusionTerms[kMomentum].grid.leftBC = BoundaryCondition::FixedValue;
    }

    setDiffusionSolverState(tNow);
    for (size_t i=0; i<nVars; i++) {
        diffusionTerms[i].splitConst = splitConstDiff.row(i);
    }

    // Production terms
    setProductionSolverState(tNow);
    for (size_t j=0; j<nPoints; j++) {
        sourceTerms[j].splitConst = splitConstProd.col(j);
    }

    // Convection terms
    setConvectionSolverState(tNow);
    dmatrix ddt = ddtConv + ddtDiff + ddtProd;
    if (options.splittingMethod == "balanced") {
        ddt += ddtCross;
    }

    // Change bottomRows to middleRows (now we also have particle moments at the bottom)
    dvec tmp = (W.inverse().matrix().transpose() * ddt.middleRows(kSpecies, nSpec).matrix()).array();
    drhodt = - rho * (ddt.row(kEnergy).transpose() / T + tmp * Wmx);

    assert(mathUtils::notnan(drhodt));
    convectionSystem.setDensityDerivative(drhodt);
    convectionSystem.setSplitConstants(splitConstConv);
    convectionSystem.updateContinuityBoundaryCondition(qDot, options.continuityBC);
    splitTimer.stop();
}

int FlameSolver::finishStep()
{
    logFile.verboseWrite("done!");

    // Burner-stabilized inlet: T(0) is prescribed (reset to Tleft in
    // setupStep), so any change the operators made to it is not a real rate.
    // Left in, it enters drhodt and continuity turns it into a spurious mass
    // source at the burner face; it also feeds the split constants and the
    // dTdt termination measure.
    if (grid.leftBC == BoundaryCondition::InletFlux) {
        ddtConv(kEnergy, 0) = 0;
        ddtDiff(kEnergy, 0) = 0;
        ddtProd(kEnergy, 0) = 0;
        ddtCross(kEnergy, 0) = 0;
    }

    // Impinging jet: the same holds for the prescribed U at both ends and
    // the wall temperature.
    if (grid.rightBC == BoundaryCondition::Wall) {
        ddtConv(kMomentum, 0) = 0;
        ddtDiff(kMomentum, 0) = 0;
        ddtProd(kMomentum, 0) = 0;
        ddtConv(kMomentum, jj) = 0;
        ddtDiff(kMomentum, jj) = 0;
        ddtProd(kMomentum, jj) = 0;
        ddtConv(kEnergy, jj) = 0;
        ddtDiff(kEnergy, jj) = 0;
        ddtProd(kEnergy, jj) = 0;
        ddtCross(kEnergy, jj) = 0;
    }

    // Prescribed temperature profile: T has no rate of change anywhere
    if (options.fixedTemperature) {
        ddtConv.row(kEnergy).setZero();
        ddtDiff.row(kEnergy).setZero();
        ddtProd.row(kEnergy).setZero();
        ddtCross.row(kEnergy).setZero();
    }

    // *** End of Strang-split integration step ***
    correctMassFractions();

    t = tNow + dt;
    tNow += dt;

    nOutput++;
    nRegrid++;
    nProfile++;
    nTerminate++;
    nCurrentState++;

    if (debugParameters::debugTimesteps) {
        int nSteps = convectionSystem.getNumSteps();
        logFile.write(format("t = %8.6f (dt = %9.3e) [C: %i]") % t % dt % nSteps);
    }

    setupTimer.resume();
    if (t + 0.5 * dt > tOutput || nOutput >= options.outputStepInterval) {
        calculateQdot();
        timeVector.push_back(t);
        heatReleaseRate.push_back(getHeatReleaseRate());
        saveTimeSeriesData("out", false);
        tOutput = t + options.outputTimeInterval;
        nOutput = 0;
    }

    // Periodic check for terminating the integration (based on steady heat
    // release rate, etc.) Quit now to skip grid adaptation on the last step
    if (t >= tEnd) {
        return 1;
    } else if (nTerminate >= options.terminateStepInterval) {
        nTerminate = 0;
        if (checkTerminationCondition()) {
            return 1;
        }
    }

    // Save the current integral and profile data in files that are
    // automatically overwritten, and save the time-series data (out.h5)
    if (nCurrentState >= options.currentStateStepInterval) {
        calculateQdot();
        nCurrentState = 0;
        saveTimeSeriesData("out", true);
        writeStateFile("profNow");
    }

    // *** Save flame profiles
    if (t + 0.5 * dt > tProfile || nProfile >= options.profileStepInterval) {
        if (options.outputProfiles) {
            writeStateFile();
        }
        tProfile = t + options.profileTimeInterval;
        nProfile = 0;
    }
    setupTimer.stop();

    if (t > tRegrid || nRegrid >= options.regridStepInterval) {
        if (debugParameters::debugAdapt || debugParameters::debugRegrid) {
            writeStateFile("preAdapt", false, false);
        }
        regridTimer.start();
        tRegrid = t + options.regridTimeInterval;
        nRegrid = 0;

        // If the left grid point moves, a new boundary value for rVzero
        // needs to be calculated from the mass flux V on the current grid points
        dvec x_prev = grid.x;
        convectionSystem.evaluate();
        dvec V_prev = convectionSystem.V;

        // dampVal sets a limit on the maximum grid size
        grid.dampVal.resize(grid.x.rows());
        for (size_t j=0; j<nPoints; j++) {
            double num = std::min(mu[j],lambda[j]/cp[j]);
            for (size_t k=0; k<nSpec; k++) {
                if (rhoD(k,j) > 0) {
                    num = std::min(num, rhoD(k,j));
                }
            }
            double den = std::max(std::abs(rho[j]*strainfunc->a(t)), 1e-100);
            grid.dampVal[j] = sqrt(num/den);
        }
        dvec dampVal_prev = grid.dampVal;

        vector<dvector> currentSolution;
        rollVectorVector(currentSolution, state);
        rollVectorVector(currentSolution, ddtConv);
        rollVectorVector(currentSolution, ddtDiff);
        rollVectorVector(currentSolution, ddtProd);

        grid.nAdapt = nVars;
        if (options.quasi2d || options.impingingJet) {
            // do not change grid extents in this case
        } else if (strainfunc->a(tNow) == 0) {
            calculateQdot();
            grid.regridUnstrained(currentSolution, qDot);
        } else {
            grid.regrid(currentSolution);
        }

        // Interpolate dampVal onto the modified grid
        grid.dampVal = mathUtils::interp1(x_prev, dampVal_prev, grid.x, false);

        grid.adapt(currentSolution);

        // Perform updates that are necessary if the grid has changed
        if (grid.updated) {
            logFile.write(format("Grid size: %i points.") % nPoints);
            resizeMappedArrays();
            unrollVectorVector(currentSolution, state, 0);
            unrollVectorVector(currentSolution, ddtConv, 1);
            unrollVectorVector(currentSolution, ddtDiff, 2);
            unrollVectorVector(currentSolution, ddtProd, 3);
            correctMassFractions();

            // Update the mass flux (including the left boundary value)
            rVzero = mathUtils::interp1(x_prev, V_prev, grid.x[0]);
            convectionSystem.utwSystem.V = mathUtils::interp1(x_prev, V_prev, grid.x);

            // Allocate the solvers and arrays for auxiliary variables
            resizeAuxiliary();
            calculateQdot();

            if (debugParameters::debugAdapt || debugParameters::debugRegrid) {
                writeStateFile("postAdapt", false, false);
            }
            grid.updated = false;
        }
        regridTimer.stop();
    }

    if (nTotal % 10 == 0) {
        printPerformanceStats();
    }
    nTotal++;
    return 0;
}

void FlameSolver::finalize()
{
    calculateQdot();
    saveTimeSeriesData("out", true);

    // *** Integration has reached the termination condition
    if (options.outputProfiles) {
        writeStateFile();
    }

    totalTimer.stop();
    printPerformanceStats();
    logFile.write(format("Runtime: %f seconds.") % totalTimer.getTime());
}

bool FlameSolver::checkTerminationCondition(void)
{
    if (tNow < options.tEndMin) {
        return false;
    }
    if (options.terminationMeasurement == "Q") {
        size_t j1 = mathUtils::findLast(timeVector < (tNow - options.terminationPeriod));
        if (j1 == npos) {
            logFile.write(format(
                    "Continuing integration: t (%8.6f) < terminationPeriod (%8.6f)") %
                    (tNow-timeVector[0]) % options.terminationPeriod);
            return false;
        }

        size_t j2 = timeVector.size()-1;
        double qMean = mathUtils::mean(heatReleaseRate,j1,j2);
        double hrrError = 0;
        for (size_t j=j1; j<=j2; j++) {
            hrrError += pow(heatReleaseRate[j]-qMean, 2);
        }
        hrrError = sqrt(hrrError) / (j2-j1+1);
        terminationCondition = hrrError / abs(qMean);

        logFile.write(format(
            "Heat release rate RMS error = %6.3f%%. absolute error: %9.4e") %
            (hrrError/qMean*100) % hrrError);

        if (hrrError/abs(qMean) < options.terminationTolerance) {
            logFile.write("Terminating integration: "
                    "Heat release RMS variation less than relative tolerance.");
            return true;
        } else if (hrrError < options.terminationAbsTol) {
            logFile.write("Terminating integration: "
                    "Heat release rate RMS variation less than absolute tolerance.");
            return true;
        }
    } else if (options.terminationMeasurement == "dTdt") {
        dvec dTdt = (ddtDiff + ddtConv + ddtProd + ddtCross).row(kEnergy);
        double value = (dTdt / T).matrix().norm() / sqrt(static_cast<double>(nPoints));
        logFile.write(format(
            "||1/T * dT/dt|| = %7.3f. Termination threshold = %7.2f") %
            value % options.termination_dTdtTol);
        if (value < options.termination_dTdtTol) {
            logFile.write("Terminating integration: "
                    "dTdt variation less than specified threshold.");
            return true;
        }
    } else if (options.terminationMeasurement == "moments") {
        // Slowest particle moment: RMS rate of change relative to the
        // moment's peak value. Not converged while a moment is still zero.
        dmatrix ddt = ddtDiff + ddtConv + ddtProd + ddtCross;
        double value = 0;
        for (size_t m=0; m<nMoments; m++) {
            double peak = moments.row(m).abs().maxCoeff();
            double rms = ddt.row(kMoments+m).matrix().norm() /
                sqrt(static_cast<double>(nPoints));
            value = (peak > 0) ? std::max(value, rms / peak) : INFINITY;
        }
        logFile.write(format(
            "max ||dM/dt|| / max|M| = %9.3e. Termination threshold = %9.3e") %
            value % options.termination_momentsTol);
        if (value < options.termination_momentsTol) {
            logFile.write("Terminating integration: "
                    "particle moment variation less than specified threshold.");
            return true;
        }
    }
    logFile.write(format(
        "Continuing integration. t = %8.6f") % (tNow-timeVector[0]));
    return false;
}

void FlameSolver::writeStateFile
(const std::string& fileNameStr, bool errorFile, bool updateDerivatives)
{
    if (stateWriter) {
        if (updateDerivatives) {
            updateChemicalProperties();
            convectionSystem.evaluate();
        }
        stateWriter->eval(fileNameStr, errorFile);
    }
}

void FlameSolver::saveTimeSeriesData(const std::string& filename, bool write)
{
    if (timeseriesWriter) {
        timeseriesWriter->eval(filename, write);
    }
}

void FlameSolver::resizeAuxiliary()
{
    size_t nPointsOld = rho.size();
    grid.setSize(T.size());

    if (nPoints == nPointsOld && !grid.updated) {
        return; // nothing to do
    }

    resizeMappedArrays();

    ddtCross.topRows(2).setZero();
    ddtCross.middleRows(kSpecies, nSpec) *= NaN;
    if (nMoments > 0) {
        ddtCross.middleRows(kMoments, nMoments).setZero();
    }

    rho.setZero(nPoints);
    drhodt.setZero(nPoints);
    Wmx.resize(nPoints);
    mu.resize(nPoints);
    lambda.setZero(nPoints);
    cp.setZero(nPoints);
    jCorr.resize(nPoints);
    sumcpj.setZero(nPoints);
    qDot.resize(nPoints);
    cpSpec.resize(nSpec, nPoints);
    rhoD.resize(nSpec, nPoints);
    Dkt.resize(nSpec, nPoints);
    particleDiameter.setZero(nPoints);
    particleTemperature = T;
    particleDiffusivity.setConstant(nPoints, options.momentDiffusivity);
    wDot.resize(nSpec, nPoints);
    hk.resize(nSpec, nPoints);
    jFick.setZero(nSpec, nPoints);
    jSoret.setZero(nSpec, nPoints);

    grid.jj = nPoints-1;
    grid.updateBoundaryIndices();

    if (nPoints > nPointsOld) {
        for (size_t j=nPointsOld; j<nPoints; j++) {
            useCVODE.push_back(options.chemistryIntegrator == "cvode");

            SourceSystem* system;
            if (useCVODE[j]) {
                system = new SourceSystemCVODE();
            } else {
                system = new SourceSystemQSS();
            }
            // initialize the new SourceSystem
            system->setGas(&gas);
            system->initialize(nSpec, nMoments);
            system->setOptions(options);
            system->setTimers(&reactionRatesTimer, &thermoTimer, &jacobianTimer);
            system->setRateMultiplierFunction(rateMultiplierFunction);
            system->setHeatLossFunction(heatLossFunction);
            system->setRhou(rhou);
            system->setPosition(j, x[j]);
            if (options.quasi2d) {
                system->setupQuasi2d(vzInterp, TInterp);
            }

            // Store the solver and system
            sourceTerms.push_back(system);
        }

    } else {
        // Delete the unwanted solvers and systems
        sourceTerms.erase(sourceTerms.begin()+nPoints, sourceTerms.end());
        useCVODE.erase(useCVODE.begin()+nPoints, useCVODE.end());
    }

    // Resize solution vector for diffusion systems / solvers
    for (size_t k=0; k<nVars; k++) {
        diffusionSolvers[k].resize(nPoints);
        diffusionTerms[k].setGrid(grid);
    }

    convectionSystem.setGrid(grid);
    convectionSystem.resize(nPoints, nSpec, nMoments, state);
    convectionSystem.setLeftBC(Tleft, Yleft, momentsLeft);

    convectionSystem.utwSystem.setStrainFunction(strainfunc);
    convectionSystem.utwSystem.fixedTemperature = options.fixedTemperature;
    convectionSystem.utwSystem.setRhou(rhou);

    if (options.quasi2d) {
        convectionSystem.setupQuasi2D(vzInterp, vrInterp);
    }

    // Resize the jCorr stabilizer
    jCorrSolver.resize(nPoints);
    jCorrSystem.setGrid(grid);

    // Set the grid position for each of the source solvers
    for (size_t j=0; j<nPoints; j++) {
        sourceTerms[j].setPosition(j, x[j]);
    }
}

void FlameSolver::resizeMappedArrays()
{
    resize(nVars, nPoints);
    remap(state, T, nPoints, kEnergy);
    remap(state, U, nPoints, kMomentum);
    remap(state, Y, nSpec, nPoints, kSpecies);
    remap(state, moments, nMoments, nPoints, kMoments);
}

void FlameSolver::updateCrossTerms()
{
    assert(mathUtils::notnan(state));
    assert(mathUtils::notnan(rhoD));

    for (size_t j=0; j<jj; j++) {
        jCorr[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            jFick(k,j) = -0.5*(rhoD(k,j)+rhoD(k,j+1)) * ((Y(k,j+1)-Y(k,j))/hh[j]);
            jSoret(k,j) = -0.5*(Dkt(k,j)/T(j) + Dkt(k,j+1)/T(j+1))
                * (T(j+1)-T(j))/hh[j];
            jCorr[j] -= jFick(k,j) + jSoret(k,j);
        }
    }
    jCorr[jj] = 0;

    assert(mathUtils::notnan(jFick));
    assert(mathUtils::notnan(jSoret));
    assert(mathUtils::notnan(lambda));
    assert(mathUtils::notnan(rho));
    assert(mathUtils::notnan(cp));

    // Add a bit of artificial diffusion to jCorr to improve stability
    jCorrSolver.y = jCorr;
    jCorrSystem.B = lambda / (rho * cp); // treat as Le = 1
    jCorrSystem.D.setConstant(nPoints, 1.0);
    jCorrSystem.splitConst.setZero(nPoints);

    double dt = options.globalTimestep;
    for (size_t j=1; j<jj; j++) {
        dt = std::min(dt, options.diffusionTimestepMultiplier*dlj[j]*dlj[j]/(jCorrSystem.B[j]));
    }

    jCorrSolver.initialize(0, dt);
    jCorrSolver.integrateToTime(options.globalTimestep);
    assert(mathUtils::notnan(jCorrSolver.y));

    jCorr = jCorrSolver.y;

    // dYdt due to gradients in other species and temperature
    Eigen::Block<dmatrix> dYdtCross = ddtCross.middleRows(kSpecies, nSpec);

    // dTdt due to gradients in species composition
    Eigen::Block<dmatrix, 1> dTdtCross = ddtCross.row(kEnergy);

    dYdtCross.col(0).setZero();
    dYdtCross.col(jj).setZero();
    for (size_t j=1; j<jj; j++) {
        sumcpj[j] = 0;
        for (size_t k=0; k<nSpec; k++) {
            dYdtCross(k,j) = -0.5 / (r[j] * rho[j] * dlj[j]) *
                (rphalf[j] * (Y(k,j) + Y(k,j+1)) * jCorr[j] -
                 rphalf[j-1] * (Y(k,j-1) + Y(k,j)) * jCorr[j-1]);
            dYdtCross(k,j) -= 1 / (r[j] * rho[j] * dlj[j]) *
                (rphalf[j] * jSoret(k,j) - rphalf[j-1] * jSoret(k,j-1));
            sumcpj[j] += 0.5*(cpSpec(k,j) + cpSpec(k,j+1)) / W[k] *
                (jFick(k,j) + jSoret(k,j) + 0.5 * (Y(k,j) + Y(k,j+1)) * jCorr[j]);
        }
        double dTdx = cfm[j] * T(j-1) + cf[j] * T(j) + cfp[j] * T(j+1);
        if (!options.quasi2d && !options.fixedTemperature) {
            dTdtCross[j] = - 0.5 * (sumcpj[j] + sumcpj[j-1]) * dTdx / (cp[j] * rho[j]);
        }
    }

    // Burner / impinging jet inlet: the control volume at j = 0 (as in the
    // InletFlux diffusion stencil) receives what leaves j = 1 across the face
    // between them, with no diffusive flux through x[0]. Without this the
    // correction flux leaving j = 1 toward the inlet is lost, which matters
    // when the flame sits on the burner (steep H2 gradient at the face).
    if (grid.leftBC == BoundaryCondition::InletFlux) {
        double c0 = (grid.alpha + 1) / (rho[0] * hh[0]);
        for (size_t k=0; k<nSpec; k++) {
            dYdtCross(k,0) = -c0 *
                (0.5 * (Y(k,0) + Y(k,1)) * jCorr[0] + jSoret(k,0));
        }
    }

    // Impinging jet wall: the half control volume at jj receives what
    // leaves j = jj-1 across the face between them (none leaves through the
    // wall), so that these fluxes are conserved like the Fickian ones.
    if (grid.rightBC == BoundaryCondition::Wall) {
        double c0 = 1 / (r[jj] * rho[jj] * 0.5 * hh[jj-1]);
        for (size_t k=0; k<nSpec; k++) {
            dYdtCross(k,jj) = c0 * rphalf[jj-1] *
                (0.5 * (Y(k,jj-1) + Y(k,jj)) * jCorr[jj-1] + jSoret(k,jj-1));
        }
    }

    assert(mathUtils::notnan(ddtCross));
    assert(mathUtils::notnan(sumcpj));
}

void FlameSolver::updateBC()
{
    BoundaryCondition::BC leftPrev = grid.leftBC;
    BoundaryCondition::BC rightPrev = grid.rightBC;

    if (options.impingingJet) {
        // Fixed domain: premixed inlet at x = 0 with the burner-style flux
        // balance (Cantera's Inlet1D), see BoundaryCondition::InletFlux
        grid.leftBC = BoundaryCondition::InletFlux;
    } else if (options.wallFlux && x[0] >= 0.0 && x[0] <= options.centerGridMin) {
        grid.leftBC = BoundaryCondition::WallFlux;
    } else if (grid.ju == 0 &&
               options.continuityBC == ContinuityBoundaryCondition::Left &&
               options.fixedLeftLoc &&
               !options.twinFlame && !options.cylindricalFlame) {
        // Burner-stabilized inlet: unburned stream enters at x = 0 with a
        // prescribed mass flux, and the user has declared x = 0 a fixed
        // physical face (fixedLeftLocation). Use a flux-balance condition so
        // the face composition can differ from the supplied composition, as
        // in Cantera's BurnerFlame -- see BoundaryCondition::InletFlux.
        grid.leftBC = BoundaryCondition::InletFlux;
    } else if (grid.ju == 0 || (grid.jb == 0 && grid.fixedBurnedVal)) {
        grid.leftBC = BoundaryCondition::FixedValue;
    } else if ((options.twinFlame || options.cylindricalFlame) &&
                x[0] >= 0.0 && x[0] <= options.centerGridMin) {
        grid.leftBC = BoundaryCondition::ControlVolume;
    } else {
        grid.leftBC = BoundaryCondition::ZeroGradient;
    }
    
    if (options.impingingJet) {
        grid.rightBC = BoundaryCondition::Wall;
    } else if (options.flameType == "premixed" && grid.jb == jj && !grid.fixedBurnedVal) {
        grid.rightBC = BoundaryCondition::Floating;
    } else {
        grid.rightBC = BoundaryCondition::FixedValue;
    }

    if (leftPrev != grid.leftBC) {
        logFile.write(format("updateBC: Left BC changed from %i to %i.") %
                      leftPrev % grid.leftBC);
    }

    if (rightPrev != grid.rightBC) {
        logFile.write(format("updateBC: Right BC changed from %i to %i.") %
                      rightPrev % grid.rightBC);
    }
}

void FlameSolver::updateChemicalProperties()
{
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nPoints,1),
                      TbbWrapper<FlameSolver>(&FlameSolver::updateChemicalProperties, this));
}

void FlameSolver::updateChemicalProperties(size_t j1, size_t j2)
{
    CanteraGas& gas = gases.local();
    if (!gas.initialized()) {
        gas.setOptions(options);
        gas.initialize();
    }

    // Calculate auxiliary data
    for (size_t j=j1; j<j2; j++) {
        // Thermodynamic properties
        thermoTimer.start();
        gas.setStateMass(&Y(0,j), T(j));
        rho[j] = gas.getDensity();
        Wmx[j] = gas.getMixtureMolecularWeight();
        cp[j] = gas.getSpecificHeatCapacity();
        gas.getSpecificHeatCapacities(&cpSpec(0,j));
        gas.getEnthalpies(&hk(0,j));
        thermoTimer.stop();

        // Transport Properties
        transportTimer.start();

        conductivityTimer.start();
        lambda[j] = gas.getThermalConductivity();
        conductivityTimer.stop();

        viscosityTimer.start();
        mu[j] = gas.getViscosity();
        viscosityTimer.stop();

        diffusivityTimer.start();
        gas.getWeightedDiffusionCoefficientsMass(&rhoD(0,j));
        gas.getThermalDiffusionCoefficients(&Dkt(0,j));
        diffusivityTimer.stop();
        transportTimer.stop();
    }
}

void FlameSolver::updateParticleDiameter()
{
    if (nMoments == 0) {
        return;
    }

    for (size_t j = 0; j < nPoints; j++) {
        double mO = (nMoments > kO) ? moments(kO, j) : 0.0;
        particleDiameter[j] = computeParticleDiameter(
            moments(kN, j), moments(kM, j), mO, options);
    }

    particleTemperature = T;
    if (nMoments > kH) {
        for (size_t j = 0; j < nPoints; j++) {
            particleTemperature[j] = ::particleTemperature(
                moments(kH, j), moments(kM, j), moments(kO, j), T(j), options);
        }
    }
}

void FlameSolver::updateParticleDiffusivity()
{
    if (nMoments == 0) {
        return;
    }

    updateParticleDiameter();

    const double kB_ = Cantera::Boltzmann;

    for (size_t j = 0; j < nPoints; j++) {
        double dp = particleDiameter[j]; // [m]

        if (!(dp > 0)) {
            // No particles present (yet) at this point -- fall back to the
            // user-specified constant so the equation stays well posed.
            particleDiffusivity[j] = options.momentDiffusivity;
            continue;
        }

        // Particles smaller than minParticleDiameter disintegrate, so a
        // smaller mean size can only come from N/mass noise where both are
        // near zero; since D_p ~ 1/dp^2 it would blow up the diffusivity.
        dp = std::max(dp, options.minParticleDiameter);

        double meanFreePath = gasMeanFreePath(mu[j], T[j], Wmx[j], options.pressure);
        double Kn = 2.0 * meanFreePath / dp;
        double Cc = cunninghamSlipCorrection(Kn);

        // Stokes-Einstein-Cunningham particle diffusivity.
        particleDiffusivity[j] = kB_ * T[j] * Cc / (3.0 * M_PI * mu[j] * dp);
    }
}

void FlameSolver::setDiffusionSolverState(double tInitial)
{
    splitTimer.resume();
    size_t k = 0;
    for (TridiagonalIntegrator& integrator : diffusionSolvers) {
        double dt = (dlj.square().tail(nPoints - 2) /
            (diffusionTerms[k].B * diffusionTerms[k].D).segment(1, nPoints-2)).minCoeff();
        dt = std::min(options.diffusionTimestepMultiplier * dt, options.globalTimestep);
        integrator.initialize(tInitial, dt);
        k++;
    }

    for (size_t i=0; i<nVars; i++) {
        diffusionSolvers[i].y = state.row(i);
    }
    splitTimer.stop();
}

void FlameSolver::setConvectionSolverState(double tInitial)
{
    splitTimer.resume();
    convectionSystem.setState(tInitial);
    splitTimer.stop();
}

void FlameSolver::setProductionSolverState(double tInitial)
{
    splitTimer.resume();
    for (size_t j=0; j<nPoints; j++) {
        sourceTerms[j].setState(tInitial, U(j), T(j), Y.col(j), moments.col(j));
    }
    splitTimer.stop();
}

void FlameSolver::integrateConvectionTerms()
{
    setConvectionSolverState(tStageStart);
    convectionTimer.start();
    try {
        convectionSystem.integrateToTime(tStageEnd);
    } catch (DebugException& e) {
        logFile.write(e.errorString);
        writeStateFile("err_convectionIntegration", true, false);
        throw;
    }
    convectionTimer.stop();

    splitTimer.resume();
    convectionSystem.unroll_y();
    splitTimer.stop();
}

void FlameSolver::integrateProductionTerms()
{
    setProductionSolverState(tStageStart);
    reactionTimer.start();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nPoints,1),
                     TbbWrapper<FlameSolver>(&FlameSolver::integrateProductionTerms, this));
    reactionTimer.stop();
}

void FlameSolver::integrateProductionTerms(size_t j1, size_t j2)
{
    CanteraGas& gas = gases.local();
    if (!gas.initialized()) {
        gas.setOptions(options);
        gas.initialize();
    }

    int err = 0;
    for (size_t j=j1; j<j2; j++) {
        if ((j == 0 && grid.leftBC == BoundaryCondition::InletFlux) ||
            (j == jj && grid.rightBC == BoundaryCondition::Wall)) {
            // Burner face: no chemistry, as in Cantera's burner boundary. T is
            // held at Tleft, and species/moments are set only by the inlet
            // flux balance. Letting chemistry run here (with T free inside
            // the stage and Y never reset) turns the face into an igniter.
            // Likewise at an inert wall (Cantera's Surface1D).
            double dtStage = tStageEnd - tStageStart;
            U(j) += splitConstProd(kMomentum, j) * dtStage;
            Y.col(j) += splitConstProd.col(j).segment(kSpecies, nSpec) * dtStage;
            if (nMoments > 0) {
                moments.col(j) += splitConstProd.col(j).segment(kMoments, nMoments) * dtStage;
            }
            continue;
        }
        SourceSystem& system = sourceTerms[j];
        logFile.verboseWrite(format("%i") % j, false);
        system.setGas(&gas);
        if (int(j) == options.debugSourcePoint &&
            tStageEnd >= options.debugSourceTime) {
            system.setDebug(true);
            std::ofstream steps("sourceTermSteps.py");
            system.writeState(steps, true);

            while (system.time() < tStageEnd - tStageStart && err >= 0) {
                err = system.integrateOneStep(tStageEnd - tStageStart);
                system.writeState(steps, false);
            }

            system.writeJacobian(steps);
            steps.close();
            std::terminate();

        } else {
            err = system.integrateToTime(tStageEnd - tStageStart);
        }
        if (err >= 0) {
            logFile.verboseWrite(format(" [%s]...") % system.getStats(), false);
            sourceTerms[j].unroll_y();
            U(j) = sourceTerms[j].U;
            T(j) = sourceTerms[j].T;
            Y.col(j) = sourceTerms[j].Y;
            if (nMoments > 0) {
                moments.col(j) = sourceTerms[j].moments;
            }
        } else {
            // Print gas mole fractions to help identify problematic reactions
            dvec X0(nSpec), X1(nSpec);
            gas.thermo->setMassFractions_NoNorm(Y.col(j).data());
            gas.getMoleFractions(X0);
            gas.thermo->setMassFractions_NoNorm(sourceTerms[j].Y.data());
            gas.getMoleFractions(X1);
            logFile.write(format("Error at j = %i. Gas state:\n") % j);
            logFile.write(" k        Species     X Initial    X Final     Delta X");
            logFile.write("----  --------------  ----------  ----------  ----------");
            logFile.write(format("      %14s  %10.4g  %10.4g  %10.4g") %
                "T" % T(j) % sourceTerms[j].T % (sourceTerms[j].T - T(j)));
            for (size_t k = 0; k < nSpec; k++) {
                if (std::abs(X0[k]) > 1e-4 ||
                    std::abs(X1[k]) > 1e-4 ||
                    std::abs(X1[k]-X0[k]) > 1e-6) {
                    logFile.write(format("%4s  %14s  %10.4g  %10.4g  %10.4g") %
                        k % gas.thermo->speciesName(k) % X0[k] % X1[k] % (X1[k]-X0[k]));
                }
            }
            if (debugParameters::veryVerbose) {
                logFile.write(format("\nT = %s") % system.T);
                logFile.write(format("U = %s") % system.U);
                logFile.write("Y = ", false);
                Eigen::IOFormat fmt(15, Eigen::DontAlignCols, ", ", ", ", "", "", "[", "]");
                logFile.write(system.Y.format(fmt));
                if (options.nThreads == 1) {
                    // This diagnostic file can only be written when running with a
                    // single thread to avoid calling Python from threads that were
                    // not initialized to use Python
                    writeStateFile((format("prod%i_error_t%.6f_j%03i") %
                            1 % tStageEnd % j).str(), true, false);
                }
            }
        }
    }
}

void FlameSolver::integrateDiffusionTerms()
{
    setDiffusionSolverState(tStageStart);
    diffusionTimer.start();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nVars, 1),
                      TbbWrapper<FlameSolver>(&FlameSolver::integrateDiffusionTerms, this));
    diffusionTimer.stop();
}

void FlameSolver::integrateDiffusionTerms(size_t k1, size_t k2)
{
    for (size_t k=k1; k<k2; k++) {
        diffusionSolvers[k].integrateToTime(tStageEnd);
        assert(mathUtils::almostEqual(diffusionSolvers[k].t, tStageEnd));
        state.row(k) = diffusionSolvers[k].y;
    }
}

void FlameSolver::rollVectorVector
(vector<dvector>& vv, const dmatrix& M) const
{
    size_t N = vv.size();
    vv.resize(N + nVars, dvector(nPoints));

    for (size_t k=0; k<nVars; k++) {
        Eigen::Map<dvec>(&vv[N+k][0], nPoints) = M.row(k);
    }
}


void FlameSolver::unrollVectorVector
(vector<dvector>& vv, dmatrix& M, size_t i) const
{
    for (size_t k=0; k<nVars; k++) {
        M.row(k) = Eigen::Map<dvec>(&vv[i*nVars+k][0], nPoints);
    }
}


void FlameSolver::update_xStag(const double t, const bool updateIntError)
{
    calculateQdot();
    xFlameActual = getFlamePosition();
    xFlameTarget = targetFlamePosition(t);
    if (updateIntError) {
        flamePosIntegralError += (xFlameTarget-xFlameActual)*(t-tFlamePrev);
        tFlamePrev = t;
    }

    // controlSignal is approximately a*xStag
    double controlSignal = options.xFlameProportionalGain *
        ((xFlameTarget - xFlameActual) +
        (flamePosIntegralError + (xFlameTarget - xFlameActual) * (t - tFlamePrev)) *
        options.xFlameIntegralGain);

    if (debugParameters::debugFlameRadiusControl) {
        logFile.write(format(
            "rFlameControl: rF=%g;  control=%g;  P=%g;  I=%g;  dt=%g") %
            xFlameActual %
            controlSignal %
            (options.xFlameProportionalGain * (xFlameTarget - xFlameActual)) %
            (options.xFlameProportionalGain * flamePosIntegralError *
            options.xFlameIntegralGain) %
            (t - tFlamePrev));
    }

    double a = strainfunc->a(t);
    if (alpha == 1) {
        rVzero = 0.5*rhoLeft*(controlSignal*abs(controlSignal)-a*x[0]*x[0]);
    } else {
        rVzero = rhoLeft*(controlSignal-a*x[0]);
    }
}


void FlameSolver::updateMassFluxControl(const double t)
{
    // Inlet mass flux from the latest evaluation of the continuity equation,
    // which is integrated from the wall (V = rV for planar and disc flames)
    double mdot = convectionSystem.utwSystem.V[0];
    double error = (options.massFluxTarget - mdot) / options.massFluxTarget;
    massFluxIntegralError += error * (t - tMassFluxPrev);
    tMassFluxPrev = t;

    // The mass flux scales roughly linearly with a, so the controller acts
    // on log(a); this also keeps a positive.
    double a = options.strainRateInitial *
        exp(options.massFluxProportionalGain * error +
            options.massFluxIntegralGain * massFluxIntegralError);
    static_cast<ControlledFunction*>(strainfunc)->value = a;

    if (debugParameters::debugFlameRadiusControl) {
        logFile.write(format("massFluxControl: mdot=%g;  a=%g;  error=%g;  I=%g") %
                      mdot % a % error % massFluxIntegralError);
    }
}

double FlameSolver::targetFlamePosition(double t)
{
    if (t <= options.xFlameT0) {
        return options.xFlameInitial;
    } else if (t >= options.xFlameT0 + options.xFlameDt) {
        return options.xFlameFinal;
    } else {
        return options.xFlameInitial + (options.xFlameFinal - options.xFlameInitial) *
            (t - options.xFlameT0) / options.xFlameDt;
    }
}


void FlameSolver::calculateQdot()
{
    reactionRatesTimer.start();
    if (rateMultiplierFunction) {
        gas.setRateMultiplier(rateMultiplierFunction->a(tNow));
    }
    for (size_t j=0; j<nPoints; j++) {
        gas.setStateMass(&Y(0,j), T(j));
        gas.getEnthalpies(&hk(0,j));
        gas.getReactionRates(&wDot(0,j));
        qDot[j] = - (wDot.col(j) * hk.col(j)).sum();
    }
    // No chemistry is integrated at the burner face (see
    // integrateProductionTerms), so don't report heat release there either
    if (grid.leftBC == BoundaryCondition::InletFlux) {
        wDot.col(0).setZero();
        qDot[0] = 0;
    }
    if (grid.rightBC == BoundaryCondition::Wall) {
        wDot.col(jj).setZero();
        qDot[jj] = 0;
    }
    reactionRatesTimer.stop();
}


void FlameSolver::correctMassFractions() {
    setupTimer.resume();
    for (size_t j=0; j<nPoints; j++) {
        // Density is an intensive property of the gas phase and must be
        // computed from a properly normalized composition, so Y is always
        // renormalized to sum to 1 here regardless of nucleation activity.
        // Species consumed by nucleation are still depleted correctly by
        // their own (negative) production term in computeNucleationRates();
        // this only removes numerical round-off drift.
        gas.setStateMass(&Y(0,j), T(j));
        gas.getMassFractions(&Y(0,j));
    }
    setupTimer.stop();
}

double FlameSolver::getHeatReleaseRate(void)
{
    return mathUtils::integrate(x, qDot);
}

double FlameSolver::getConsumptionSpeed(void)
{
    double QoverCp = mathUtils::integrate(x, qDot / cp);
    double rhouDeltaT = rhou*(T(grid.jb)-T(grid.ju));
    return QoverCp/rhouDeltaT;
}

double FlameSolver::getFlamePosition(void)
{
    return mathUtils::trapz(x, (x*qDot).eval())/mathUtils::trapz(x,qDot);
}

void FlameSolver::loadProfile(void)
{
    grid.unburnedLeft = options.unburnedLeft;

    // Read initial condition specified in the configuration file
    x = options.x_initial;
    nPoints = x.size();
    resizeMappedArrays();

    U = options.U_initial;
    T = options.T_initial;
    if (options.Y_initial.rows() == static_cast<dmatrix::Index>(x.size())) {
        Y = options.Y_initial.transpose();
    } else {
        Y = options.Y_initial;
    }
    if (nMoments > 0) {
        moments.setZero();
    }
    convectionSystem.V = options.V_initial;
    convectionSystem.utwSystem.V = options.V_initial;
    rVzero = convectionSystem.utwSystem.V[0];

    grid.setSize(x.size());
    grid.updateValues();
    grid.updateBoundaryIndices();

    tNow = (options.haveTStart) ? options.tStart : 0.0;

    if (options.flameType == "premixed") {
        gas.setStateMass(&Y(0,grid.ju), T(grid.ju));
        rhou = gas.getDensity();
        gas.setStateMass(&Y(0,grid.jb), T(grid.ju));
        rhob = gas.getDensity();

        if (options.unburnedLeft) {
            rhoLeft = rhou;
            Tleft = T(grid.ju);
            Yleft = Y.col(grid.ju);
            rhoRight = rhob;
            Tright = T(grid.jb);
            Yright = Y.col(grid.jb);
        } else {
            rhoLeft = rhob;
            Tleft = T(grid.jb);
            Yleft = Y.col(grid.jb);
            rhoRight = rhou;
            Tright = T(grid.ju);
            Yright = Y.col(grid.ju);
        }

    } else if (options.flameType == "diffusion") {
        // Fuel composition
        size_t jFuel = (options.fuelLeft) ? 0 : jj;
        gas.thermo->setState_TPY(T(jFuel), options.pressure, &Y(0,jFuel));
        double rhoFuel = gas.getDensity();
        dvec Yfuel(nSpec);
        gas.getMassFractions(Yfuel);
        double Tfuel = gas.thermo->temperature();

        // Oxidizer composition
        size_t jOxidizer = (options.fuelLeft) ? jj : 0;
        gas.thermo->setState_TPY(T(jOxidizer),
                                options.pressure,
                                &Y(0,jOxidizer));
        double rhoOxidizer = gas.getDensity();
        dvec Yoxidizer(nSpec);
        gas.getMassFractions(Yoxidizer);
        double Toxidizer = gas.thermo->temperature();

        if (options.fuelLeft) {
            rhoLeft = rhoFuel;
            Tleft = Tfuel;
            Yleft = Yfuel;
            rhoRight = rhoOxidizer;
            Tright = Toxidizer;
            Yright = Yoxidizer;
        } else {
            rhoLeft = rhoOxidizer;
            Tleft = Toxidizer;
            Yleft = Yoxidizer;
            rhoRight = rhoFuel;
            Tright = Tfuel;
            Yright = Yfuel;
        }

        rhou = rhoOxidizer;

    } else if (options.flameType == "quasi2d") {
        gas.thermo->setState_TPY(T(0), options.pressure, &Y(0,0));
        rhoLeft = gas.thermo->density();
        Tleft = T(0);
        Yleft.resize(nSpec);
        gas.getMassFractions(Yleft);
        gas.thermo->setState_TPY(T(jj), options.pressure, &Y(0,jj));
        rhoRight = gas.thermo->density();
        Tright = T(jj);
        Yright.resize(nSpec);
        gas.getMassFractions(Yright);

        rhou = rhoRight;

    } else {
        throw DebugException("Invalid flameType: " + options.flameType);
    }

    if (options.impingingJet) {
        Tright = options.wallTemperature;
    }

    updateBC();

    if (grid.leftBC == BoundaryCondition::ControlVolume &&
        options.xFlameControl)
    {
        double controlSignal;
        if (alpha == 0) {
            controlSignal = rVcenter/rhoLeft;
        } else {
            double tmp = pow(x[0],2) + 2*rVcenter/rhoLeft;
            controlSignal = mathUtils::sign(tmp)*sqrt(abs(tmp));
        }
        flamePosIntegralError = controlSignal /
            (options.xFlameProportionalGain * options.xFlameIntegralGain);
    }
}

void FlameSolver::printPerformanceStats(void)
{
    std::string filename = options.outputDir + "/stats";
    std::ofstream stats(filename.c_str(), std::ios::trunc | std::ios::out);
    totalTimer.stop();
    totalTimer.resume();
    stats << "\n   *** Performance Stats ***       time   ( call count )\n";
    printPerfString(stats, "                General Setup: ", setupTimer);
    printPerfString(stats, "             Split Term Setup: ", splitTimer);
    printPerfString(stats, "              Grid Adaptation: ", regridTimer);
    printPerfString(stats, "    Reaction Term Integration: ", reactionTimer);
    printPerfString(stats, "   Diffusion Term Integration: ", diffusionTimer);
    printPerfString(stats, "  Convection Term Integration: ", convectionTimer);
    printPerfString(stats, "                        Total: ", totalTimer);
    stats << "\n Subcomponents:\n";
    printPerfString(stats, "               Reaction Rates: ", reactionRatesTimer);
    printPerfString(stats, "         Transport Properties: ", transportTimer);
    printPerfString(stats, "          - thermal cond.    : ", conductivityTimer);
    printPerfString(stats, "          - viscosity        : ", viscosityTimer);
    printPerfString(stats, "          - diffusion coeff. : ", diffusivityTimer);
    printPerfString(stats, "     Thermodynamic Properties: ", thermoTimer);
    printPerfString(stats, "   Source Jacobian Evaluation: ", jacobianTimer);
    printPerfString(stats, "   UTW Convection Integration: ", convectionSystem.utwTimer);
    printPerfString(stats, "    Yk Convection Integration: ", convectionSystem.speciesTimer);
}

void FlameSolver::printPerfString(std::ostream& stats, const std::string& label,
                                  const PerfTimer& T)
{
    if (T.getTime() > 0.0) {
        stats << format("%s %9.3f (%12i)\n") % label % T.getTime() % T.getCallCount();
    }
}
