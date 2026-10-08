#cython: embedsignature=True, language_level=3
#distutils: language = c++

import numpy as np
cimport numpy as np
import os
import sys
import cantera

from cython.operator cimport dereference as deref
from ._ember cimport *

cdef int _pythonMajorVersion = sys.version_info[0]

cdef string stringify(x) except *:
    """ Converts Python strings to std::string. """
    # This method works with both Python 2.x and 3.x.
    if isinstance(x, bytes):
        return string(<bytes>x)
    else:
        tmp = bytes(x.encode())
        return string(tmp)

cdef pystr(string x):
    cdef bytes s = x.c_str()
    if _pythonMajorVersion == 2:
        # Python 2.x
        return s
    else:
        # Python 3.x
        return s.decode()

cdef np.ndarray[np.double_t, ndim=1] getArray_vector(vector[double]& vec):
    cdef np.ndarray[np.double_t, ndim=1] v = np.empty(vec.size())
    cdef int i
    for i in range(len(v)):
        v[i] = vec[i]
    return v

cdef np.ndarray[np.double_t, ndim=1] getArray_Vec(CxxEigenVec& vec):
    cdef np.ndarray[np.double_t, ndim=1] v = np.empty(vec.size())
    cdef int i
    for i in range(len(v)):
        v[i] = vec.item(i)
    return v

cdef np.ndarray[np.double_t, ndim=2] getArray_Matrix(CxxEigenMatrix& M):
    cdef np.ndarray[np.double_t, ndim=2] v = np.empty((M.rows(), M.cols()))
    cdef int i, j
    for i in range(M.rows()):
        for j in range(M.cols()):
            v[i,j] = M.item(i,j)
    return v

cdef np.ndarray[np.double_t, ndim=1] getArray_MatrixRow(CxxEigenMatrix& M, int i):
    cdef np.ndarray[np.double_t, ndim=1] v = np.empty(M.cols())
    cdef int j
    for j in range(M.cols()):
        v[j] = M.item(i,j)
    return v

cdef np.ndarray[np.double_t, ndim=1] getArray_VecMap(CxxEigenVecMap& vec):
    cdef np.ndarray[np.double_t, ndim=1] v = np.empty(vec.size())
    cdef int i
    for i in range(len(v)):
        v[i] = vec.item(i)
    return v

cdef np.ndarray[np.double_t, ndim=2] getArray_MatrixMap(CxxEigenMatrixMap& M):
    cdef np.ndarray[np.double_t, ndim=2] v = np.empty((M.rows(), M.cols()))
    cdef int i, j
    for i in range(M.rows()):
        for j in range(M.cols()):
            v[i,j] = M.item(i,j)
    return v


cdef class LoggerCallback:
    """
    A wrapper for functions that are used to extract information from a
    FlameSolver object at specific points during the integration, e.g. to save
    periodic output files.
    """
    def __cinit__(self, func):
        """
        :param func:
            A callable object that takes the default output file name as an
            argument.
        """
        self.callback = new CxxLoggerCallback(logger_func_callback, <void*>self)
        self.func = func
        self.exception = None

    def __dealloc__(self):
        del self.callback

    def eval(self, name, flag):
        self.func(name, flag)


cdef void logger_func_callback(const string& name, int flag,
                               void* obj, void** err) noexcept nogil:
    """
    This function is called from C/C++ to evaluate a `LoggerCallback` object
    *obj*. If an exception occurs while evaluating the function, the Python
    exception info is saved in the two-element array *err*.
    """
    with gil:
        try:
            (<LoggerCallback>obj).eval(pystr(name), flag)
        except BaseException as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()

            # Stash the exception info to prevent it from being garbage collected
            (<LoggerCallback>obj).exception = exc_type, exc_value, exc_traceback
            err[0] = <void*>exc_type
            err[1] = <void*>exc_value
            err[2] = <void*>exc_traceback


cdef class IntegratorCallback:
    """
    A wrapper for functions that are evaluate a function of the form
    z = f(x, t, U, T, Y) where x, t, U, T and z are scalars, and Y is a vector.
    """
    def __cinit__(self, func):
        """
        :param func:
            A callable object with the signature
            `float = f(float, float, ndarray)`
        """
        self.callback = new CxxIntegratorCallback(integrator_func_callback,
                                                  <void*>self)
        self.func = func
        self.exception = None

    def __dealloc__(self):
        del self.callback

    def eval(self, x, t, U, T, Y):
        return self.func(x, t, U, T, Y)


cdef double integrator_func_callback(double x, double t, double U, double T,
                                     CxxEigenVec& y, void* obj, void** err) noexcept nogil:
    """
    This function is called from C/C++ to evaluate a `IntegratorCallback` object
    *obj*. If an exception occurs while evaluating the function, the Python
    exception info is saved in the two-element array *err*.
    """
    with gil:
        try:
            return (<IntegratorCallback>obj).eval(x, t, U, T, getArray_Vec(y))
        except BaseException as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()

            # Stash the exception info to prevent it from being garbage collected
            (<IntegratorCallback>obj).exception = exc_type, exc_value, exc_traceback
            err[0] = <void*>exc_type
            err[1] = <void*>exc_value
            err[2] = <void*>exc_traceback
            return np.nan

def addCanteraDirectory(dirname):
    CxxAddCanteraDirectory(stringify(dirname))


def _speciesCollisionProperties(gas, name):
    """(index, collision diameter [m], molecular mass [kg/molecule]) for *name*."""
    k = gas.species_index(name)
    transport = gas.species(k).transport
    if transport is None or transport.diameter is None:
        raise ValueError("Species '%s' has no transport data (collision"
                         " diameter) in the mechanism." % name)
    d = transport.diameter
    m = gas.molecular_weights[k] / cantera.avogadro
    return k, d, m


def writelog(text):
    CxxSingletonLogfile.write(stringify(text))


cdef class ConfigOptions:
    def __cinit__(self, *args, **kwargs):
        self.opts = new CxxConfigOptions()

    def __dealloc__(self):
        del self.opts

    def apply_options(self):
        cdef CxxConfigOptions* opts = self.opts

        def get(option, default):
            if option is None:
                return default, False
            else:
                return option, True

        # Paths
        opts.outputDir = stringify(self.paths.outputDir)
        if self.paths.logFile is not None:
            CxxSingletonLogfile.open(stringify(self.paths.logFile))

        # General
        G = self.general
        opts.fixedBurnedVal = G.fixedBurnedVal
        opts.fixedLeftLoc = G.fixedLeftLocation
        opts.unburnedLeft = G.unburnedLeft
        opts.fuelLeft = G.fuelLeft
        opts.nThreads = G.nThreads
        opts.cylindricalFlame = G.flameGeometry == 'cylindrical'
        opts.discFlame = G.flameGeometry == 'disc'
        opts.twinFlame = G.twinFlame
        opts.chemistryIntegrator = stringify(G.chemistryIntegrator)
        opts.splittingMethod = stringify(G.splittingMethod)
        opts.setContinuityBC(stringify(G.continuityBC))
        opts.errorStopCount = G.errorStopCount
        opts.stopIfError = G.errorStopCount > 0

        # Chemistry
        opts.gasMechanismFile = stringify(self.chemistry.mechanismFile)
        opts.gasPhaseID = stringify(self.chemistry.phaseID)
        opts.transportModel = stringify(self.chemistry.transportModel)
        opts.kineticsModel = stringify(self.chemistry.kineticsModel)
        opts.transportThreshold = self.chemistry.threshold
        if self.chemistry.rateMultiplierFunction is not None:
            opts.rateMultiplierFunctionType = stringify('chebyshev')

        # Initial condition
        IC = self.initialCondition
        cdef np.ndarray[np.double_t, ndim=1] data
        cdef np.ndarray[np.double_t, ndim=2] Y

        data = np.ascontiguousarray(IC.x)
        opts.x_initial = map_vector(&data[0], len(data), 1)

        data = np.ascontiguousarray(IC.T)
        opts.T_initial = map_vector(&data[0], len(data), 1)

        data = np.ascontiguousarray(IC.U)
        opts.U_initial = map_vector(&data[0], len(data), 1)

        data = np.ascontiguousarray(IC.V)
        opts.V_initial = map_vector(&data[0], len(data), 1)

        # Numpy defaults to row major; Eigen is column major
        Y = np.ascontiguousarray(IC.Y.T)
        w, h = Y.shape[0], Y.shape[1]
        opts.Y_initial = map_matrix(&Y[0,0], h, w, h, 1)

        opts.flameType = stringify(IC.flameType)
        opts.pressure = IC.pressure
        opts.quasi2d = IC.flameType == 'quasi2d'

        # Wall flux boundary condition
        if self.wallFlux:
            opts.wallFlux = True
            opts.Tinf = self.wallFlux.Tinf
            opts.Kwall = self.wallFlux.Kwall
        else:
            opts.wallFlux = 0
            opts.Tinf = IC.Tu
            opts.Kwall = 0

        # ignition parameters
        opts.ignition_tStart = self.ignition.tStart
        opts.ignition_duration = self.ignition.duration
        opts.ignition_energy = self.ignition.energy
        opts.ignition_center = self.ignition.center
        opts.ignition_stddev = self.ignition.stddev

        # external heat flux
        opts.alwaysUpdateHeatFlux = self.externalHeatFlux.alwaysUpdate

        # strain rate parameters
        if self.strainParameters.function is not None:
            opts.strainFunctionType = stringify('chebyshev')
        else:
            opts.strainFunctionType = stringify('linear')
            opts.strainRateInitial = self.strainParameters.initial
            opts.strainRateFinal = self.strainParameters.final
            opts.strainRateT0 = self.strainParameters.tStart
            opts.strainRateDt = self.strainParameters.dt

        # Flame position control
        if self.positionControl is not None:
            PC = self.positionControl
            opts.xStagControl = True
            opts.xFlameControl = True
            opts.xFlameInitial = PC.xInitial
            opts.xFlameFinal = PC.xFinal
            opts.xFlameT0 = PC.tStart
            opts.xFlameDt = PC.dt
            opts.xFlameIntegralGain = PC.integralGain
            opts.xFlameProportionalGain = PC.proportionalGain
        else:
            opts.xStagControl = False
            opts.xFlameControl = False

        # Impinging jet: inlet at x = 0, wall at the right end of a fixed domain
        if self.impingingJet is not None:
            IJ = self.impingingJet
            opts.impingingJet = True
            opts.wallTemperature = IJ.wallTemperature
            opts.fixedLeftLoc = True
            opts.setContinuityBC(stringify('wall'))
            opts.massFluxControl = IJ.massFlux is not None
            opts.massFluxTarget = IJ.massFlux or 0.0
            opts.massFluxProportionalGain = IJ.proportionalGain
            opts.massFluxIntegralGain = IJ.integralGain
            opts.fixedTemperature = IJ.temperatureProfile is not None
            if opts.fixedTemperature:
                data = np.ascontiguousarray(IJ.temperatureProfile[0], dtype=float)
                opts.Tfixed_x = map_vector(&data[0], len(data), 1)
                data = np.ascontiguousarray(IJ.temperatureProfile[1], dtype=float)
                opts.Tfixed_T = map_vector(&data[0], len(data), 1)
        else:
            opts.impingingJet = False
            opts.massFluxControl = False
            opts.fixedTemperature = False

        # Grid
        opts.centerGridMin = self.grid.centerGridMin
        opts.vtol = self.grid.vtol
        opts.dvtol = self.grid.dvtol
        opts.rmTol = self.grid.rmTol
        opts.dampConst = self.grid.dampConst
        opts.gridMax = self.grid.gridMax
        opts.gridMin = self.grid.gridMin
        opts.uniformityTol = self.grid.uniformityTol
        opts.absvtol = self.grid.absvtol
        opts.boundaryTol = self.grid.boundaryTol
        opts.boundaryTolRm = self.grid.boundaryTolRm
        opts.addPointCount = self.grid.addPointCount
        opts.unstrainedDownstreamWidth = self.grid.unstrainedDownstreamWidth

        # Particles (passive scalar moment transport)
        opts.nMoments = self.particles.nMoments
        opts.momentDiffusivity = self.particles.momentDiffusivity
        opts.momentBCLeft = self.particles.momentBCLeft
        opts.momentBCRight = self.particles.momentBCRight

        # Particle material: metal + oxygen, densities of its phases
        metal = self.particles.metal
        if self.particles.nMoments > 0 and metal not in self.gas.element_names:
            raise ValueError("particles.metal: element '%s' is not in the"
                             " mechanism." % metal)
        phases = [(float(x), float(rho)) for x, rho in self.particles.phases]
        if (not phases or phases[0][0] != 0.0 or
                any(b[0] <= a[0] for a, b in zip(phases, phases[1:])) or
                any(rho <= 0 for _, rho in phases)):
            raise ValueError("particles.phases must be (O/metal ratio, density)"
                             " pairs with ascending ratios starting at 0 (pure"
                             " metal) and positive densities; got %r." % (phases,))
        opts.phaseRatio = [x for x, _ in phases]
        opts.phaseDensity = [rho for _, rho in phases]
        opts.metalWeight = cantera.Element(metal).weight
        opts.oxygenWeight = cantera.Element('O').weight
        opts.coagulation = self.particles.coagulation
        opts.minParticleDiameter = self.particles.minParticleDiameter

        surfReactant, surfProduct = [], []
        surfA, surfTa = [], []
        surfDeltaM, surfDeltaO = [], []
        for reaction in self.particles.surfaceReactions:
            dM, dO = reaction.particleChange(self.gas, metal)
            if dO != 0 and self.particles.nMoments < 3:
                raise ValueError(
                    "Particles.surfaceReactions: '%s -> %s' changes the"
                    " particle's oxygen content, which requires nMoments = 3."
                    % (reaction.reactant.value, reaction.product.value))
            surfReactant.append(self.gas.species_index(reaction.reactant.value))
            if reaction.product.value is None:
                surfProduct.append(-1)
            else:
                surfProduct.append(self.gas.species_index(reaction.product.value))
            surfA.append(reaction.A.value)
            surfTa.append(reaction.Ta.value)
            surfDeltaM.append(dM)
            surfDeltaO.append(dO)
        opts.surfReactant = surfReactant
        opts.surfProduct = surfProduct
        opts.surfA = surfA
        opts.surfTa = surfTa
        opts.surfDeltaM = surfDeltaM
        opts.surfDeltaO = surfDeltaO

        condSpecies, condDiameter = [], []
        condDeltaM, condDeltaO = [], []
        for name in self.particles.condensationSpecies:
            composition = self.gas.species(name).composition
            dM = int(composition.get(metal, 0))
            dO = int(composition.get('O', 0))
            if set(composition) - {metal, 'O'} or dM == 0:
                raise ValueError(
                    "Particles.condensationSpecies: '%s' must contain the"
                    " particle metal '%s' and otherwise only O; found %r." %
                    (name, metal, composition))
            if dO != 0 and self.particles.nMoments < 3:
                raise ValueError(
                    "Particles.condensationSpecies: '%s' contains oxygen,"
                    " which requires nMoments = 3." % name)
            k, d, _ = _speciesCollisionProperties(self.gas, name)
            condSpecies.append(k)
            condDiameter.append(d)
            condDeltaM.append(dM)
            condDeltaO.append(dO)
        opts.condSpecies = condSpecies
        opts.condDiameter = condDiameter
        opts.condDeltaM = condDeltaM
        opts.condDeltaO = condDeltaO

        nucleation = self.particles.nucleation
        if nucleation is not None:
            nucleation.validate(self.gas, metal)

            nucSpeciesA, nucSpeciesB = [], []
            nucStoichA, nucStoichB = [], []
            nucMassPerEvent = []
            nucCollisionPrefactor = []
            stoichA = nucleation.stoichA.value
            stoichB = nucleation.stoichB.value
            for nameA, nameB in zip(nucleation.collisionSpeciesA.value,
                                    nucleation.collisionSpeciesB.value):
                kA, dA, mA = _speciesCollisionProperties(self.gas, nameA)
                kB, dB, mB = _speciesCollisionProperties(self.gas, nameB)
                WA = self.gas.molecular_weights[kA]
                WB = self.gas.molecular_weights[kB]

                dAB = 0.5 * (dA + dB)
                muAB = mA * mB / (mA + mB)
                # Convert kmol/m^3 concentrations to molecule/m^3 and the
                # resulting event rate back to kmol/m^3/s.
                prefactor = (np.pi * dAB**2 *
                             np.sqrt(8.0 * cantera.boltzmann / (np.pi * muAB)) *
                             cantera.avogadro)

                nucSpeciesA.append(kA)
                nucSpeciesB.append(kB)
                nucStoichA.append(stoichA)
                nucStoichB.append(stoichB)
                nucCollisionPrefactor.append(prefactor)
                # Always derived, never user-specified: this is what preserves
                # the exact invariant sum(Y) + mM + mO == 1.
                nucMassPerEvent.append(stoichA * WA + stoichB * WB)

            opts.nucSpeciesA = nucSpeciesA
            opts.nucSpeciesB = nucSpeciesB
            opts.nucStoichA = nucStoichA
            opts.nucStoichB = nucStoichB
            opts.nucMassPerEvent = nucMassPerEvent
            opts.nucCollisionPrefactor = nucCollisionPrefactor

            nucPrecursorSpecies = []
            nucPrecursorAtomCount = []
            nucMonomerIndex = -1
            for name in nucleation.precursorSpecies.value:
                k = self.gas.species_index(name)
                atomCount = int(next(iter(self.gas.species(k).composition.values())))
                if atomCount == 1:
                    nucMonomerIndex = len(nucPrecursorSpecies)
                nucPrecursorSpecies.append(k)
                nucPrecursorAtomCount.append(atomCount)

            opts.nucPrecursorSpecies = nucPrecursorSpecies
            opts.nucPrecursorAtomCount = nucPrecursorAtomCount
            opts.nucMonomerIndex = nucMonomerIndex

            A, B, C = nucleation.antoineCoeffs.value
            opts.nucAntoineA, opts.nucAntoineB, opts.nucAntoineC = A, B, C
            sA, sB = nucleation.surfaceTensionCoeffs.value
            opts.nucSurfaceTensionA, opts.nucSurfaceTensionB = sA, sB


        # Times
        opts.tStart, opts.haveTStart = get(self.times.tStart, 0.0)
        opts.regridTimeInterval = self.times.regridTimeInterval
        opts.regridStepInterval = self.times.regridStepInterval
        opts.outputTimeInterval = self.times.outputTimeInterval
        opts.outputStepInterval = self.times.outputStepInterval
        opts.profileTimeInterval = self.times.profileTimeInterval
        opts.profileStepInterval = self.times.profileStepInterval
        opts.currentStateStepInterval = self.times.currentStateStepInterval
        opts.terminateStepInterval = self.times.terminateStepInterval
        opts.globalTimestep = self.times.globalTimestep
        opts.diffusionTimestepMultiplier = self.times.diffusionTimestepMultiplier

        # Debugging options
        CxxSetDebugParameters(self.debug.adaptation,
                              self.debug.regridding,
                              self.debug.timesteps,
                              self.debug.flameRadiusControl,
                              self.debug.veryVerbose)

        opts.debugSourcePoint = self.debug.sourcePoint
        opts.debugSourceTime = self.debug.sourceTime
        opts.debugStartTime = self.debug.startTime
        opts.debugStopTime = self.debug.stopTime

        # CVODE options
        CVODE = self.cvodeTolerances
        opts.integratorRelTol = CVODE.relativeTolerance
        opts.integratorMomentumAbsTol = CVODE.momentumAbsTol
        opts.integratorEnergyAbsTol = CVODE.energyAbsTol
        opts.integratorSpeciesAbsTol = CVODE.speciesAbsTol
        opts.integratorMinTimestep = CVODE.minimumTimestep

        # QSS integrator options
        QSS = self.qssTolerances
        opts.qss_epsmax = QSS.epsmax
        opts.qss_epsmin = QSS.epsmin
        opts.qss_dtmin = QSS.dtmin
        opts.qss_dtmax = QSS.dtmax
        opts.qss_iterationCount = QSS.iterationCount
        opts.qss_abstol = QSS.abstol
        opts.qss_minval = QSS.minval
        opts.qss_stabilityCheck = QSS.stabilityCheck

        # Output files
        opts.outputProfiles = self.outputFiles.saveProfiles
        opts.outputDebugIntegratorStages = self.outputFiles.debugIntegratorStages

        # Termination Conditions
        TC = self.terminationCondition
        opts.tEnd = TC.tEnd
        opts.tEndMin = TC.tMin
        opts.terminationMeasurement = stringify(TC.measurement if TC.measurement is not None else '')
        if TC.measurement:
            opts.terminationPeriod = TC.steadyPeriod
            opts.terminationTolerance = TC.tolerance
            opts.terminationAbsTol = TC.abstol
            opts.termination_dTdtTol = TC.dTdtTol
            opts.termination_momentsTol = TC.momentsTol


cdef np.ndarray[np.double_t, ndim=1] chebyshev1(double x, int N):
    """
    Return the values of the Chebyshev polynomials of the first kind
    evaluated at *x* up to order *N*.
    """
    cdef np.ndarray[np.double_t, ndim=1] T = np.empty(N+1)
    T[0] = 1
    if N == 1:
        return T

    T[1] = x
    if N == 2:
        return T

    cdef int i
    for i in range(1, N-1):
        T[i+1] = 2 * x * T[i] - T[i-1]

    return T


cdef np.ndarray[np.double_t, ndim=1] getChebyshevCoeffs(f, int N, double t0, double t1):
    cdef np.ndarray[np.double_t, ndim=1] a = np.zeros(N)
    cdef np.ndarray[np.double_t, ndim=1] T = np.zeros(N)
    cdef int k, n
    cdef np.ndarray[np.double_t, ndim=1] X = np.array([np.cos(np.pi * (k + 0.5) / N) for k in range(N)])
    cdef double fx, t
    for k in range(N):
        T = chebyshev1(X[k], N)
        t = 0.5 * (X[k] * (t1-t0) + t1 + t0)
        fx = f(t)
        a[0] += 1.0 / N * T[0] * fx
        for n in range(1, N):
            a[n] += 2.0 / N * T[n] * fx
    return a


cdef class FlameSolver:
    def __cinit__(self, *args, **kwargs):
        self.solver = new CxxFlameSolver()

    def __init__(self, options):
        if not isinstance(options, ConfigOptions):
            options = options.evaluate()
        self.options = <ConfigOptions?>options # keep to prevent garbage collection
        self.solver.setOptions(deref(self.options.opts))

        self.strainFunction = options.strainParameters.function
        self.rateMultiplierFunction = options.chemistry.rateMultiplierFunction
        self._updateStrainFunction()

        self.heatLossFunction = options.externalHeatFlux.heatLoss

        self.timeseriesWriter = options.outputFiles.timeSeriesWriter
        self.stateWriter = options.outputFiles.stateWriter

        if self.options.general.interpFile is not None:
            self._setup_interp_data(self.options.initialCondition.interpData)

    def __dealloc__(self):
        del self.solver

    def initialize(self):
        self.solver.initialize()

    def finalize(self):
        self.solver.finalize()

    def step(self):
        if self.strainFunction is not None:
            self._updateStrainFunction()
        if self.rateMultiplierFunction is not None:
            self._updateRateMultiplierFunction()

        try:
            with nogil:
                done = self.solver.step()
        except Exception as e:
            writelog(str(e))
            raise

        return done

    def _setup_interp_data(self, data):
        cdef np.ndarray[np.double_t, ndim=1] r = np.ascontiguousarray(data.r)
        cdef np.ndarray[np.double_t, ndim=1] z = np.ascontiguousarray(data.z)

        cdef np.ndarray[np.double_t, ndim=2] T = np.ascontiguousarray(data.T)
        cdef np.ndarray[np.double_t, ndim=2] vr = np.ascontiguousarray(data.vr)
        cdef np.ndarray[np.double_t, ndim=2] vz = np.ascontiguousarray(data.vz)

        nz, nr = vz.shape[0], vz.shape[1]
        assert T.shape[0] == vr.shape[0] == vz.shape[0] == len(z)
        assert T.shape[1] == vr.shape[1] == vz.shape[1] == len(r)

        self.solver.vzInterp.setup(map_matrix(&vz[0,0], nr, nz, nr, 1),
                                   map_vector(&r[0], nr, 1),
                                   map_vector(&z[0], nz, 1))
        self.solver.vrInterp.setup(map_matrix(&vr[0,0], nr, nz, nr, 1),
                                   map_vector(&r[0], nr, 1),
                                   map_vector(&z[0], nz, 1))
        self.solver.TInterp.setup(map_matrix(&T[0,0], nr, nz, nr, 1),
                                  map_vector(&r[0], nr, 1),
                                  map_vector(&z[0], nz, 1))

    property heatLossFunction:
        def __set__(self, func):
            if func is not None:
                self._heatLossFunction = IntegratorCallback(func)
                self.solver.heatLossFunction = self._heatLossFunction.callback
            else:
                self._heatLossFunction = None
                self.solver.heatLossFunction = NULL
        def __get__(self):
            return self._heatLossFunction.func

    property timeseriesWriter:
        def __set__(self, writer):
            if writer is not None:
                self._timeseriesWriter = LoggerCallback(writer(self, self.options))
                self.solver.timeseriesWriter = self._timeseriesWriter.callback
            else:
                self._timeseriesWriter = None
                self.solver.timeseriesWriter = NULL
        def __get__(self):
            return self._timeseriesWriter.func

    property stateWriter:
        def __set__(self, writer):
            if writer is not None:
                self._stateWriter = LoggerCallback(writer(self, self.options))
                self.solver.stateWriter = self._stateWriter.callback
            else:
                self._stateWriter = None
                self.solver.stateWriter = NULL
        def __get__(self):
            return self._stateWriter.func

    def writeStateFile(self, filename):
        """
        Create an HDF5 output file named *filename* containing the spatial
        profiles of flame at the current time.
        """
        self.stateWriter(filename)

    def writeTimeseriesFile(self, filename):
        """
        Create an HDF5 output file named *filename* containing integral flame
        properties (flame speed, total heat release, flame position, etc.) as
        a function of time, from the start of the simulation up to the current
        time.
        """
        self.timeseriesWriter(filename)

    def _updateStrainFunction(self):
        if self.solver.strainfunc == NULL:
            return
        cdef int N = self.options.general.chebyshevOrder
        cdef np.ndarray[np.double_t, ndim=1] coeffs = np.empty(N+2)
        cdef double t0 = self.solver.tNow
        cdef double t1 = self.solver.tNow + self.solver.dt
        coeffs[0] = t0
        coeffs[1] = t1
        coeffs[2:] = getChebyshevCoeffs(self.strainFunction, N, t0, t1)
        self.solver.strainfunc.setCoefficients(N+2, &coeffs[0])

    def _updateRateMultiplierFunction(self):
        if self.solver.rateMultiplierFunction == NULL:
            return
        cdef int N = self.options.general.chebyshevOrder
        cdef np.ndarray[np.double_t, ndim=1] coeffs = np.empty(N+2)
        cdef double t0 = self.solver.tNow
        cdef double t1 = self.solver.tNow + self.solver.dt
        coeffs[0] = t0
        coeffs[1] = t1
        coeffs[2:] = getChebyshevCoeffs(self.rateMultiplierFunction, N, t0, t1)
        self.solver.rateMultiplierFunction.setCoefficients(N+2, &coeffs[0])

    property tNow:
        def __get__(self):
            return self.solver.tNow

    property dt:
        def __get__(self):
            return self.solver.dt

    property heatReleaseRate:
        def __get__(self):
            return self.solver.getHeatReleaseRate()

    property consumptionSpeed:
        def __get__(self):
            return self.solver.getConsumptionSpeed()

    property flamePosition:
        def __get__(self):
            return self.solver.getFlamePosition()

    property x:
        def __get__(self):
            return getArray_Vec(self.solver.grid.x)

    property gridAlpha:
        def __get__(self):
            return self.solver.grid.alpha

    property gridBeta:
        def __get__(self):
            return self.solver.grid.beta

    property T:
        def __get__(self):
            return getArray_VecMap(self.solver.T)

    property U:
        def __get__(self):
            return getArray_VecMap(self.solver.U)

    property Y:
        def __get__(self):
            return getArray_MatrixMap(self.solver.Y)

    property moments:
        def __get__(self):
            return getArray_MatrixMap(self.solver.moments)

    property numberDensity:
        """
        Particle number density [particles/m^3], derived from the N moment
        (``moments[0]``, transported internally as kmol of particles per kg
        of gas -- consistent with the kmol convention used for gas species)
        and the local gas density. This is purely a unit conversion of
        already-computed state; the underlying transport/source terms are
        unchanged.
        """
        def __get__(self):
            if self.nMoments == 0:
                raise AttributeError(
                    "numberDensity is unavailable: particles.nMoments == 0")
            return self.moments[0, :] * self.rho * cantera.avogadro

    property nSpec:
        def __get__(self):
            return self.solver.nSpec

    property nMoments:
        def __get__(self):
            return self.solver.nMoments

    property V:
        def __get__(self):
            return getArray_Vec(self.solver.convectionSystem.V)

    property a:
        def __get__(self):
            return self.solver.strainfunc.a(self.solver.tNow)

    property dadt:
        def __get__(self):
            return self.solver.strainfunc.dadt(self.solver.tNow)

    property terminationCondition:
        def __get__(self):
            return self.solver.terminationCondition

    property qDot:
        def __get__(self):
            return getArray_Vec(self.solver.qDot)

    property rho:
        def __get__(self):
            return getArray_Vec(self.solver.rho)

    property particleDiameter:
        """Mean single-particle diameter [m] implied by the N, mM, mO moments."""
        def __get__(self):
            return getArray_Vec(self.solver.particleDiameter)

    property particleOxygenRatio:
        """
        Particle-phase O/metal atomic ratio, from the mM (``moments[1]``) and
        mO (``moments[2]``) moments. 0 where no metal is present, or for
        pure-metal particles (``nMoments`` < 3).
        """
        def __get__(self):
            if self.nMoments < 2:
                raise AttributeError(
                    "particleOxygenRatio is unavailable: particles.nMoments < 2")
            mM = self.moments[1, :]
            if self.nMoments < 3:
                return np.zeros_like(mM)
            nM = mM / self.options.opts.metalWeight
            nO = self.moments[2, :] / self.options.opts.oxygenWeight
            return np.where(mM > 0, nO / np.where(mM > 0, nM, 1.0), 0.0)

    property particleDiffusivity:
        """Particle Brownian diffusivity [m^2/s] (Stokes-Einstein-Cunningham)."""
        def __get__(self):
            return getArray_Vec(self.solver.particleDiffusivity)

    property splitConstConv:
        def __get__(self):
            return getArray_Matrix(self.solver.splitConstConv)

    property splitConstDiff:
        def __get__(self):
            return getArray_Matrix(self.solver.splitConstDiff)

    property splitConstProd:
        def __get__(self):
            return getArray_Matrix(self.solver.splitConstProd)

    property dUdtDiff:
        def __get__(self):
            return getArray_MatrixRow(self.solver.ddtDiff, 0)

    property dUdtConv:
        def __get__(self):
            return getArray_MatrixRow(self.solver.ddtConv, 0)

    property dUdtProd:
        def __get__(self):
            return getArray_MatrixRow(self.solver.ddtProd, 0)

    property dTdtDiff:
        def __get__(self):
            return getArray_MatrixRow(self.solver.ddtDiff, 1)

    property dTdtConv:
        def __get__(self):
            return getArray_MatrixRow(self.solver.ddtConv, 1)

    property dTdtProd:
        def __get__(self):
            return getArray_MatrixRow(self.solver.ddtProd, 1)

    property dTdtCross:
        def __get__(self):
            return getArray_MatrixRow(self.solver.ddtCross, 1)

    property dYdtDiff:
        def __get__(self):
            return getArray_Matrix(self.solver.ddtDiff)[2:2+self.solver.nSpec]

    property dYdtConv:
        def __get__(self):
            return getArray_Matrix(self.solver.ddtConv)[2:2+self.solver.nSpec]

    property dYdtProd:
        def __get__(self):
            return getArray_Matrix(self.solver.ddtProd)[2:2+self.solver.nSpec]

    property dYdtCross:
        def __get__(self):
            return getArray_Matrix(self.solver.ddtCross)[2:2+self.solver.nSpec]

    property dWdt:
        def __get__(self):
            return getArray_Vec(self.solver.convectionSystem.dWdt)

    property drhodt:
        def __get__(self):
            return getArray_Vec(self.solver.drhodt)

    property sumcpj:
        def __get__(self):
            return getArray_Vec(self.solver.sumcpj)

    property Tleft:
        def __get__(self):
            return self.solver.Tleft

    property Yleft:
        def __get__(self):
            return getArray_Vec(self.solver.Yleft)

    property dWdx:
        def __get__(self):
            return getArray_Vec(self.solver.convectionSystem.utwSystem.dWdx)

    property dTdx:
        def __get__(self):
            return getArray_Vec(self.solver.convectionSystem.utwSystem.dTdx)

    property wdot:
        def __get__(self):
            return getArray_Matrix(self.solver.wDot)

    property rhoD:
        def __get__(self):
            return getArray_Matrix(self.solver.rhoD)

    property cp:
        def __get__(self):
            return getArray_Vec(self.solver.cp)

    property mu:
        def __get__(self):
            return getArray_Vec(self.solver.mu)

    property k:
        def __get__(self):
            return getArray_Vec(self.solver.k)

    property W:
        def __get__(self):
            return getArray_Vec(self.solver.W)

    property Wmx:
        def __get__(self):
            return getArray_Vec(self.solver.Wmx)

    property cfp:
        def __get__(self):
            return getArray_Vec(self.solver.grid.cfp)

    property cf:
        def __get__(self):
            return getArray_Vec(self.solver.grid.cf)

    property cfm:
        def __get__(self):
            return getArray_Vec(self.solver.grid.cfm)

    property hh:
        def __get__(self):
            return getArray_Vec(self.solver.grid.hh)

    property rphalf:
        def __get__(self):
            return getArray_Vec(self.solver.grid.rphalf)

    property jFick:
        def __get__(self):
            return getArray_Matrix(self.solver.jFick)

    property jSoret:
        def __get__(self):
            return getArray_Matrix(self.solver.jSoret)

    property jCorr:
        def __get__(self):
            return getArray_Vec(self.solver.jCorr)
