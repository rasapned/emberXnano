#!/usr/bin/python
"""
As single_Igor_impingingJet_fixedT_particles.py, plus the particle enthalpy
moment (Particles energy=True), so that the particles have their own temperature: they
take up the enthalpy of the nucleating/condensing gas molecules (latent heat,
consistent with the mechanism and the condensed-phase NASA thermo), the
surface-reaction heats (deltaH below), and lose heat by conduction to the gas
(thermal accommodation) and Rayleigh-limit radiation. The surface-reaction
Arrhenius rates use the particle temperature. The gas temperature is the
measured profile; the gas energy is not coupled back.

Original description of the particle model:
Igor's burner (lean, 500 ppm FEC5O5, 3000 Pa) as an impinging jet with the
measured temperature profile, as single_Igor_impingingJet_fixedT.py, plus
iron particles with the same phenomena as burnerFla_particles_oxid.py:
nucleation from Fe-cluster collisions, coagulation, particle-phase oxygen
(nMoments=3), the OpenFOAM gas-particle surface reactions (oxidation by
O2/O/H2O, reduction by H2/H, etching by OH) and heterogeneous condensation of
Fe, the Fe clusters and FeO.

Cantera's ImpingingJet (gas phase only, same fixed T profile) gives the
initial gas profiles and the equivalent strain rate; the particle moments
start from zero and build up from nucleation. Particles do not cross the
plate (zero flux, no deposition or thermophoresis).

Termination: with T prescribed, dT/dt = 0 everywhere, so the 'dTdt' check
would stop at tMin, and the heat release ('Q') barely sees the particles.
The 'moments' measurement checks the particle fields, which settle last.
"""

from ember import *
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct

output = 'run/single_Igor_impingingJet_fixedT_particlesEnergy_Hrates'
mech = 'Iron_elte_Syngas-newTransp.yaml'

# Lean case, 500 ppm FEC5O5, volume flow rates in sccm
X_FEC5O5 = 0.0005
VS_total = (400 + 400 + 600) / (1 - X_FEC5O5)
X_H2 = 400 / VS_total
X_O2 = 400 / VS_total
X_AR = 1.0 - (X_H2 + X_O2 + X_FEC5O5)
react = f"H2:{X_H2:.8e}, AR:{X_AR:.8e}, FEC5O5:{X_FEC5O5:.8e}, O2:{X_O2:.8e}"

pressure = 3000.0  # [Pa]

# Measured temperature profile: nozzle at x = 0, plate at x = L
xProfile, TProfile = np.genfromtxt('T_of_x-L500.csv', delimiter=',', comments='#').T
L = xProfile[-1]
Tburner = TProfile[0]
Tplate = TProfile[-1]

# Inlet mass flux: standard (1 bar, 273.15 K) volume flow through the nozzle
gas = ct.Solution(mech, transport_model='mixture-averaged')
gas.TPX = 273.15, 1e5, react
mdot = gas.density * VS_total / 60e6 / (0.25 * np.pi * 0.036**2)
gas.TPX = Tburner, pressure, react
rhou = gas.density

# *** Cantera (gas phase only), fixed temperature profile ***
jet = ct.ImpingingJet(gas=gas, width=L)
jet.inlet.mdot = mdot
jet.surface.T = Tplate
jet.set_refine_criteria(ratio=3, slope=0.05, curve=0.1, prune=0.02)
jet.flame.set_fixed_temp_profile(xProfile / L, TProfile)
jet.energy_enabled = False
jet.solve(loglevel=0, refine_grid=True)
aCantera = 2.0 * np.sqrt(-jet.L[0] / rhou)  # beta = 2 (disc)
print('Cantera: {} points, equivalent strain rate a = {:.3f} 1/s'.format(
    len(jet.grid), aCantera))

# ember takes the inlet stream composition from the first point of the given
# profile, so that point must hold the supplied mixture, not Cantera's inlet
# node (depleted in H2, which diffuses upstream against the inlet flux)
Y0 = jet.Y.copy()
gas.TPX = Tburner, pressure, react
Y0[:, 0] = gas.Y

# *** Particle model, as burnerFla_particles_oxid.py ***
# Collisional nucleation: the mechanism tracks gas-phase Fe clusters up to
# FE7; any collision whose combined size reaches FE8 or larger forms a
# particle nucleus.
maxClusterSize = 7
clusterSpecies = ['FE'] + [f'FE{n}' for n in range(2, maxClusterSize + 1)]
collisionSpeciesA, collisionSpeciesB = [], []
for i in range(1, maxClusterSize + 1):
    for j in range(i, maxClusterSize + 1):
        if i + j >= 8:
            collisionSpeciesA.append(clusterSpecies[i-1])
            collisionSpeciesB.append(clusterSpecies[j-1])

# Surface reactions translated from the OpenFOAM setup (see
# burnerFla_particles_oxid.py): A_ember [m/s] = |A_OF| / N_A[1/kmol] * 1000,
# Ta = Ea / R; the sign of A_OF (+ adds to, - removes from the particle) is
# checked against the particle change ember derives from the compositions.
R_JmolK = 8.314
openfoamReactions = [
    # reactant, product, A_OF,     Ea [J/mol]
    ('O2',  'O',    2.324e26,  8.314e3),   # oxidation
    ('O',   None,   2.629e26,  0.000),     # oxidation (OpenFOAM product: AR)
    ('H2O', 'H2',   2.714e26,  1.081e4),   # oxidation
    ('H2',  'H2O', -8.040e26,  1.663e4),   # reduction
    ('H',   'OH',  -1.079e27,  4.157e3),   # reduction
    ('OH',  'FEOH', -2.568e27, 1.829e4),   # etching (removes Fe)
]
# Reaction enthalpies deltaH [J/mol] of reactant(g) + particle -> product(g)
# + particle': only -deltaH heats the particle (its composition change itself
# is booked at the particle temperature). Currently 0 for all reactions (no
# reaction heat; the particles are heated only by the latent heat of
# nucleation/condensation): set the OpenFOAM values in deltaH below. For
# reference, the standard reaction enthalpies at 298 K for Fe/FeO particles
# (mechanism gas thermo + Cantera's condensed-phase data) are printed.
condensed = {s.name: s for s in ct.Species.list_from_file('nasa_condensed.yaml')}
hFe = condensed['Fe(a)'].thermo.h(298.15)       # [J/kmol]
hFeO = condensed['FeO(s)'].thermo.h(298.15)

def reactionEnthalpy298(reactant, product, dM, dO):
    """Standard reaction enthalpy [J/mol]; particle change (dM, dO) taken as Fe -> FeO."""
    h = -gas.species(reactant).thermo.h(298.15)
    if product is not None:
        h += gas.species(product).thermo.h(298.15)
    h += hFe * dM + (hFeO - hFe) * dO
    return h / 1000.0

# (reactant, product) -> deltaH [J/mol]; reactions not listed: 0. OpenFOAM
# values are the heat released into the particle, so deltaH = -value: every
# surface reaction heats the particle.
deltaH = {
    ('O2',  'O'):    -80000.0,
    ('O',   None):   -578000.0,
    ('H2O', 'H2'):    -87000.0,
    ('H2',  'H2O'):   -52000.0,
    ('H',   'OH'):    -53000.0,
    ('OH',  'FEOH'): -157000.0,
}

surfaceReactions = []
for reactant, product, A_OF, Ea in openfoamReactions:
    reaction = SurfaceReaction(reactant=reactant, product=product,
                               A=abs(A_OF) / ct.avogadro*1000, Ta=Ea / R_JmolK)
    dM, dO = reaction.particleChange(gas, 'Fe')
    dH = deltaH.get((reactant, product), 0.0)
    reaction.deltaH.value = dH
    print('Surface reaction %-4s -> %-5s deltaH = %8.1f kJ/mol (298 K standard: %8.1f)' % (
        reactant, product, dH / 1e3, reactionEnthalpy298(reactant, product, dM, dO) / 1e3))
    if (dM + dO > 0) != (A_OF > 0):
        raise ValueError("Surface reaction %s -> %s: sign of A (%g) does not"
                         " match the particle change (dFe=%d, dO=%d)" %
                         (reactant, product, A_OF, dM, dO))
    surfaceReactions.append(reaction)

conf = Config(
    Paths(outputDir=output),
    Chemistry(mechanismFile=mech,
              transportModel='Mix'),
    General(flameGeometry='disc',
            nThreads=16,
            chemistryIntegrator='cvode'),
    InitialCondition(reactants=react,
                     Tu=Tburner,
                     pressure=pressure,
                     xRight=L,
                     haveProfiles=True,
                     x=jet.grid,
                     T=jet.T,
                     U=jet.spread_rate,
                     V=jet.density * jet.velocity,
                     Y=Y0),
    ImpingingJet(wallTemperature=Tplate,
                 massFlux=mdot,
                 temperatureProfile=(xProfile, TProfile)),
    StrainParameters(initial=aCantera,
                     final=aCantera),
    CvodeTolerances(relativeTolerance=1e-7,
                    speciesAbsTol=1e-12,
                    energyAbsTol=1e-8,
                    minimumTimestep=1e-16),
    # Stop when every particle moment changes by less than momentsTol (RMS
    # over the grid, relative to its peak) per second. tMin lets nucleation
    # build the particle fields up from zero first.
    TerminationCondition(measurement='moments',
                         momentsTol=0.1,
                         tMin=0.02,
                         tEnd=0.5),
    Particles(
        nMoments=3,                 # N, Fe mass, O mass
        energy=True,                # + particle enthalpy -> particle temperature
        momentBCLeft=0.0,
        metal='Fe',
        # (O/metal atomic ratio, density [kg/m3]): Fe, FeO, Fe3O4, Fe2O3
        # condensed species from Cantera's nasa_condensed.yaml (lowest Gibbs
        # energy at each T: polymorphs and melting included)
        phases=[(0.0, 7874.0, ['Fe(a)', 'Fe(c)', 'Fe(d)', 'Fe(L)']),
                (1.0, 5745.0, ['FeO(s)', 'FeO(L)']),
                (4/3, 5170.0, ['Fe3O4(s)']),
                (1.5, 5240.0, ['Fe2O3(s)'])],
        thermalAccommodation=0.3,   # typical for Fe / iron oxide nanoparticles
        radiationAbsorption=0.3,    # E(m), Rayleigh-limit radiation
        radiationTsurr=300.0,       # surroundings [K]
        nucleation=NucleationChannel(
            collisionSpeciesA=collisionSpeciesA,
            collisionSpeciesB=collisionSpeciesB,
        ),
        coagulation=True,
        # Particles shrinking below this diameter [m] disintegrate (N removed
        # with the lost mass); also the floor for the particle-diffusivity size.
        minParticleDiameter=0.5e-9,
        surfaceReactions=surfaceReactions,
        condensationSpecies=clusterSpecies + ['FEO'],
    ),
)

if __name__ == '__main__':
    conf.run()

    struct = utils.load(output + '/profNow.h5')
    print('Stopped at t = {:.4f} s'.format(struct.t))
    print('Inlet mass flux: Cantera {:.5f}, ember {:.5f} kg/m^2/s'.format(
        mdot, struct.V[0]))
    print('Equivalent strain rate: Cantera (no particles) {:.3f}, ember {:.3f} 1/s'.format(
        aCantera, struct.a))

    def gasFe(Y):
        kFe = gas.element_index('Fe')
        nFe = np.array([gas.n_atoms(k, kFe) for k in range(gas.n_species)])
        return (nFe * gas.atomic_weights[kFe] / gas.molecular_weights) @ Y

    # --- Gas phase vs Cantera (which has no particle sink) ---
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    ax = ax.ravel()
    ax[0].plot(jet.grid, jet.velocity, 'k--', label='Cantera (gas only)')
    ax[0].plot(struct.x, struct.V / struct.rho, label='ember')
    ax[0].set_ylabel('Axial velocity [m/s]')
    ax[1].plot(jet.grid, jet.spread_rate, 'k--')
    ax[1].plot(struct.x, struct.U)
    ax[1].set_ylabel('Radial velocity / r [1/s]')
    for i, species in enumerate(['H2O', 'FE', 'OH']):
        k = gas.species_index(species)
        ax[2+i].plot(jet.grid, jet.Y[k], 'k--')
        ax[2+i].plot(struct.x, struct.Y[k])
        ax[2+i].set_ylabel('Y ' + species)
    ax[5].plot(jet.grid, gasFe(jet.Y), 'k--')
    ax[5].plot(struct.x, gasFe(struct.Y))
    ax[5].set_ylabel('Fe in gas phase (element mass fraction)')
    for a in ax:
        a.set_xlabel('Position [m]')
    ax[0].legend()
    plt.suptitle('Impinging jet, fixed T profile, with particles: gas phase')
    plt.tight_layout()
    plt.savefig(output + '/gas_comparison.png')
    plt.close()

    # --- Particle fields ---
    # Size and O/Fe ratio are ratios of moments; where there are practically
    # no particles they are 0/0 noise, so they are only shown where the
    # particle mass is above 1e-3 of its peak
    particleMass = struct.moments[1:3].sum(axis=0)
    present = particleMass > 1e-3 * particleMass.max()
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    ax = ax.ravel()
    ax[0].plot(struct.x, struct.T, 'r-')
    ax[0].set_ylabel('Temperature [K]')
    ax[1].plot(struct.x, struct.numberDensity)
    ax[1].set_ylabel('Number density [1/m³]')
    ax[2].plot(struct.x, np.where(present, struct.particleDiameter * 1e9, np.nan))
    ax[2].set_ylabel('Particle diameter [nm]')
    ax[3].plot(struct.x, particleMass)
    ax[3].set_ylabel('Particle mass fraction [-]')
    ax[4].plot(struct.x, np.where(present, struct.particleOxygenRatio, np.nan))
    ax[4].axhline(1.5, color='k', ls='--', lw=0.8, label='Fe2O3')
    ax[4].axhline(4/3, color='0.5', ls=':', lw=0.8, label='Fe3O4')
    ax[4].axhline(1.0, color='0.7', ls=':', lw=0.8, label='FeO')
    ax[4].set_ylabel('Particle O/Fe atomic ratio')
    ax[4].legend(fontsize=8)
    ax[5].plot(struct.x, gasFe(struct.Y), label='gas phase')
    ax[5].plot(struct.x, struct.moments[1], label='particles')
    ax[5].plot(struct.x, gasFe(struct.Y) + struct.moments[1], 'k:', label='total')
    ax[5].set_ylabel('Fe mass per mass of gas [-]')
    ax[5].legend(fontsize=8)
    for a in ax:
        a.set_xlabel('Position [m]')
    plt.suptitle('Impinging jet, fixed T profile: particles')
    plt.tight_layout()
    plt.savefig(output + '/particles.png')
    plt.close()

    # --- Overview, as burnerFla_particles_oxid.py: T, particle fields and
    # mole fractions of the Fe-O-H and radical species. Cantera has no
    # particle phase, so N, Y_p and d_p are ember only.
    Xember = struct.Y / gas.molecular_weights[:, None]
    Xember /= Xember.sum(axis=0)
    fields = [
        ('T', 'T [K]', struct.T, jet.T),
        ('N', 'N [1/m³]', struct.numberDensity, None),
        ('Y_p', 'Particle mass fraction [-]', particleMass, None),
        ('d_p', 'd_p [nm]', np.where(present, struct.particleDiameter * 1e9, np.nan), None),
    ]
    for name in ('FE', 'FEO', 'FEO2', 'FEOH', 'FEO2H2', 'O', 'H', 'OH'):
        k = gas.species_index(name)
        fields.append((f'x_{name}', 'Mole fraction', Xember[k], jet.X[k]))

    fig, axes = plt.subplots(3, 4, figsize=(18, 10), sharex=True)
    for ax, (title, ylabel, yE, yC) in zip(axes.flat, fields):
        ax.plot(struct.x, yE, lw=2, label='Ember')
        if yC is not None:
            ax.plot(jet.grid, yC, '--', lw=2, label='Cantera (gas only)')
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.legend(loc='best', fontsize=8)
    for ax in axes[-1, :]:
        ax.set_xlabel('Position [m]')
    fig.suptitle('Impinging jet, fixed T profile: ember (with particles) vs Cantera ImpingingJet')
    fig.tight_layout()
    fig.savefig(output + '/Particles_overview.png')
    plt.close(fig)

    # --- Particle temperature and heat losses ---
    # Losses per unit volume [W/m^3], recomputed from the profiles with the
    # same formulas as SourceSystem::computeParticleHeatTransfer()
    P = conf.particles
    kB, NA = ct.boltzmann, ct.avogadro
    nP = struct.numberDensity
    qCond = np.zeros_like(struct.x)
    qRad = np.zeros_like(struct.x)
    for j in range(len(struct.x)):
        if not present[j]:
            continue
        gas.TPY = struct.T[j], pressure, struct.Y[:, j]
        Tg, Tp, dp = struct.T[j], struct.particleTemperature[j], struct.particleDiameter[j]
        cBar = np.sqrt(8 * kB * Tg / (np.pi * gas.mean_molecular_weight / NA))
        gamma = gas.cp_mass / gas.cv_mass
        qFM = (P.thermalAccommodation.value * np.pi * dp**2 * pressure / 8 * cBar *
               (gamma + 1) / (gamma - 1) * (Tp / Tg - 1))
        qC = 2 * np.pi * dp * gas.thermal_conductivity * (Tp - Tg)
        qCond[j] = nP[j] * qFM * qC / (qFM + qC)
        qRad[j] = nP[j] * (8 * 24 * 1.0369277551433699 * np.pi**3 * dp**3 *
                           P.radiationAbsorption.value * kB**5 /
                           (ct.planck**4 * ct.light_speed**3) *
                           (Tp**5 - P.radiationTsurr.value**5))
    dT = np.where(present, struct.particleTemperature - struct.T, np.nan)
    jMax = np.nanargmax(dT)
    print('Particle overheating: max Tp - Tg = {:.1f} K at x = {:.2f} mm '
          '(Tg = {:.0f} K, d_p = {:.2f} nm)'.format(
              dT[jMax], 1e3 * struct.x[jMax], struct.T[jMax],
              1e9 * struct.particleDiameter[jMax]))

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    ax[0].plot(struct.x, struct.T, 'r-', label='gas (prescribed)')
    ax[0].plot(struct.x, np.where(present, struct.particleTemperature, np.nan), 'k-',
               label='particles')
    ax[0].set_ylabel('Temperature [K]')
    ax[0].legend()
    ax[1].plot(struct.x, dT)
    ax[1].set_ylabel('T_p - T_gas [K]')
    ax[2].plot(struct.x, qCond, label='conduction to gas')
    ax[2].plot(struct.x, qRad, label='radiation')
    ax[2].set_ylabel('Particle heat loss [W/m³]')
    ax[2].legend()
    for a in ax:
        a.set_xlabel('Position [m]')
        a.set_xlim(0, struct.x[present].max() * 1.2)
    plt.suptitle('Impinging jet, fixed gas T: particle temperature')
    plt.tight_layout()
    plt.savefig(output + '/particle_temperature.png')
    plt.close()
