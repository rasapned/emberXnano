#!/usr/bin/python
"""
Similar to as is done experimentally (https://doi.org/10.1016/j.combustflame.2009.06.011),
a single disc flame opposing a cold inert is established at a given strain rate. Infinite
burner separation distance assumed. The converged axial velocity profile is plotted.

NEW CASE BASED ON EXAMPLE_SINGLE, recreates a setup closer to igor's burner.
"""

from ember import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cantera as ct

output = 'run/single_Igor_nucl'

# Set strain rate (1/s) for the counterflow flame
a = 150.0
# Set the mole fraction of FEC5O5
X_FEC5O5 = 0.0005

# Volume flow rates in sccm for the LEAN 500 ppm FEC5O5 case
VS_total = 400 + 400 + 600
VS_total /= (1 - X_FEC5O5)
X_H2 = 400 / VS_total
X_O2 = 400 / VS_total
X_AR = 1.0 - (X_H2 + X_O2 + X_FEC5O5)

# Write the reactants
react = f"H2:{X_H2:.8e}, AR:{X_AR:.8e}, FEC5O5:{X_FEC5O5:.8e}, O2:{X_O2:.8e}"
print("Reactants lean:", react)

# Alternative: set fuel and oxidiser separately (it doesn't work xd) + equivalence ratio
fuel = f"H2:{X_H2:.8e}, FEC5O5:{X_FEC5O5:.8e}"
oxidizer = f"O2:{X_O2:.8e}, AR:{X_AR:.8e}"
equivalenceRatio = 0.5

# Collisional nucleation: the mechanism only tracks gas-phase Fe clusters up
# to FE7 (FE, FE2, ..., FE7); any collision whose combined size reaches FE8
# or larger is treated as leaving the gas phase to form a particle nucleus.
maxClusterSize = 7
clusterSpecies = ['FE'] + [f'FE{n}' for n in range(2, maxClusterSize + 1)]
collisionSpeciesA, collisionSpeciesB = [], []
for i in range(1, maxClusterSize + 1):
    for j in range(i, maxClusterSize + 1):
        if i + j >= 8:
            collisionSpeciesA.append(clusterSpecies[i-1])
            collisionSpeciesB.append(clusterSpecies[j-1])
print("Nucleation channels:",
      list(zip(collisionSpeciesA, collisionSpeciesB)))

# Configure the ember flame
conf = Config(
    Paths(outputDir=output),
    Chemistry(mechanismFile='Iron_elte_Syngas-newTransp.yaml',
              transportModel='Mix'),
    General(twinFlame=False,
            flameGeometry= 'planar',
            nThreads=1,
            chemistryIntegrator='cvode',
            splittingMethod='strang',
            continuityBC='stagnationPoint'
    ),
    InitialCondition(reactants=react,
                     pressure=3000.0,
                     counterflow = 'N2:1.0',
                     #Tu=300.0,
                     Tcounterflow = 300.0,
                     xLeft=-0.0,
                     xRight=0.065,
                     centerWidth=0.01,
                     slopeWidth=0.005,
                     ),
    # Grid(
    #     vtol=0.14,           # coarser interior refinement (default 0.12)
    #     dvtol=0.25,          # coarser gradient refinement (default 0.2)
    #     gridMax=3e-4,       # larger max spacing → fewer points (default 2e-4)
    #     boundaryTol=1e-4,   # less aggressive domain extension (default 5e-5)
    #     boundaryTolRm=5e-5, # remove unneeded boundary points (default 1e-5)
    #     addPointCount=2,    # add 1 point per extension event, not 3 (default 3)
    # ),
    StrainParameters(initial=a/5,
                     final=a,
                     dt=0.005,),
    CvodeTolerances(
        relativeTolerance=1e-7,
        speciesAbsTol=1e-12,
        energyAbsTol=1e-8,
        minimumTimestep=1e-16,
    ),
    TerminationCondition(
        #steadyPeriod=0.002,   # check over longer window
        tolerance=5e-4,        # relative tolerance (default 1e-4)
        tMin=0.008,             # don't check until flame is established
        tEnd=0.015
    ),
    Particles(
        nMoments=2,
        momentBCLeft=0.0,
        particleDensity=7874.0,
        nucleation=NucleationChannel(
            collisionSpeciesA=collisionSpeciesA,
            collisionSpeciesB=collisionSpeciesB,
        ),
        coagulation=True,
    )
)

# Run the simulation and plot the results
if __name__ == '__main__':
    conf.run()

    struct = utils.load(output + '/profNow.h5')

    print(struct.keys())    
    
    # Plot the number particles and volume of the particles
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(struct.x, struct.numberDensity, 'b-', label='Number of particles')
    ax2.plot(struct.x, struct.moments[1,:], 'r-', label='Particle diameter')
    ax1.set_xlabel('Position [m]')
    ax1.set_ylabel('Particle number density / m⁻³', color='b')
    ax2.set_ylabel('Particle diameter / m', color='r')
    ax1.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    #plt.title(f'Ember planar, p=3000 Pa, a={a} s⁻¹')
    plt.tight_layout()
    plt.savefig(output + '/Particles_nucl.png')
    plt.close()

    # Plot the temperature of the flame
    fig, ax1 = plt.subplots()
    ax1.plot(struct.x, struct.T, 'r-', label='Temperature')
    ax1.set_xlabel('Position [m]')
    ax1.set_ylabel('Temperature / K', color='r')
    ax1.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    #plt.title(f'Ember planar, p=3000 Pa, a={a} s⁻¹')
    plt.tight_layout()
    plt.savefig(output + '/Temperature.png')
    plt.close()


    # # Plot the mass flux profile
    # fig2, ax3 = plt.subplots()
    # ax3.plot(struct.x, struct.V, 'g-', label='Mass flux')
    # ax3.set_xlabel('Position [m]')
    # ax3.set_ylabel('Mass flux [kg/m²·s]', color='g')
    # plt.title(f'Ember disc, p=3000 Pa, a={a} s⁻¹')
    # plt.tight_layout()
    # plt.savefig(output + '/massflow.png')
    # plt.close()
    
    # # Plot the progress variable and mixture fraction
    # fig4, (ax6, ax8) = plt.subplots(1,2, figsize=(12, 5))
    # ax7 = ax6.twinx()
    # # Progress variable   
    # gas = ct.Solution('Iron_elte_Syngas-newTransp.yaml')
    # i_H2O = gas.species_index('H2O')
    # i_H2  = gas.species_index('H2')
    # i_N2  = gas.species_index('N2')
    # # Idea 1: add water and N2 and normalise
    # #progress = struct.Y[i_H2O, :] + struct.Y[i_N2, :] 
    # #progress /= progress.max()  # Normalize progress variable to [0, 1]
    # # Idea 2: make a separate mixture fraction and progress variable
    # Y_H2_u = struct.Y[i_H2, 0]  # Unburned H2 mass fraction
    # Y_H2O_max = struct.Y[i_H2O, :].max()  # Maximum H2O mass fraction
    # Z = 1.0 - struct.Y[i_N2, :]  # Mixture fraction based on N2
    # progress_H2 = 1.0 - struct.Y[i_H2, :] / (Z * Y_H2_u + 1e-12)  # Progress-mixture variable
    # progress_H2 = np.clip(progress_H2, 0.0, 1.0)  # Ensure progress variable is within [0, 1]
    # progress_H2O = struct.Y[i_H2O, :] / (Z * Y_H2O_max + 1e-12)  # Progress-mixture variable
    # progress_H2O = np.clip(progress_H2O, 0.0, 1.0)  # Ensure progress variable is within [0, 1]
    # # Clip the progress var. until it reaches local maximum
    # i = 0
    # while progress_H2O[i] <= progress_H2O[i+1]:
    #     i += 1
    # j = 0
    # while progress_H2[j] <= progress_H2[j+1]:
    #     j += 1
    # ax6.plot(struct.x[:i+1], progress_H2O[:i+1], 'y-', label='Ember c_H2O')
    # #ax6.plot(struct.x[:j+1], progress_H2[:j+1], 'c-', label='Ember c_H2')
    # ax7.plot(struct.x, Z, 'r-', label='Mixture fraction')
    # ax8.plot(progress_H2O[:i+1], struct.T[:i+1], 'y-', label='Ember c_H2O')
    # #ax8.plot(progress_H2[:j+1], struct.T[:j+1], 'c-', label='Ember c_H2')
    # ax6.set_xlabel('Position [m]')
    # ax6.set_ylabel('Progress variable', color='b')
    # ax7.set_ylabel('Mixture fraction', color='r')
    # ax8.set_xlabel('Progress variable')
    # ax8.set_ylabel('Temperature [K]')
    # # Read cantera data and plot for comparison
    # path_cant = output + '/cant_progress_s150.txt'
    # print("Reading cantera data from:", path_cant)
    # cant_x, cant_T, cant_c_H2O, cant_c_H2 = np.genfromtxt(path_cant, unpack=True, skip_header=1)
    # i = 0
    # j = 0
    # while cant_c_H2O[i] <= cant_c_H2O[i+1]:
    #     i += 1
    # while cant_c_H2[j] <= cant_c_H2[j+1]:
    #     j += 1
    # ax6.plot(cant_x[:i+1], cant_c_H2O[:i+1], 'y--', label='Cantera c_H2O')
    # #ax6.plot(cant_x[:j+1], cant_c_H2[:j+1], 'c--', label='Cantera c_H2')
    # ax8.plot(cant_c_H2O[:i+1], cant_T[:i+1], 'y--', label='Cantera c_H2O')
    # #ax8.plot(cant_c_H2[:j+1], cant_T[:j+1], 'c--', label='Cantera c_H2')
    # ax6.legend(loc='best')
    # ax7.legend(loc='best') 
    # ax8.legend(loc='best')
    # plt.title(f'Ember planar, p=3000 Pa, a={a} s⁻¹')
    # plt.tight_layout()
    # plt.savefig(output + '/ProgressV.png')
    # plt.close()
    # plt.close()

   
    # for t in range(17):
    #     fig3, ax4 = plt.subplots()
    #     ax4.set_xlabel('Position [m]')
    #     ax4.set_ylabel('Moment', color='b')
    #     ax4.tick_params(axis='y', labelcolor='b')
       
    #     currFile = output + '/prof{:06d}.h5'.format(t)
    #     struct_t = utils.load(currFile) 
    #     ax4.plot(struct_t.x, struct_t.moments[0], 'g-',  label=f'time = {t*0.001:.3f} s')

    #     ax4.legend(loc='best')
    #     plt.title(f'Test')
    #     plt.tight_layout()
    #     plt.savefig(output + f'/moment_test-{t:06d}.png')
    #     plt.close()
    
