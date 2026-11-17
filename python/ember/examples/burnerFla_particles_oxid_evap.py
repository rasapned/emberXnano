#!/usr/bin/python
"""
burnerFla_particles_oxid.py with reversible condensation: Particles(evaporation=True)
(see doc/EvaporationModel.md). Everything else -- mechanism, OpenFOAM grid,
Cantera initial profile, surface reactions, tolerances, termination -- is
taken unchanged from the reference script by importing its Config. Particle
energy is off, so evaporation runs at T_p = T_gas.

Compare with run/burnerFla/burnerFla_particles_oxid.
"""
import burnerFla_particles_oxid as case

conf = case.conf
conf.paths.outputDir.value = 'run/burnerFla/burnerFla_particles_oxid_evap'
conf.particles.evaporation.value = True
# Same (ratio, density) as the reference; evaporation also needs each phase's
# condensed species (saturation pressures from nasa_condensed.yaml)
conf.particles.phases.value = [(0.0, 7874.0, ['Fe(a)', 'Fe(c)', 'Fe(d)', 'Fe(L)']),
                               (1.0, 5745.0, ['FeO(s)', 'FeO(L)']),
                               (4/3, 5170.0, ['Fe3O4(s)']),
                               (1.5, 5240.0, ['Fe2O3(s)'])]

if __name__ == '__main__':
    conf.run()
