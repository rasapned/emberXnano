#!/usr/bin/python
"""
burnerFla_particles_cond.py with reversible condensation: Particles(evaporation=True)
(Kelvin evaporation of Fe, Fe2..Fe7 from the Fe phase and FeO from the FeO
phase; see doc/EvaporationModel.md). Everything else -- mechanism, OpenFOAM
grid, Cantera initial profile, tolerances, termination -- is taken unchanged
from the reference script by importing its Config. Particle energy is off, so
evaporation runs at T_p = T_gas.

Compare with run/burnerFla/burnerFla_particles_cond: run/burnerFla/OFcomparison/plot_cond_evap_ember_OF.py
"""
import burnerFla_particles_cond as case

conf = case.conf
conf.paths.outputDir.value = 'run/burnerFla/burnerFla_particles_cond_evap'
conf.particles.evaporation.value = True
# Same (ratio, density) as the reference; evaporation also needs each phase's
# condensed species (saturation pressures from nasa_condensed.yaml)
conf.particles.phases.value = [(0.0, 7874.0, ['Fe(a)', 'Fe(c)', 'Fe(d)', 'Fe(L)']),
                               (1.0, 5745.0, ['FeO(s)', 'FeO(L)']),
                               (4/3, 5170.0, ['Fe3O4(s)']),
                               (1.5, 5240.0, ['Fe2O3(s)'])]
conf.general.nThreads.value = 8

if __name__ == '__main__':
    conf.run()
