#!/usr/bin/python
"""
single_Igor_impingingJet_fixedT_particles.py (no particle energy, surface
reactions on) with reversible condensation and, optionally, Kelvin-limited
nucleation. Everything else is taken unchanged from the reference script by
importing its Config.

    python single_Igor_impingingJet_fixedT_particles_evap.py              # no evaporation (current code)
    python single_Igor_impingingJet_fixedT_particles_evap.py evap         # Particles(evaporation=True)
    python single_Igor_impingingJet_fixedT_particles_evap.py evap kelvin  # + NucleationChannel(kelvinLimit=True)

T is prescribed, so evaporation runs at T_p = T_gas = the measured profile.
"""
import sys
evaporation = 'evap' in sys.argv[1:]
kelvin = 'kelvin' in sys.argv[1:]
sys.argv = sys.argv[:1]
import single_Igor_impingingJet_fixedT_particles as case

conf = case.conf
conf.paths.outputDir.value = ('run/single_Igor/single_Igor_impingingJet_fixedT_particles_'
                              + ('evap' if evaporation else 'noevap')
                              + ('_kelvin' if kelvin else ''))
conf.general.nThreads.value = 5
conf.particles.evaporation.value = evaporation
conf.particles.nucleation.value.kelvinLimit.value = kelvin
# Same (ratio, density) as the reference; evaporation and the Kelvin limit
# need each phase's condensed species (nasa_condensed.yaml)
conf.particles.phases.value = [(0.0, 7874.0, ['Fe(a)', 'Fe(c)', 'Fe(d)', 'Fe(L)']),
                               (1.0, 5745.0, ['FeO(s)', 'FeO(L)']),
                               (4/3, 5170.0, ['Fe3O4(s)']),
                               (1.5, 5240.0, ['Fe2O3(s)'])]

if __name__ == '__main__':
    conf.run()
