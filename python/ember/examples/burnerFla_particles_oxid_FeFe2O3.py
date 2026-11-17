#!/usr/bin/python
"""
burnerFla_particles_oxid.py with only two particle phases, Fe and Fe2O3 (no
FeO / Fe3O4), as OpenFOAM, which knows only the Fe and the fully oxidized
density. Everything else is taken unchanged from the reference script by
importing its Config. The density is then linear in O/Fe between 7874 and
5240 kg/m3 (OpenFOAM adds the Fe and O volumes, bulkRho 7874 and
volatileRho 2946, which gives the same end points).

    python burnerFla_particles_oxid_FeFe2O3.py          # no evaporation
    python burnerFla_particles_oxid_FeFe2O3.py evap     # with evaporation

With evaporation only the Fe clusters evaporate (lever-rule Fe activity
1 - x_O/1.5); FeO has no matching phase and only condenses.
"""
import sys
evaporation = 'evap' in sys.argv[1:]
sys.argv = sys.argv[:1]
import burnerFla_particles_oxid as case

conf = case.conf
conf.paths.outputDir.value = 'run/burnerFla/burnerFla_particles_oxid_FeFe2O3' + ('_evap' if evaporation else '')
conf.general.nThreads.value = 8
conf.particles.evaporation.value = evaporation
conf.particles.phases.value = [(0.0, 7874.0, ['Fe(a)', 'Fe(c)', 'Fe(d)', 'Fe(L)']),
                               (1.5, 5240.0, ['Fe2O3(s)'])]
# 0.5 nm in the reference script; 0.4 nm as OpenFOAM
conf.particles.minParticleDiameter.value = 0.4e-9

if __name__ == '__main__':
    conf.run()
