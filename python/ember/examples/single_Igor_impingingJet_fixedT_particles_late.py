#!/usr/bin/python
"""
single_Igor_impingingJet_fixedT_particlesEnergy.py (early particles: Fe
clusters, single-mode moments with particle energy) plus an independent
population of late iron oxide particles formed from the FE2O3 precursor,
with the bi-modal model of Sellmann et al. (Proc. Combust. Inst. 37, 2019;
Jeong & Choi 2003): every FE2O3 molecule becomes a mode-1 particle
(monomer conversion), which coagulate (Fuchs kernel, fractal aggregates,
Df = 1.8) into mode-2 aggregates; FE2O3 also condenses onto them.
Everything else is taken unchanged from the reference script by importing
its Config. The two populations interact only through the gas.

    python single_Igor_impingingJet_fixedT_particles_late.py

Output: the early particles under the keys *_early (moments_early,
numberDensity_early, ...), the late ones under *_late (moments_late = N1,
N2, M2, A2; numberDensity1_late, massMeanDiameter_late, collisionDiameter_late,
...).
"""
import single_Igor_impingingJet_fixedT_particlesEnergy as case
from ember import Particles, NucleationChannel

conf = case.conf
conf.paths.outputDir.value = 'run/single_Igor/single_Igor_impingingJet_fixedT_particles_late'
conf.general.nThreads.value = 6
conf.particles.name.value = 'early'

conf.particleSets.append(Particles(
    name='late',
    model='bimodal',
    metal='Fe',
    # the particle material: Fe2O3 (O/Fe = 1.5) [kg/m3]
    phases=[(1.5, 5240.0)],
    nucleation=NucleationChannel(
        conversionSpecies=['FE2O3(s)'],   # gas-phase Fe2O3 of the mechanism
        conversionRate=1e5,               # [1/s]: fast, "every molecule a monomer"
    ),
    coagulation=True,
    fractalDimension=1.8,                 # Sellmann et al. (after Poliak et al.)
    condensationSpecies=['FE2O3(s)'],
    # Sintering (A [s/(K m^4)], Ta [K]) of tau = A T d_pp^4 exp(Ta/T) has no
    # default: set the values for Fe2O3 to enable it. Without, the aggregates
    # keep their surface.
    sintering=None,
))

if __name__ == '__main__':
    conf.run()
