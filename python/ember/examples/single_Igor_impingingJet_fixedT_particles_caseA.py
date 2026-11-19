#!/usr/bin/python
"""
Candidate (a) of the jet sensitivity grid (doc/ParticleTemperatureConclusions.md,
"Within the collision limit"), optionally with the late iron oxide particles:

    python single_Igor_impingingJet_fixedT_particles_caseA.py          # early particles only
    python single_Igor_impingingJet_fixedT_particles_caseA.py late     # + late particles
    python single_Igor_impingingJet_fixedT_particles_caseA.py none     # no particles (gas Fe reference)

Early particles (case (a)): single_Igor_impingingJet_fixedT_particlesEnergy.py
with particle energy and evaporation, the OpenFOAM surface heats as prescribed
heats, the H2O -> H2 / H2 -> H2O pair heat-neutral (H2 -> H2O = +87 kJ/mol),
OH etching rate x0.1, H -> OH rate x0.55, collision cap on surface reactions,
momentsFloor = 1e-3. Reported for it: T_p - T = 607 K (median over 2-8 mm),
particles gone at 8.6 mm, gas Fe / no-particle gas Fe = 0.28 (min) / 0.62
(mean up to there).

Late particles: the model of Sellmann et al. (Proc. Combust. Inst. 37, 2019),
bi-modal (Jeong & Choi 2003): every gas-phase Fe2O3 molecule becomes a mode-1
particle (monomer conversion), coagulation within and between the modes with
the Fuchs kernel and the fractal collision diameter (Df = 1.8), size-dependent
diffusion. No sintering (the paper gives no sintering time) and no
condensation (not in the paper).
"""
import sys
mode = sys.argv[1] if len(sys.argv) > 1 else 'early'
sys.argv = sys.argv[:1]
import single_Igor_impingingJet_fixedT_particlesEnergy as case
from ember import Particles, NucleationChannel

conf = case.conf
conf.paths.outputDir.value = 'run/single_Igor/caseA_' + mode
conf.general.nThreads.value = 5
conf.terminationCondition.momentsFloor.value = 1e-3

early = conf.particles
early.name.value = 'early'
early.evaporation.value = True
early.surfaceReactionHeat.value = 'prescribed'
early.surfaceCollisionLimit.value = True
# OpenFOAM heats [J/mol] (negative heats the particle), with the H2O/H2 pair
# heat-neutral; rate factors of case (a)
deltaH = {('O2', 'O'): -80e3, ('O', None): -578e3, ('H2O', 'H2'): -87e3,
          ('H2', 'H2O'): +87e3, ('H', 'OH'): -53e3, ('OH', 'FEOH'): -157e3}
rateFactor = {('OH', 'FEOH'): 0.1, ('H', 'OH'): 0.55}
for reaction in early.surfaceReactions.value:
    key = (reaction.reactant.value, reaction.product.value)
    reaction.deltaH.value = deltaH[key]
    reaction.A.value *= rateFactor.get(key, 1.0)

if mode == 'late':
    conf.particleSets.append(Particles(
        name='late',
        model='bimodal',
        metal='Fe',
        phases=[(1.5, 5240.0)],               # Fe2O3 [kg/m3]
        nucleation=NucleationChannel(
            conversionSpecies=['FE2O3(s)'],   # gas-phase Fe2O3 of the mechanism
            conversionRate=1e5,               # [1/s]: every molecule a monomer
        ),
        coagulation=True,
        fractalDimension=1.8,
    ))
elif mode == 'none':
    conf.particleSets[:] = [Particles()]
    conf.particles = conf.particleSets[0]
    conf.terminationCondition.measurement.value = None
    conf.terminationCondition.tEnd.value = 0.05
elif mode != 'early':
    raise ValueError("mode must be 'early', 'late' or 'none'")

if __name__ == '__main__':
    conf.run()
