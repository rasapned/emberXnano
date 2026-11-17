#!/usr/bin/python
"""
single_Igor_impingingJet_fixedT_particlesEnergy.py (particle energy on,
surface reactions on) without / with reversible condensation. Everything
else is taken unchanged from the reference script by importing its Config.
Every SurfaceReaction.deltaH is set to 0. With the default
surfaceReactionHeat = 'thermo' the surface-reaction heat is the
thermochemical one (reactant enthalpy at T_gas in, product enthalpy at T_p
out); with 'prescribed' it is -deltaH, so here no surface-reaction heat.

    python single_Igor_impingingJet_fixedT_particlesEnergy_evap.py [evap] [prescribed | ofdH] [HOH=<kJ/mol>] [OHx<f>]

evap: Particles(evaporation=True); prescribed: surfaceReactionHeat = 'prescribed';
ofdH: 'prescribed' with the OpenFOAM deltaH values (all exothermic);
HOH=<kJ/mol>: with ofdH, replace the H -> OH heat (OpenFOAM -53);
OHx<f>: multiply the OH etching rate A by f.

Evaporation runs at the particle temperature T_p.
"""
import sys
evaporation = 'evap' in sys.argv[1:]
ofdH = 'ofdH' in sys.argv[1:]
HOH = [float(a[4:]) for a in sys.argv[1:] if a.startswith('HOH=')]
OHx = [float(a[3:]) for a in sys.argv[1:] if a.startswith('OHx')]
prescribed = 'prescribed' in sys.argv[1:] or ofdH
sys.argv = sys.argv[:1]
import single_Igor_impingingJet_fixedT_particlesEnergy as case

conf = case.conf
conf.paths.outputDir.value = ('run/single_Igor/single_Igor_impingingJet_fixedT_particlesEnergy_'
                              + ('evap' if evaporation else 'noevap')
                              + ('_prescribedOF' if ofdH else '_prescribed0' if prescribed else '')
                              + (f'_HOH{HOH[0]:+g}' if HOH else '')
                              + (f'_OHx{OHx[0]:g}' if OHx else ''))
conf.general.nThreads.value = 4
# Leave the near-zero, flickering moments at the evaporation front out of the
# 'moments' termination check
conf.terminationCondition.momentsFloor.value = 1e-3
if prescribed:
    conf.particles.surfaceReactionHeat.value = 'prescribed'
conf.particles.evaporation.value = evaporation
# OpenFOAM setup [J/mol], (reactant, product) -> deltaH; negative heats the particle
deltaHOpenFOAM = {('O2', 'O'): -80e3, ('O', None): -578e3, ('H2O', 'H2'): -87e3,
                  ('H2', 'H2O'): -52e3, ('H', 'OH'): -53e3, ('OH', 'FEOH'): -157e3}
if HOH:
    deltaHOpenFOAM[('H', 'OH')] = HOH[0] * 1e3
for reaction in conf.particles.surfaceReactions.value:
    key = (reaction.reactant.value, reaction.product.value)
    reaction.deltaH.value = deltaHOpenFOAM[key] if ofdH else 0.0
    if OHx and key == ('OH', 'FEOH'):
        reaction.A.value *= OHx[0]

if __name__ == '__main__':
    conf.run()
