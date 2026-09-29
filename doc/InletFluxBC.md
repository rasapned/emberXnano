# Burner-face boundary condition (`BoundaryCondition::InletFlux`)

The burner face fixes the **total flux** of each species entering the domain, not its
**concentration**. Light species (H, H2, OH) can therefore diffuse back toward the burner,
and the face composition settles wherever the flux balance puts it. This is the same
condition Cantera's `BurnerFlame` uses (`Inlet1D`).

It is selected automatically in `FlameSolver::updateBC` for a burner-stabilized case:
unburned stream on the left, `fixedLeftLocation=True`, not twin/cylindrical.

## Why `FixedValue` was wrong

`FixedValue` forces the face composition to the supplied mixture, `Y_k(0) = Y_k,u`. The mass
entering the domain is then more than what the burner supplies:

```
actual inflow = mdot*Y_k(0) + j_k(0) = mdot*Y_k,u + j_k(0)
```

Here `j_k` is the diffusive flux. In a low-pressure H2 flame, diffusion near the face is large
and runs in both directions:

- **H2** diffuses downstream toward the flame, where it is consumed. The fixed face acts as an
  extra H2 source, so more H2 enters than the burner supplies.
- **H, OH and H2O** diffuse upstream toward the burner. The face held them at zero, so it acted
  as a sink that destroyed them.

So element fluxes were not conserved at the inlet.

## The flux-balance condition

The physical requirement is that what leaves the burner equals what enters the flame:

```
mdot*Y_k,u  =  mdot*Y_k(0)  +  j_k(0)
  supply        convection     diffusion at the face
```

`Y_k(0)` is left free to satisfy this equation:

- **H2** diffuses into the domain (`j > 0`), so `Y_H2(0) < Y_H2,u`. The H2 value at the face
  sits below the feed value.
- **H** is absent from the feed but diffuses upstream (`j < 0`), so `mdot*Y_H(0) = -j_H`. H is
  present at the face even though the burner supplies none.

## How ember implements it

Ember treats the face as a small control volume from `x0 = 0` to `x1`. Each split operator
contributes one piece of the balance.

1. **Convection** (`convectionSystem.cpp`, left-BC block of `f()`) is upwind inflow from the
   supply stream:

   ```
   dY0/dt = -mdot*(Y0 - Yleft)/x1
   ```

   `Yleft` must be the supplied composition. It is taken from the first point of the initial
   profile (`loadProfile`, `Yleft = Y.col(grid.ju)`), so an initial profile must have
   `T[0] = Tu`, `Y[0] = reactants`. Otherwise, a Cantera-based start feeds its face composition
   in as the feed.

2. **Diffusion** (`diffusionSystem.cpp`, `get_A`) is the flux across the `x0`-`x1` face:

   ```
   dY0/dt = rhoD_half*(Y1 - Y0)/x1^2
   ```

   `c2[0]` makes point 1 receive exactly what point 0 loses. Without it, the flux appeared from
   nothing or vanished at point 1, which drained about a third of the inlet H.

3. **At steady state** the storage term is zero, and the two terms combine to

   ```
   mdot*(Yu - Y0) = -rhoD*(Y1 - Y0)/x1 = j(0)
   ```

   That is exactly the balance above.

## What the face does not do

- **Temperature stays fixed at `Tleft`**, so the burner acts as a heat sink, as in Cantera
  (`FlameSolver::setupStep`). Because `T(0)` is prescribed:
  - the energy equation uses the fixed-value stencil (`prepareIntegrators`);
  - the T rates at `j = 0` are zeroed (`finishStep`).

  Without the zeroing, those rates entered `drhodt`, continuity turned them into a false mass
  source (mass flux 71% too high), and the dTdt termination check never converged. The heat
  loss to the burner comes naturally from the conduction gradient between points 0 and 1.

- **No chemistry at the face** (`integrateProductionTerms`, `calculateQdot`). Cantera's inlet
  boundary has no reactions either. With reactions allowed, the diffused H and OH made the face
  act as an igniter.

- **The mass flux is fixed:** `continuityBC='fixedLeft'` holds `V(0) = mdot`.

## Comparison with Cantera

Cantera imposes the balance on a zero-width boundary node, while ember uses a finite cell of
width `x1`. At steady state the storage term vanishes, so both reduce to the same equation.
They differ only through how each code computes `j_k`:

- ember uses a mass-fraction gradient: `-rho*D_k*grad(Y_k)`;
- Cantera uses a mole-fraction gradient: `-rho*(W_k/W_mix)*D_k*grad(X_k)`.

That is the source of the small H2 and H difference at the face that remains.

## Related notes

- **Strang splitting:** `drhodt` only includes `ddtCross` when `splittingMethod == "balanced"`.
  Under Strang this leaks mass (V about +7.6%), so use `splittingMethod='balanced'` for
  burner-stabilized flames.
- **Validation case:** `python/ember/examples/burnerFla_test.py`. `initMode='cantera'` starts
  from the Cantera profile; `initMode='shifted'` starts from the same profile moved 4 mm
  downstream, to test independence from the initial condition. At t = 0.0248 s the shifted run
  had its flame at `x(1000 K)` = 4.05 mm (Cantera: 4.07 mm), with the outlet temperature still
  relaxing (2222 K, heading toward about 2178 K).
