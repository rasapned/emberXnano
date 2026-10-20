---
name: physics-changelog
description: Write a technical physics/math summary of an implementation change in ember (equations, discretisation, boundary conditions, particle moments/source terms) into doc/PhysicsChangelog.md and keep doc/PhysicsAndNumerics.md current. Use whenever a feat/fix commit touching src/ or python/ember/*.py|.pyx is being prepared or was just made, or when the user asks to document what a commit/range changed physically or numerically.
---

# Physics changelog for ember

Goal: for every implementation commit, a reader who knows CFD/combustion should see **what
physics or numerics changed, as equations**, without reading the diff.

## Inputs

- No argument: document the **staged** change (`git diff --cached`) before it is committed, so
  the doc update goes into the same commit. If nothing is staged, document `HEAD`.
- Argument = a commit hash or range (`abc123`, `abc123..HEAD`): document each implementation
  commit in it, oldest first.

## Steps

1. **Collect the change.** Run `git diff --cached --stat` (or `git show --stat <hash>`), then
   read the full diff of `src/`, `python/ember/_ember.pyx`, `python/ember/input.py`, and
   `python/ember/output.py`. For case scripts in `python/ember/examples/`, only note new
   options they use or new validation results.
2. **Classify.** If the commit is pure `ref`/`test`/`docs` with no change to equations,
   discretisation, BCs, or the options' meaning, add a one-line entry ("no physics/numerics
   change") and stop.
3. **Read enough context** around each changed function to state the equation it now solves:
   the caller, units of the inputs, and which operator it belongs to (diffusion, convection,
   production, cross term, or BC). `doc/PhysicsAndNumerics.md` gives the notation and the
   current model.
4. **Write the entry** at the top of `doc/PhysicsChangelog.md`, below the header, using the
   template below.
5. **Update `doc/PhysicsAndNumerics.md`:**
   - edit the affected section so it describes the new state;
   - add a row to the table in §7 "Recent implementations";
   - add or remove items in §8 "Known limitations".
   Keep the notation of that file.
6. Report to the user, in a few lines, what was documented and any inconsistency you found.

## Entry template

```markdown
## <YYYY-MM-DD> — <commit subject> (<short hash, or "staged">)

**Area:** gas transport | chemistry | particles: <process> | BC: <name> | numerics: <operator> | I/O

**What changed (physics):** 1–3 sentences: the phenomenon, the model choice, and the
assumptions.

**Equations:** the continuous form actually implemented, in LaTeX ($$…$$), with every symbol
and its unit defined once. Give the old form too if it was replaced.

**Discretisation / numerics:** stencil, BC row, operator it lives in, integrator, tolerances,
explicit vs implicit, and what is frozen over the split step. Only if relevant.

**Conservation & units:** which invariant holds (e.g. $\sum Y_k + M_M + M_O = 1$, element
fluxes), and the unit convention (kmol-based $N_A$, mass-specific moments).

**Code:** `file.cpp::function`, the option names in `input.py`.

**Verification:** test/case and numbers, if any were run; otherwise "not verified".

**Limits / open points:** what is neglected, and when the model is invalid.
```

## Rules

- **Derive equations from the code, not from the commit message or from literature.** If the
  code differs from the textbook form (e.g. a missing factor, or a different limit), write what
  the code does and flag the difference under "Limits / open points".
- **Units:** check them for every new term against ember's conventions: Cantera kmol units,
  $N_A = 6.022\times10^{26}$/kmol, and moments per kg of gas. Report mismatches; do not
  silently fix code.
- **Brevity:** at most about 40 lines per entry. One equation per process; no derivations.
- **Transient vs steady:** say whether the new term affects transients only (e.g. controllers,
  split constants) or the steady solution as well.
- **Scope:** do not touch code files; this skill only writes the two docs. Do not commit unless
  the user asked to commit.
