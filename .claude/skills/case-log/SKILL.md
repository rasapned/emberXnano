---
name: case-log
description: Dated log of parameter changes tested in ember example cases (python/ember/examples/*.py, run/ folders), which are local-only and not in git. Use when a case parameter, mechanism, grid, tolerance or BC is changed for a test, when a test's outcome is known (kept / reverted), when the user asks whether something was already tried, and after starting or finishing ember runs. Also invoked as /case-log [note | query].
---

# Case parameter log

Files (all in `python/ember/examples/caseLog/`, git-ignored, local only):
- `CASE_LOG.md` — human log, one entry per tested change, **newest first**.
- `configs/<run folder>/<YYYY-MM-DD_HHMMSS>/` — automatic snapshots: raw `config`,
  normalised `config.norm`, and the case script as it was at snapshot time.
- `snapshot_configs.py` — makes the snapshots; run from `python/ember/examples/`
  with the ember-build python (`/home/piotr/miniforge3/envs/ember-build/bin/python`).

## Always first

Run `python caseLog/snapshot_configs.py`. It snapshots every `run/*/config` whose
normalised content changed since the last snapshot and prints the diffs. Run it
also **before** starting a run into an existing output folder (the run overwrites
`config`).

## Recording (a change was made, a run finished, or an outcome is known)

1. Take the parameter change from the snapshot diff, `git diff` of the case script
   (if tracked), or the session. `config` misses script-level variables (e.g.
   `u_in`, initMode, mechanism-file edits) and the content of nucleation/surface
   reaction objects — read the saved script copy for those.
2. Append an entry at the top of `CASE_LOG.md` (below the template), using the
   template exactly:
   ```
   ## YYYY-MM-DD — <case script> (<run folder>)
   - **Change:** <parameter>: <old> → <new>
   - **Why:** <reason / hypothesis>
   - **Result:** <numbers, comparison, figure path; "pending" if not run yet>
   - **Status:** kept | reverted | superseded | open
   ```
   Use absolute dates (today's date from the session). One entry per logical change;
   group parameters changed together for one test.
3. When an outcome becomes known later, **update the existing entry** (Result,
   Status) rather than adding a duplicate; if a change is undone, set it to
   `reverted` and say when. If the user's reason or outcome is unknown, ask in one
   line instead of inventing it.

## Querying ("was X tried?", "what did we do with Y?")

1. `grep -n -i` the term in `CASE_LOG.md`.
2. `python caseLog/snapshot_configs.py --grep <option or value>` — every snapshot
   (config or script copy) containing it.
3. `python caseLog/snapshot_configs.py --history <run folder>` — all snapshots of a
   case with diffs between them.

Answer with dates, values tried, status and result; say if something appears in
snapshots but has no log entry (tested but never annotated).
