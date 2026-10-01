# GasFreeSinkCells

## What this tests

The `HomogeneousBox` example never produces a cell that holds sink
particles but zero gas particles: a debug probe for exactly that
condition fired 0 times across full runs. This example builds ICs that
*reliably* produce such cells, to exercise two suspected bugs in the sink
task graph:

1. **Sink-sink merger marking in gas-free cells.** `DOSELF1_SINKS` /
   `DOPAIR1_SINKS_NAIVE` gate gas-swallow interactions on
   `cell->hydro.count != 0`, but the sink-sink swallow loop must not be
   gated the same way. This was fixed in `af185ff5e` (exempt sink-sink
   swallow marking from the gas-count gates). The debug probe in
   `src/runner_doiact_functions_sinks.h` ("sink-sink marking reached
   gas-free cell ...") should fire while running this example.

2. **A task-creation gap in `scheduler_splittasks.c`.** When a hydro
   density self/pair task recurses into progeny, the progeny filter only
   checks `hydro.count` (and `stars.count`/`feedback`), never
   `sinks.count`. A progeny cell that has sinks but no gas and no stars
   is silently dropped -- it gets no density task, and therefore no
   `do_sink_swallow` task either (cloned 1:1 off the density task). The
   debug probes "DROPPED SELF DENSITY TASK for sink-only progeny" and
   "DROPPED INTRA-PARENT PAIR DENSITY TASKS for sink-only progeny" in
   `src/scheduler_splittasks.c` should fire.

## Geometry

The box reproduces the `HomogeneousBox` uniform gas background (same
internal unit system, gas density `0.1 atom/cm3`, gas particle mass
`50 Msun`) with one exception: a cubic corner region, one quarter of the
box side length, contains **zero gas particles**. That corner cube is
sized to coincide exactly with the `[0,0,0]` top-level cell (see
`Scheduler:max_top_level_cells: 4` in `params.yml`, `1 / 4 ==`
`--corner-frac 0.25` in `makeIC.py`).

Inside that corner cube, `makeIC.py` places **2000 sink particles**
in the ICs (gravity particles are matched automatically). `params.yml` raises
`Scheduler:cell_split_size` from its default (400,
`src/space.h:space_splitsize_default`) to **800**, above the ~500 gas
particles the bulk of the box has per top-level cell, so gas-only
top-level cells stay leaves. The 2000 sinks (2.5x that raised threshold)
still force `src/space_split.c` to recursively split the gas-free corner
cell purely on sink `gcount` (`with_self_gravity && gcount >
space_splitsize` does not check `hydro.count`/`stars.count`). Every
descendant of that corner cell inherits zero gas and zero stars, so the
recursive split manufactures deep "sink-only" progeny cells at several
tree depths, confined to that one corner instead of the whole box.
Splitting stops once `cell_can_split_self_hydro_task`'s smoothing-length
guard binds (`src/cell.h`): with the default `cut_off_radius: 5e-3`, that
guard caps the recursion at depth ~4 below the corner cell -- if you
raise `--n-sinks` or shrink `--cluster-frac` enough to want deeper
splitting, shrink `cut_off_radius`/`sink_h` accordingly.

`params.yml` also forces `Scheduler:cell_sub_size_self_hydro` and
`cell_sub_size_self_stars` to **0**. By default (32000,
`src/space.h:space_subsize_self_hydro_default`) a 0-gas-count cell always
satisfies `hydro.count < cell_sub_size_self_hydro`, so
`scheduler_splittasks.c` takes the lazily-recursed "sub-task" branch and
never reaches the explicit per-progeny split code where the "DROPPED
..." probes (bug 2 below) live. Forcing the threshold to 0 makes every
split-eligible cell take the explicit branch instead; since
`cell_split_size` above confines actual splitting (`c->split`) to the
corner hierarchy, this only affects that hierarchy, not the rest of the
box.

The 2000 sinks are packed into a small cube (`--cluster-frac 0.12` of
the corner side, ~26 pc) offset away from the corner cell's midplanes
(`--cluster-offset-frac 0.3`) so the cluster stays inside a single
octant at every recursive split -- this drives the splitting
consistently down one branch instead of spreading thinly across all 8
children. The resulting mean nearest-neighbour separation inside the
cluster (~2 pc) is well inside `GEARSink:cut_off_radius` (5e-3 kpc = 5
pc), so sink-sink mergers are actually attempted once the simulation
runs, not just theoretically possible.

`GEARSink:disable_sink_formation` is set to `1`: the sinks already exist
in the ICs, and we want to isolate sink-sink/task-graph behaviour from
gas -> sink formation happening elsewhere in the box.

## Expected diagnostics

Built with `--enable-debugging-checks --enable-debug` (so
`SWIFT_DEBUG_CHECKS` is defined), the run should print:

- `sink-sink marking reached gas-free cell (...)` from
  `runner_doiact_functions_sinks.h` -- confirms the gas-count gate
  exemption lets sink-sink marking proceed.
- `DROPPED SELF DENSITY TASK for sink-only progeny: ...` and/or
  `DROPPED INTRA-PARENT PAIR DENSITY TASKS for sink-only progeny: ...`
  from `scheduler_splittasks.c` -- confirms the progeny filter gap. This
  fires via *intra-cell* recursion of the corner cell's own self task
  (all 2000 sinks are in one cluster, so no top-level *pair* task has
  sinks on both sides -- to also exercise the cross-cell pair path,
  scatter a few dozen sinks into the neighbouring gas-bearing top-level
  cell).

If the second bug gets fixed, the "DROPPED ..." messages should stop
appearing (progeny with sinks now get their density/swallow tasks), while
the first probe should keep firing (its outcome is independent).

## Files

- `makeIC.py` -- generates the ICs (uniform gas background + gas-free
  sink cluster corner). Run `python3 makeIC.py --help` for all options.
- `params.yml` -- adapted from `HomogeneousBox/params.yml`: same units
  and physics, `max_top_level_cells: 4`, `cell_split_size: 800` and
  `cell_sub_size_self_hydro/stars: 0` to force the corner hierarchy
  through the explicit task-split code path (see Geometry above),
  shortened `time_end` (5 Myr, this test only needs a tree build and a
  few active steps), sink formation disabled.
- `run.sh` -- env-var interface (`n_ranks`, `n_threads`, `level`) similar to
  `HomogeneousBox/run.sh`. This example is sinks-only: `with_star_formation`
  is not supported (GEAR star formation combined with sinks has never been
  tested). ICs carry sink particles directly, and `--sinks` must always be
  passed to ensure they are read.
- `getGrackleCoolingTable.sh`, `getChemistryTable.sh` -- copied
  unmodified from `HomogeneousBox` (same cooling/feedback tables).

## Running

```
n_ranks=8 n_threads=1 examples/SinkParticles/GasFreeSinkCells/run.sh
```

Configure exactly as documented in the repository `CLAUDE.md`:

```
./configure --with-chemistry=GEAR_10 --with-feedback=GEAR --with-cooling=grackle_0 \
            --with-stars=GEAR --with-sink=GEAR --with-star-formation=GEAR \
            --with-kernel=wendland-C2 --with-grackle=$GRACKLE_ROOT \
            --enable-debugging-checks --enable-debug
make -j$(nproc)
```
