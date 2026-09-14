# GasFreeSinkCells

## What this tests

The `HomogeneousBox` example never produces a cell that holds sink
particles but zero gas particles. This example builds ICs that reliably
produce such cells, at several tree depths, to exercise the parts of the
sink task graph that only see that topology:

1. **Sink-sink merging in gas-free cells.** The sink density and
   gas-swallow loops are gated on `cell->hydro.count != 0`; the sink-sink
   swallow (merger) loop must run even when a cell has no gas at all.

2. **Explicit task splitting of sink-only progeny.** When a hydro
   density self/pair task is split over a cell's progeny
   (`scheduler_splittasks.c`), the progeny filters decide which children
   keep a task. A child with sinks but no gas and no stars must keep its
   density task, since every sink task is cloned off it.

3. **Star formation from sinks in a gas-free region.** A dense sink
   cluster spawns many stars in one top-level cell; the `cell_extra_*`
   pools in `params.yml` are sized for that (see the comments there).

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
guard caps the recursion at depth ~4 below the corner cell. If you
raise `--n-sinks` or shrink `--cluster-frac` enough to want deeper
splitting, shrink `cut_off_radius`/`sink_h` accordingly.

`params.yml` also forces `Scheduler:cell_sub_size_self_hydro` and
`cell_sub_size_self_stars` to **0**. By default (32000,
`src/space.h:space_subsize_self_hydro_default`) a 0-gas-count cell always
satisfies `hydro.count < cell_sub_size_self_hydro`, so
`scheduler_splittasks.c` takes the lazily-recursed "sub-task" branch and
never reaches the explicit per-progeny split filters (item 2 above).
Forcing the threshold to 0 makes every split-eligible cell take the
explicit branch instead; since `cell_split_size` above confines actual
splitting (`c->split`) to the corner hierarchy, this only affects that
hierarchy, not the rest of the box.

The 2000 sinks are packed into a small cube (`--cluster-frac 0.12` of
the corner side, ~26 pc) offset away from the corner cell's midplanes
(`--cluster-offset-frac 0.3`) so the cluster stays inside a single
octant at every recursive split. This drives the splitting
consistently down one branch instead of spreading thinly across all 8
children. The resulting mean nearest-neighbour separation inside the
cluster (~2 pc) is well inside `GEARSink:cut_off_radius` (5e-3 kpc = 5
pc), so sink-sink mergers are actually attempted once the simulation
runs, not just theoretically possible.

`GEARSink:disable_sink_formation` is set to `1`: the sinks already exist
in the ICs, and we want to isolate sink-sink/task-graph behaviour from
gas -> sink formation happening elsewhere in the box.

## Files

- `makeIC.py`: generates the ICs (uniform gas background + gas-free
  sink cluster corner). Run `python3 makeIC.py --help` for all options.
- `params.yml`: adapted from `HomogeneousBox/params.yml`: same units
  and physics, `max_top_level_cells: 4`, `cell_split_size: 800` and
  `cell_sub_size_self_hydro/stars: 0` to force the corner hierarchy
  through the explicit task-split code path (see Geometry above),
  larger `cell_extra_gparts`/`cell_extra_sparts` pools for the stars the
  sink cluster spawns, shortened `time_end` (5 Myr, this test only needs
  a tree build and a few active steps), sink formation disabled.
- `run.sh`: env-var interface (`n_ranks`, `n_threads`, `level`) similar to
  `HomogeneousBox/run.sh`. This example is sinks-only: `with_star_formation`
  is not supported (GEAR star formation combined with sinks has never been
  tested). ICs carry sink particles directly, and `--sinks` must always be
  passed to ensure they are read.
- `getGrackleCoolingTable.sh`, `getChemistryTable.sh`: copied
  unmodified from `HomogeneousBox` (same cooling/feedback tables).

## Running

```
n_ranks=8 n_threads=1 examples/SinkParticles/GasFreeSinkCells/run.sh
```

Configure with the GEAR sink model, e.g.

```
./configure --with-chemistry=GEAR_10 --with-feedback=GEAR --with-cooling=grackle_0 \
            --with-stars=GEAR --with-sink=GEAR --with-star-formation=GEAR \
            --with-kernel=wendland-C2 --with-grackle=$GRACKLE_ROOT \
            --enable-debugging-checks --enable-debug
make -j$(nproc)
```
