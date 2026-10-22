# GasFreeSinkCells

A uniform gas box with one corner that holds sink particles but no gas.
`HomogeneousBox` never makes cells like this.

## What this tests

1. **Sink-sink mergers in gas-free cells.** The sink density and gas-swallow
   loops are skipped when `cell->hydro.count == 0`. The merger loop must
   still run.
2. **Task splitting of sink-only cells.** A child cell with sinks but no gas
   and no stars must keep its density task, because all sink tasks are
   cloned from it (`scheduler_splittasks.c`).
3. **Star formation from a dense sink cluster.** Many stars spawn in one
   top-level cell. The `cell_extra_*` pools in `params.yml` are sized for
   this.

## Setup

- The gas background is the same as in `HomogeneousBox`: density
  `0.1 atom/cm3`, particle mass `50 Msun`.
- A corner cube of 1/4 of the box side has no gas. It matches the
  `[0,0,0]` top-level cell (`max_top_level_cells: 4`, `--corner-frac 0.25`).
- `makeIC.py` puts 2000 sinks in a small cube inside that corner
  (`--cluster-frac 0.12`, offset by `--cluster-offset-frac 0.3`). The mean
  neighbour distance (~2 pc) is below `GEARSink:cut_off_radius` (5 pc), so
  mergers do happen.
- `cell_split_size: 800` keeps gas-only cells as leaves. The 2000 sinks
  still force the corner cell to split, down to depth ~4. For more sinks or
  a tighter cluster, reduce `cut_off_radius` and `sink_h`.
- `cell_sub_size_self_hydro` and `cell_sub_size_self_stars` are 0. With the
  defaults, a cell without gas would skip the explicit split code that item
  2 needs.
- `GEARSink:disable_sink_formation: 1`. The sinks come from the ICs, so the
  test only sees sink-sink and task-graph behaviour.

## Files

- `makeIC.py`: makes the ICs. `python3 makeIC.py --help` lists the options.
- `params.yml`: based on `HomogeneousBox/params.yml`, with the changes above,
  larger `cell_extra_gparts`/`cell_extra_sparts` pools, and `time_end` of
  5 Myr.
- `run.sh`: options `n_ranks`, `n_threads`, `level`. It always passes
  `--sinks`. Star formation is not supported here, because GEAR star
  formation with sinks is untested.

`run.sh` downloads the cooling and feedback tables with the scripts in
`examples/GEAR_ICs_and_SCRIPTS/`.

## Running

```
n_ranks=8 n_threads=1 examples/SinkParticles/GasFreeSinkCells/run.sh
```

Configure with the GEAR sink model:

```
./configure --with-chemistry=GEAR_10 --with-feedback=GEAR --with-cooling=grackle_0 \
            --with-stars=GEAR --with-sink=GEAR --with-star-formation=GEAR \
            --with-kernel=wendland-C2 --with-grackle=$GRACKLE_ROOT \
            --enable-debugging-checks --enable-debug
make -j$(nproc)
```
