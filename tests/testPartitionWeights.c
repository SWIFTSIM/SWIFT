/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Will J. Roper (w.roper@sussex.ac.uk).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/**
 * @file testPartitionWeights.c
 * @brief Regression test for uniform and zoom partition weight indexing.
 *
 * Calls the production part, gpart, and spart accumulation mappers directly,
 * before MPI reduction, gravity percentile clipping, or integer-range scaling.
 * Expected byte weights come from known occupancy, not the coordinate-to-cell
 * helper under test. This catches both collapsing zoom particles onto the
 * background grid and omitting the background grid's cell-array offset.
 *
 * The fixture is a periodic unit box with a 4^3 background grid. Background
 * cell (1, 1, 1), index 21, contains a 2^3 zoom grid covering [0.25, 0.5)^3.
 * In zoom mode the eight zoom cells precede the 64 background cells, so the
 * void cell is at combined index 8 + 21. In uniform mode there is no offset.
 *
 * Built and registered with Automake only when MPI is enabled. Run from the
 * repository root after building SWIFT:
 * @code
 * make -C tests check TESTS=testPartitionWeights
 * @endcode
 * No launcher or IC files are needed. MPI is initialised for SWIFT's error
 * reporting; under mpirun each rank repeats the same local checks
 * independently. This does not test collective reduction, partition selection,
 * or gravity.
 */

/* Config parameters. */
#include <config.h>

/* Standard headers. */
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

/**
 * @brief Context consumed by the production accumulation mappers.
 *
 * Must match struct counts_mapper_data in src/partition.c. Keeping this test
 * declaration local avoids exposing the private mapper context as a public API.
 */
struct counts_mapper_data {
  /*! Per-cell accumulated byte weights. */
  double *counts;
  /*! Size in bytes of the particle type being accumulated. */
  size_t size;
  /*! Space describing the uniform or nested zoom cell grids. */
  struct space *s;
};

/* Private mapper entry points defined in src/partition.c. */
void partition_accumulate_sizes_mapper_part(void *map_data, int num_elements,
                                            void *extra_data);
void partition_accumulate_sizes_mapper_gpart(void *map_data, int num_elements,
                                             void *extra_data);
void partition_accumulate_sizes_mapper_spart(void *map_data, int num_elements,
                                             void *extra_data);

/* Eight zoom-cell centres, two first-background-cell particles (one requiring
 * periodic wrapping), and one last-background-cell particle. */
static const double positions[11][3] = {
    {0.3125, 0.3125, 0.3125}, {0.3125, 0.3125, 0.4375},
    {0.3125, 0.4375, 0.3125}, {0.3125, 0.4375, 0.4375},
    {0.4375, 0.3125, 0.3125}, {0.4375, 0.3125, 0.4375},
    {0.4375, 0.4375, 0.3125}, {0.4375, 0.4375, 0.4375},
    {0.125, 0.125, 0.125},    {1.125, -0.875, 1.125},
    {0.875, 0.875, 0.875}};

/**
 * @brief Check every raw cell weight against the analytical fixture occupancy.
 *
 * Zoom cells 0--7 each receive one particle; background cells 8 and 71 receive
 * two and one respectively. On the uniform grid the eight interior particles
 * all occupy cell 21, and background particles occupy cells 0 and 63.
 * All other entries, including the unused tail in uniform mode, must stay zero.
 * Byte counts are small integers exactly representable as doubles, so no
 * floating-point tolerance is needed.
 *
 * @param s The mock #space defining the cell grids.
 * @param particles Eleven particles with the fixture coordinates. The mapper
 * normalises their coordinates in place for periodic boundaries.
 * @param size The size of one particle in bytes.
 * @param mapper The production accumulation mapper for this particle type.
 */
static void check_weights(struct space *s, void *particles, const size_t size,
                          void (*mapper)(void *, int, void *)) {
  double counts[72] = {0};
  double expected[72] = {0};
  struct counts_mapper_data data = {counts, size, s};

  /* Independently known occupancy: eight fine cell centres in the zoom,
   * two particles in the first background cell (one periodically wrapped),
   * and one particle in the last background cell. */
  if (s->with_zoom_region) {
    for (int i = 0; i < 8; i++) expected[i] = size;
    expected[8] = 2 * size;
    expected[71] = size;
  } else {
    expected[21] = 8 * size;
    expected[0] = 2 * size;
    expected[63] = size;
  }

  mapper(particles, 11, &data);
  for (int i = 0; i < 72; i++)
    if (counts[i] != expected[i])
      error("Cell %d weight is %g, expected %g (zoom=%d, particle size=%zu)", i,
            counts[i], expected[i], s->with_zoom_region, size);
}

/* Use fresh particles for each check because the mapper wraps them in place. */
#define CHECK_TYPE(TYPE)                                          \
  do {                                                            \
    struct TYPE particles[11] = {0};                              \
    for (int i = 0; i < 11; i++)                                  \
      memcpy(particles[i].x, positions[i], sizeof(positions[i])); \
    check_weights(&s, particles, sizeof(struct TYPE),             \
                  partition_accumulate_sizes_mapper_##TYPE);      \
  } while (0)

/**
 * @brief Check gas, gravity, and star particle weights on both grid layouts.
 */
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &engine_rank);

  /* Initialize CPU frequency, this also starts time for error messages. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Only geometry and the void-cell marker are needed by the index lookup. */
  struct space s = {0};
  struct zoom_region_properties zoom = {0};
  if (posix_memalign((void **)&s.cells_top, cell_align,
                     72 * sizeof(struct cell)) != 0)
    error("Failed to allocate test cells");
  bzero(s.cells_top, 72 * sizeof(struct cell));
  s.periodic = 1;
  s.zoom_props = &zoom;
  zoom.bkg_cell_offset = 8;
  for (int i = 0; i < 3; i++) {
    s.dim[i] = 1.;
    s.cdim[i] = 4;
    s.iwidth[i] = 4.;
    zoom.cdim[i] = 2;
    zoom.iwidth[i] = 8.;
    zoom.region_lower_bounds[i] = 0.25;
  }
  s.cells_top[8 + 21].subtype = cell_subtype_void;

  for (int with_zoom = 0; with_zoom < 2; with_zoom++) {
    s.with_zoom_region = with_zoom;
    s.nr_cells = with_zoom ? 72 : 64;
    CHECK_TYPE(part);
    CHECK_TYPE(gpart);
    CHECK_TYPE(spart);
  }

  free(s.cells_top);
  if (engine_rank == 0)
    message("Uniform and zoom pre-clipping partition weights passed.");
  MPI_Finalize();
  return 0;
}
