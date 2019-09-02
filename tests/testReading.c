/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

/* Some standard headers. */
#include "../config.h"

/* Some standard headers. */
#include <stdlib.h>

/* Includes. */
#include "swift.h"

int main(int argc, char *argv[]) {

  size_t Ngas = 0, Ngpart = 0, Ngpart_background = 0, Nspart = 0, Nbpart = 0;
  int flag_entropy_ICs = -1;
  int i, j, k;
  double dim[3];
  struct part *parts = NULL;
  struct gpart *gparts = NULL;
  struct spart *sparts = NULL;
  struct bpart *bparts = NULL;

  /* Default unit system */
  struct unit_system us;
  units_init_cgs(&us);

  /* Properties of the ICs */
  const double boxSize = 1.;
  const size_t L = 4;
  const double rho = 2.;
#ifdef CHEMISTRY_GRACKLE
  const float he_density = rho * 0.24;
#endif

  /* Read data */
  read_ic_single("input.hdf5", &us, dim, &parts, &gparts, &sparts, &bparts,
                 &Ngas, &Ngpart, &Ngpart_background, &Nspart, &Nbpart,
                 &flag_entropy_ICs,
                 /*with_hydro=*/1,
                 /*with_gravity=*/1,
                 /*with_stars=*/0,
                 /*with_black_holes=*/0,
                 /*with_cosmology=*/0,
                 /*cleanup_h=*/0,
                 /*cleanup_sqrt_a=*/0,
                 /*h=*/1., /*a=*/1., /*n_threads=*/1, /*dry_run=*/0);

  /* Check global properties read are correct */
  assert(dim[0] == boxSize);
  assert(dim[1] == boxSize);
  assert(dim[2] == boxSize);
  assert(Ngas == L * L * L);
  assert(Ngpart == L * L * L);

  /* Check particles */
  for (size_t n = 0; n < Ngas; ++n) {

    /* Check that indices are in a reasonable range */
    unsigned long long index = parts[n].id;
    assert(index < Ngas);

    /* Check masses */
    float mass = hydro_get_mass(&parts[n]);
    float correct_mass = boxSize * boxSize * boxSize * rho / Ngas;
    assert(mass == correct_mass);

    /* Check smoothing length */
    float h = parts[n].h;
    float correct_h = 2.251 * boxSize / L;
    assert(h == correct_h);

    /* Check velocity */
    assert(parts[n].v[0] == 0.);
    assert(parts[n].v[1] == 0.);
    assert(parts[n].v[2] == 0.);

    /* Check positions */
    k = index % 4;
    j = ((index - k) / 4) % 4;
    i = (index - k - 4 * j) / 16;
    double correct_x = i * boxSize / L + boxSize / (2 * L);
    double correct_y = j * boxSize / L + boxSize / (2 * L);
    double correct_z = k * boxSize / L + boxSize / (2 * L);
    assert(parts[n].x[0] == correct_x);
    assert(parts[n].x[1] == correct_y);
    assert(parts[n].x[2] == correct_z);

    /* Check accelerations */
    assert(parts[n].a_hydro[0] == 0.);
    assert(parts[n].a_hydro[1] == 0.);
    assert(parts[n].a_hydro[2] == 0.);

#ifdef CHEMISTRY_GRACKLE
    assert(parts[n].chemistry_data.he_density == he_density);
#endif
  }

  /* Clean-up */
  free(parts);
  free(gparts);

  return 0;
}
