/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
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

#ifndef SWIFT_DEFAULT_NEUTRINO_MESH_DELTA_F_H
#define SWIFT_DEFAULT_NEUTRINO_MESH_DELTA_F_H

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

#include "cosmology.h"
#include "neutrino_properties.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Shared information for delta-f neutrino weighting of a cell.
 */
struct neutrino_data {
  char use_mesh_delta_f;
  double *m_eV_array;
  int N_nu;
  double fac;
  long long neutrino_seed;
};

void gather_neutrino_data(const struct space *s, struct neutrino_data *nu_data);
void gpart_neutrino_weight(const struct gpart *gp,
                           const struct neutrino_data *nu_data, double *weight);

#endif /* SWIFT_DEFAULT_NEUTRINO_MESH_DELTA_F_H */
