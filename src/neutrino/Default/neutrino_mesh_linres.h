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

#ifndef SWIFT_DEFAULT_NEUTRINO_MESH_LINRES_H
#define SWIFT_DEFAULT_NEUTRINO_MESH_LINRES_H

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

#include "cosmology.h"
#include "neutrino_properties.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Structure for handling the linear neutrino response on the mesh
 */
struct neutrino_mesh {

  /*! Logarithm of minimum scale factor for which we have transfer functions */
  double log_a_min;

  /*! Logarithm of maximum scale factor for which we have transfer functions */
  double log_a_max;

  /*! Logarithm of minimum wavenumber for which we have transfer functions */
  double log_k_min;

  /*! Logarithm of maximum wavenumber for which we have transfer functions */
  double log_k_max;

  /*! Spacing of log scale factor for transfer function interpolation */
  double delta_log_a;

  /*! Spacing of log wavenumber for transfer function interpolation */
  double delta_log_k;

  /*! Table of the neutrino overdensity to cdm & baryon overdensity ratio */
  double *ncdm_over_cb;
};

void neutrino_mesh_init(struct neutrino_mesh *numesh,
                        struct swift_params *params,
                        const struct unit_system *us, const double dim[3],
                        const struct cosmology *c,
                        const struct neutrino_props *np,
                        const struct gravity_props *gp, int verbose);
void neutrino_mesh_clean(struct neutrino_mesh *numesh);
void neutrino_mesh_compute(const struct space *s, struct pm_mesh *mesh,
                           struct threadpool *tp, fftw_complex *frho,
                           const int slice_offset, const int slice_width,
                           int verbose);
void neutrino_mesh_struct_dump(int enabled, const struct neutrino_mesh *numesh,
                               FILE *stream);
void neutrino_mesh_struct_restore(int enabled, struct neutrino_mesh *numesh,
                                  FILE *stream);

#endif /* SWIFT_DEFAULT_NEUTRINO_MESH_LINRES_H */
