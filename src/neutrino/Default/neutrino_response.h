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

#ifndef SWIFT_DEFAULT_NEUTRINO_RESPONSE_H
#define SWIFT_DEFAULT_NEUTRINO_RESPONSE_H

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
struct neutrino_response {

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
  double *pt_density_ratio;

  /*! Size of the transfer function ratio table */
  hsize_t tf_size;

  /*! Whether to use a fixed present-day background density */
  char fixed_bg_density;
};

void neutrino_response_init(struct neutrino_response *numesh,
                            struct swift_params *params,
                            const struct unit_system *us, const double dim[3],
                            const struct cosmology *c,
                            const struct neutrino_props *np,
                            const struct gravity_props *gp, int rank,
                            int verbose);
void neutrino_response_clean(struct neutrino_response *numesh);
void neutrino_response_compute(const struct space *s, struct pm_mesh *mesh,
                               struct threadpool *tp, fftw_complex *frho,
                               const int slice_offset, const int slice_width,
                               int verbose);
void neutrino_response_struct_dump(const struct neutrino_response *numesh,
                                   FILE *stream);
void neutrino_response_struct_restore(struct neutrino_response *numesh,
                                      FILE *stream);

#endif /* SWIFT_DEFAULT_NEUTRINO_RESPONSE_H */
