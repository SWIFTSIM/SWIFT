/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_GRAVITY_PROPERTIES
#define SWIFT_GRAVITY_PROPERTIES

/* Config parameters. */
#include "../config.h"

#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/* Local includes. */
#include "parser.h"

/**
 * @brief Contains all the constants and parameters of the self-gravity scheme
 */
struct gravity_props {

  /*! Mesh smoothing scale in units of top-level cell size */
  float a_smooth;

  /*! Distance below which the truncated mesh force is Newtonian in units of
   * a_smooth */
  float r_cut_min;

  /*! Distance above which the truncated mesh force is negligible in units of
   * a_smooth */
  float r_cut_max;

  /*! Time integration dimensionless multiplier */
  float eta;

  /*! Tree opening angle (Multipole acceptance criterion) */
  double theta_crit;

  /*! Inverse of opening angle */
  double theta_crit_inv;

  /*! Softening length */
  double epsilon;

  /*! Square of softening length */
  double epsilon2;

  /*! Inverse of softening length */
  double epsilon_inv;
};

void gravity_props_print(const struct gravity_props *p);
void gravity_props_init(struct gravity_props *p,
                        const struct swift_params *params);

#if defined(HAVE_HDF5)
void gravity_props_print_snapshot(hid_t h_grpsph,
                                  const struct gravity_props *p);
#endif

#endif /* SWIFT_GRAVITY_PROPERTIES */
