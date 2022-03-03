/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

#ifndef SWIFT_PROJECTED_KERNEL_H
#define SWIFT_PROJECTED_KERNEL_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "error.h"
#include "inline.h"

#define PROJECTED_KERNEL_NTAB 1000

struct projected_kernel_table {
  int n;
  double du;
  double inv_du;
  double u_max;
  double *value;
};

/**
 * @brief Computes 2D projection of the 3D kernel function.
 *
 * This version interpolates the value from the supplied
 * look up table.
 *
 * @param u The ratio of the (2D) distance to the smoothing length
 */
__attribute__((always_inline)) INLINE static double projected_kernel_eval(
    struct projected_kernel_table *tab, double u) {

  /* Check u is in range */
  if (u >= tab->u_max) return 0.0;
  if (u < 0.0) error("Negative u in projected kernel!");

  /* Determine which interval we're in */
  int i = u * tab->inv_du;

  /* Find where we are in the interval */
  double f = (u - i * tab->du) * tab->inv_du;

  /* Linear interpolation */
  return (1.0 - f) * tab->value[i] + f * tab->value[i + 1];
}

void projected_kernel_init(struct projected_kernel_table *tab);
void projected_kernel_clean(struct projected_kernel_table *tab);
void projected_kernel_dump(void);

#endif
