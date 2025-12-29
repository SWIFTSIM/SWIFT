/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_STARS_H
#define SWIFT_STARS_H

/* Config parameters. */
#include <config.h>

/* Select the correct star model */
#if defined(STARS_NONE)
#include "./stars/None/stars.h"
#include "./stars/None/stars_iact.h"
#elif defined(STARS_BASIC)
#include "./stars/Basic/stars.h"
#include "./stars/Basic/stars_iact.h"
#elif defined(STARS_EAGLE)
#include "./stars/EAGLE/stars.h"
#include "./stars/EAGLE/stars_iact.h"
#elif defined(STARS_GEAR)
#include "./stars/GEAR/stars.h"
#include "./stars/GEAR/stars_iact.h"
#include "./stars/GEAR/stars_stellar_type.h"
#else
#error "Invalid choice of star model"
#endif

struct engine;
struct space;

void stars_exact_density_compute(struct space *s, const struct engine *e);
void stars_exact_density_check(struct space *s, const struct engine *e,
                               const double rel_tol);

/**
 * @brief Accumulate the displacement of a spart since the last tree rebuild
 * and cell sort.
 *
 * We need to update these offsets everytime we spawn a new spart.
 *
 * @param sp The #spart to update.
 * @param displacement The 3D displacement vector to add.
 */
__attribute__((always_inline))
INLINE static void spart_add_displacement(struct spart *restrict sp,
					  const float displacement[3]) {

  sp->x_diff[0] += displacement[0];
  sp->x_diff[1] += displacement[1];
  sp->x_diff[2] += displacement[2];

  sp->x_diff_sort[0] += displacement[0];
  sp->x_diff_sort[1] += displacement[1];
  sp->x_diff_sort[2] += displacement[2];
}

#endif /* SWIFT_STARS_H */
