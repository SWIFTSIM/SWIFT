/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Will Roper (w.roper@sussex.ac.uk)
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

/* Standard includes. */
#include <math.h>

/* Local includes. */
#include "cell.h"
#include "engine.h"
#include "proxy.h"
#include "space.h"
#include "zoom.h"

/**
 * @brief Calculate the distance at which we can truncate the parent volume.
 *
 * This uses a simple geometric argument (based on the (L/R)^3 tidal criterion)
 * to find the distance from the zoom region at which the contributions from the
 * background drop below the desired accuracy.
 *
 * @param zoom_dim The zoom region dimensions.
 * @param tidal_factor The tidal factor accounting for anisotropies in the
 *     background (>1, higher means more background preserved, i.e. more
 *     accurate).
 * @param epsilon The desired accuracy.
 * @return The truncation distance.
 */
double zoom_compute_bkg_truncate_dist(const double zoom_dim[3],
                                      const double tidal_factor,
                                      const double epsilon) {

  /* Get the maximum zoom dimension in case we are not cubic */
  double dim = fmax(zoom_dim[0], fmax(zoom_dim[1], zoom_dim[2]));

  return tidal_factor * dim / pow(epsilon, 1.0 / 3.0);
}

/**
 * @brief Truncate the simulation volume to remove distant background.
 *
 * This removes all cells that are further away from the zoom region than the
 * truncation distance computed with zoom_compute_bkg_truncate_dist.
 *
 * @param s The #space.
 * @param verbose Whether to be verbose or not.
 */
void zoom_truncate_background(struct space *s, const int verbose) {

  /* Extract some useful pointers and information. */
  struct cell *cells = s->cells_top;
  double *zoom_dim = s->zoom_props->zoom_dim;
  double tidal_factor = s->zoom_props->tidal_factor;
  double epsilon = s->zoom_props->truncate_epsilon;

  /* Compute the truncation distance. */
  const double r_trunc =
      zoom_compute_bkg_truncate_dist(zoom_dim, tidal_factor, epsilon);
  if (verbose)
    message(
        "Computed a truncation distance of %.2e (with %.2f x %.2e * "
        "(%.1e)^(1/3))",
        r_trunc, tidal_factor,
        fmax(zoom_dim[0], fmax(zoom_dim[1], zoom_dim[2])), epsilon);

  /* If the truncation distance exceeds the box size we can't truncate. */
  if (r_trunc * 2.0 >= fmin(s->dim[0], fmin(s->dim[1], s->dim[2]))) {
    error(
        "Truncation distance (%.2e) exceeds box size (%.2e), cannot truncate. "
        "You probably don't need truncation in this case, turn off "
        "ZoomRegion:truncate_background.",
        r_trunc, fmin(s->dim[0], fmin(s->dim[1], s->dim[2])));
    return;
  }
}
