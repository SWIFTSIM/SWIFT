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

/* Config */
#include <config.h>

/* Standard includes */
#include <math.h>

/* Local includes */
#include "engine.h"
#include "minmax.h"
#include "parser.h"
#include "part.h"
#include "periodic.h"
#include "space.h"
#include "timers.h"
#include "zoom.h"

/**
 * @brief Calculate the half-width of the box after truncation.
 *
 * Truncation removes all background beyond a cube of half-width R centred on
 * the high-resolution region and adopts the truncated cube as the new periodic
 * volume. The error this introduces in the tidal field across the protected
 * high-resolution region is bounded as described below assuming a uniform
 * distribution of background mass.
 *
 * The mean background density exerts no tidal field across the region (shell
 * theorem), so removing distant, statistically uniform background loses
 * nothing at zeroth order; masses therefore never enter the criterion. The
 * dominant residual error is the spurious tidal field from the periodic
 * images of the retained volume (whose internal structure is dominated by the
 * high-resolution region itself, mass M) at the new truncated box size D = 2 R.
 *
 *   - The nearest image pair at +/- D produces a differential (tidal)
 *     acceleration across the protected half-extent l of
 *         delta_a = 4 G M l / D^3.
 *   - Comparing with the region's self-gravity at its edge,
 *         a_self = G M / l^2,
 *     and requiring delta_a / a_self <= epsilon (M cancels) gives
 *         D >= l (4 / epsilon)^(1/3),
 *     i.e.
 *         R = D / 2 = l (2 epsilon)^(-1/3).
 *
 * The discreteness error from a single unbalanced removed background lump of
 * mass m at distance R is smaller than the image term by ~ m / M << 1; it,
 * the more distant diagonal images, and any anisotropy are absorbed into the
 * user-tunable safety factor (>= 1).
 *
 * Note that truncation also (deliberately) discards the real large-scale
 * tidal field sourced beyond R. That loss is a physical choice inherent to
 * truncating and is *not* controlled by epsilon due to the consequences of the
 * shell theorem. Note that these assumptions breakdown if the uniformity
 * assumption is violated.
 *
 * @param protected_half_extent The half-extent l of the protected region.
 * @param tidal_factor Safety factor (>=1) for anisotropies and the more
 *     distant periodic images. Higher means more background is kept.
 * @param epsilon The maximum tolerated fractional tidal error across the
 *     protected region, relative to its self-gravity.
 * @return The truncation half-width R.
 */
static double zoom_compute_bkg_truncate_half_width(
    const double protected_half_extent, const double tidal_factor,
    const double epsilon) {

  return tidal_factor * protected_half_extent * cbrt(0.5 / epsilon);
}

/**
 * @brief Truncate the simulation volume to remove distant background.
 *
 * This removes background particles outside a cube centred on the
 * high-resolution particle distribution and adopts that cube as the new
 * (smaller) periodic simulation volume. The half-width of the retained cube
 * is computed with zoom_compute_bkg_truncate_half_width().
 *
 * Note that this truncation process requires the assumption that the background
 * is statistically uniform and isotropic, so that the tidal field across the
 * high-resolution region is dominated by the nearest periodic images of the
 * retained volume.
 *
 * This must be called at the end of space_init(), before anything downstream
 * consumes the box dimensions (the gravity properties, the PM mesh, the
 * neutrino response, line-of-sight properties, etc. are all derived from
 * s->dim at start-up). Since space_init() only runs on a fresh start, the
 * box is never re-truncated on a restart or during regridding.
 *
 * @param params The swift parameters (needed to refresh the mesh-derived
 *     neighbour distance for the truncated volume).
 * @param s The #space.
 * @param verbose Whether to be verbose or not.
 */
void zoom_truncate_bkg(struct swift_params *params, struct space *s,
                       const int verbose) {

  const ticks tic = getticks();

  /* Nothing to do if truncation is off. */
  if (!s->zoom_props->truncate_background) return;

  /* Truncation redefines the periodic volume; it cannot be applied to a
   * non-periodic box. */
  if (!s->periodic) {
    error("ZoomRegion:truncate_background requires a periodic box.");
  }

  /* Extract some useful information. */
  const double tidal_factor = s->zoom_props->tidal_factor;
  const double epsilon = s->zoom_props->truncate_epsilon;
  const double old_dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  /* Measure the high-resolution particle extent and compute the shift that
   * centres it in the (current) box. */
  zoom_get_region_dim_and_shift(s, verbose);

  /* Protect twice the measured high-resolution extent, i.e. a protected
   * half-extent equal to the measured extent. This gives factor-2 headroom
   * for the high-resolution region to grow during the run. */
  const double high_res_dim =
      max3(s->zoom_props->part_dim[0], s->zoom_props->part_dim[1],
           s->zoom_props->part_dim[2]);
  const double protected_dim = 2.0 * high_res_dim;
  const double protected_half_extent = high_res_dim;

  /* Remember what we protected so we can warn later if the high-resolution
   * region outgrows it. */
  s->zoom_props->truncate_protected_dim = protected_dim;

  /* Compute the retained half-box width. */
  const double retained_half_width = zoom_compute_bkg_truncate_half_width(
      protected_half_extent, tidal_factor, epsilon);

  if (verbose) {
    message(
        "Computed a truncation half-width of %.2f internal units for a "
        "protected high-resolution extent of %.2f (with %.2f x %.2f x "
        "(2 x %.1e)^(-1/3))",
        retained_half_width, protected_dim, tidal_factor, protected_half_extent,
        epsilon);
  }

  /* The retained box must comfortably contain the protected region. */
  if (2.0 * retained_half_width <= protected_dim) {
    error(
        "The truncated box (%.2e) does not contain the protected "
        "high-resolution region (extent %.2e). Lower "
        "ZoomRegion:truncate_epsilon.",
        2.0 * retained_half_width, protected_dim);
  }

  /* If the retained box exceeds the original box we can't truncate. */
  if (retained_half_width * 2.0 >=
      fmin(old_dim[0], fmin(old_dim[1], old_dim[2]))) {
    error(
        "The truncated box (%.2e) exceeds the parent box (%.2e), cannot "
        "truncate. You probably don't need truncation in this case, turn off "
        "ZoomRegion:truncate_background.",
        retained_half_width * 2.0,
        fmin(old_dim[0], fmin(old_dim[1], old_dim[2])));
  }

  /* Add the shift that moves the high-resolution centre from the middle of
   * the parent box to the middle of the retained box. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->zoom_shift[i] += retained_half_width - old_dim[i] / 2.0;
  }

  /* Apply the combined (centring + truncation) shift to all particles. */
  zoom_apply_zoom_shift_to_particles(s, verbose);

  /* Wrap the shifted positions back into the parent box. In the shifted
   * frame the retained cube is exactly [0, 2R) along each axis, so after
   * wrapping a particle is retained if and only if all its coordinates are
   * below 2R. */
  for (size_t k = 0; k < s->nr_parts; k++) {
    for (int i = 0; i < 3; i++) {
      s->parts[k].x[i] = box_wrap_multiple(s->parts[k].x[i], 0.0, old_dim[i]);
    }
  }
  for (size_t k = 0; k < s->nr_gparts; k++) {
    for (int i = 0; i < 3; i++) {
      s->gparts[k].x[i] = box_wrap_multiple(s->gparts[k].x[i], 0.0, old_dim[i]);
    }
  }
  for (size_t k = 0; k < s->nr_sparts; k++) {
    for (int i = 0; i < 3; i++) {
      s->sparts[k].x[i] = box_wrap_multiple(s->sparts[k].x[i], 0.0, old_dim[i]);
    }
  }
  for (size_t k = 0; k < s->nr_bparts; k++) {
    for (int i = 0; i < 3; i++) {
      s->bparts[k].x[i] = box_wrap_multiple(s->bparts[k].x[i], 0.0, old_dim[i]);
    }
  }
  for (size_t k = 0; k < s->nr_sinks; k++) {
    for (int i = 0; i < 3; i++) {
      s->sinks[k].x[i] = box_wrap_multiple(s->sinks[k].x[i], 0.0, old_dim[i]);
    }
  }

  /* Loop over all the gparts and inhibit background particles outside the
   * retained cube. */
  size_t ntrunc = 0;
  size_t nbkg = 0;
  for (size_t k = 0; k < s->nr_gparts; k++) {

    /* Skip non-background particles. */
    if (s->gparts[k].type != swift_type_dark_matter_background) {
      continue;
    }
    nbkg++;

    /* Inhibit background particles that are too far away. */
    for (int i = 0; i < 3; i++) {
      if (s->gparts[k].x[i] >= 2.0 * retained_half_width) {
        s->gparts[k].time_bin = time_bin_inhibited;
        s->nr_inhibited_gparts++;
        ntrunc++;
        break;
      }
    }
  }

  /* Adopt the truncated volume as the simulation box. */
  for (int i = 0; i < 3; i++) {
    s->dim[i] = 2.0 * retained_half_width;
  }

  /* The gravity mesh scale is derived from the box size, so refresh the
   * neighbour distance for the truncated volume (same calculation as in
   * zoom_props_init, which ran before the box was resized). */
  const int mesh_size =
      parser_get_param_int(params, "Gravity:mesh_side_length");
  const float a_smooth =
      parser_get_opt_param_float(params, "Gravity:a_smooth", 1.25);
  const float r_cut_max_ratio =
      parser_get_opt_param_float(params, "Gravity:r_cut_max", 4.5);
  const float r_s = a_smooth * s->dim[0] / mesh_size;
  s->zoom_props->neighbour_distance = r_s * r_cut_max_ratio;

  if (verbose) {
    message("Removing %zu background particles out of %zu.", ntrunc, nbkg);
    message("Truncated box dimensions are [%.2f, %.2f, %.2f]", s->dim[0],
            s->dim[1], s->dim[2]);
    message("Truncating the box took %f %s",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
  }
}
