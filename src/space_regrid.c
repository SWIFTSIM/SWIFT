/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* Config parameters. */
#include <config.h>

/* This object's header. */
#include "space.h"

/* Local headers. */
#include "cell.h"
#include "engine.h"
#include "scheduler.h"

/**
 * @brief Upper bound applied to a star's radiation term (the larger of its
 * HII search radius and its ordinary smoothing length) before it
 * contributes to the top-level grid sizing.
 *
 * Exposed (rather than kept local to space_regrid()) so that the rebuild
 * criterion's own regression test can ask the real function what grid a
 * given radiation term produces, instead of re-deriving the expression and
 * drifting from it.
 *
 * @param cell_max_width The coarsest top-level cell width the run allows.
 */
float space_regrid_radiation_cap_for(const double cell_max_width) {

  return 0.99f * (float)cell_max_width /
         (kernel_gamma * space_stretch * radiation_search_radius_factor);
}

/**
 * @brief A star's contribution to the top-level grid's radiation term:
 * the larger of its HII search radius and its own smoothing length,
 * capped by space_regrid_radiation_cap_for().
 *
 * Exposed so that every accumulation site in space_regrid() below, and
 * the rebuild criterion's own regression test, share this one routing
 * decision instead of each re-deriving max(h_hii_max, h_max) locally and
 * risking one site silently drifting from the others. This function feeds
 * ONLY h_max_radiation (for cell_need_rebuild_for_radiation_pair()); the
 * same raw h_max is separately and additionally folded into h_max_no_hii
 * by the caller (for cell_need_rebuild_for_stars_pair(), which bounds a
 * star's reach without this function's cap or factor) -- the two
 * accumulators serve two different criteria, so a star's h_max belongs in
 * both, not routed through only one of them.
 *
 * @param h_hii_max A star's (or a cell's max over its stars') HII search
 * radius.
 * @param h_max A star's (or a cell's max over its stars') ordinary
 * smoothing length.
 * @param cell_max_width The coarsest top-level cell width the run allows.
 */
float space_regrid_star_radiation_term_for(const float h_hii_max,
                                           const float h_max,
                                           const double cell_max_width) {

  const float star_reach = fmaxf(h_hii_max, h_max);
  return fminf(star_reach, space_regrid_radiation_cap_for(cell_max_width));
}

/**
 * @brief A star's contribution to the top-level grid's UNCAPPED,
 * unfactored h_max_no_hii accumulator: its own smoothing length, with no
 * cap and no radiation_search_radius_factor, folded into the running
 * maximum so far.
 *
 * Exposed alongside space_regrid_star_radiation_term_for() so the
 * regression test can exercise both of a star's routing decisions
 * directly instead of only the capped/factored one; this is the
 * counterpart that keeps cell_need_rebuild_for_stars_pair() (unfactored,
 * unexempted) correctly served.
 *
 * @param h_max_no_hii_so_far The running accumulator before this star.
 * @param h_max A star's (or a cell's max over its stars') ordinary
 * smoothing length.
 */
float space_regrid_star_h_max_no_hii_term_for(const float h_max_no_hii_so_far,
                                              const float h_max) {

  return fmaxf(h_max_no_hii_so_far, h_max);
}

/**
 * @brief Search radius the top-level grid must be sized for.
 *
 * @param h_max_no_hii Largest ordinary hydro/black-hole/sink smoothing
 * length (already passed through space::h_max_no_hii_hwm's ratchet).
 * @param h_max_radiation Largest star radiation term (max of h_hii and the
 * star's own smoothing length, per star, then maxed over stars), already
 * capped by space_regrid_radiation_cap_for().
 */
double space_regrid_search_radius_for(const float h_max_no_hii,
                                      const float h_max_radiation) {

  /* The radiation runner (runner_do_stars_hii_ionization_feedback()), the
     rebuild criterion (cell_need_rebuild_for_radiation_pair()) and the split
     predicates (cell_can_split_{pair,self}_radiation_subgrid_task()) all
     size a star's radiation reach as radiation_search_radius_factor *
     kernel_gamma * max(h_hii_max, h_max) -- the split predicates express
     this as two independent bounds on h_hii_max and h_max rather than a
     literal max(), but bounding both terms against the same threshold is
     equivalent to bounding their max. The top-level grid must be sized the
     same way. Sizing it on the unfactored kernel_gamma * h_max_radiation
     instead let cell_need_rebuild_for_radiation_pair() demand a coarsening
     the regrid would not perform, which engine_rebuild() reports as
     "engine_unskip failed after a rebuild!". The factor applies ONLY to the
     star radiation term: hydro/black-hole/sink smoothing lengths keep their
     historical sizing. */
  const float radius_no_hii = h_max_no_hii * kernel_gamma * space_stretch;
  const float radius_radiation = h_max_radiation * kernel_gamma *
                                 space_stretch * radiation_search_radius_factor;
  return fmaxf(radius_no_hii, radius_radiation);
}

/**
 * @brief Effective top-level cell width for a given pair of maxima.
 *
 * @param h_max_no_hii As in space_regrid_search_radius_for().
 * @param h_max_radiation As in space_regrid_search_radius_for().
 * @param cell_min The finest top-level cell width the run allows.
 * @param cell_max_width The coarsest top-level cell width the run allows.
 */
double space_regrid_cell_width_for(const float h_max_no_hii,
                                   const float h_max_radiation,
                                   const double cell_min,
                                   const double cell_max_width) {

  const double search_radius =
      space_regrid_search_radius_for(h_max_no_hii, h_max_radiation);
  return fmin(fmax(search_radius, cell_min), cell_max_width);
}

/**
 * @brief Re-build the top-level cell grid.
 *
 * @param s The #space.
 * @param verbose Print messages to stdout or not.
 */
void space_regrid(struct space *s, int verbose) {

  const size_t nr_parts = s->nr_parts;
  const size_t nr_sparts = s->nr_sparts;
  const size_t nr_bparts = s->nr_bparts;
  const size_t nr_sinks = s->nr_sinks;
  const ticks tic = getticks();
  const integertime_t ti_current = (s->e != NULL) ? s->e->ti_current : 0;

  /* Run through the cells and get the current h_max, tracked separately for
     ordinary hydro/black-hole/sink smoothing lengths (h_max_no_hii) and for
     the star radiation term -- max(h_hii, h) per star, then maxed over
     stars (h_max_radiation). They are combined further down via
     space::h_max_no_hii_hwm's ratchet -- see the comment there for why the
     star radiation term must stay independent of the ordinary quantities
     instead of folding into one combined running value.

     A star's ordinary smoothing length h feeds BOTH accumulators, not
     either/or: h_max_no_hii, alongside hydro/black_holes/sinks, because
     cell_need_rebuild_for_stars_pair() bounds a star-hydro pair by the
     unfactored kernel_gamma * max(stars.h_max, hydro.h_max) with no
     exemption, so a star's h must keep driving the same uncapped
     ratchet/abort those quantities always have; and h_max_radiation, via
     space_regrid_star_radiation_term_for(), because
     cell_need_rebuild_for_radiation_pair() separately bounds a star's
     reach by radiation_search_radius_factor * kernel_gamma *
     max(h_hii_max, h_max), not h_hii_max alone -- a star resting at its
     own h with no active HII region is the same growing-search-radius
     hazard h_hii was fixed for, in a sparse top-level cell where h binds.
     Routing stars.h_max through only the capped/factored term (dropping it
     from h_max_no_hii) silently caps the grid's response to
     cell_need_rebuild_for_stars_pair's own uncapped, unexempted bound,
     reproducing the "engine_unskip failed after a rebuild!" deadlock this
     fix line exists to close, via a criterion it was never meant to
     touch. */
  // tic = getticks();
  const float h_max_floor = s->cell_min / kernel_gamma / space_stretch;
  float h_max_no_hii = h_max_floor;
  /* Starts at zero, not at h_max_floor: the radiation term carries an extra
     radiation_search_radius_factor (see space_regrid_search_radius_for()),
     so seeding it with the floor would size a star-free grid as
     radiation_search_radius_factor * cell_min and silently coarsen every
     run in the tree. The floor is preserved by the cell_min clamp below,
     which is where it belongs. */
  float h_max_radiation = 0.f;

  /* A star's radiation term (max(h_hii, h)) is allowed to grow far beyond
     the box's characteristic cell size. Once cdim is already at the
     coarsest grid the run allows (Scheduler:min_top_level_cells, at least 3
     when periodic), the 27-cell radiation stencil wraps around and already
     covers the entire box, so the radiation term cannot miss any gas no
     matter how large it grows from there. Cap its own contribution to
     h_max at the value that would size cells down to exactly that minimum:
     it can still widen cells toward full coverage at finer configurations,
     just never push cdim to the floor for a quantity that buys nothing once
     coverage is already total. Ordinary hydro/black-hole/sink h_max
     growth, and a star's own h_max as it separately feeds h_max_no_hii
     below, is NOT capped -- that remains a genuine problem worth the "too
     few cells" error below. Only the copy of stars.h_max routed through
     space_regrid_star_radiation_term_for() into h_max_radiation is
     capped; the same value's other copy in h_max_no_hii is not. */
  if (nr_parts > 0) {

    /* Can we use the list of local non-empty top-level cells? */
    if (s->local_cells_with_particles_top != NULL) {
      for (int k = 0; k < s->nr_local_cells_with_particles; ++k) {
        const struct cell *c =
            &s->cells_top[s->local_cells_with_particles_top[k]];
        if (c->hydro.h_max > h_max_no_hii) {
          h_max_no_hii = c->hydro.h_max;
        }
        h_max_no_hii = space_regrid_star_h_max_no_hii_term_for(h_max_no_hii,
                                                               c->stars.h_max);
        const float star_radiation_term = space_regrid_star_radiation_term_for(
            c->stars.h_hii_max, c->stars.h_max, s->cell_max_width);
        if (star_radiation_term > h_max_radiation) {
          h_max_radiation = star_radiation_term;
        }
        if (c->black_holes.h_max > h_max_no_hii) {
          h_max_no_hii = c->black_holes.h_max;
        }
        if (c->sinks.h_max > h_max_no_hii) {
          h_max_no_hii = c->sinks.h_max;
        }
      }

      /* Can we instead use all the top-level cells? */
    } else if (s->cells_top != NULL) {
      for (int k = 0; k < s->nr_cells; k++) {
        const struct cell *c = &s->cells_top[k];
        if (c->nodeID == engine_rank && c->hydro.h_max > h_max_no_hii) {
          h_max_no_hii = c->hydro.h_max;
        }
        if (c->nodeID == engine_rank) {
          h_max_no_hii = space_regrid_star_h_max_no_hii_term_for(
              h_max_no_hii, c->stars.h_max);
          const float star_radiation_term =
              space_regrid_star_radiation_term_for(
                  c->stars.h_hii_max, c->stars.h_max, s->cell_max_width);
          if (star_radiation_term > h_max_radiation) {
            h_max_radiation = star_radiation_term;
          }
        }
        if (c->nodeID == engine_rank && c->black_holes.h_max > h_max_no_hii) {
          h_max_no_hii = c->black_holes.h_max;
        }
        if (c->nodeID == engine_rank && c->sinks.h_max > h_max_no_hii) {
          h_max_no_hii = c->sinks.h_max;
        }
      }

      /* Last option: run through the particles */
    } else {
      for (size_t k = 0; k < nr_parts; k++) {
        if (s->parts[k].h > h_max_no_hii) h_max_no_hii = s->parts[k].h;
      }
      for (size_t k = 0; k < nr_sparts; k++) {
        h_max_no_hii = space_regrid_star_h_max_no_hii_term_for(h_max_no_hii,
                                                               s->sparts[k].h);
        const float star_radiation_term = space_regrid_star_radiation_term_for(
            s->sparts[k].h_hii, s->sparts[k].h, s->cell_max_width);
        if (star_radiation_term > h_max_radiation)
          h_max_radiation = star_radiation_term;
      }
      for (size_t k = 0; k < nr_bparts; k++) {
        if (s->bparts[k].h > h_max_no_hii) h_max_no_hii = s->bparts[k].h;
      }
      for (size_t k = 0; k < nr_sinks; k++) {
        if (s->sinks[k].h > h_max_no_hii) h_max_no_hii = s->sinks[k].h;
      }
    }
  }

/* If we are running in parallel, make sure everybody agrees on
   how large the largest cell should be. Reduce both quantities
   independently -- max() distributes over independent reductions, so this
   gives the same global h_max_no_hii/h_max_radiation as reducing one
   combined value would, while keeping them separable for the ratchet
   below. */
#ifdef WITH_MPI
  {
    float buff_in[2] = {h_max_no_hii, h_max_radiation};
    float buff_out[2];
    if (MPI_Allreduce(buff_in, buff_out, 2, MPI_FLOAT, MPI_MAX,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate the rebuild flag across nodes.");
    h_max_no_hii = buff_out[0];
    h_max_radiation = buff_out[1];
  }
#endif

  /* space::h_max_no_hii_hwm is a one-way ratchet that exactly reproduces
     this function's historical behaviour for ordinary hydro/black-hole/
     sink/star smoothing lengths: once the grid has coarsened for one of
     these (including a star's own h_max, folded into this accumulator
     above alongside hydro/black_holes/sinks), it stays coarse even if
     that h later shrinks (previously an implicit effect of the
     cdim-can-only-shrink trigger below; now made explicit so it can be
     kept separate from the star radiation term). This is what keeps
     cell_need_rebuild_for_stars_pair() correctly served: that criterion's
     own bound on stars.h_max is unfactored and unexempted, so the grid
     must never un-coarsen for it, same as for hydro/black_holes/sinks.

     The star radiation term (h_max_radiation) is a SEPARATE accumulator,
     partly built from the same stars.h_max values, and is deliberately
     NOT folded into this ratchet: h_hii is a transient, per-star quantity
     (reset to 0 when a star dies or ages out of HII eligibility, see
     feedback_common.c) that should stop requiring a coarse grid the
     moment no star's search radius needs it any more, and
     cell_need_rebuild_for_radiation_pair() -- the only criterion
     h_max_radiation serves -- is itself exempted once
     space_radiation_top_stencil_covers_box() holds, so un-coarsening this
     term cannot reopen a demand that criterion still enforces. A star's
     own h_max is therefore ratcheted once (via h_max_no_hii, for
     cell_need_rebuild_for_stars_pair) and separately re-evaluated fresh
     every regrid (via h_max_radiation, for
     cell_need_rebuild_for_radiation_pair): two accumulators serving two
     different criteria, not one unified quantity. */
  s->h_max_no_hii_hwm = fmaxf(s->h_max_no_hii_hwm, h_max_no_hii);
  const float h_max = fmaxf(s->h_max_no_hii_hwm, h_max_radiation);

  if (verbose) message("h_max is %.3e (cell_min=%.3e).", h_max, s->cell_min);

  /* Get the new putative cell dimensions. The effective cell width is
     clamped between cell_min (Scheduler:max_top_level_cells' complement,
     existing) and cell_max_width (Scheduler:min_top_level_cells' own
     complement): once h_max would push the width past cell_max_width, it
     is capped there instead, so cdim never drops below
     Scheduler:min_top_level_cells regardless of what is driving the
     coarsening. */
  const double search_radius =
      space_regrid_search_radius_for(s->h_max_no_hii_hwm, h_max_radiation);
  const double cell_width = space_regrid_cell_width_for(
      s->h_max_no_hii_hwm, h_max_radiation, s->cell_min, s->cell_max_width);
  const int cdim[3] = {(int)floor(s->dim[0] / cell_width),
                       (int)floor(s->dim[1] / cell_width),
                       (int)floor(s->dim[2] / cell_width)};

  /* check that we have at least 1 cell in each dimension */
  if (cdim[0] == 0 || cdim[1] == 0 || cdim[2] == 0) {
    error(
        "Top level cell dimension of size 0 detected (cdim = [%i %i "
        "%i])!\nThis usually indicates a problem with the initial smoothing "
        "lengths of the particles, e.g. a smoothing length that is comparable "
        "in size to the box size.",
        cdim[0], cdim[1], cdim[2]);
  }

  /* Check if we have enough cells for periodicity. The clamp above bounds
     the grid at cell_max_width, so cdim can no longer report this: test the
     search radius directly instead. The h_max > h_max_floor term keeps the
     test on real particle reach. h_max is floored at cell_min, which on a
     non-cubic box can exceed cell_max_width (cell_min derives from dmax,
     cell_max_width from dmin), and that floor alone must not be reported as
     a smoothing length too large for the box. Comparing against the float
     that seeded h_max is exact, where comparing the reconstructed width
     against cell_min is not. */
  if (s->periodic && h_max > h_max_floor && search_radius > s->cell_max_width)
    error(
        "Must have at least Scheduler:min_top_level_cells cells in each "
        "spatial dimension when periodicity is switched on (h_max = %g, from "
        "an ordinary smoothing length of %g and a star radiation term "
        "(max of h_hii and a star's own h) of %g, gives a search radius of "
        "%g, but the coarsest top-level cell allowed is %g).\nThe bound is "
        "the shortest box side divided by Scheduler:min_top_level_cells, so "
        "raising Scheduler:max_top_level_cells does not help. This error is "
        "often caused by any of the followings:\n"
        " - too few particles to generate a sensible grid,\n"
        " - 'Scheduler:min_top_level_cells' is too large (it cannot go below "
        "3 when periodicity is switched on, in which case the box itself is "
        "too small for these smoothing lengths),\n"
        " - the (minimal) time-step is too large leading to particles with "
        "predicted smoothing lengths too large for the box size,\n"
        " - particles with velocities so large that they move by more than two "
        "box sizes per time-step.\n",
        h_max, s->h_max_no_hii_hwm, h_max_radiation, search_radius,
        s->cell_max_width);

/* In MPI-Land, changing the top-level cell size requires that the
 * global partition is recomputed and the particles redistributed.
 * Be prepared to do that. */
#ifdef WITH_MPI
  double oldwidth[3] = {0., 0., 0.};
  double oldcdim[3] = {0., 0., 0.};
  int *oldnodeIDs = NULL;
  if (s->cells_top != NULL && (cdim[0] != s->cdim[0] || cdim[1] != s->cdim[1] ||
                               cdim[2] != s->cdim[2])) {

    /* Capture state of current space. h_max_no_hii_hwm's ratchet (see
       above) means cdim can now also come out larger (finer) than before
       -- only ever via the star radiation term (h_hii or a star's own h)
       shrinking, never via ordinary hydro/black-hole/sink physics -- so
       this capture must trigger on any change, not just a decrease. It
       still requires an existing grid: on the first call (from
       space_init) and on restart, s->cells_top is NULL while s->cdim may
       already differ, and both the capture loop and the repartition it
       triggers need live cells and a live engine. Those cases go through
       the no_regrid path below instead. */
    oldcdim[0] = s->cdim[0];
    oldcdim[1] = s->cdim[1];
    oldcdim[2] = s->cdim[2];
    oldwidth[0] = s->width[0];
    oldwidth[1] = s->width[1];
    oldwidth[2] = s->width[2];

    if ((oldnodeIDs =
             (int *)swift_malloc("nodeIDs", sizeof(int) * s->nr_cells)) == NULL)
      error("Failed to allocate temporary nodeIDs.");

    int cid = 0;
    for (int i = 0; i < s->cdim[0]; i++) {
      for (int j = 0; j < s->cdim[1]; j++) {
        for (int k = 0; k < s->cdim[2]; k++) {
          cid = cell_getid(oldcdim, i, j, k);
          oldnodeIDs[cid] = s->cells_top[cid].nodeID;
        }
      }
    }
  }

  /* Are we about to allocate new top level cells without a regrid?
   * Can happen when restarting the application. */
  const int no_regrid = (s->cells_top == NULL && oldnodeIDs == NULL);
#endif

  /* Do we need to re-build the upper-level cells? Any change now
     triggers a regrid, not just a decrease: h_max_no_hii_hwm's ratchet
     (above) guarantees cdim can only come out larger (finer) than
     s->cdim because the star radiation term's own contribution shrank,
     never because ordinary hydro/black-hole/sink h_max did -- so allowing
     that direction here does not reopen the general one-way-coarsening
     behaviour those quantities still rely on. */
  // tic = getticks();
  if (s->cells_top == NULL || cdim[0] != s->cdim[0] || cdim[1] != s->cdim[1] ||
      cdim[2] != s->cdim[2]) {

/* Be verbose about this. */
#ifdef SWIFT_DEBUG_CHECKS
    message("(re)griding space cdim=(%d %d %d)", cdim[0], cdim[1], cdim[2]);
    fflush(stdout);
#endif

    /* Free the old cells, if they were allocated. */
    if (s->cells_top != NULL) {
      space_free_cells(s);
      swift_free("local_cells_with_tasks_top", s->local_cells_with_tasks_top);
      swift_free("local_cells_top", s->local_cells_top);
      swift_free("cells_with_particles_top", s->cells_with_particles_top);
      swift_free("local_cells_with_particles_top",
                 s->local_cells_with_particles_top);
      swift_free("cells_top", s->cells_top);
      swift_free("cells_top_updated", s->cells_top_updated);
      swift_free("multipoles_top", s->multipoles_top);
    }

    /* Also free the task arrays, these will be regenerated and we can use the
     * memory while copying the particle arrays. */
    if (s->e != NULL) scheduler_free_tasks(&s->e->sched);

    /* Set the new cell dimensions. */
    for (int k = 0; k < 3; k++) {
      s->cdim[k] = cdim[k];
      s->width[k] = s->dim[k] / cdim[k];
      s->iwidth[k] = 1.0 / s->width[k];
    }
    const float dmin = min3(s->width[0], s->width[1], s->width[2]);

    /* Allocate the highest level of cells. */
    s->tot_cells = s->nr_cells = cdim[0] * cdim[1] * cdim[2];

    if (swift_memalign("cells_top", (void **)&s->cells_top, cell_align,
                       s->nr_cells * sizeof(struct cell)) != 0)
      error("Failed to allocate top-level cells.");
    bzero(s->cells_top, s->nr_cells * sizeof(struct cell));

    /* Allocate the multipoles for the top-level cells. */
    if (s->with_self_gravity) {
      if (swift_memalign("multipoles_top", (void **)&s->multipoles_top,
                         multipole_align,
                         s->nr_cells * sizeof(struct gravity_tensors)) != 0)
        error("Failed to allocate top-level multipoles.");
      bzero(s->multipoles_top, s->nr_cells * sizeof(struct gravity_tensors));
    }

    if (swift_memalign("cells_top_updated", (void **)&s->cells_top_updated,
                       cell_align, s->nr_cells * sizeof(char)) != 0)
      error("Failed to allocate top-level cells.");
    bzero(s->cells_top_updated, s->nr_cells * sizeof(char));

    /* Allocate the indices of local cells */
    if (swift_memalign("local_cells_top", (void **)&s->local_cells_top,
                       SWIFT_STRUCT_ALIGNMENT, s->nr_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of local top-level cells.");
    bzero(s->local_cells_top, s->nr_cells * sizeof(int));

    /* Allocate the indices of local cells with tasks */
    if (swift_memalign("local_cells_with_tasks_top",
                       (void **)&s->local_cells_with_tasks_top,
                       SWIFT_STRUCT_ALIGNMENT, s->nr_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of local top-level cells with tasks.");
    bzero(s->local_cells_with_tasks_top, s->nr_cells * sizeof(int));

    /* Allocate the indices of cells with particles */
    if (swift_memalign("cells_with_particles_top",
                       (void **)&s->cells_with_particles_top,
                       SWIFT_STRUCT_ALIGNMENT, s->nr_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of top-level cells with particles.");
    bzero(s->cells_with_particles_top, s->nr_cells * sizeof(int));

    /* Allocate the indices of local cells with particles */
    if (swift_memalign("local_cells_with_particles_top",
                       (void **)&s->local_cells_with_particles_top,
                       SWIFT_STRUCT_ALIGNMENT, s->nr_cells * sizeof(int)) != 0)
      error(
          "Failed to allocate indices of local top-level cells with "
          "particles.");
    bzero(s->local_cells_with_particles_top, s->nr_cells * sizeof(int));

    /* Set the cells' locks */
    for (int k = 0; k < s->nr_cells; k++) {
      if (lock_init(&s->cells_top[k].hydro.lock) != 0)
        error("Failed to init spinlock for hydro.");
      if (lock_init(&s->cells_top[k].hydro.extra_sort_lock) != 0)
        error("Failed to init spinlock for hydro extra sort.");
      if (lock_init(&s->cells_top[k].grav.plock) != 0)
        error("Failed to init spinlock for gravity.");
      if (lock_init(&s->cells_top[k].grav.mlock) != 0)
        error("Failed to init spinlock for multipoles.");
      if (lock_init(&s->cells_top[k].grav.star_formation_lock) != 0)
        error("Failed to init spinlock for star formation (gpart).");
      if (lock_init(&s->cells_top[k].stars.lock) != 0)
        error("Failed to init spinlock for stars.");
      if (lock_init(&s->cells_top[k].sinks.lock) != 0)
        error("Failed to init spinlock for sinks.");
      if (lock_init(&s->cells_top[k].sinks.sink_formation_lock) != 0)
        error("Failed to init spinlock for sink formation.");
      if (lock_init(&s->cells_top[k].black_holes.lock) != 0)
        error("Failed to init spinlock for black holes.");
      if (lock_init(&s->cells_top[k].stars.star_formation_lock) != 0)
        error("Failed to init spinlock for star formation (spart).");
    }

    /* Set the cell location and sizes. */
    for (int i = 0; i < cdim[0]; i++)
      for (int j = 0; j < cdim[1]; j++)
        for (int k = 0; k < cdim[2]; k++) {
          const size_t cid = cell_getid(cdim, i, j, k);
          struct cell *restrict c = &s->cells_top[cid];
          c->loc[0] = i * s->width[0];
          c->loc[1] = j * s->width[1];
          c->loc[2] = k * s->width[2];
          c->width[0] = s->width[0];
          c->width[1] = s->width[1];
          c->width[2] = s->width[2];
          c->dmin = dmin;
          c->h_min_allowed = c->dmin * 0.5 * (1. / kernel_gamma);
          c->h_max_allowed = c->dmin * (1. / kernel_gamma);
          c->depth = 0;
          c->split = 0;
          c->hydro.count = 0;
          c->grav.count = 0;
          c->stars.count = 0;
          c->sinks.count = 0;
          c->top = c;
          c->super = c;
          c->hydro.super = c;
          c->grav.super = c;
          c->hydro.ti_old_part = ti_current;
          c->grav.ti_old_part = ti_current;
          c->stars.ti_old_part = ti_current;
          c->sinks.ti_old_part = ti_current;
          c->black_holes.ti_old_part = ti_current;
          c->grav.ti_old_multipole = ti_current;
#ifdef WITH_MPI
          c->mpi.tag = -1;
          c->mpi.recv = NULL;
          c->mpi.send = NULL;
#endif  // WITH_MPI
          if (s->with_self_gravity) c->grav.multipole = &s->multipoles_top[cid];
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
          cell_assign_top_level_cell_index(c, s->cdim, s->dim, s->iwidth);
#endif
        }

    /* Be verbose about the change. */
    if (verbose)
      message("set cell dimensions to [ %i %i %i ].", cdim[0], cdim[1],
              cdim[2]);

#ifdef WITH_MPI
    if (oldnodeIDs != NULL) {
      /* We have changed the top-level cell dimension, so need to redistribute
       * cells around the nodes. We repartition using the old space node
       * positions as a grid to resample. */
      if (s->e->nodeID == 0)
        message(
            "basic cell dimensions have increased - recalculating the "
            "global partition.");

      if (!partition_space_to_space(oldwidth, oldcdim, oldnodeIDs, s)) {

        /* Failed, try another technique that requires no settings. */
        message("Failed to get a new partition, trying less optimal method");
        struct partition initial_partition;
#if defined(HAVE_PARMETIS) || defined(HAVE_METIS)
        initial_partition.type = INITPART_METIS_NOWEIGHT;
#else
        initial_partition.type = INITPART_VECTORIZE;
#endif
        partition_initial_partition(&initial_partition, s->e->nodeID,
                                    s->e->nr_nodes, s);
      }

      /* Re-distribute the particles to their new nodes. */
      engine_redistribute(s->e);

      /* Make the proxies. */
      engine_makeproxies(s->e);

      /* Finished with these. */
      swift_free("nodeIDs", oldnodeIDs);

    } else if (no_regrid && s->e != NULL) {
      /* If we have created the top-levels cells and not done an initial
       * partition (can happen when restarting), then the top-level cells
       * are not assigned to a node, we must do that and then associate the
       * particles with the cells. Note requires that
       * partition_store_celllist() was called once before, or just before
       * dumping the restart files.*/
      partition_restore_celllist(s, s->e->reparttype);

      /* Now re-distribute the particles, should just add to cells? */
      engine_redistribute(s->e);

      /* Make the proxies. */
      engine_makeproxies(s->e);
    }
#endif /* WITH_MPI */

    // message( "rebuilding upper-level cells took %.3f %s." ,
    // clocks_from_ticks(double)(getticks() - tic), clocks_getunit());

  } /* re-build upper-level cells? */
  else { /* Otherwise, just clean up the cells. */

    /* Free the old cells, if they were allocated. */
    space_free_cells(s);
  }

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
