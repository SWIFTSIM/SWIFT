/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Will Roper (w.roper@sussex.ac.uk)
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

/* Some standard headers. */
#include <string.h>

/* Local headers. */
#include "active.h"
#include "atomic.h"
#include "cell.h"
#include "engine.h"
#include "gravity_properties.h"
#include "runner.h"
#include "runner_doiact_grav.h"
#include "space.h"
#include "timers.h"

/**
 * @brief Performs M-M interactions between a given top-level cell and
 *        all other top level cells not interacted with via pair tasks.
 *
 * This is the non-periodic case where there is no mesh so all cells not
 * handled by a pair task are interacted with here in this long range
 * gravity function.
 *
 * This function is used when running a non-periodic uniform box (i.e. non-zoom
 * simulation) and just loops over every other cell with particles and interacts
 * with them.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_uniform_non_periodic(struct runner *r,
                                                    struct cell *ci,
                                                    struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;

  /* Get the multipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* Recover the list of top-level cells */
  struct cell *cells = e->s->cells_top;
  int *cells_with_particles = e->s->cells_with_particles_top;
  const int nr_cells_with_particles = e->s->nr_cells_with_particles;

  /* Loop over all the top-level cells and go for a M-M interaction if
   * well-separated */
  for (int n = 0; n < nr_cells_with_particles; ++n) {

    /* Handle on the top-level cell and it's gravity business*/
    struct cell *cj = &cells[cells_with_particles[n]];
    struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Avoid self contributions */
    if (top == cj) continue;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

    if (cell_can_use_pair_mm(top, cj, e, e->s, /*use_rebuild_data=*/1,
                             /*is_tree_walk=*/0,
                             /*periodic boundaries*/ s->periodic,
                             /*use_mesh*/ s->periodic)) {

      /* Call the PM interaction function on the active sub-cells of ci */
      runner_dopair_grav_mm_nonsym(r, ci, cj);
      // runner_dopair_recursive_grav_pm(r, ci, cj);

      /* Record that this multipole received a contribution */
      multi_i->pot.interacted = 1;

    } /* We are in charge of this pair */
  } /* Loop over top-level cells */
}

/**
 * @brief Performs M-M interactions between a given top-level cell and
 *        all other top level cells not interacted with via pair tasks.
 *
 * This is the non-periodic case where there is no mesh so all cells not
 * handled by a pair task are interacted with here in this long range
 * gravity function.
 *
 * This function is used when running a non-periodic zoom simulation and
 * will loop over all non-zoom cells but use the void cell hierarchy to
 * interact with all zoom cells.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_zoom_non_periodic(struct runner *r,
                                                 struct cell *ci,
                                                 struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;

  /* Get the multipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* Recover the list of top-level cells */
  struct cell *bkg_cells = e->s->zoom_props->bkg_cells_top;

#ifdef SWIFT_DEBUG_CHECKS
  /* Define counters used to count gparts. */
  size_t tested_gparts = 0;
#endif

  /* Since the zoom cells will be handled by the void cell hierarchy we can
   * just loop over all other cells which are not zoom cells. This is
   * trivial since the zoom cells are first in cells_top. */
  for (int cjd = 0; cjd < s->zoom_props->nr_bkg_cells; cjd++) {

    /* Handle on the top-level cell and it's gravity business*/
    struct cell *cj = &bkg_cells[cjd];
    struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

#ifdef SWIFT_DEBUG_CHECKS
    tested_gparts += multi_j->m_pole.num_gpart;
#endif

    /* Avoid self contributions */
    if (top == cj) continue;

    if (cell_can_use_pair_mm(top, cj, e, e->s, /*use_rebuild_data=*/1,
                             /*is_tree_walk=*/0,
                             /*periodic boundaries*/ s->periodic,
                             /*use_mesh*/ s->periodic)) {
      /* Call the PM interaction function on the active sub-cells of ci */
      runner_dopair_grav_mm_nonsym(r, ci, cj);
      // runner_dopair_recursive_grav_pm(r, ci, cj);

      /* Record that this multipole received a contribution */
      multi_i->pot.interacted = 1;

    } /* We are in charge of this pair */
  } /* Loop over top-level cells */

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure we at tested against all possible gparts. */
  if (tested_gparts != e->s->nr_gparts) {
    error(
        "Not all gparts were tested in long range gravity task! (tested: %ld, "
        "total: %ld)",
        tested_gparts, e->s->nr_gparts);
  }
#endif
}

/**
 * @brief Performs M-M interactions between a given top-level cell and all other
 * top level cells not interacted with via pair tasks or the mesh.
 *
 * This is used when running a zoom simulation, the space is periodic and there
 * is a mesh, therefore we only interact with cells that are closer than the
 * mesh interaction distance but further than the direct interaction distance.
 *
 * This function will be "clever" and only loop over the parts of the
 * top-level grid that are not covered by the mesh. Add some buffer for
 * safety.
 *
 * This function handles the following long range interactions:
 * - zoom -> bkg
 * - bkg -> bkg
 *
 * Since we define the original pair tasks at the void level we only need to
 * consider combinations of background cells. However, If a zoom cell has a long
 * range task it means there were no tasks in the void level but we still only
 * need to consider interactions between the zoom cell and the background cells.
 * Note that the top cell pointer here is already going to be the top level void
 * parent cell if ci is a zoom cell.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_zoom_periodic(struct runner *r, struct cell *ci,
                                             struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *bkg_cells = s->zoom_props->bkg_cells_top;

  /* Get the maximum distance at which we can have a non-mesh interaction. */
  const double max_distance = e->mesh->r_cut_max;

  /* Get the multipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* Get the (i,j,k) location of the top-level cell in the grid. */
  int top_i = top->loc[0] * s->iwidth[0];
  int top_j = top->loc[1] * s->iwidth[1];
  int top_k = top->loc[2] * s->iwidth[2];

  /* Maximal distance any interaction can take place before the mesh kicks in,
   * rounded up to the next integer */
  int d =
      ceil(max_distance * max3(s->iwidth[0], s->iwidth[1], s->iwidth[2])) + 1;

  /* Ensure we don't go out of bounds */
  if (d > s->cdim[0] / 2) d = s->cdim[0] / 2;
  if (d > s->cdim[1] / 2) d = s->cdim[1] / 2;
  if (d > s->cdim[2] / 2) d = s->cdim[2] / 2;

  /* Loop over plausibly useful cells */
  for (int ii = top_i - d; ii <= top_i + d; ++ii) {
    for (int jj = top_j - d; jj <= top_j + d; ++jj) {
      for (int kk = top_k - d; kk <= top_k + d; ++kk) {

        /* Box wrap */
        const int iii = (ii + s->cdim[0]) % s->cdim[0];
        const int jjj = (jj + s->cdim[1]) % s->cdim[1];
        const int kkk = (kk + s->cdim[2]) % s->cdim[2];

        /* Get the cell */
        const int cell_index = cell_getid(s->cdim, iii, jjj, kkk);

        /* Handle on the top-level cell */
        struct cell *cj = &bkg_cells[cell_index];

        /* Avoid self contributions  */
        if (top == cj) continue;

        /* Handle on the top-level cell's gravity business*/
        const struct gravity_tensors *multi_j = cj->grav.multipole;

        /* Skip empty cells */
        if (multi_j->m_pole.M_000 == 0.f) continue;

        /* Are we beyond the distance where the truncated forces are 0 ?*/
        if (cell_can_use_mesh(e, top, cj)) {

          /* Record that this multipole received a contribution */
          multi_i->pot.interacted = 1;

          /* We are done here. */
          continue;
        }

        /* Shall we interact with this cell? */
        if (cell_can_use_pair_mm(top, cj, e, e->s, /*use_rebuild_data=*/1,
                                 /*is_tree_walk=*/0,
                                 /*periodic boundaries*/ s->periodic,
                                 /*use_mesh*/ s->periodic)) {

          /* Call the PM interaction function on the active sub-cells of ci */
          runner_dopair_grav_mm_nonsym(r, ci, cj);

          /* Record that this multipole received a contribution */
          multi_i->pot.interacted = 1;

        } /* We can interact with this cell */
      } /* Loop over relevant top-level cells (k) */
    } /* Loop over relevant top-level cells (j) */
  } /* Loop over relevant top-level cells (i) */
}

/**
 * @brief Performs M-M interactions between a given top-level cell and all other
 * top level cells not interacted with via pair tasks or the mesh.
 *
 * This is used when running a uniform box, the space is periodic and there is a
 * mesh, therefore we only interact with cells that are closer than the mesh
 * interaction distance but further than the direct interaction distance.
 *
 * This function will be "clever" and only loop over the parts of the
 * top-level grid that are not covered by the mesh. Add some buffer for
 * safety.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param top The top-level parent of the #cell of interest.
 */
void runner_do_grav_long_range_uniform_periodic(struct runner *r,
                                                struct cell *ci,
                                                struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;

  /* Get the maximum distance at which we can have a non-mesh interaction. */
  const double max_distance = e->mesh->r_cut_max;

  /* Get the multipole of the cell we are interacting. */
  struct gravity_tensors *const multi_i = ci->grav.multipole;

  /* Get the (i,j,k) location of the top-level cell in the grid. */
  int top_i = top->loc[0] * s->iwidth[0];
  int top_j = top->loc[1] * s->iwidth[1];
  int top_k = top->loc[2] * s->iwidth[2];

  /* Maximal distance any interaction can take place before the mesh kicks in,
   * rounded up to the next integer */
  int d =
      ceil(max_distance * max3(s->iwidth[0], s->iwidth[1], s->iwidth[2])) + 1;

  /* Ensure we don't go out of bounds */
  if (d > s->cdim[0] / 2) d = s->cdim[0] / 2;
  if (d > s->cdim[1] / 2) d = s->cdim[1] / 2;
  if (d > s->cdim[2] / 2) d = s->cdim[2] / 2;

  /* Loop over plausibly useful cells */
  for (int ii = top_i - d; ii <= top_i + d; ++ii) {
    for (int jj = top_j - d; jj <= top_j + d; ++jj) {
      for (int kk = top_k - d; kk <= top_k + d; ++kk) {

        /* Box wrap */
        const int iii = (ii + s->cdim[0]) % s->cdim[0];
        const int jjj = (jj + s->cdim[1]) % s->cdim[1];
        const int kkk = (kk + s->cdim[2]) % s->cdim[2];

        /* Get the cell */
        const int cell_index = cell_getid(s->cdim, iii, jjj, kkk);

        /* Handle on the top-level cell */
        struct cell *cj = &cells[cell_index];

        /* Avoid self contributions  */
        if (top == cj) continue;

        /* Handle on the top-level cell's gravity business*/
        const struct gravity_tensors *multi_j = cj->grav.multipole;

        /* Skip empty cells */
        if (multi_j->m_pole.M_000 == 0.f) continue;

        /* Are we beyond the distance where the truncated forces are 0 ?*/
        if (cell_can_use_mesh(e, top, cj)) {

          /* Record that this multipole received a contribution */
          multi_i->pot.interacted = 1;

          /* We are done here. */
          continue;
        }

        /* Shall we interact with this cell? */
        if (cell_can_use_pair_mm(top, cj, e, e->s, /*use_rebuild_data=*/1,
                                 /*is_tree_walk=*/0,
                                 /*periodic boundaries*/ s->periodic,
                                 /*use_mesh*/ s->periodic)) {
          /* Call the PM interaction function on the active sub-cells of ci
           */
          runner_dopair_grav_mm_nonsym(r, ci, cj);
          // runner_dopair_recursive_grav_pm(r, ci, cj);

          /* Record that this multipole received a contribution */
          multi_i->pot.interacted = 1;

        } /* We can interact with this cell */
      } /* Loop over relevant top-level cells (k) */
    } /* Loop over relevant top-level cells (j) */
  } /* Loop over relevant top-level cells (i) */
}

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)

static void runner_accumulate_interaction(
    struct gravity_tensors *restrict multi_i,
    struct gravity_tensors *restrict multi_j);

static void runner_count_mesh_interactions_uniform(struct runner *r,
                                                   struct cell *ci,
                                                   struct cell *top);
static void runner_count_mesh_interactions_zoom(struct runner *r,
                                                struct cell *ci,
                                                struct cell *top);

#ifdef SWIFT_DEBUG_CHECKS
void runner_debug_get_top_level_methods_by_type(const struct engine *e,
                                                const struct cell *ci,
                                                long long mesh_counts[4],
                                                long long mm_counts[4],
                                                long long p2p_counts[4],
                                                int mesh_nr_cells[4],
                                                int mm_nr_cells[4],
                                                int p2p_nr_cells[4]);
#endif

#ifdef SWIFT_DEBUG_CHECKS
struct runner_debug_mesh_overlap_entry {
  unsigned long long recipient_ptr_key;
  unsigned long long source_ptr_key;
};

struct runner_debug_recursive_mesh_budget_entry {
  unsigned long long recipient_ptr_key;
  unsigned long long source_top_ptr_key;
  const char *attachment_case;
  long long counted_gparts;
};

struct runner_debug_zoom_pair_visit_entry {
  unsigned long long recipient_ptr_key;
  unsigned long long ci_ptr_key;
  unsigned long long cj_ptr_key;
  int visit_count;
  int mesh_stop_count;
  int recurse_count;
};

static int runner_debug_mesh_overlap_report_count = 0;
static int runner_debug_mesh_overlap_limit_reported = 0;
static int runner_debug_recursive_mesh_bkg_neigh_report_count = 0;
static int runner_debug_recursive_mesh_bkg_neigh_limit_reported = 0;
static int runner_debug_mesh_pair_overlap_report_count = 0;
static int runner_debug_mesh_pair_overlap_limit_reported = 0;
static int runner_debug_zoom_mesh_stop_overlap_report_count = 0;
static int runner_debug_zoom_mesh_stop_overlap_limit_reported = 0;
static int runner_debug_pair_recursive_bkg_neigh_attach_report_count = 0;
static int runner_debug_pair_recursive_bkg_neigh_attach_limit_reported = 0;

#define runner_debug_mesh_overlap_report_max 64
#define runner_debug_recursive_mesh_bkg_neigh_report_max 128
#define runner_debug_mesh_pair_overlap_report_max 128
#define runner_debug_recursive_mesh_budget_report_max 64
#define runner_debug_zoom_pair_revisit_report_max 64
#define runner_debug_zoom_mesh_stop_overlap_report_max 128
#define runner_debug_pair_recursive_bkg_neigh_attach_report_max 128
#define runner_debug_mesh_overlap_table_size (1 << 14)
#define runner_debug_mesh_overlap_probe_limit 16
#define runner_debug_recursive_mesh_budget_table_size (1 << 15)
#define runner_debug_recursive_mesh_budget_probe_limit 32
#define runner_debug_zoom_pair_visit_table_size (1 << 15)
#define runner_debug_zoom_pair_visit_probe_limit 32

static struct runner_debug_mesh_overlap_entry
    runner_debug_mesh_overlap_table[runner_debug_mesh_overlap_table_size];
static struct runner_debug_recursive_mesh_budget_entry
    runner_debug_recursive_mesh_budget_table
        [runner_debug_recursive_mesh_budget_table_size];
static struct runner_debug_zoom_pair_visit_entry
    runner_debug_zoom_pair_visit_table[runner_debug_zoom_pair_visit_table_size];

static INLINE const struct cell *runner_debug_get_recursive_mesh_source_top(
    const struct cell *source) {

  const struct cell *source_top = source->top;
  if (source_top->void_parent != NULL) source_top = source_top->void_parent->top;
  return source_top;
}

static void runner_debug_log_recursive_mesh_bkg_neigh(
    const struct cell *super, const struct cell *ci, const struct cell *cj,
    const struct cell *recipient, const struct cell *source,
    const char *attachment_case) {

  if (runner_debug_get_source_class(source) != runner_debug_source_class_bkg_neigh)
    return;

  const int report_index =
      atomic_inc(&runner_debug_recursive_mesh_bkg_neigh_report_count);

  if (report_index < runner_debug_recursive_mesh_bkg_neigh_report_max) {
    message(
        "mesh-recursive-bkg-neigh: super=%llu (%s/%s depth=%d top=%p) "
        "ci=%llu (%s/%s depth=%d top=%p) cj=%llu (%s/%s depth=%d top=%p) "
        "recipient=%llu (%s/%s depth=%d top=%p) "
        "source=%llu (%s/%s depth=%d top=%p gparts=%lld) case=%s",
        super->cellID, cellID_names[super->type], subcellID_names[super->subtype],
        super->depth, (void *)super->top, ci->cellID, cellID_names[ci->type],
        subcellID_names[ci->subtype], ci->depth, (void *)ci->top, cj->cellID,
        cellID_names[cj->type], subcellID_names[cj->subtype], cj->depth,
        (void *)cj->top, recipient->cellID, cellID_names[recipient->type],
        subcellID_names[recipient->subtype], recipient->depth,
        (void *)recipient->top, source->cellID, cellID_names[source->type],
        subcellID_names[source->subtype], source->depth, (void *)source->top,
        source->grav.multipole->m_pole.num_gpart, attachment_case);
  } else if (report_index == runner_debug_recursive_mesh_bkg_neigh_report_max &&
             atomic_cas(&runner_debug_recursive_mesh_bkg_neigh_limit_reported, 0,
                        1) == 0) {
    message("mesh-recursive-bkg-neigh reporting capped at %d lines",
            runner_debug_recursive_mesh_bkg_neigh_report_max);
  }
}

static void runner_debug_check_mesh_overlap(const struct cell *recipient,
                                            const struct cell *source) {

  const unsigned long long recipient_ptr_key =
      (unsigned long long)(uintptr_t)recipient;
  const unsigned long long source_ptr_key =
      (unsigned long long)(uintptr_t)source;
  const unsigned long long hash =
      recipient_ptr_key * 11400714819323198485ull ^
      source_ptr_key * 14029467366897019727ull;

  for (int probe = 0; probe < runner_debug_mesh_overlap_probe_limit; probe++) {
    struct runner_debug_mesh_overlap_entry *slot =
        &runner_debug_mesh_overlap_table
             [(hash + (unsigned long long)probe) &
              (runner_debug_mesh_overlap_table_size - 1)];

    if (slot->recipient_ptr_key == 0ULL) {
      if (atomic_cas(&slot->recipient_ptr_key, 0ULL, recipient_ptr_key) == 0ULL) {
        slot->source_ptr_key = source_ptr_key;
        return;
      }
    }

    if (slot->recipient_ptr_key != recipient_ptr_key) continue;

    if (slot->source_ptr_key == source_ptr_key) return;

    const struct cell *other_source =
        (const struct cell *)(uintptr_t)slot->source_ptr_key;

    if (cell_contains_progeny(other_source, source) ||
        cell_contains_progeny(source, other_source)) {
      const int report_index = atomic_inc(&runner_debug_mesh_overlap_report_count);

      if (report_index < runner_debug_mesh_overlap_report_max) {
        message(
            "mesh-overlap: recipient=%p cellID=%llu (%s/%s depth=%d) "
            "first_source=%p cellID=%llu (%s/%s depth=%d top=%p) "
            "second_source=%p cellID=%llu (%s/%s depth=%d top=%p)",
            (void *)recipient, recipient->cellID, cellID_names[recipient->type],
            subcellID_names[recipient->subtype], recipient->depth,
            (void *)other_source, other_source->cellID,
            cellID_names[other_source->type], subcellID_names[other_source->subtype],
            other_source->depth, (void *)other_source->top, (void *)source,
            source->cellID, cellID_names[source->type],
            subcellID_names[source->subtype], source->depth, (void *)source->top);
      } else if (report_index == runner_debug_mesh_overlap_report_max &&
                 atomic_cas(&runner_debug_mesh_overlap_limit_reported, 0, 1) == 0) {
        message("mesh-overlap reporting capped at %d lines",
                runner_debug_mesh_overlap_report_max);
      }
    }

    return;
  }
}

static void runner_debug_check_mesh_pair_overlap(const struct cell *recipient,
                                                 const struct cell *source,
                                                 const int origin,
                                                 const char *attachment_case) {

  if (origin == 0) return;

  for (const struct link *l = recipient->grav.grav; l != NULL; l = l->next) {

    const struct task *t = l->t;

    if (t->type != task_type_pair || t->subtype != task_subtype_grav) continue;

    const struct cell *other = (t->ci == recipient) ? t->cj : t->ci;
    if (other == NULL) continue;

    if (other != source && !cell_contains_progeny(other, source) &&
        !cell_contains_progeny(source, other))
      continue;

    const int report_index = atomic_inc(&runner_debug_mesh_pair_overlap_report_count);

    if (report_index < runner_debug_mesh_pair_overlap_report_max) {
      message(
          "mesh-pair-overlap: recipient=%llu (%s/%s depth=%d top=%p) "
          "source=%llu (%s/%s depth=%d top=%p gparts=%lld) "
          "pair_other=%llu (%s/%s depth=%d top=%p gparts=%lld) "
          "task_cells=[%llu,%llu] origin=%d case=%s",
          recipient->cellID, cellID_names[recipient->type],
          subcellID_names[recipient->subtype], recipient->depth,
          (void *)recipient->top, source->cellID, cellID_names[source->type],
          subcellID_names[source->subtype], source->depth, (void *)source->top,
          source->grav.multipole->m_pole.num_gpart, other->cellID,
          cellID_names[other->type], subcellID_names[other->subtype],
          other->depth, (void *)other->top, other->grav.multipole->m_pole.num_gpart,
          t->ci != NULL ? t->ci->cellID : 0, t->cj != NULL ? t->cj->cellID : 0,
          origin, attachment_case);
    } else if (report_index == runner_debug_mesh_pair_overlap_report_max &&
               atomic_cas(&runner_debug_mesh_pair_overlap_limit_reported, 0,
                          1) == 0) {
      message("mesh-pair-overlap reporting capped at %d lines",
              runner_debug_mesh_pair_overlap_report_max);
    }

    return;
  }
}

static void runner_debug_check_zoom_mesh_stop_pair_overlap(
    const struct cell *recipient, const struct cell *ci, const struct cell *cj) {

  const struct cell *branches[2] = {ci, cj};
  const struct cell *partners[2] = {cj, ci};

  for (int side = 0; side < 2; side++) {
    const struct cell *branch = branches[side];
    const struct cell *partner = partners[side];

    for (const struct link *l = branch->grav.grav; l != NULL; l = l->next) {
      const struct task *t = l->t;

      if (t->type != task_type_pair || t->subtype != task_subtype_grav) continue;

      const struct cell *other = (t->ci == branch) ? t->cj : t->ci;
      if (other == NULL) continue;

      if (other != partner && !cell_contains_progeny(other, partner) &&
          !cell_contains_progeny(partner, other))
        continue;

      const int report_index =
          atomic_inc(&runner_debug_zoom_mesh_stop_overlap_report_count);

      if (report_index < runner_debug_zoom_mesh_stop_overlap_report_max) {
        message(
            "zoom-mesh-stop-pair-overlap: recipient=%llu (%s/%s depth=%d top=%p) "
            "branch=%llu (%s/%s depth=%d top=%p gparts=%lld) "
            "partner=%llu (%s/%s depth=%d top=%p gparts=%lld) "
            "pair_other=%llu (%s/%s depth=%d top=%p gparts=%lld) "
            "task_cells=[%llu,%llu]",
            recipient->cellID, cellID_names[recipient->type],
            subcellID_names[recipient->subtype], recipient->depth,
            (void *)recipient->top, branch->cellID, cellID_names[branch->type],
            subcellID_names[branch->subtype], branch->depth, (void *)branch->top,
            branch->grav.multipole->m_pole.num_gpart, partner->cellID,
            cellID_names[partner->type], subcellID_names[partner->subtype],
            partner->depth, (void *)partner->top,
            partner->grav.multipole->m_pole.num_gpart, other->cellID,
            cellID_names[other->type], subcellID_names[other->subtype],
            other->depth, (void *)other->top, other->grav.multipole->m_pole.num_gpart,
            t->ci != NULL ? t->ci->cellID : 0, t->cj != NULL ? t->cj->cellID : 0);
      } else if (report_index == runner_debug_zoom_mesh_stop_overlap_report_max &&
                 atomic_cas(&runner_debug_zoom_mesh_stop_overlap_limit_reported, 0,
                            1) == 0) {
        message("zoom-mesh-stop-pair-overlap reporting capped at %d lines",
                runner_debug_zoom_mesh_stop_overlap_report_max);
      }

      return;
    }
  }
}

static void runner_debug_log_pair_recursive_bkg_neigh_attachment(
    const struct cell *super, const struct cell *ci, const struct cell *cj,
    const struct cell *recipient, const struct cell *source,
    const char *attachment_case) {

  if (runner_debug_get_source_class(source) != runner_debug_source_class_bkg_neigh)
    return;

  const int report_index =
      atomic_inc(&runner_debug_pair_recursive_bkg_neigh_attach_report_count);

  if (report_index < runner_debug_pair_recursive_bkg_neigh_attach_report_max) {
    message(
        "pair-recursive-bkg-neigh-attach: super=%llu (%s/%s depth=%d top=%p) "
        "ci=%llu (%s/%s depth=%d top=%p) cj=%llu (%s/%s depth=%d top=%p) "
        "recipient=%llu (%s/%s depth=%d top=%p parent=%p super=%llu) "
        "source=%llu (%s/%s depth=%d top=%p gparts=%lld) case=%s",
        super->cellID, cellID_names[super->type], subcellID_names[super->subtype],
        super->depth, (void *)super->top, ci->cellID, cellID_names[ci->type],
        subcellID_names[ci->subtype], ci->depth, (void *)ci->top, cj->cellID,
        cellID_names[cj->type], subcellID_names[cj->subtype], cj->depth,
        (void *)cj->top, recipient->cellID, cellID_names[recipient->type],
        subcellID_names[recipient->subtype], recipient->depth,
        (void *)recipient->top, (void *)recipient->parent,
        recipient->grav.super != NULL ? recipient->grav.super->cellID : 0ULL,
        source->cellID, cellID_names[source->type],
        subcellID_names[source->subtype], source->depth, (void *)source->top,
        source->grav.multipole->m_pole.num_gpart, attachment_case);
  } else if (report_index ==
                 runner_debug_pair_recursive_bkg_neigh_attach_report_max &&
             atomic_cas(&runner_debug_pair_recursive_bkg_neigh_attach_limit_reported,
                        0, 1) == 0) {
    message("pair-recursive-bkg-neigh-attach reporting capped at %d lines",
            runner_debug_pair_recursive_bkg_neigh_attach_report_max);
  }
}

static void runner_debug_record_recursive_mesh_budget(
    const struct cell *recipient, const struct cell *source,
    const char *attachment_case) {

  const struct cell *source_top = runner_debug_get_recursive_mesh_source_top(source);
  const unsigned long long recipient_ptr_key =
      (unsigned long long)(uintptr_t)recipient;
  const unsigned long long source_top_ptr_key =
      (unsigned long long)(uintptr_t)source_top;
  const unsigned long long case_key = (unsigned long long)(uintptr_t)attachment_case;
  const unsigned long long hash =
      recipient_ptr_key * 11400714819323198485ull ^
      source_top_ptr_key * 14029467366897019727ull ^
      case_key * 1609587929392839161ull;

  for (int probe = 0; probe < runner_debug_recursive_mesh_budget_probe_limit;
       probe++) {
    struct runner_debug_recursive_mesh_budget_entry *slot =
        &runner_debug_recursive_mesh_budget_table
             [(hash + (unsigned long long)probe) &
              (runner_debug_recursive_mesh_budget_table_size - 1)];

    if (slot->recipient_ptr_key == 0ULL) {
      if (atomic_cas(&slot->recipient_ptr_key, 0ULL, recipient_ptr_key) == 0ULL) {
        slot->source_top_ptr_key = source_top_ptr_key;
        slot->attachment_case = attachment_case;
        slot->counted_gparts = source->grav.multipole->m_pole.num_gpart;
        return;
      }
    }

    if (slot->recipient_ptr_key != recipient_ptr_key) continue;
    if (slot->source_top_ptr_key != source_top_ptr_key) continue;
    if (slot->attachment_case != attachment_case) continue;

    accumulate_add_ll(&slot->counted_gparts, source->grav.multipole->m_pole.num_gpart);
    return;
  }
}

int runner_debug_dump_recursive_mesh_budget_overlaps(const struct cell *ci) {

  int report_count = 0;
  struct runner_debug_recursive_mesh_budget_entry path_totals
      [runner_debug_recursive_mesh_budget_report_max];
  int path_totals_count = 0;

  bzero(path_totals, sizeof(path_totals));

  for (const struct cell *recipient = ci; recipient != NULL;
       recipient = (recipient == recipient->grav.super) ? NULL : recipient->parent) {
    const unsigned long long recipient_ptr_key =
        (unsigned long long)(uintptr_t)recipient;

    for (int i = 0; i < runner_debug_recursive_mesh_budget_table_size; i++) {
      const struct runner_debug_recursive_mesh_budget_entry *slot =
          &runner_debug_recursive_mesh_budget_table[i];

      if (slot->recipient_ptr_key != recipient_ptr_key) continue;
      if (slot->source_top_ptr_key == 0ULL) continue;
      if (slot->attachment_case == NULL) continue;

      const struct cell *source_top =
          (const struct cell *)(uintptr_t)slot->source_top_ptr_key;
      const long long budget = source_top->grav.multipole->m_pole.num_gpart;
      const long long counted = slot->counted_gparts;

      int found_path_total = 0;
      for (int j = 0; j < path_totals_count; j++) {
        if (path_totals[j].source_top_ptr_key != slot->source_top_ptr_key) continue;
        if (path_totals[j].attachment_case != slot->attachment_case) continue;

        path_totals[j].counted_gparts += counted;
        found_path_total = 1;
        break;
      }

      if (!found_path_total &&
          path_totals_count < runner_debug_recursive_mesh_budget_report_max) {
        path_totals[path_totals_count].source_top_ptr_key = slot->source_top_ptr_key;
        path_totals[path_totals_count].attachment_case = slot->attachment_case;
        path_totals[path_totals_count].counted_gparts = counted;
        path_totals_count++;
      }

      if (counted <= budget) continue;

      if (report_count < runner_debug_recursive_mesh_budget_report_max) {
        message(
            "recursive-mesh-budget-overlap: recipient=%llu (%s/%s depth=%d super=%llu) "
            "source_top=%llu (%s/%s depth=%d top=%p budget=%lld) "
            "case=%s counted=%lld excess=%lld",
            recipient->cellID, cellID_names[recipient->type],
            subcellID_names[recipient->subtype], recipient->depth,
            recipient->grav.super->cellID, source_top->cellID,
            cellID_names[source_top->type], subcellID_names[source_top->subtype],
            source_top->depth, (void *)source_top->top, budget,
            slot->attachment_case, counted, counted - budget);
      } else if (report_count == runner_debug_recursive_mesh_budget_report_max) {
        message("recursive-mesh-budget-overlap reporting capped at %d lines",
                runner_debug_recursive_mesh_budget_report_max);
      }

      report_count++;
    }
  }

  for (int i = 0; i < path_totals_count; i++) {
    const struct cell *source_top =
        (const struct cell *)(uintptr_t)path_totals[i].source_top_ptr_key;
    const long long budget = source_top->grav.multipole->m_pole.num_gpart;
    const long long counted = path_totals[i].counted_gparts;

    if (counted <= budget) continue;

    if (report_count < runner_debug_recursive_mesh_budget_report_max) {
      message(
          "recursive-mesh-budget-path-overlap: leaf=%llu (%s/%s depth=%d super=%llu) "
          "source_top=%llu (%s/%s depth=%d top=%p budget=%lld) "
          "case=%s counted=%lld excess=%lld",
          ci->cellID, cellID_names[ci->type], subcellID_names[ci->subtype], ci->depth,
          ci->grav.super->cellID, source_top->cellID,
          cellID_names[source_top->type], subcellID_names[source_top->subtype],
          source_top->depth, (void *)source_top->top, budget,
          path_totals[i].attachment_case, counted, counted - budget);
    } else if (report_count == runner_debug_recursive_mesh_budget_report_max) {
      message("recursive-mesh-budget-overlap reporting capped at %d lines",
              runner_debug_recursive_mesh_budget_report_max);
    }

    report_count++;
  }

  return report_count;
}

static void runner_debug_record_zoom_pair_visit(const struct cell *recipient,
                                                const struct cell *ci,
                                                const struct cell *cj,
                                                const int mesh_stop) {

  const unsigned long long recipient_ptr_key =
      (unsigned long long)(uintptr_t)recipient;
  unsigned long long ci_ptr_key = (unsigned long long)(uintptr_t)ci;
  unsigned long long cj_ptr_key = (unsigned long long)(uintptr_t)cj;

  if (cj_ptr_key < ci_ptr_key) {
    const unsigned long long tmp = ci_ptr_key;
    ci_ptr_key = cj_ptr_key;
    cj_ptr_key = tmp;
  }

  const unsigned long long hash =
      recipient_ptr_key * 11400714819323198485ull ^
      ci_ptr_key * 14029467366897019727ull ^ cj_ptr_key * 1609587929392839161ull;

  for (int probe = 0; probe < runner_debug_zoom_pair_visit_probe_limit; probe++) {
    struct runner_debug_zoom_pair_visit_entry *slot =
        &runner_debug_zoom_pair_visit_table
             [(hash + (unsigned long long)probe) &
              (runner_debug_zoom_pair_visit_table_size - 1)];

    if (slot->recipient_ptr_key == 0ULL) {
      if (atomic_cas(&slot->recipient_ptr_key, 0ULL, recipient_ptr_key) == 0ULL) {
        slot->ci_ptr_key = ci_ptr_key;
        slot->cj_ptr_key = cj_ptr_key;
        slot->visit_count = 1;
        slot->mesh_stop_count = mesh_stop;
        slot->recurse_count = !mesh_stop;
        return;
      }
    }

    if (slot->recipient_ptr_key != recipient_ptr_key) continue;
    if (slot->ci_ptr_key != ci_ptr_key) continue;
    if (slot->cj_ptr_key != cj_ptr_key) continue;

    slot->visit_count += 1;
    slot->mesh_stop_count += mesh_stop;
    slot->recurse_count += !mesh_stop;
    return;
  }
}

int runner_debug_dump_zoom_pair_recursive_revisits(const struct cell *ci) {

  int report_count = 0;

  for (const struct cell *recipient = ci; recipient != NULL;
       recipient = (recipient == recipient->grav.super) ? NULL : recipient->parent) {
    const unsigned long long recipient_ptr_key =
        (unsigned long long)(uintptr_t)recipient;

    for (int i = 0; i < runner_debug_zoom_pair_visit_table_size; i++) {
      const struct runner_debug_zoom_pair_visit_entry *slot =
          &runner_debug_zoom_pair_visit_table[i];

      if (slot->recipient_ptr_key != recipient_ptr_key) continue;
      if (slot->visit_count <= 1) continue;

      const struct cell *ci_branch = (const struct cell *)(uintptr_t)slot->ci_ptr_key;
      const struct cell *cj_branch = (const struct cell *)(uintptr_t)slot->cj_ptr_key;

      if (report_count < runner_debug_zoom_pair_revisit_report_max) {
        message(
            "zoom-pair-revisit: recipient=%llu (%s/%s depth=%d super=%llu) "
            "ci=%llu (%s/%s depth=%d top=%p) "
            "cj=%llu (%s/%s depth=%d top=%p) visits=%d mesh_stops=%d recurses=%d",
            recipient->cellID, cellID_names[recipient->type],
            subcellID_names[recipient->subtype], recipient->depth,
            recipient->grav.super->cellID, ci_branch->cellID,
            cellID_names[ci_branch->type], subcellID_names[ci_branch->subtype],
            ci_branch->depth, (void *)ci_branch->top, cj_branch->cellID,
            cellID_names[cj_branch->type], subcellID_names[cj_branch->subtype],
            cj_branch->depth, (void *)cj_branch->top, slot->visit_count,
            slot->mesh_stop_count, slot->recurse_count);
      } else if (report_count == runner_debug_zoom_pair_revisit_report_max) {
        message("zoom-pair-revisit reporting capped at %d lines",
                runner_debug_zoom_pair_revisit_report_max);
      }

      report_count++;
    }
  }

  return report_count;
}

int runner_debug_dump_zoom_pair_handoff_overlaps(const struct cell *ci) {

  struct {
    unsigned long long source_top_ptr_key;
    long long void_stage_count;
    long long zoom_stage_count;
  } stage_totals[runner_debug_recursive_mesh_budget_report_max];
  int stage_totals_count = 0;
  int report_count = 0;

  bzero(stage_totals, sizeof(stage_totals));

  for (const struct cell *recipient = ci; recipient != NULL;) {
    int stage = -1;

    if (recipient->type == cell_type_bkg &&
        recipient->subtype == cell_subtype_void) {
      stage = 0;
    } else if (recipient->type == cell_type_zoom) {
      stage = 1;
    }

    if (stage != -1) {
      const unsigned long long recipient_ptr_key =
          (unsigned long long)(uintptr_t)recipient;

      for (int i = 0; i < runner_debug_recursive_mesh_budget_table_size; i++) {
        const struct runner_debug_recursive_mesh_budget_entry *slot =
            &runner_debug_recursive_mesh_budget_table[i];

        if (slot->recipient_ptr_key != recipient_ptr_key) continue;
        if (slot->source_top_ptr_key == 0ULL) continue;
        if (slot->attachment_case == NULL) continue;
        if (strcmp(slot->attachment_case, "zoom_pair_recursive") != 0) continue;

        const struct cell *source_top =
            (const struct cell *)(uintptr_t)slot->source_top_ptr_key;
        if (runner_debug_get_source_class(source_top) !=
            runner_debug_source_class_bkg_neigh)
          continue;

        int found = 0;
        for (int j = 0; j < stage_totals_count; j++) {
          if (stage_totals[j].source_top_ptr_key != slot->source_top_ptr_key)
            continue;

          if (stage == 0)
            stage_totals[j].void_stage_count += slot->counted_gparts;
          else
            stage_totals[j].zoom_stage_count += slot->counted_gparts;

          found = 1;
          break;
        }

        if (!found &&
            stage_totals_count < runner_debug_recursive_mesh_budget_report_max) {
          stage_totals[stage_totals_count].source_top_ptr_key =
              slot->source_top_ptr_key;
          if (stage == 0)
            stage_totals[stage_totals_count].void_stage_count =
                slot->counted_gparts;
          else
            stage_totals[stage_totals_count].zoom_stage_count =
                slot->counted_gparts;
          stage_totals_count++;
        }
      }
    }

    if (recipient->void_parent != NULL)
      recipient = recipient->void_parent;
    else
      recipient = recipient->parent;
  }

  for (int i = 0; i < stage_totals_count; i++) {
    const struct cell *source_top =
        (const struct cell *)(uintptr_t)stage_totals[i].source_top_ptr_key;
    const long long void_stage = stage_totals[i].void_stage_count;
    const long long zoom_stage = stage_totals[i].zoom_stage_count;

    if (void_stage == 0 || zoom_stage == 0) continue;

    if (report_count < runner_debug_recursive_mesh_budget_report_max) {
      message(
          "zoom-pair-handoff-overlap: leaf=%llu (%s/%s depth=%d super=%llu) "
          "source_top=%llu (%s/%s depth=%d top=%p budget=%lld) "
          "void_stage=%lld zoom_stage=%lld total=%lld excess=%lld",
          ci->cellID, cellID_names[ci->type], subcellID_names[ci->subtype],
          ci->depth, ci->grav.super->cellID, source_top->cellID,
          cellID_names[source_top->type], subcellID_names[source_top->subtype],
          source_top->depth, (void *)source_top->top,
          source_top->grav.multipole->m_pole.num_gpart, void_stage, zoom_stage,
          void_stage + zoom_stage,
          void_stage + zoom_stage - source_top->grav.multipole->m_pole.num_gpart);
    } else if (report_count == runner_debug_recursive_mesh_budget_report_max) {
      message("zoom-pair-handoff-overlap reporting capped at %d lines",
              runner_debug_recursive_mesh_budget_report_max);
    }

    report_count++;
  }

  return report_count;
}

int runner_debug_dump_path_pair_recursive_sources(const struct cell *leaf) {

  const struct cell *path[64];
  int depth = 0;
  int report_count = 0;

  for (const struct cell *c = leaf; c != NULL;) {
    if (depth == 64) error("Cell path too deep for debug dump.");
    path[depth++] = c;

    if (c->void_parent != NULL)
      c = c->void_parent;
    else
      c = c->parent;
  }

  for (int p = depth - 1; p >= 0; p--) {
    const struct cell *recipient = path[p];
    const unsigned long long recipient_ptr_key =
        (unsigned long long)(uintptr_t)recipient;
    const int path_level = depth - 1 - p;

    for (int i = 0; i < runner_debug_recursive_mesh_budget_table_size; i++) {
      const struct runner_debug_recursive_mesh_budget_entry *slot =
          &runner_debug_recursive_mesh_budget_table[i];

      if (slot->recipient_ptr_key != recipient_ptr_key) continue;
      if (slot->source_top_ptr_key == 0ULL) continue;
      if (slot->attachment_case == NULL) continue;
      if (strcmp(slot->attachment_case, "pair_recursive") != 0 &&
          strcmp(slot->attachment_case, "zoom_pair_recursive") != 0)
        continue;

      const struct cell *source_top =
          (const struct cell *)(uintptr_t)slot->source_top_ptr_key;
      const enum runner_debug_source_class source_class =
          runner_debug_get_source_class(source_top);

      if (report_count < runner_debug_recursive_mesh_budget_report_max) {
        message(
            "path-pair-recursive-source: level=%d/%d recipient=%llu (%s/%s depth=%d super=%llu) "
            "case=%s source_top=%llu (%s/%s depth=%d top=%p class=%d budget=%lld) counted=%lld",
            path_level, depth, recipient->cellID, cellID_names[recipient->type],
            subcellID_names[recipient->subtype], recipient->depth,
            recipient->grav.super != NULL ? recipient->grav.super->cellID : 0ULL,
            slot->attachment_case, source_top->cellID,
            cellID_names[source_top->type], subcellID_names[source_top->subtype],
            source_top->depth, (void *)source_top->top, (int)source_class,
            source_top->grav.multipole->m_pole.num_gpart, slot->counted_gparts);
      } else if (report_count == runner_debug_recursive_mesh_budget_report_max) {
        message("path-pair-recursive-source reporting capped at %d lines",
                runner_debug_recursive_mesh_budget_report_max);
      }

      report_count++;
    }
  }

  return report_count;
}
#endif

static INLINE void runner_record_mesh_attachment(
    struct cell *super, struct cell *ci, struct cell *cj,
    struct cell *recipient, struct cell *source, const int origin,
    const char *attachment_case) {

  (void)super;
  (void)ci;
  (void)cj;
  (void)origin;
  (void)attachment_case;

  runner_accumulate_interaction(recipient->grav.multipole, source->grav.multipole);

#ifdef SWIFT_DEBUG_CHECKS
  runner_debug_check_mesh_overlap(recipient, source);
  runner_debug_check_mesh_pair_overlap(recipient, source, origin,
                                       attachment_case);
  if (origin != 0)
    runner_debug_record_recursive_mesh_budget(recipient, source, attachment_case);
  if (origin == 2)
    runner_debug_log_pair_recursive_bkg_neigh_attachment(
        super, ci, cj, recipient, source, attachment_case);
  if (origin != 0)
    runner_debug_log_recursive_mesh_bkg_neigh(super, ci, cj, recipient, source,
                                              attachment_case);
  runner_debug_add_tensor_interactions_by_type(
      recipient->grav.multipole->pot.num_interacted_pm_by_type, source,
      source->grav.multipole->m_pole.num_gpart);
  runner_debug_add_tensor_interactions_by_type(
      recipient->grav.multipole->pot.num_interacted_pm_long_range_by_type,
      source, source->grav.multipole->m_pole.num_gpart);
  if (origin == 0) {
    runner_debug_add_tensor_interactions_by_type(
        recipient->grav.multipole->pot.num_interacted_pm_long_range_direct_by_type,
        source, source->grav.multipole->m_pole.num_gpart);
  } else if (origin == 1) {
    runner_debug_add_tensor_interactions_by_type(
        recipient->grav.multipole
            ->pot.num_interacted_pm_long_range_self_recursive_by_type,
        source, source->grav.multipole->m_pole.num_gpart);
  } else if (origin == 2) {
    runner_debug_add_tensor_interactions_by_type(
        recipient->grav.multipole
            ->pot.num_interacted_pm_long_range_pair_recursive_by_type,
        source, source->grav.multipole->m_pole.num_gpart);
  }
#endif
}

/**
 * @brief Increment the mesh interaction counters.
 *
 * This is a helper function for incrementing the mesh interaction counters
 * for debugging purposes.
 *
 * @param multi_i The multipole receiving the interaction.
 * @param multi_j The multipole giving the interaction.
 */
static void runner_accumulate_interaction(
    struct gravity_tensors *restrict multi_i,
    struct gravity_tensors *restrict multi_j) {

  /* Ensure we aren't self-interacting */
  if (multi_i == multi_j) {
    error("Self interactions should not be handled in this function!");
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Need to account for the mesh interactions we missed */
  accumulate_add_ll(&multi_i->pot.num_interacted, multi_j->m_pole.num_gpart);
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Need to account for the mesh interactions we missed */
  accumulate_add_ll(&multi_i->pot.num_interacted_pm, multi_j->m_pole.num_gpart);
#endif
  /* Record that this multipole received a contribution */
  multi_i->pot.interacted = 1;
}

/**
 * @brief Count a mesh interaction between two related cells.
 *
 * Since the counts are accumulated downwards from the super level in the
 * grav/down task we need to update different cells based on the cells we have
 * been passed. This function will select what cell should be updated based on
 * the nesting:
 *   - If the super cell is nested within ci, then the super cell is updated.
 *   - If the super cell is ci, then they are both the same and it does not
 *     matter which is updated.
 *   - If ci is nested within the super cell, then ci is updated.
 *   - If we have a self interaction (pair nested within the same top cell):
 *     - If the super cell is nested within cj, then the super cell is updated.
 *     - If the super cell is cj, then they are both the same and it does not
 *       matter which is updated.
 *     - If cj is nested within the super cell, then cj is updated.
 *
 * @param super The super-cell being updated.
 * @param ci The first #cell in the pair.
 * @param cj The second #cell in the pair.
 */
static void runner_count_mesh_interaction(struct cell *super, struct cell *ci,
                                          struct cell *cj, const int origin,
                                          const char *attachment_case) {

  /* Get the correct top level cells for self-interaction check */
  struct cell *top_i = ci->top;
  if (top_i->void_parent != NULL) top_i = top_i->void_parent->top;
  struct cell *top_j = cj->top;
  if (top_j->void_parent != NULL) top_j = top_j->void_parent->top;

  /* Do we share the same top level cell? i.e. are we self-interacting? */
  int is_self = top_i == top_j;

  /* Decide which cell we are updating. */
  if (super == ci) {
    runner_record_mesh_attachment(super, ci, cj, super, cj, origin,
                                  attachment_case);
  } else if (cell_contains_progeny(ci, super)) {
    runner_record_mesh_attachment(super, ci, cj, super, cj, origin,
                                  attachment_case);
  } else if (cell_contains_progeny(super, ci)) {
    runner_record_mesh_attachment(super, ci, cj, ci, cj, origin,
                                  attachment_case);
  }

  /* Handle the symmetric case for self interactions */
  if (is_self) {
    if (super == cj) {
      runner_record_mesh_attachment(super, ci, cj, super, ci, origin,
                                    attachment_case);
    } else if (cell_contains_progeny(cj, super)) {
      runner_record_mesh_attachment(super, ci, cj, super, ci, origin,
                                    attachment_case);
    } else if (cell_contains_progeny(super, cj)) {
      runner_record_mesh_attachment(super, ci, cj, cj, ci, origin,
                                    attachment_case);
    }
  }
}

/**
 * @brief Recursively accumulate mesh interactions for pair interactions.
 *
 * This function mirrors the logic in scheduler_splittask_gravity for pair
 * tasks, recursing down the cell hierarchy and counting mesh interactions.
 *
 * @param ci The #cell of interest (active cell receiving interactions).
 * @param cpi The current #cell from ci's hierarchy being processed.
 * @param cpj The current #cell from cj's hierarchy being processed.
 * @param s The #space.
 */
static void runner_count_mesh_interactions_pair_recursive(struct cell *c,
                                                          struct cell *ci,
                                                          struct cell *cj,
                                                          struct space *s) {

  /* No self interactions here */
  if (ci == cj) {
    return;
  }

  struct engine *e = s->e;

  /* Maintain the invariant that ci is the branch containing c whenever one
   * side of the pair does. */
  const int ci_contains_c = (ci == c || cell_contains_progeny(ci, c));
  const int cj_contains_c = (cj == c || cell_contains_progeny(cj, c));
  if (!ci_contains_c && cj_contains_c) {
    struct cell *tmp = ci;
    ci = cj;
    cj = tmp;
  }

  /* Foreign pair? This mirrors scheduler_splittask_gravity. */
  if (ci->nodeID != engine_rank && cj->nodeID != engine_rank) {
    return;
  }

  /* Should this pair be split? */
  if (cell_can_split_pair_gravity_task(ci, cj)) {

    /* Check particle count threshold - mirrors scheduler_splittask_gravity */
    const long long gcount_i = ci->grav.count;
    const long long gcount_j = cj->grav.count;
    if (gcount_i * gcount_j < ((long long)space_subsize_pair_grav)) {
      return;
    }

    /* Recurse on all progeny pairs */
    for (int i = 0; i < 8; i++) {
      if (ci->progeny[i] == NULL) continue;
      struct cell *cpi = ci->progeny[i];
      for (int j = 0; j < 8; j++) {
        if (cj->progeny[j] == NULL) continue;
        struct cell *cpj = cj->progeny[j];

        /* Can we use the mesh for this pair? */
        if (cell_can_use_mesh(e, cpi, cpj)) {
          /* Record the mesh interaction */
          runner_count_mesh_interaction(c, cpi, cpj, 2, "pair_recursive");
          continue;
        }

        const int cpi_contains_c = (cpi == c || cell_contains_progeny(cpi, c));
        const int cpj_contains_c = (cpj == c || cell_contains_progeny(cpj, c));

        /* We would create real tasks, so recurse to find mesh interactions */
        if (!cpi_contains_c && cpj_contains_c) {
          runner_count_mesh_interactions_pair_recursive(c, cpj, cpi, s);
        } else {
          runner_count_mesh_interactions_pair_recursive(c, cpi, cpj, s);
        }
      }
    }
  }
  /* else: We have a real task that doesn't split further, no mesh
   * interactions to count */
}

/**
 * @brief Recursively accumulate mesh interactions for self interactions.
 *
 * This function mirrors the logic in scheduler_splittask_gravity for self
 * tasks, recursing down the cell hierarchy and counting mesh interactions.
 *
 * @param c The #cell of interest (active cell receiving interactions).
 * @param ci The current #cell from c's hierarchy being processed.
 * @param s The #space.
 */
static void runner_count_mesh_interactions_self_recursive(struct cell *c,
                                                          struct cell *ci,
                                                          struct space *s) {

  struct engine *e = s->e;

  /* Foreign task? This mirrors scheduler_splittask_gravity. */
  if (ci->nodeID != engine_rank) {
    return;
  }

  /* Should this self task be split? */
  if (cell_can_split_self_gravity_task(ci)) {

    /* Check particle count threshold - mirrors scheduler_splittask_gravity */
    if (ci->grav.count < space_subsize_self_grav) {
      return;
    }

    /* Recurse on self interactions for each progeny */
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] == NULL) continue;
      runner_count_mesh_interactions_self_recursive(c, ci->progeny[k], s);
    }

    /* Now handle pair interactions between progeny */
    for (int j = 0; j < 8; j++) {
      if (ci->progeny[j] == NULL) continue;
      struct cell *cpj = ci->progeny[j];
      for (int k = j + 1; k < 8; k++) {
        if (ci->progeny[k] == NULL) continue;
        struct cell *cpk = ci->progeny[k];

        /* Can we use the mesh for this pair? */
        if (cell_can_use_mesh(e, cpj, cpk)) {
          /* Record the mesh interaction */
          runner_count_mesh_interaction(c, cpj, cpk, 1, "self_recursive");
          continue;
        }

        /* Check which cell contains c to decide how to recurse. */
        const int cpj_contains_c = (cpj == c || cell_contains_progeny(cpj, c));
        const int cpk_contains_c = (cpk == c || cell_contains_progeny(cpk, c));

        /* Recurse as a pair interaction, ensuring we always pass the cell
         * containing c as the first argument after c. */
        if (cpj_contains_c) {
          runner_count_mesh_interactions_pair_recursive(c, cpj, cpk, s);
        } else if (cpk_contains_c) {
          runner_count_mesh_interactions_pair_recursive(c, cpk, cpj, s);
        } else {
          runner_count_mesh_interactions_pair_recursive(c, cpj, cpk, s);
        }
      }
    }
  }
  /* else: We have a real task that doesn't split further, no mesh
   * interactions to count */
}

/**
 * @brief Accumulate the number of particle mesh interactions for debugging
 * purposes.
 *
 * This function mirrors the task creation and splitting logic to count
 * mesh interactions that ci would receive from all top-level cells.
 *
 * NOTE: This will recurse over cells that are not directly realted to c (the
 * super cell of ci). It will not add their contribution though so is "safe".
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest (active cell receiving interactions).
 * @param top The current top-level cell (ci's top-level parent).
 */
static void runner_count_mesh_interactions_uniform(struct runner *r,
                                                   struct cell *ci,
                                                   struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;

  /* First, handle self interactions from the top-level cell.
   * This mirrors the self task created at the top level. */
  runner_count_mesh_interactions_self_recursive(ci, top, s);

  /* Now loop over all other top-level cells for pair interactions.
   * This mirrors the pair tasks created between top-level cells. */
  for (int n = 0; n < s->nr_cells; n++) {

    /* Handle on the top-level cell and its gravity business */
    struct cell *cj = &cells[n];
    struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Avoid self contributions (already handled above) */
    if (top == cj) continue;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

    /* Can we use the mesh for this top-level pair? */
    if (cell_can_use_mesh(e, top, cj)) {

      /* If so, record the mesh interaction */
      runner_count_mesh_interaction(ci, top, cj, 0, "top_level_direct");
      continue;
    }

    /* Can we use M-M for this top-level pair? */
    if (cell_can_use_pair_mm(top, cj, e, s, /*use_rebuild_data=*/1,
                             /*is_tree_walk=*/0,
                             /*periodic boundaries*/ s->periodic,
                             /*use_mesh*/ s->periodic)) {

      /* M-M task handles this, nothing to count */
      continue;
    }

    /* We would create a pair task here, so recurse to count mesh interactions
     * that arise from task splitting */
    runner_count_mesh_interactions_pair_recursive(ci, top, cj, s);
  }
}

/**
 * @brief Recursively accumulate mesh interactions for zoom pair interactions.
 *
 * This function mirrors the logic in
 * zoom_scheduler_splittask_gravity_void_pair, recursing through the void cell
 * hierarchy and then using the normal pair recursive function once we reach
 * non-void cells.
 *
 * @param c The #cell of interest (active cell receiving interactions).
 * @param ci The current #cell from c's hierarchy being processed.
 * @param cj The current #cell being paired with.
 * @param s The #space.
 */
static void runner_count_mesh_interactions_zoom_pair_recursive(
    struct cell *c, struct cell *ci, struct cell *cj, struct space *s) {

  struct engine *e = s->e;

  /* Maintain the invariant that ci is the branch containing c whenever one
   * side of the pair does. */
  const int ci_contains_c = (ci == c || cell_contains_progeny(ci, c));
  const int cj_contains_c = (cj == c || cell_contains_progeny(cj, c));
  if (!ci_contains_c && cj_contains_c) {
    struct cell *tmp = ci;
    ci = cj;
    cj = tmp;
  }

  /* If neither cell is a void cell, use the normal pair recursive function */
  if (ci->subtype != cell_subtype_void && cj->subtype != cell_subtype_void) {
    runner_count_mesh_interactions_pair_recursive(c, ci, cj, s);
    return;
  }

#ifdef SWIFT_DEBUG_CHECKS
  runner_debug_record_zoom_pair_visit(c, ci, cj, 0);
#endif

  /* Loop over progeny pairs, mirroring
   * zoom_scheduler_splittask_gravity_void_pair */
  for (int i = 0; i < 8; i++) {
    struct cell *cpi = ci->progeny[i];

    /* Skip NULL progeny */
    if (cpi == NULL) continue;

    /* Skip void progeny that do not contain any zoom cells. */
    if (cpi->subtype == cell_subtype_void && !cpi->contains_zoom_cells)
      continue;

    /* Skip any empty progeny of a void cell (void cells themselves always
     * have 0 particles but are never "empty"). */
    if (cell_is_empty_mpole(cpi)) continue;

    for (int j = 0; j < 8; j++) {
      struct cell *cpj = cj->progeny[j];

      /* Skip NULL progeny */
      if (cpj == NULL) continue;

      /* Skip void progeny that do not contain any zoom cells. */
      if (cpj->subtype == cell_subtype_void && !cpj->contains_zoom_cells)
        continue;

      /* Skip leaf neighbours interacting with void cells. */
      if (!ci->split && cpj->subtype == cell_subtype_void) continue;

      /* Skip any empty progeny of a void cell (void cells themselves always
       * have 0 particles but are never "empty"). */
      if (cell_is_empty_mpole(cpj)) continue;

      /* Skip entirely foreign pairs. */
      if (cpi->nodeID != engine_rank && cpj->nodeID != engine_rank) continue;

      /* Can we use the mesh for this pair? */
      if (cell_can_use_mesh(e, cpi, cpj)) {
#ifdef SWIFT_DEBUG_CHECKS
        runner_debug_record_zoom_pair_visit(c, cpi, cpj, 1);
        runner_debug_check_zoom_mesh_stop_pair_overlap(c, cpi, cpj);
#endif
        /* Record the mesh interaction */
        runner_count_mesh_interaction(c, cpi, cpj, 2, "zoom_pair_recursive");
        continue;
      }

      const int cpi_contains_c = (cpi == c || cell_contains_progeny(cpi, c));
      const int cpj_contains_c = (cpj == c || cell_contains_progeny(cpj, c));

      /* Recurse to find more mesh interactions */
      if (!cpi_contains_c && cpj_contains_c) {
        runner_count_mesh_interactions_zoom_pair_recursive(c, cpj, cpi, s);
      } else {
        runner_count_mesh_interactions_zoom_pair_recursive(c, cpi, cpj, s);
      }
    }
  }
}

/**
 * @brief Recursively accumulate mesh interactions for zoom self interactions.
 *
 * This function mirrors the logic in
 * zoom_scheduler_splittask_gravity_void_self, recursing through the void cell
 * hierarchy and then using the normal self recursive function once we reach
 * non-void cells.
 *
 * @param c The #cell of interest (active cell receiving interactions).
 * @param ci The current #cell from c's hierarchy being processed.
 * @param s The #space.
 */
static void runner_count_mesh_interactions_zoom_self_recursive(
    struct cell *c, struct cell *ci, struct space *s) {

  struct engine *e = s->e;

  /* If not a void cell, use the normal self recursive function */
  if (ci->subtype != cell_subtype_void) {
    runner_count_mesh_interactions_self_recursive(c, ci, s);
    return;
  }

  /* Loop over progeny for self interactions */
  for (int k = 0; k < 8; k++) {
    if (ci->progeny[k] == NULL) continue;

    /* Skip void progeny that do not contain any zoom cells. */
    if (ci->progeny[k]->subtype == cell_subtype_void &&
        !ci->progeny[k]->contains_zoom_cells)
      continue;

    /* Skip empty progeny (void cells are never empty). */
    if (cell_is_empty_mpole(ci->progeny[k])) continue;

    /* Skip foreign progeny (no such thing as a foreign self task). */
    if (ci->progeny[k]->type == cell_type_zoom &&
        ci->progeny[k]->nodeID != engine_rank)
      continue;

    runner_count_mesh_interactions_zoom_self_recursive(c, ci->progeny[k], s);
  }

  /* Now handle pair interactions between progeny */
  for (int j = 0; j < 8; j++) {
    if (ci->progeny[j] == NULL) continue;

    /* Skip void progeny that do not contain any zoom cells. */
    if (ci->progeny[j]->subtype == cell_subtype_void &&
        !ci->progeny[j]->contains_zoom_cells)
      continue;

    /* Skip empty non-void progeny. */
    if (ci->progeny[j]->subtype != cell_subtype_void &&
        ci->progeny[j]->grav.count == 0)
      continue;

    /* Skip empty progeny. */
    if (cell_is_empty_mpole(ci->progeny[j])) continue;

    struct cell *cpj = ci->progeny[j];

    for (int k = j + 1; k < 8; k++) {
      if (ci->progeny[k] == NULL) continue;

      /* Skip void progeny that do not contain any zoom cells. */
      if (ci->progeny[k]->subtype == cell_subtype_void &&
          !ci->progeny[k]->contains_zoom_cells)
        continue;

      /* Skip empty progeny. */
      if (cell_is_empty_mpole(ci->progeny[k])) continue;

      struct cell *cpk = ci->progeny[k];

      /* Skip entirely foreign pairs. */
      if ((cpj->type == cell_type_zoom && cpk->type == cell_type_zoom) &&
          cpj->nodeID != engine_rank && cpk->nodeID != engine_rank)
        continue;

      /* Can we use the mesh for this pair? */
      if (cell_can_use_mesh(e, cpj, cpk)) {
        /* Record the mesh interaction */
        runner_count_mesh_interaction(c, cpj, cpk, 1, "zoom_self_recursive");
        continue;
      }

      /* Check which cell contains c to decide how to recurse. */
      const int cpj_contains_c = (cpj == c || cell_contains_progeny(cpj, c));
      const int cpk_contains_c = (cpk == c || cell_contains_progeny(cpk, c));

      /* Recurse as a pair interaction, ensuring we always pass the cell
       * containing c as the first argument after c. */
      if (cpj_contains_c) {
        runner_count_mesh_interactions_zoom_pair_recursive(c, cpj, cpk, s);
      } else if (cpk_contains_c) {
        runner_count_mesh_interactions_zoom_pair_recursive(c, cpk, cpj, s);
      } else {
        runner_count_mesh_interactions_zoom_pair_recursive(c, cpj, cpk, s);
      }
    }
  }
}

/**
 * @brief Accumulate the number of particle mesh interactions for debugging
 * purposes in zoom simulations.
 *
 * This function mirrors the task creation and splitting logic for zoom
 * simulations to count mesh interactions that ci would receive from all
 * top-level cells.
 *
 * NOTE: This will recurse over cells that are not directly related to c (the
 * super cell of ci). It will not add their contribution though so is "safe".
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest (active cell receiving interactions).
 * @param top The current top-level cell (ci's top-level parent).
 */
static void runner_count_mesh_interactions_zoom(struct runner *r,
                                                struct cell *ci,
                                                struct cell *top) {

  struct engine *e = r->e;
  struct space *s = e->s;
  struct cell *bkg_cells = s->zoom_props->bkg_cells_top;

  /* First, handle self interactions from the top-level cell.
   * This mirrors the self task created at the void level. */
  runner_count_mesh_interactions_zoom_self_recursive(ci, top, s);

  /* Now loop over all other background/void top-level cells for pair
   * interactions. This mirrors the pair tasks created between void cells. */
  for (int n = 0; n < s->zoom_props->nr_bkg_cells; n++) {

    /* Handle on the top-level cell and its gravity business */
    struct cell *cj = &bkg_cells[n];
    struct gravity_tensors *const multi_j = cj->grav.multipole;

    /* Avoid self contributions (already handled above) */
    if (top == cj) continue;

    /* Skip empty cells */
    if (multi_j->m_pole.M_000 == 0.f) continue;

    /* Can we use the mesh for this top-level pair? */
    if (cell_can_use_mesh(e, top, cj)) {

      /* If so, record the mesh interaction */
      runner_count_mesh_interaction(ci, top, cj, 0, "zoom_top_level_direct");
      continue;
    }

    /* Can we use M-M for this top-level pair? */
    if (cell_can_use_pair_mm(top, cj, e, s, /*use_rebuild_data=*/1,
                             /*is_tree_walk=*/0,
                             /*periodic boundaries*/ s->periodic,
                             /*use_mesh*/ s->periodic)) {

      /* M-M task handles this, nothing to count */
      continue;
    }

    /* We would create a pair task here, so recurse to count mesh interactions
     * that arise from task splitting through the void hierarchy */
    runner_count_mesh_interactions_zoom_pair_recursive(ci, top, cj, s);
  }
}

#ifdef SWIFT_DEBUG_CHECKS
void runner_debug_get_top_level_methods_by_type(const struct engine *e,
                                                const struct cell *ci,
                                                long long mesh_counts[4],
                                                long long mm_counts[4],
                                                long long p2p_counts[4],
                                                int mesh_nr_cells[4],
                                                int mm_nr_cells[4],
                                                int p2p_nr_cells[4]) {

  struct cell *top = ci->top;
  while (top->void_parent != NULL) top = top->void_parent->top;

  for (int i = 0; i < 4; i++) {
    mesh_counts[i] = 0;
    mm_counts[i] = 0;
    p2p_counts[i] = 0;
    mesh_nr_cells[i] = 0;
    mm_nr_cells[i] = 0;
    p2p_nr_cells[i] = 0;
  }

  if (!e->s->periodic) return;

  struct space *s = e->s;
  struct cell *cells = s->zoom_props->bkg_cells_top;

  for (int n = 0; n < s->zoom_props->nr_bkg_cells; n++) {
    struct cell *source = &cells[n];
    const struct gravity_tensors *multi = source->grav.multipole;
    const enum runner_debug_source_class source_class =
        runner_debug_get_source_class(source);

    if (source == top) continue;
    if (multi->m_pole.M_000 == 0.f) continue;
    if (cell_can_use_mesh((struct engine *)e, top, source)) {
      mesh_counts[source_class] += multi->m_pole.num_gpart;
      mesh_nr_cells[source_class] += 1;
    } else if (cell_can_use_pair_mm(top, source, (struct engine *)e, s,
                                    /*use_rebuild_data=*/1,
                                    /*is_tree_walk=*/0,
                                    /*periodic boundaries=*/s->periodic,
                                    /*use_mesh=*/s->periodic)) {
      mm_counts[source_class] += multi->m_pole.num_gpart;
      mm_nr_cells[source_class] += 1;
    } else {
      p2p_counts[source_class] += multi->m_pole.num_gpart;
      p2p_nr_cells[source_class] += 1;
    }
  }
}
#endif

#endif /* SWIFT_DEBUG_CHECKS || SWIFT_GRAVITY_FORCE_CHECKS */

/**
 * @brief Performs all M-M interactions between a given top-level cell and
 * all the other top-levels that are far enough.
 *
 * @param r The thread #runner.
 * @param ci The #cell of interest.
 * @param timer Are we timing this ?
 */
void runner_do_grav_long_range(struct runner *r, struct cell *ci,
                               const int timer) {

  TIMER_TIC;

  struct space *s = r->e->s;

  /* Is the space periodic? */
  const int periodic = s->periodic;

  /* Anything to do here? */
  if (!cell_is_active_gravity(ci, r->e)) return;

  if (ci->nodeID != engine_rank)
    error("Non-local cell in long-range gravity task!");

  /* Check multipole has been drifted */
  if (ci->grav.ti_old_multipole < r->e->ti_current)
    cell_drift_multipole(ci, r->e);

  /* Find this cell's top-level (great-)parent */
  struct cell *top = ci;
  while (top->parent != NULL) top = top->parent;

  /* If we have a nested cell the true top level cell where we defined the
   * interactions is the background void parent cell. */
  /* TODO: This is a bit of a hack, we should actually have the top pointers set
   * properly during void tree splitting. */
  while (top->void_parent != NULL) top = top->void_parent->top;

  /* Call the appropriate interaction function based on the type of the
   * cell in question. */
  if (periodic) {
    switch (ci->type) {

      case cell_type_regular:
        runner_do_grav_long_range_uniform_periodic(r, ci, top);
        break;
      case cell_type_zoom:
        runner_do_grav_long_range_zoom_periodic(r, ci, top);
        break;
      case cell_type_bkg:
        runner_do_grav_long_range_zoom_periodic(r, ci, top);
        break;
      default:
        error("Unknown cell type in long-range gravity task!");
    }
  } else {

    switch (ci->type) {

      case cell_type_regular:
        runner_do_grav_long_range_uniform_non_periodic(r, ci, top);
        break;
      case cell_type_zoom:
        runner_do_grav_long_range_zoom_non_periodic(r, ci, top);
        break;
      case cell_type_bkg:
        runner_do_grav_long_range_zoom_non_periodic(r, ci, top);
        break;
      default:
        error("Unknown cell type in long-range gravity task!");
    }
  }

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
  /* Count the number of mesh interactions if using the mesh. */
  if (periodic && s->with_zoom_region) {
    runner_count_mesh_interactions_zoom(r, ci, top);
  } else if (periodic) {
    runner_count_mesh_interactions_uniform(r, ci, top);
  }
#endif

  if (timer) TIMER_TOC(timer_dograv_long_range);
}
