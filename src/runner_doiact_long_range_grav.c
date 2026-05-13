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

#ifdef SWIFT_DEBUG_CHECKS
enum runner_debug_mesh_count_origin {
  runner_debug_mesh_count_origin_uniform_pair = 0,
  runner_debug_mesh_count_origin_uniform_self = 1,
  runner_debug_mesh_count_origin_zoom_pair = 2,
  runner_debug_mesh_count_origin_zoom_self = 3,
};

struct runner_debug_mesh_attachment_entry {
  unsigned long long top_ptr_key;
  unsigned long long recipient_ptr_key;
  unsigned long long source_ptr_key;
  unsigned long long first_super_ptr_key;
  unsigned long long first_ci_ptr_key;
  unsigned long long first_cj_ptr_key;
  int first_origin;
  int seen_count;
  int duplicate_reported;
  const char *first_attachment_case;
};

static int runner_debug_mesh_zoom_top_repeat_report_count = 0;
static int runner_debug_mesh_zoom_top_repeat_limit_reported = 0;

#define runner_debug_mesh_zoom_top_repeat_report_max 32

struct runner_debug_mesh_zoom_top_entry {
  unsigned long long top_ptr_key;
  int seen_count;
};

#define runner_debug_mesh_zoom_top_table_size (1 << 14)
#define runner_debug_mesh_zoom_top_probe_limit 16

static struct runner_debug_mesh_zoom_top_entry
    runner_debug_mesh_zoom_top_table[runner_debug_mesh_zoom_top_table_size];

static int runner_debug_mesh_attachment_report_count = 0;
static int runner_debug_mesh_attachment_limit_reported = 0;

#define runner_debug_mesh_attachment_report_max 32
#define runner_debug_mesh_attachment_table_size (1 << 14)
#define runner_debug_mesh_attachment_probe_limit 16

static struct runner_debug_mesh_attachment_entry
    runner_debug_mesh_attachment_table[runner_debug_mesh_attachment_table_size];

static int runner_debug_mesh_overlap_report_count = 0;
static int runner_debug_mesh_overlap_limit_reported = 0;

#define runner_debug_mesh_overlap_report_max 32

static int runner_debug_mesh_payload_report_count = 0;
static int runner_debug_mesh_payload_limit_reported = 0;

#define runner_debug_mesh_payload_report_max 64

static void runner_debug_record_zoom_top_invocation(const struct cell *top,
                                                    const struct cell *ci) {

  if (ci != top) return;
  if (top->type != cell_type_bkg || top->subtype != cell_subtype_neighbour)
    return;

  const unsigned long long top_ptr_key = (unsigned long long)(uintptr_t)top;
  const unsigned long long hash = top_ptr_key * 11400714819323198485ull;

  for (int probe = 0; probe < runner_debug_mesh_zoom_top_probe_limit; probe++) {
    struct runner_debug_mesh_zoom_top_entry *slot =
        &runner_debug_mesh_zoom_top_table
            [(hash + (unsigned long long)probe) &
             (runner_debug_mesh_zoom_top_table_size - 1)];

    if (slot->top_ptr_key == 0ULL) {
      if (atomic_cas(&slot->top_ptr_key, 0ULL, top_ptr_key) == 0ULL) {
        slot->seen_count = 1;
        return;
      }
    }

    if (slot->top_ptr_key == top_ptr_key) {
      const int previous_seen_count = atomic_inc(&slot->seen_count);

      if (previous_seen_count >= 1) {
        const int report_index =
            atomic_inc(&runner_debug_mesh_zoom_top_repeat_report_count);

        if (report_index < runner_debug_mesh_zoom_top_repeat_report_max) {
          message(
              "zoom-mesh repeated top invocation: top_ptr=%p top_cellID=%llu "
              "top_type=%s/%s top_loc=(%.3f,%.3f,%.3f) ci_ptr=%p count=%d",
              (void *)top, top->cellID, cellID_names[top->type],
              subcellID_names[top->subtype], top->loc[0], top->loc[1],
              top->loc[2], (void *)ci, previous_seen_count + 1);
        } else if (report_index == runner_debug_mesh_zoom_top_repeat_report_max &&
                   atomic_cas(&runner_debug_mesh_zoom_top_repeat_limit_reported,
                              0, 1) == 0) {
          message("zoom-mesh repeated top invocation reporting capped at %d lines",
                  runner_debug_mesh_zoom_top_repeat_report_max);
        }
      }

      return;
    }
  }
}

static void runner_debug_record_mesh_attachment(
    const struct cell *super, const struct cell *ci, const struct cell *cj,
    const struct cell *recipient, const struct cell *source,
    const enum runner_debug_mesh_count_origin origin,
    const char *attachment_case) {

  if (recipient->top->type != cell_type_bkg ||
      recipient->top->subtype != cell_subtype_neighbour)
    return;

  const int payload_report_index = atomic_inc(&runner_debug_mesh_payload_report_count);
  if (payload_report_index < runner_debug_mesh_payload_report_max) {
    message(
        "mesh-attachment payload: top_ptr=%p top_cellID=%llu recipient=%p (%s/%s depth=%d count=%d mpole_gparts=%lld) "
        "source=%p (%s/%s depth=%d count=%d mpole_gparts=%lld contains_zoom=%d top=%p top_type=%s/%s) "
        "super=%p ci=%p cj=%p attach=%s origin=%d",
        (void *)recipient->top, recipient->top->cellID, (void *)recipient,
        cellID_names[recipient->type], subcellID_names[recipient->subtype],
        recipient->depth, recipient->grav.count,
        recipient->grav.multipole->m_pole.num_gpart, (void *)source,
        cellID_names[source->type], subcellID_names[source->subtype],
        source->depth, source->grav.count, source->grav.multipole->m_pole.num_gpart,
        source->contains_zoom_cells, (void *)source->top,
        cellID_names[source->top->type], subcellID_names[source->top->subtype],
        (void *)super, (void *)ci, (void *)cj, attachment_case, (int)origin);
  } else if (payload_report_index == runner_debug_mesh_payload_report_max &&
             atomic_cas(&runner_debug_mesh_payload_limit_reported, 0, 1) == 0) {
    message("mesh-attachment payload reporting capped at %d lines",
            runner_debug_mesh_payload_report_max);
  }

  const unsigned long long top_ptr_key =
      (unsigned long long)(uintptr_t)recipient->top;
  const unsigned long long recipient_ptr_key =
      (unsigned long long)(uintptr_t)recipient;
  const unsigned long long source_ptr_key =
      (unsigned long long)(uintptr_t)source;
  const unsigned long long hash =
      top_ptr_key * 11400714819323198485ull ^
      recipient_ptr_key * 14029467366897019727ull ^
      source_ptr_key * 1609587929392839161ull;

  for (int probe = 0; probe < runner_debug_mesh_attachment_probe_limit; probe++) {
    struct runner_debug_mesh_attachment_entry *slot =
        &runner_debug_mesh_attachment_table
            [(hash + (unsigned long long)probe) &
             (runner_debug_mesh_attachment_table_size - 1)];

    if (slot->top_ptr_key == 0ULL) {
      if (atomic_cas(&slot->top_ptr_key, 0ULL, top_ptr_key) == 0ULL) {
        slot->recipient_ptr_key = recipient_ptr_key;
        slot->source_ptr_key = source_ptr_key;
        slot->first_super_ptr_key = (unsigned long long)(uintptr_t)super;
        slot->first_ci_ptr_key = (unsigned long long)(uintptr_t)ci;
        slot->first_cj_ptr_key = (unsigned long long)(uintptr_t)cj;
        slot->first_origin = (int)origin;
        slot->seen_count = 1;
        slot->first_attachment_case = attachment_case;
        return;
      }
    }

    if (slot->top_ptr_key == top_ptr_key &&
        slot->recipient_ptr_key == recipient_ptr_key &&
        slot->source_ptr_key != source_ptr_key) {

      const struct cell *other_source = (const struct cell *)(uintptr_t)slot->source_ptr_key;

      if (cell_contains_progeny(other_source, source) ||
          cell_contains_progeny(source, other_source)) {

        const int report_index = atomic_inc(&runner_debug_mesh_overlap_report_count);

        if (report_index < runner_debug_mesh_overlap_report_max) {
          message(
              "mesh-attachment overlap(ptr): top_ptr=%p top_cellID=%llu recipient=%p (%s/%s depth=%d) "
              "first[source=%p (%s/%s depth=%d) super=%p ci=%p cj=%p attach=%s origin=%d] "
              "second[source=%p (%s/%s depth=%d) super=%p ci=%p cj=%p attach=%s origin=%d]",
              (void *)recipient->top, recipient->top->cellID, (void *)recipient,
              cellID_names[recipient->type], subcellID_names[recipient->subtype],
              recipient->depth, (void *)other_source,
              cellID_names[other_source->type], subcellID_names[other_source->subtype],
              other_source->depth, (void *)slot->first_super_ptr_key,
              (void *)slot->first_ci_ptr_key, (void *)slot->first_cj_ptr_key,
              slot->first_attachment_case, slot->first_origin, (void *)source,
              cellID_names[source->type], subcellID_names[source->subtype],
              source->depth, (void *)super, (void *)ci, (void *)cj,
              attachment_case, (int)origin);
        } else if (report_index == runner_debug_mesh_overlap_report_max &&
                   atomic_cas(&runner_debug_mesh_overlap_limit_reported, 0, 1) == 0) {
          message("mesh-attachment overlap(ptr) reporting capped at %d lines",
                  runner_debug_mesh_overlap_report_max);
        }
      }
    }

    if (slot->top_ptr_key == top_ptr_key &&
        slot->recipient_ptr_key == recipient_ptr_key &&
        slot->source_ptr_key == source_ptr_key) {
      const int previous_seen_count = atomic_inc(&slot->seen_count);

      if (previous_seen_count == 1 &&
          atomic_cas(&slot->duplicate_reported, 0, 1) == 0) {
        const int report_index =
            atomic_inc(&runner_debug_mesh_attachment_report_count);

        if (report_index < runner_debug_mesh_attachment_report_max) {
          message(
              "mesh-attachment duplicate(ptr): top_ptr=%p top_cellID=%llu "
              "first[super=%p ci=%p cj=%p attach=%s origin=%d] "
              "duplicate[super=%p ci=%p cj=%p attach=%s origin=%d] "
              "recipient=%p (%s/%s depth=%d) source=%p (%s/%s depth=%d)",
              (void *)recipient->top, recipient->top->cellID,
              (void *)slot->first_super_ptr_key, (void *)slot->first_ci_ptr_key,
              (void *)slot->first_cj_ptr_key, slot->first_attachment_case,
              slot->first_origin, (void *)super, (void *)ci, (void *)cj,
              attachment_case, (int)origin, (void *)recipient,
              cellID_names[recipient->type], subcellID_names[recipient->subtype],
              recipient->depth, (void *)source, cellID_names[source->type],
              subcellID_names[source->subtype], source->depth);
        } else if (report_index == runner_debug_mesh_attachment_report_max &&
                   atomic_cas(&runner_debug_mesh_attachment_limit_reported, 0,
                              1) == 0) {
          message("mesh-attachment duplicate(ptr) reporting capped at %d lines",
                  runner_debug_mesh_attachment_report_max);
        }
      }

      return;
    }
  }
}

void runner_debug_dump_mesh_attachments_for_top(const struct cell *top) {

  if (top == NULL) return;

  const unsigned long long top_ptr_key = (unsigned long long)(uintptr_t)top;
  int dumped = 0;

  for (int i = 0; i < runner_debug_mesh_attachment_table_size; i++) {
    const struct runner_debug_mesh_attachment_entry *slot =
        &runner_debug_mesh_attachment_table[i];

    if (slot->top_ptr_key != top_ptr_key) continue;

    const struct cell *recipient =
        (const struct cell *)(uintptr_t)slot->recipient_ptr_key;
    const struct cell *source = (const struct cell *)(uintptr_t)slot->source_ptr_key;
    const struct cell *super =
        (const struct cell *)(uintptr_t)slot->first_super_ptr_key;
    const struct cell *ci = (const struct cell *)(uintptr_t)slot->first_ci_ptr_key;
    const struct cell *cj = (const struct cell *)(uintptr_t)slot->first_cj_ptr_key;

    message(
        "mesh-attachment dump: top_ptr=%p top_cellID=%llu recipient=%p (%s/%s depth=%d count=%d mpole_gparts=%lld) "
        "source=%p (%s/%s depth=%d count=%d mpole_gparts=%lld contains_zoom=%d top=%p top_type=%s/%s) "
        "first[super=%p ci=%p cj=%p attach=%s origin=%d seen=%d]",
        (void *)top, top->cellID, (void *)recipient, cellID_names[recipient->type],
        subcellID_names[recipient->subtype], recipient->depth,
        recipient->grav.count, recipient->grav.multipole->m_pole.num_gpart,
        (void *)source, cellID_names[source->type],
        subcellID_names[source->subtype], source->depth, source->grav.count,
        source->grav.multipole->m_pole.num_gpart, source->contains_zoom_cells,
        (void *)source->top, cellID_names[source->top->type],
        subcellID_names[source->top->subtype], (void *)super, (void *)ci,
        (void *)cj, slot->first_attachment_case, slot->first_origin,
        slot->seen_count);
    dumped++;
  }

  message("mesh-attachment dump summary: top_ptr=%p top_cellID=%llu entries=%d",
          (void *)top, top->cellID, dumped);
}

static INLINE void runner_record_mesh_attachment(
    struct cell *super, struct cell *ci, struct cell *cj,
    struct cell *recipient, struct cell *source,
    const enum runner_debug_mesh_count_origin origin,
    const char *attachment_case) {

  runner_accumulate_interaction(recipient->grav.multipole, source->grav.multipole);
  runner_debug_add_tensor_origin_count(&recipient->grav.multipole->pot, source,
                                       source->grav.multipole->m_pole.num_gpart,
                                       runner_debug_tensor_origin_mesh);
  runner_debug_record_mesh_attachment(super, ci, cj, recipient, source, origin,
                                      attachment_case);
}
#else
static INLINE void runner_record_mesh_attachment(
    struct cell *super, struct cell *ci, struct cell *cj,
    struct cell *recipient, struct cell *source, const int origin,
    const char *attachment_case) {

  runner_accumulate_interaction(recipient->grav.multipole, source->grav.multipole);
}
#endif

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
                                          struct cell *cj,
                                          const enum runner_debug_mesh_count_origin
                                              origin) {

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
                                  "super_eq_ci");
  } else if (cell_contains_progeny(ci, super)) {
    runner_record_mesh_attachment(super, ci, cj, super, cj, origin,
                                  "super_in_ci");
  } else if (cell_contains_progeny(super, ci)) {
    runner_record_mesh_attachment(super, ci, cj, ci, cj, origin,
                                  "ci_in_super");
  }

  /* Handle the symmetric case for self interactions */
  if (is_self) {
    if (super == cj) {
      runner_record_mesh_attachment(super, ci, cj, super, ci, origin,
                                    "self_super_eq_cj");
    } else if (cell_contains_progeny(cj, super)) {
      runner_record_mesh_attachment(super, ci, cj, super, ci, origin,
                                    "self_super_in_cj");
    } else if (cell_contains_progeny(super, cj)) {
      runner_record_mesh_attachment(super, ci, cj, cj, ci, origin,
                                    "self_cj_in_super");
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
                                                          struct space *s,
                                                          const enum runner_debug_mesh_count_origin
                                                              origin) {

  /* No self interactions here */
  if (ci == cj) {
    return;
  }

  struct engine *e = s->e;

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
          runner_count_mesh_interaction(c, cpi, cpj, origin);
          continue;
        }

        /* Can we use M-M for this pair? */
        if (cell_can_use_pair_mm(cpi, cpj, e, s, /*use_rebuild_data=*/1,
                                 /*is_tree_walk=*/1,
                                 /*periodic boundaries*/ s->periodic,
                                 /*use_mesh*/ s->periodic)) {
          /* This would be handled by a M-M task, nothing to count */
          continue;
        }

        /* We would create real tasks, so recurse to find mesh interactions */
        runner_count_mesh_interactions_pair_recursive(c, cpi, cpj, s, origin);
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
                                                          struct space *s,
                                                          const enum runner_debug_mesh_count_origin
                                                              origin) {

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
      runner_count_mesh_interactions_self_recursive(c, ci->progeny[k], s,
                                                    origin);
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
          runner_count_mesh_interaction(c, cpj, cpk, origin);
          continue;
        }

        /* Check which cell contains c to decide how to recurse. */
        const int cpj_contains_c = (cpj == c || cell_contains_progeny(cpj, c));
        const int cpk_contains_c = (cpk == c || cell_contains_progeny(cpk, c));

        /* Recurse as a pair interaction, ensuring we always pass the cell
         * containing c as the first argument after c. */
        if (cpj_contains_c) {
          runner_count_mesh_interactions_pair_recursive(c, cpj, cpk, s, origin);
        } else if (cpk_contains_c) {
          runner_count_mesh_interactions_pair_recursive(c, cpk, cpj, s, origin);
        } else {
          runner_count_mesh_interactions_pair_recursive(c, cpj, cpk, s, origin);
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
  runner_count_mesh_interactions_self_recursive(
      ci, top, s, runner_debug_mesh_count_origin_uniform_self);

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
      runner_count_mesh_interaction(ci, top, cj,
                                    runner_debug_mesh_count_origin_uniform_pair);
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
    runner_count_mesh_interactions_pair_recursive(
        ci, top, cj, s, runner_debug_mesh_count_origin_uniform_pair);
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
    struct cell *c, struct cell *ci, struct cell *cj, struct space *s,
    const enum runner_debug_mesh_count_origin origin) {

  struct engine *e = s->e;

  /* If neither cell is a void cell, use the normal pair recursive function */
  if (ci->subtype != cell_subtype_void && cj->subtype != cell_subtype_void) {
    runner_count_mesh_interactions_pair_recursive(c, ci, cj, s, origin);
    return;
  }

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
        /* Record the mesh interaction */
        runner_count_mesh_interaction(c, cpi, cpj, origin);
        continue;
      }

      /* Can we use M-M for this pair? */
      if (cell_can_use_pair_mm(cpi, cpj, e, s, /*use_rebuild_data=*/1,
                               /*is_tree_walk=*/1,
                               /*periodic boundaries*/ s->periodic,
                               /*use_mesh*/ s->periodic)) {
        /* M-M task handles this, nothing to count */
        continue;
      }

      /* Recurse to find more mesh interactions */
      runner_count_mesh_interactions_zoom_pair_recursive(c, cpi, cpj, s,
                                                         origin);
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
    struct cell *c, struct cell *ci, struct space *s,
    const enum runner_debug_mesh_count_origin origin) {

  struct engine *e = s->e;

  /* If not a void cell, use the normal self recursive function */
  if (ci->subtype != cell_subtype_void) {
    runner_count_mesh_interactions_self_recursive(c, ci, s, origin);
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

    runner_count_mesh_interactions_zoom_self_recursive(c, ci->progeny[k], s,
                                                       origin);
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
        runner_count_mesh_interaction(c, cpj, cpk, origin);
        continue;
      }

      /* Check which cell contains c to decide how to recurse. */
      const int cpj_contains_c = (cpj == c || cell_contains_progeny(cpj, c));
      const int cpk_contains_c = (cpk == c || cell_contains_progeny(cpk, c));

      /* Recurse as a pair interaction, ensuring we always pass the cell
       * containing c as the first argument after c. */
      if (cpj_contains_c) {
        runner_count_mesh_interactions_zoom_pair_recursive(c, cpj, cpk, s,
                                                           origin);
      } else if (cpk_contains_c) {
        runner_count_mesh_interactions_zoom_pair_recursive(c, cpk, cpj, s,
                                                           origin);
      } else {
        runner_count_mesh_interactions_zoom_pair_recursive(c, cpj, cpk, s,
                                                           origin);
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

#ifdef SWIFT_DEBUG_CHECKS
  runner_debug_record_zoom_top_invocation(top, ci);
#endif

  /* First, handle self interactions from the top-level cell.
   * This mirrors the self task created at the void level. */
  runner_count_mesh_interactions_zoom_self_recursive(
      ci, top, s, runner_debug_mesh_count_origin_zoom_self);

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
      runner_count_mesh_interaction(ci, top, cj,
                                    runner_debug_mesh_count_origin_zoom_pair);
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
    runner_count_mesh_interactions_zoom_pair_recursive(
        ci, top, cj, s, runner_debug_mesh_count_origin_zoom_pair);
  }
}

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
