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
#include "cell.h"

/* Local headers. */
#include "active.h"
#include "engine.h"
#include "hydro.h"
#include "sink_properties.h"

/**
 * @brief Recursively update the pointer and counter for #spart after the
 * addition of a new particle.
 *
 * @param c The cell we are working on.
 * @param progeny_list The list of the progeny index at each level for the
 * leaf-cell where the particle was added.
 * @param main_branch Are we in a cell directly above the leaf where the new
 * particle was added?
 */
void cell_recursively_shift_sparts(struct cell *c,
                                   const int progeny_list[space_cell_maxdepth],
                                   const int main_branch) {
  if (c->split) {
    /* No need to recurse in progenies located before the insestion point */
    const int first_progeny = main_branch ? progeny_list[(int)c->depth] : 0;

    for (int k = first_progeny; k < 8; ++k) {
      if (c->progeny[k] != NULL)
        cell_recursively_shift_sparts(c->progeny[k], progeny_list,
                                      main_branch && (k == first_progeny));
    }
  }

  /* When directly above the leaf with the new particle: increase the particle
   * count */
  /* When after the leaf with the new particle: shift by one position */
  if (main_branch) {
    c->stars.count++;

    /* Indicate that the cell is not sorted and cancel the pointer sorting
     * arrays. */
    c->stars.sorted = 0;
    cell_free_stars_sorts(c);

  } else {
    c->stars.parts++;
  }
}

/**
 * @brief Recursively update the pointer and counter for #sink after the
 * addition of a new particle.
 *
 * @param c The cell we are working on.
 * @param progeny_list The list of the progeny index at each level for the
 * leaf-cell where the particle was added.
 * @param main_branch Are we in a cell directly above the leaf where the new
 * particle was added?
 */
void cell_recursively_shift_sinks(struct cell *c,
                                  const int progeny_list[space_cell_maxdepth],
                                  const int main_branch) {
  if (c->split) {
    /* No need to recurse in progenies located before the insestion point */
    const int first_progeny = main_branch ? progeny_list[(int)c->depth] : 0;

    for (int k = first_progeny; k < 8; ++k) {
      if (c->progeny[k] != NULL)
        cell_recursively_shift_sinks(c->progeny[k], progeny_list,
                                     main_branch && (k == first_progeny));
    }
  }

  /* When directly above the leaf with the new particle: increase the particle
   * count */
  /* When after the leaf with the new particle: shift by one position */
  if (main_branch) {
    c->sinks.count++;
  } else {
    c->sinks.parts++;
  }
}

/**
 * @brief Recursively update the pointer and counter for #gpart after the
 * addition of a new particle.
 *
 * @param c The cell we are working on.
 * @param progeny_list The list of the progeny index at each level for the
 * leaf-cell where the particle was added.
 * @param main_branch Are we in a cell directly above the leaf where the new
 * particle was added?
 */
void cell_recursively_shift_gparts(struct cell *c,
                                   const int progeny_list[space_cell_maxdepth],
                                   const int main_branch) {
  if (c->split) {
    /* No need to recurse in progenies located before the insestion point */
    const int first_progeny = main_branch ? progeny_list[(int)c->depth] : 0;

    for (int k = first_progeny; k < 8; ++k) {
      if (c->progeny[k] != NULL)
        cell_recursively_shift_gparts(c->progeny[k], progeny_list,
                                      main_branch && (k == first_progeny));
    }
  }

  /* When directly above the leaf with the new particle: increase the particle
   * count */
  /* When after the leaf with the new particle: shift by one position */
  if (main_branch) {
    c->grav.count++;
  } else {
    c->grav.parts++;
  }
}

/**
 * @brief "Add" a #spart in a given #cell.
 *
 * This function will add a #spart at the start of the current cell's array by
 * shifting all the #spart in the top-level cell by one position. All the
 * pointers and cell counts are updated accordingly.
 *
 * @param e The #engine.
 * @param c The leaf-cell in which to add the #spart.
 *
 * @return A pointer to the newly added #spart. The spart has a been zeroed
 * and given a position within the cell as well as set to the minimal active
 * time bin.
 */
struct spart *cell_add_spart(struct engine *e, struct cell *const c) {
  /* Perform some basic consitency checks */
  if (c->nodeID != engine_rank) error("Adding spart on a foreign node");
  if (c->stars.ti_old_part != e->ti_current) error("Undrifted cell!");
  if (c->split) error("Addition of spart performed above the leaf level");

  /* Progeny number at each level */
  int progeny[space_cell_maxdepth];
#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < space_cell_maxdepth; ++i) progeny[i] = -1;
#endif

  /* Get the top-level this leaf cell is in and compute the progeny indices at
     each level */
  struct cell *top = c;
  while (top->parent != NULL) {
    /* What is the progeny index of the cell? */
    for (int k = 0; k < 8; ++k) {
      if (top->parent->progeny[k] == top) {
        progeny[(int)top->parent->depth] = k;
      }
    }

    /* Check that the cell was indeed drifted to this point to avoid future
     * issues */
#ifdef SWIFT_DEBUG_CHECKS
    if (top->hydro.super != NULL && top->stars.count > 0 &&
        top->stars.ti_old_part != e->ti_current) {
      error("Cell had not been correctly drifted before star formation");
    }
#endif

    /* Climb up */
    top = top->parent;
  }

  /* Lock the top-level cell as we are going to operate on it */
  lock_lock(&top->stars.star_formation_lock);

  /* Are there any extra particles left? */
  if (top->stars.count == top->stars.count_total) {

    message("We ran out of free star particles!");

    /* Release the local lock before exiting. */
    if (lock_unlock(&top->stars.star_formation_lock) != 0)
      error("Failed to unlock the top-level cell.");

    atomic_inc(&e->forcerebuild);
    return NULL;
  }

  /* Number of particles to shift in order to get a free space. */
  const size_t n_copy = &top->stars.parts[top->stars.count] - c->stars.parts;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->stars.parts + n_copy > top->stars.parts + top->stars.count)
    error("Copying beyond the allowed range");
#endif

  if (n_copy > 0) {
    // MATTHIEU: This can be improved. We don't need to copy everything, just
    // need to swap a few particles.
    memmove(&c->stars.parts[1], &c->stars.parts[0],
            n_copy * sizeof(struct spart));

    /* Update the spart->gpart links (shift by 1) */
    for (size_t i = 0; i < n_copy; ++i) {

#ifdef SWIFT_DEBUG_CHECKS
      if (c->stars.parts[i + 1].gpart == NULL) {
        error("Incorrectly linked spart!");
      }
#endif
      c->stars.parts[i + 1].gpart->id_or_neg_offset--;
    }
  }

  /* Recursively shift all the stars to get a free spot at the start of the
   * current cell*/
  cell_recursively_shift_sparts(top, progeny, /* main_branch=*/1);

  /* Make sure the gravity will be recomputed for this particle in the next
   * step
   */
  struct cell *top2 = c;
  while (top2->parent != NULL) {
    top2->stars.ti_old_part = e->ti_current;
    top2 = top2->parent;
  }
  top2->stars.ti_old_part = e->ti_current;

  /* Release the lock */
  if (lock_unlock(&top->stars.star_formation_lock) != 0)
    error("Failed to unlock the top-level cell.");

  /* We now have an empty spart as the first particle in that cell */
  struct spart *sp = &c->stars.parts[0];
  bzero(sp, sizeof(struct spart));

  /* Give it a decent position */
  sp->x[0] = c->loc[0] + 0.5 * c->width[0];
  sp->x[1] = c->loc[1] + 0.5 * c->width[1];
  sp->x[2] = c->loc[2] + 0.5 * c->width[2];

  /* Set it to the current time-bin */
  sp->time_bin = e->min_active_bin;

#ifdef SWIFT_DEBUG_CHECKS
  /* Specify it was drifted to this point */
  sp->ti_drift = e->ti_current;
#endif

  /* Register that we used one of the free slots. */
  const size_t one = 1;
  atomic_sub(&e->s->nr_extra_sparts, one);

  return sp;
}

/**
 * @brief "Add" a #sink in a given #cell.
 *
 * This function will add a #sink at the start of the current cell's array by
 * shifting all the #sink in the top-level cell by one position. All the
 * pointers and cell counts are updated accordingly.
 *
 * @param e The #engine.
 * @param c The leaf-cell in which to add the #sink.
 *
 * @return A pointer to the newly added #sink. The sink has a been zeroed
 * and given a position within the cell as well as set to the minimal active
 * time bin.
 */
struct sink *cell_add_sink(struct engine *e, struct cell *const c) {
  /* Perform some basic consitency checks */
  if (c->nodeID != engine_rank) error("Adding sink on a foreign node");
  if (c->sinks.ti_old_part != e->ti_current) error("Undrifted cell!");
  if (c->split) error("Addition of sink performed above the leaf level");

  /* Progeny number at each level */
  int progeny[space_cell_maxdepth];
#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < space_cell_maxdepth; ++i) progeny[i] = -1;
#endif

  /* Get the top-level this leaf cell is in and compute the progeny indices at
     each level */
  struct cell *top = c;
  while (top->parent != NULL) {
    /* What is the progeny index of the cell? */
    for (int k = 0; k < 8; ++k) {
      if (top->parent->progeny[k] == top) {
        progeny[(int)top->parent->depth] = k;
      }
    }

    /* Check that the cell was indeed drifted to this point to avoid future
     * issues */
#ifdef SWIFT_DEBUG_CHECKS
    if (top->hydro.super != NULL && top->sinks.count > 0 &&
        top->sinks.ti_old_part != e->ti_current) {
      error("Cell had not been correctly drifted before sink formation");
    }
#endif

    /* Climb up */
    top = top->parent;
  }

  /* Lock the top-level cell as we are going to operate on it */
  lock_lock(&top->sinks.sink_formation_lock);

  /* Are there any extra particles left? */
  if (top->sinks.count == top->sinks.count_total) {

    error("We ran out of free sink particles!");

    /* Release the local lock before exiting. */
    if (lock_unlock(&top->sinks.sink_formation_lock) != 0)
      error("Failed to unlock the top-level cell.");

    atomic_inc(&e->forcerebuild);
    return NULL;
  }

  /* Number of particles to shift in order to get a free space. */
  const size_t n_copy = &top->sinks.parts[top->sinks.count] - c->sinks.parts;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->sinks.parts + n_copy > top->sinks.parts + top->sinks.count)
    error("Copying beyond the allowed range");
#endif

  if (n_copy > 0) {
    // MATTHIEU: This can be improved. We don't need to copy everything, just
    // need to swap a few particles.
    memmove(&c->sinks.parts[1], &c->sinks.parts[0],
            n_copy * sizeof(struct sink));

    /* Update the sink->gpart links (shift by 1) */
    for (size_t i = 0; i < n_copy; ++i) {

      // TODO: Matthieu figure out whether this is strictly needed
      /* Skip inhibited (swallowed) sink particles */
      if (sink_is_inhibited(&c->sinks.parts[i + 1], e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (c->sinks.parts[i + 1].gpart == NULL) {
        error("Incorrectly linked sink!");
      }
#endif
      c->sinks.parts[i + 1].gpart->id_or_neg_offset--;
    }
  }

  /* Recursively shift all the sinks to get a free spot at the start of the
   * current cell*/
  cell_recursively_shift_sinks(top, progeny, /* main_branch=*/1);

  /* Make sure the gravity will be recomputed for this particle in the next
   * step
   */
  struct cell *top2 = c;
  while (top2->parent != NULL) {
    top2->sinks.ti_old_part = e->ti_current;
    top2 = top2->parent;
  }
  top2->sinks.ti_old_part = e->ti_current;

  /* Release the lock */
  if (lock_unlock(&top->sinks.sink_formation_lock) != 0)
    error("Failed to unlock the top-level cell.");

  /* We now have an empty spart as the first particle in that cell */
  struct sink *sp = &c->sinks.parts[0];
  bzero(sp, sizeof(struct sink));

  /* Give it a decent position */
  sp->x[0] = c->loc[0] + 0.5 * c->width[0];
  sp->x[1] = c->loc[1] + 0.5 * c->width[1];
  sp->x[2] = c->loc[2] + 0.5 * c->width[2];

  /* Set it to the current time-bin */
  sp->time_bin = e->min_active_bin;

#ifdef SWIFT_DEBUG_CHECKS
  /* Specify it was drifted to this point */
  sp->ti_drift = e->ti_current;
#endif

  /* Register that we used one of the free slots. */
  const size_t one = 1;
  atomic_sub(&e->s->nr_extra_sinks, one);

  return sp;
}

/**
 * @brief "Add" a #gpart in a given #cell.
 *
 * This function will add a #gpart at the start of the current cell's array by
 * shifting all the #gpart in the top-level cell by one position. All the
 * pointers and cell counts are updated accordingly.
 *
 * @param e The #engine.
 * @param c The leaf-cell in which to add the #gpart.
 *
 * @return A pointer to the newly added #gpart. The gpart has a been zeroed
 * and given a position within the cell as well as set to the minimal active
 * time bin.
 */
struct gpart *cell_add_gpart(struct engine *e, struct cell *c) {
  /* Perform some basic consitency checks */
  if (c->nodeID != engine_rank) error("Adding gpart on a foreign node");
  if (c->grav.ti_old_part != e->ti_current) error("Undrifted cell!");
  if (c->split) error("Addition of gpart performed above the leaf level");

  struct space *s = e->s;

  /* Progeny number at each level */
  int progeny[space_cell_maxdepth];
#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < space_cell_maxdepth; ++i) progeny[i] = -1;
#endif

  /* Get the top-level this leaf cell is in and compute the progeny indices at
     each level */
  struct cell *top = c;
  while (top->parent != NULL) {
    /* What is the progeny index of the cell? */
    for (int k = 0; k < 8; ++k) {
      if (top->parent->progeny[k] == top) {
        progeny[(int)top->parent->depth] = k;
      }
    }

    /* Check that the cell was indeed drifted to this point to avoid future
     * issues */
#ifdef SWIFT_DEBUG_CHECKS
    if (top->grav.super != NULL && top->grav.count > 0 &&
        top->grav.ti_old_part != e->ti_current) {
      error("Cell had not been correctly drifted before adding a gpart");
    }
#endif

    /* Climb up */
    top = top->parent;
  }

  /* Lock the top-level cell as we are going to operate on it */
  lock_lock(&top->grav.star_formation_lock);

  /* Are there any extra particles left? */
  if (top->grav.count == top->grav.count_total) {

    message("We ran out of free gravity particles!");

    /* Release the local lock before exiting. */
    if (lock_unlock(&top->grav.star_formation_lock) != 0)
      error("Failed to unlock the top-level cell.");

    atomic_inc(&e->forcerebuild);
    return NULL;
  }

  /* Number of particles to shift in order to get a free space. */
  const size_t n_copy = &top->grav.parts[top->grav.count] - c->grav.parts;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->grav.parts + n_copy > top->grav.parts + top->grav.count)
    error("Copying beyond the allowed range");
#endif

  if (n_copy > 0) {
    // MATTHIEU: This can be improved. We don't need to copy everything, just
    // need to swap a few particles.
    memmove(&c->grav.parts[1], &c->grav.parts[0],
            n_copy * sizeof(struct gpart));

    /* Update the gpart->spart links (shift by 1) */
    struct gpart *gparts = c->grav.parts;
    for (size_t i = 0; i < n_copy; ++i) {

      /* Skip inhibited particles */
      if (gpart_is_inhibited(&c->grav.parts[i + 1], e)) continue;

      if (gparts[i + 1].type == swift_type_gas) {
        s->parts[-gparts[i + 1].id_or_neg_offset].gpart++;
      } else if (gparts[i + 1].type == swift_type_stars) {
        s->sparts[-gparts[i + 1].id_or_neg_offset].gpart++;
      } else if (gparts[i + 1].type == swift_type_sink) {
        s->sinks[-gparts[i + 1].id_or_neg_offset].gpart++;
      } else if (gparts[i + 1].type == swift_type_black_hole) {
        s->bparts[-gparts[i + 1].id_or_neg_offset].gpart++;
      }
    }
  }

  /* Recursively shift all the gpart to get a free spot at the start of the
   * current cell*/
  cell_recursively_shift_gparts(top, progeny, /* main_branch=*/1);

  /* Make sure the gravity will be recomputed for this particle in the next
   * step
   */
  struct cell *top2 = c;
  while (top2->parent != NULL) {
    top2->grav.ti_old_part = e->ti_current;
    top2 = top2->parent;
  }
  top2->grav.ti_old_part = e->ti_current;

  /* Release the lock */
  if (lock_unlock(&top->grav.star_formation_lock) != 0)
    error("Failed to unlock the top-level cell.");

  /* We now have an empty gpart as the first particle in that cell */
  struct gpart *gp = &c->grav.parts[0];
  bzero(gp, sizeof(struct gpart));

  /* Give it a decent position */
  gp->x[0] = c->loc[0] + 0.5 * c->width[0];
  gp->x[1] = c->loc[1] + 0.5 * c->width[1];
  gp->x[2] = c->loc[2] + 0.5 * c->width[2];

  /* Set it to the current time-bin */
  gp->time_bin = e->min_active_bin;

#ifdef SWIFT_DEBUG_CHECKS
  /* Specify it was drifted to this point */
  gp->ti_drift = e->ti_current;
#endif

  /* Register that we used one of the free slots. */
  const size_t one = 1;
  atomic_sub(&e->s->nr_extra_gparts, one);

  return gp;
}

/**
 * @brief "Remove" a gas particle from the calculation.
 *
 * The particle is inhibited and will officially be removed at the next
 * rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param p The #part to remove.
 * @param xp The extended data of the particle to remove.
 */
void cell_remove_part(const struct engine *e, struct cell *c, struct part *p,
                      struct xpart *xp) {
  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  /* Don't remove a particle twice */
  if (p->time_bin == time_bin_inhibited) return;

  /* Mark the particle as inhibited */
  p->time_bin = time_bin_inhibited;
  /* Mark the RT time bin as inhibited as well,
   * so part_is_rt_active() checks work as intended */
  p->rt_time_data.time_bin = time_bin_inhibited;

  /* Mark the gpart as inhibited and stand-alone */
  if (p->gpart) {
    p->gpart->time_bin = time_bin_inhibited;
    p->gpart->id_or_neg_offset = 1;
    p->gpart->type = swift_type_dark_matter;
  }

  /* Update the space-wide counters */
  const size_t one = 1;
  atomic_add(&e->s->nr_inhibited_parts, one);
  if (p->gpart) {
    atomic_add(&e->s->nr_inhibited_gparts, one);
  }

  /* Un-link the part */
  p->gpart = NULL;
}

/**
 * @brief "Remove" a gravity particle from the calculation.
 *
 * The particle is inhibited and will officially be removed at the next
 * rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param gp The #gpart to remove.
 */
void cell_remove_gpart(const struct engine *e, struct cell *c,
                       struct gpart *gp) {

  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  /* Don't remove a particle twice */
  if (gp->time_bin == time_bin_inhibited) return;

  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  if (gp->type == swift_type_dark_matter_background)
    error("Can't remove a DM background particle!");

  /* Mark the particle as inhibited */
  gp->time_bin = time_bin_inhibited;

  /* Update the space-wide counters */
  const size_t one = 1;
  atomic_add(&e->s->nr_inhibited_gparts, one);
}

/**
 * @brief "Remove" a star particle from the calculation.
 *
 * The particle is inhibited and will officially be removed at the next
 * rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param sp The #spart to remove.
 */
void cell_remove_spart(const struct engine *e, struct cell *c,
                       struct spart *sp) {
  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  /* Don't remove a particle twice */
  if (sp->time_bin == time_bin_inhibited) return;

  /* Mark the particle as inhibited and stand-alone */
  sp->time_bin = time_bin_inhibited;
  if (sp->gpart) {
    sp->gpart->time_bin = time_bin_inhibited;
    sp->gpart->id_or_neg_offset = 1;
    sp->gpart->type = swift_type_dark_matter;
  }

  /* Update the space-wide counters */
  const size_t one = 1;
  atomic_add(&e->s->nr_inhibited_sparts, one);
  if (sp->gpart) {
    atomic_add(&e->s->nr_inhibited_gparts, one);
  }

  /* Un-link the spart */
  sp->gpart = NULL;
}

/**
 * @brief "Remove" a black hole particle from the calculation.
 *
 * The particle is inhibited and will officially be removed at the next
 * rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param bp The #bpart to remove.
 */
void cell_remove_bpart(const struct engine *e, struct cell *c,
                       struct bpart *bp) {

  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  /* Don't remove a particle twice */
  if (bp->time_bin == time_bin_inhibited) return;

  /* Mark the particle as inhibited and stand-alone */
  bp->time_bin = time_bin_inhibited;
  if (bp->gpart) {
    bp->gpart->time_bin = time_bin_inhibited;
    bp->gpart->id_or_neg_offset = 1;
    bp->gpart->type = swift_type_dark_matter;
  }

  /* Update the space-wide counters */
  const size_t one = 1;
  atomic_add(&e->s->nr_inhibited_bparts, one);
  if (bp->gpart) {
    atomic_add(&e->s->nr_inhibited_gparts, one);
  }

  /* Un-link the bpart */
  bp->gpart = NULL;
}

/**
 * @brief "Remove" a sink particle from the calculation.
 *
 * The particle is inhibited and will officially be removed at the next
 * rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param sp The #sink to remove.
 */
void cell_remove_sink(const struct engine *e, struct cell *c,
                      struct sink *sink) {
  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  /* Don't remove a particle twice */
  if (sink->time_bin == time_bin_inhibited) return;

  /* Mark the particle as inhibited and stand-alone */
  sink->time_bin = time_bin_inhibited;
  if (sink->gpart) {
    sink->gpart->time_bin = time_bin_inhibited;
    sink->gpart->id_or_neg_offset = 1;
    sink->gpart->type = swift_type_dark_matter;
  }

  /* Update the space-wide counters */
  const size_t one = 1;
  atomic_add(&e->s->nr_inhibited_sinks, one);
  if (sink->gpart) {
    atomic_add(&e->s->nr_inhibited_gparts, one);
  }

  /* Un-link the sink */
  sink->gpart = NULL;
}

/**
 * @brief "Remove" a gas particle from the calculation and convert its gpart
 * friend to a dark matter particle.
 *
 * Note that the #part is not destroyed. The pointer is still valid
 * after this call and the properties of the #part are not altered
 * apart from the time-bin and #gpart pointer.
 * The particle is inhibited and will officially be removed at the next
 * rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param p The #part to remove.
 * @param xp The extended data of the particle to remove.
 *
 * @return Pointer to the #gpart the #part has become. It carries the
 * ID of the #part and has a dark matter type.
 */
struct gpart *cell_convert_part_to_gpart(const struct engine *e, struct cell *c,
                                         struct part *p, struct xpart *xp) {
  /* Quick cross-checks */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  if (p->gpart == NULL)
    error("Trying to convert part without gpart friend to dark matter!");

  /* Get a handle */
  struct gpart *gp = p->gpart;

  /* Mark the particle as inhibited */
  p->time_bin = time_bin_inhibited;

  /* Un-link the part */
  p->gpart = NULL;

  /* Mark the gpart as dark matter */
  gp->type = swift_type_dark_matter;
  gp->id_or_neg_offset = p->id;

#ifdef SWIFT_DEBUG_CHECKS
  gp->ti_kick = p->ti_kick;
#endif

  /* Update the space-wide counters */
  atomic_inc(&e->s->nr_inhibited_parts);

  return gp;
}

/**
 * @brief "Remove" a spart particle from the calculation and convert its gpart
 * friend to a dark matter particle.
 *
 * Note that the #spart is not destroyed. The pointer is still valid
 * after this call and the properties of the #spart are not altered
 * apart from the time-bin and #gpart pointer.
 * The particle is inhibited and will officially be removed at the next
 * rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param sp The #spart to remove.
 *
 * @return Pointer to the #gpart the #spart has become. It carries the
 * ID of the #spart and has a dark matter type.
 */
struct gpart *cell_convert_spart_to_gpart(const struct engine *e,
                                          struct cell *c, struct spart *sp) {
  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  if (sp->gpart == NULL)
    error("Trying to convert spart without gpart friend to dark matter!");

  /* Get a handle */
  struct gpart *gp = sp->gpart;

  /* Mark the particle as inhibited */
  sp->time_bin = time_bin_inhibited;

  /* Un-link the spart */
  sp->gpart = NULL;

  /* Mark the gpart as dark matter */
  gp->type = swift_type_dark_matter;
  gp->id_or_neg_offset = sp->id;

#ifdef SWIFT_DEBUG_CHECKS
  gp->ti_kick = sp->ti_kick;
#endif

  /* Update the space-wide counters */
  atomic_inc(&e->s->nr_inhibited_sparts);

  return gp;
}

/**
 * @brief "Remove" a #part from a #cell and replace it with a #spart
 * connected to the same #gpart.
 *
 * Note that the #part is not destroyed. The pointer is still valid
 * after this call and the properties of the #part are not altered
 * apart from the time-bin and #gpart pointer.
 * The particle is inhibited and will officially be removed at the next
 * rebuild.
 *
 * @param e The #engine.
 * @param c The #cell from which to remove the #part.
 * @param p The #part to remove (must be inside c).
 * @param xp The extended data of the #part.
 *
 * @return A fresh #spart with the same ID, position, velocity and
 * time-bin as the original #part.
 */
struct spart *cell_convert_part_to_spart(struct engine *e, struct cell *c,
                                         struct part *p, struct xpart *xp) {
  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  if (p->gpart == NULL)
    error("Trying to convert part without gpart friend to star!");

  /* Create a fresh (empty) spart */
  struct spart *sp = cell_add_spart(e, c);

  /* Did we run out of free spart slots? */
  if (sp == NULL) return NULL;

  /* Copy over the distance since rebuild */
  sp->x_diff[0] = xp->x_diff[0];
  sp->x_diff[1] = xp->x_diff[1];
  sp->x_diff[2] = xp->x_diff[2];

  /* Destroy the gas particle and get it's gpart friend */
  struct gpart *gp = cell_convert_part_to_gpart(e, c, p, xp);

  /* Assign the ID back */
  sp->id = gp->id_or_neg_offset;
  gp->type = swift_type_stars;

  /* Re-link things */
  sp->gpart = gp;
  gp->id_or_neg_offset = -(sp - e->s->sparts);

  /* Synchronize clocks */
  gp->time_bin = sp->time_bin;

  /* Synchronize masses, positions and velocities */
  sp->mass = gp->mass;
  sp->x[0] = gp->x[0];
  sp->x[1] = gp->x[1];
  sp->x[2] = gp->x[2];
  sp->v[0] = gp->v_full[0];
  sp->v[1] = gp->v_full[1];
  sp->v[2] = gp->v_full[2];

#ifdef SWIFT_DEBUG_CHECKS
  sp->ti_kick = gp->ti_kick;
  gp->ti_drift = sp->ti_drift;
#endif

  /* Set a smoothing length */
  sp->h = p->h;

  /* Here comes the Sun! */
  return sp;
}

/**
 * @brief Add a new #spart based on a #part and link it to a new #gpart.
 * The part and xpart are not changed.
 *
 * @param e The #engine.
 * @param c The #cell from which to remove the #part.
 * @param p The #part to remove (must be inside c).
 * @param xp The extended data of the #part.
 *
 * @return A fresh #spart with a different ID, but same position,
 * velocity and time-bin as the original #part.
 */
struct spart *cell_spawn_new_spart_from_part(struct engine *e, struct cell *c,
                                             const struct part *p,
                                             const struct xpart *xp) {
  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't spawn a particle in a foreign cell.");

  if (p->gpart == NULL)
    error("Trying to create a new spart from a part without gpart friend!");

  /* Create a fresh (empty) spart */
  struct spart *sp = cell_add_spart(e, c);

  /* Did we run out of free spart slots? */
  if (sp == NULL) return NULL;

  /* Copy over the distance since rebuild */
  sp->x_diff[0] = xp->x_diff[0];
  sp->x_diff[1] = xp->x_diff[1];
  sp->x_diff[2] = xp->x_diff[2];

  /* Create a new gpart */
  struct gpart *gp = cell_add_gpart(e, c);

  /* Did we run out of free gpart slots? */
  if (gp == NULL) {
    /* Remove the particle created */
    cell_remove_spart(e, c, sp);
    return NULL;
  }

  /* Copy the gpart */
  *gp = *p->gpart;

  /* Assign the ID. */
  sp->id = space_get_new_unique_id(e->s);
  gp->type = swift_type_stars;

  /* Re-link things */
  sp->gpart = gp;
  gp->id_or_neg_offset = -(sp - e->s->sparts);

  /* Synchronize clocks */
  gp->time_bin = sp->time_bin;

  /* Synchronize masses, positions and velocities */
  sp->mass = hydro_get_mass(p);
  sp->x[0] = p->x[0];
  sp->x[1] = p->x[1];
  sp->x[2] = p->x[2];
  sp->v[0] = xp->v_full[0];
  sp->v[1] = xp->v_full[1];
  sp->v[2] = xp->v_full[2];

#ifdef SWIFT_DEBUG_CHECKS
  sp->ti_kick = p->ti_kick;
  sp->ti_drift = p->ti_drift;
#endif

  /* Set a smoothing length */
  sp->h = p->h;

  /* Here comes the Sun! */
  return sp;
}

/**
 * @brief "Remove" a #part from a #cell and replace it with a #sink
 * connected to the same #gpart.
 *
 * Note that the #part is not destroyed. The pointer is still valid
 * after this call and the properties of the #part are not altered
 * apart from the time-bin and #gpart pointer.
 * The particle is inhibited and will officially be removed at the next
 * rebuild.
 *
 * @param e The #engine.
 * @param c The #cell from which to remove the #part.
 * @param p The #part to remove (must be inside c).
 * @param xp The extended data of the #part.
 *
 * @return A fresh #sink with the same ID, position, velocity and
 * time-bin as the original #part.
 */
struct sink *cell_convert_part_to_sink(struct engine *e, struct cell *c,
                                       struct part *p, struct xpart *xp) {
  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  if (p->gpart == NULL)
    error("Trying to convert part without gpart friend to sink!");

  /* Create a fresh (empty) sink */
  struct sink *sp = cell_add_sink(e, c);

  /* Did we run out of free sink slots? */
  if (sp == NULL) return NULL;

  /* Copy over the distance since rebuild */
  sp->x_diff[0] = xp->x_diff[0];
  sp->x_diff[1] = xp->x_diff[1];
  sp->x_diff[2] = xp->x_diff[2];

  /* Destroy the gas particle and get it's gpart friend */
  struct gpart *gp = cell_convert_part_to_gpart(e, c, p, xp);

  /* Assign the ID back */
  sp->id = p->id;
  gp->type = swift_type_sink;

  /* Re-link things */
  sp->gpart = gp;
  gp->id_or_neg_offset = -(sp - e->s->sinks);

  /* Synchronize clocks */
  gp->time_bin = sp->time_bin;

  /* Synchronize masses, positions and velocities */
  sp->mass = gp->mass;
  sp->x[0] = gp->x[0];
  sp->x[1] = gp->x[1];
  sp->x[2] = gp->x[2];
  sp->v[0] = gp->v_full[0];
  sp->v[1] = gp->v_full[1];
  sp->v[2] = gp->v_full[2];

#ifdef SWIFT_DEBUG_CHECKS
  sp->ti_kick = gp->ti_kick;
  gp->ti_drift = sp->ti_drift;

  message("A new sink (%lld) is born !", sp->id);
#endif

  /* Here comes the Sink! */
  return sp;
}

/**
 * @brief Spawn a new #spart along with its related #gpart from a #sink.
 * The sink is not changed.
 *
 * @param e The #engine.
 * @param c The #cell from which to remove the #part.
 * @param s The #sink particle that spawns the #spart.
 *
 * @return A fresh #spart with a different ID, but same position,
 * velocity and time-bin as the original #sink.
 */
struct spart *cell_spawn_new_spart_from_sink(struct engine *e, struct cell *c,
                                             const struct sink *s) {
  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't spawn a particle in a foreign cell.");

  if (s->gpart == NULL)
    error("Trying to create a new spart from a part without gpart friend!");

  /* Create a fresh (empty) spart */
  struct spart *sp = cell_add_spart(e, c);

  /* Did we run out of free spart slots? */
  if (sp == NULL) return NULL;

  /* Copy over the distance since rebuild */
  sp->x_diff[0] = s->x_diff[0];
  sp->x_diff[1] = s->x_diff[1];
  sp->x_diff[2] = s->x_diff[2];

  /* Create a new gpart */
  struct gpart *gp = cell_add_gpart(e, c);

  /* Did we run out of free gpart slots? */
  if (gp == NULL) {
    /* Remove the particle created */
    cell_remove_spart(e, c, sp);
    return NULL;
  }

  /* Copy the gpart */
  *gp = *s->gpart;

  /* Assign the ID. */
  sp->id = space_get_new_unique_id(e->s);
  gp->type = swift_type_stars;

  /* Re-link things */
  sp->gpart = gp;
  gp->id_or_neg_offset = -(sp - e->s->sparts);

  /* Synchronize clocks */
  gp->time_bin = sp->time_bin;

  /* Synchronize masses, positions and velocities */
  sp->mass = s->mass;
  sp->x[0] = s->x[0];
  sp->x[1] = s->x[1];
  sp->x[2] = s->x[2];
  sp->v[0] = s->v[0];
  sp->v[1] = s->v[1];
  sp->v[2] = s->v[2];

#ifdef SWIFT_DEBUG_CHECKS
  sp->ti_kick = s->ti_kick;
  sp->ti_drift = s->ti_drift;
#endif

  /* Set a smoothing length */
  sp->h = s->r_cut;

  /* Here comes the Sun! */
  return sp;
}
