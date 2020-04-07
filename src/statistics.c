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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "statistics.h"

/* Local headers. */
#include "cooling.h"
#include "engine.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "threadpool.h"

/**
 * @brief Information required to compute the statistics in the mapper
 */
struct index_data {
  /*! The space we play with */
  const struct space *s;

  /*! The #statistics aggregator to fill */
  struct statistics *stats;
};

/**
 * @brief Adds the content of one #statistics aggregator to another one.
 *
 * Performs a += b;
 *
 * @param a The #statistics structure to update.
 * @param b The #statistics structure to add to a.
 */
void stats_add(struct statistics *a, const struct statistics *b) {

  /* Add everything */
  a->E_kin += b->E_kin;
  a->E_int += b->E_int;
  a->E_pot_self += b->E_pot_self;
  a->E_pot_ext += b->E_pot_ext;
  a->E_rad += b->E_rad;
  a->entropy += b->entropy;
  a->mass += b->mass;
  a->mom[0] += b->mom[0];
  a->mom[1] += b->mom[1];
  a->mom[2] += b->mom[2];
  a->ang_mom[0] += b->ang_mom[0];
  a->ang_mom[1] += b->ang_mom[1];
  a->ang_mom[2] += b->ang_mom[2];
  a->centre_of_mass[0] += b->centre_of_mass[0];
  a->centre_of_mass[1] += b->centre_of_mass[1];
  a->centre_of_mass[2] += b->centre_of_mass[2];
}

/**
 * @brief Initialises a statistics aggregator to a valid state.
 *
 * @param s The #statistics aggregator to initialise
 */
void stats_init(struct statistics *s) {

  /* Zero everything */
  bzero(s, sizeof(struct statistics));

  /* Set the lock */
  lock_init(&s->lock);
}

/**
 * @brief The #threadpool mapper function used to collect statistics for #part.
 *
 * @param map_data Pointer to the particles.
 * @param nr_parts The number of particles in this chunk
 * @param extra_data The #statistics aggregator.
 */
void stats_collect_part_mapper(void *map_data, int nr_parts, void *extra_data) {

  /* Unpack the data */
  const struct index_data *data = (struct index_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = s->e;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int with_ext_grav = (e->policy & engine_policy_external_gravity);
  const int with_self_grav = (e->policy & engine_policy_self_gravity);
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;
  const double time = e->time;
  const struct part *restrict parts = (struct part *)map_data;
  const struct xpart *restrict xparts =
      s->xparts + (ptrdiff_t)(parts - s->parts);
  struct statistics *const global_stats = data->stats;

  /* Some information about the physical model */
  const struct external_potential *potential = e->external_potential;
  const struct phys_const *phys_const = e->physical_constants;
  const struct cosmology *cosmo = e->cosmology;

  /* Some constants from cosmology */
  const float a_inv = cosmo->a_inv;
  const float a_inv2 = a_inv * a_inv;

  /* Local accumulator */
  struct statistics stats;
  stats_init(&stats);

  /* Loop over particles */
  for (int k = 0; k < nr_parts; k++) {

    /* Get the particle */
    const struct part *p = &parts[k];
    const struct xpart *xp = &xparts[k];
    const struct gpart *gp = p->gpart;

    /* Get useful time variables */
    const integertime_t ti_beg =
        get_integer_time_begin(ti_current, p->time_bin);
    const integertime_t ti_end = get_integer_time_end(ti_current, p->time_bin);

    /* Get time-step since the last kick */
    float dt_kick_grav, dt_kick_hydro, dt_therm;
    if (with_cosmology) {
      dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_beg, ti_current);
      dt_kick_grav -=
          cosmology_get_grav_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
      dt_kick_hydro =
          cosmology_get_hydro_kick_factor(cosmo, ti_beg, ti_current);
      dt_kick_hydro -=
          cosmology_get_hydro_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
      dt_therm = cosmology_get_therm_kick_factor(cosmo, ti_beg, ti_current);
      dt_therm -=
          cosmology_get_therm_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
    } else {
      dt_kick_grav = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
      dt_kick_hydro = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
      dt_therm = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
    }

    float v[3];
    hydro_get_drifted_velocities(p, xp, dt_kick_hydro, dt_kick_grav, v);
    const double x[3] = {p->x[0], p->x[1], p->x[2]};
    const float m = hydro_get_mass(p);
    const float entropy = hydro_get_drifted_physical_entropy(p, cosmo);
    const float u_inter = hydro_get_drifted_physical_internal_energy(p, cosmo);

    /* Collect mass */
    stats.mass += m;

    /* Collect centre of mass */
    stats.centre_of_mass[0] += m * x[0];
    stats.centre_of_mass[1] += m * x[1];
    stats.centre_of_mass[2] += m * x[2];

    /* Collect momentum */
    stats.mom[0] += m * v[0];
    stats.mom[1] += m * v[1];
    stats.mom[2] += m * v[2];

    /* Collect angular momentum */
    stats.ang_mom[0] += m * (x[1] * v[2] - x[2] * v[1]);
    stats.ang_mom[1] += m * (x[2] * v[0] - x[0] * v[2]);
    stats.ang_mom[2] += m * (x[0] * v[1] - x[1] * v[0]);

    /* Collect energies. */
    stats.E_kin += 0.5f * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) *
                   a_inv2; /* 1/2 m a^2 \dot{r}^2 */
    stats.E_int += m * u_inter;
    stats.E_rad += cooling_get_radiated_energy(xp);
    if (gp != NULL && with_self_grav)
      stats.E_pot_self += 0.5f * m * gravity_get_physical_potential(gp, cosmo);
    if (gp != NULL && with_ext_grav)
      stats.E_pot_ext += m * external_gravity_get_potential_energy(
                                 time, potential, phys_const, gp);

    /* Collect entropy */
    stats.entropy += m * entropy;
  }

  /* Now write back to memory */
  if (lock_lock(&global_stats->lock) == 0) stats_add(global_stats, &stats);
  if (lock_unlock(&global_stats->lock) != 0) error("Failed to unlock stats.");
}

/**
 * @brief The #threadpool mapper function used to collect statistics for #gpart.
 *
 * @param map_data Pointer to the g-particles.
 * @param nr_gparts The number of g-particles in this chunk
 * @param extra_data The #statistics aggregator.
 */
void stats_collect_gpart_mapper(void *map_data, int nr_gparts,
                                void *extra_data) {

  /* Unpack the data */
  const struct index_data *data = (struct index_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = s->e;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int with_ext_grav = (e->policy & engine_policy_external_gravity);
  const int with_self_grav = (e->policy & engine_policy_self_gravity);
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;
  const double time = e->time;
  const struct gpart *restrict gparts = (struct gpart *)map_data;
  struct statistics *const global_stats = data->stats;

  /* Some information about the physical model */
  const struct external_potential *potential = e->external_potential;
  const struct phys_const *phys_const = e->physical_constants;
  const struct cosmology *cosmo = e->cosmology;

  /* Some constants from cosmology */
  const float a_inv = cosmo->a_inv;
  const float a_inv2 = a_inv * a_inv;

  /* Local accumulator */
  struct statistics stats;
  stats_init(&stats);

  /* Loop over particles */
  for (int k = 0; k < nr_gparts; k++) {

    /* Get the particle */
    const struct gpart *gp = &gparts[k];

    /* If the g-particle has a counterpart, ignore it */
    if (gp->id_or_neg_offset < 0) continue;

    /* Get useful variables */
    const integertime_t ti_beg =
        get_integer_time_begin(ti_current, gp->time_bin);
    const integertime_t ti_end = get_integer_time_end(ti_current, gp->time_bin);

    /* Get time-step since the last kick */
    float dt_kick_grav;
    if (with_cosmology) {
      dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_beg, ti_current);
      dt_kick_grav -=
          cosmology_get_grav_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
    } else {
      dt_kick_grav = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
    }

    /* Extrapolate velocities */
    const float v[3] = {gp->v_full[0] + gp->a_grav[0] * dt_kick_grav,
                        gp->v_full[1] + gp->a_grav[1] * dt_kick_grav,
                        gp->v_full[2] + gp->a_grav[2] * dt_kick_grav};

    const float m = gravity_get_mass(gp);
    const double x[3] = {gp->x[0], gp->x[1], gp->x[2]};

    /* Collect mass */
    stats.mass += m;

    /* Collect centre of mass */
    stats.centre_of_mass[0] += m * x[0];
    stats.centre_of_mass[1] += m * x[1];
    stats.centre_of_mass[2] += m * x[2];

    /* Collect momentum */
    stats.mom[0] += m * v[0];
    stats.mom[1] += m * v[1];
    stats.mom[2] += m * v[2];

    /* Collect angular momentum */
    stats.ang_mom[0] += m * (x[1] * v[2] - x[2] * v[1]);
    stats.ang_mom[1] += m * (x[2] * v[0] - x[0] * v[2]);
    stats.ang_mom[2] += m * (x[0] * v[1] - x[1] * v[0]);

    /* Collect energies. */
    stats.E_kin += 0.5f * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) *
                   a_inv2; /* 1/2 m a^2 \dot{r}^2 */
    if (with_self_grav)
      stats.E_pot_self += 0.5f * m * gravity_get_physical_potential(gp, cosmo);
    if (with_ext_grav)
      stats.E_pot_ext += m * external_gravity_get_potential_energy(
                                 time, potential, phys_const, gp);
  }

  /* Now write back to memory */
  if (lock_lock(&global_stats->lock) == 0) stats_add(global_stats, &stats);
  if (lock_unlock(&global_stats->lock) != 0) error("Failed to unlock stats.");
}

/**
 * @brief Collect physical statistics over all particles in a #space.
 *
 * @param s The #space to collect from.
 * @param stats The #statistics aggregator to fill.
 */
void stats_collect(const struct space *s, struct statistics *stats) {

  /* Prepare the data */
  struct index_data extra_data;
  extra_data.s = s;
  extra_data.stats = stats;

  /* Run parallel collection of statistics for parts */
  if (s->nr_parts > 0)
    threadpool_map(&s->e->threadpool, stats_collect_part_mapper, s->parts,
                   s->nr_parts, sizeof(struct part), threadpool_auto_chunk_size,
                   &extra_data);

  /* Run parallel collection of statistics for gparts */
  if (s->nr_gparts > 0)
    threadpool_map(&s->e->threadpool, stats_collect_gpart_mapper, s->gparts,
                   s->nr_gparts, sizeof(struct gpart),
                   threadpool_auto_chunk_size, &extra_data);
}

/**
 * @brief Apply final opetations on the #statistics.
 *
 * @param stats The #statistics to work on.
 */
void stats_finalize(struct statistics *stats) {

  if (stats->mass > 0.) {
    stats->centre_of_mass[0] /= stats->mass;
    stats->centre_of_mass[1] /= stats->mass;
    stats->centre_of_mass[2] /= stats->mass;
  }
}

/**
 * @brief Prints the content of a #statistics aggregator to a file
 *
 * @param file File to write to.
 * @param stats The #statistics object to write to the file
 * @param time The current physical time.
 */
void stats_print_to_file(FILE *file, const struct statistics *stats,
                         double time) {

  const double E_pot = stats->E_pot_self + stats->E_pot_ext;
  const double E_tot = stats->E_kin + stats->E_int + E_pot;

  fprintf(file,
          " %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e "
          "%14e %14e %14e %14e %14e %14e\n",
          time, stats->mass, E_tot, stats->E_kin, stats->E_int, E_pot,
          stats->E_pot_self, stats->E_pot_ext, stats->E_rad, stats->entropy,
          stats->mom[0], stats->mom[1], stats->mom[2], stats->ang_mom[0],
          stats->ang_mom[1], stats->ang_mom[2], stats->centre_of_mass[0],
          stats->centre_of_mass[1], stats->centre_of_mass[2]);
  fflush(file);
}

/* Extra stuff in MPI-land */
#ifdef WITH_MPI

/**
 * @brief MPI datatype corresponding to the #statistics structure.
 */
MPI_Datatype statistics_mpi_type;

/**
 * @brief MPI operator used for the reduction of #statistics structure.
 */
MPI_Op statistics_mpi_reduce_op;

/**
 * @brief MPI reduce operator for #statistics structures.
 */
void stats_add_mpi(void *in, void *inout, int *len, MPI_Datatype *datatype) {

  for (int i = 0; i < *len; ++i)
    stats_add(&((struct statistics *)inout)[0],
              &((const struct statistics *)in)[i]);
}

/**
 * @brief Registers MPI #statistics type and reduction function.
 */
void stats_create_mpi_type(void) {

  /* This is not the recommended way of doing this.
     One should define the structure field by field
     But as long as we don't do serialization via MPI-IO
     we don't really care.
     Also we would have to modify this function everytime something
     is added to the statistics structure. */
  if (MPI_Type_contiguous(sizeof(struct statistics) / sizeof(unsigned char),
                          MPI_BYTE, &statistics_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&statistics_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for statistics.");
  }

  /* Create the reduction operation */
  MPI_Op_create(stats_add_mpi, 1, &statistics_mpi_reduce_op);
}

void stats_free_mpi_type(void) {
  MPI_Type_free(&statistics_mpi_type);
  MPI_Op_free(&statistics_mpi_reduce_op);
}
#endif
