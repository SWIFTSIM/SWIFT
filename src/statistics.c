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

/* Config parameters. */
#include <config.h>

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
#include "black_holes.h"
#include "black_holes_io.h"
#include "chemistry.h"
#include "cooling.h"
#include "engine.h"
#include "error.h"
#include "gravity_io.h"
#include "hydro_io.h"
#include "mhd_io.h"
#include "potential.h"
#include "sink_io.h"
#include "stars_io.h"
#include "threadpool.h"

/**
 * @brief Information required to compute the statistics in the mapper
 */
struct space_index_data {
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
  a->gas_mass += b->gas_mass;
  a->dm_mass += b->dm_mass;
  a->sink_mass += b->sink_mass;
  a->star_mass += b->star_mass;
  a->bh_mass += b->bh_mass;
  a->bh_subgrid_mass += b->bh_subgrid_mass;
  a->gas_Z_mass += b->gas_Z_mass;
  a->star_Z_mass += b->star_Z_mass;
  a->bh_Z_mass += b->bh_Z_mass;
  a->bh_accretion_rate += b->bh_accretion_rate;
  a->bh_accreted_mass += b->bh_accreted_mass;
  a->bh_bolometric_luminosity += b->bh_bolometric_luminosity;
  a->bh_jet_power += b->bh_jet_power;
  a->mom[0] += b->mom[0];
  a->mom[1] += b->mom[1];
  a->mom[2] += b->mom[2];
  a->ang_mom[0] += b->ang_mom[0];
  a->ang_mom[1] += b->ang_mom[1];
  a->ang_mom[2] += b->ang_mom[2];
  a->centre_of_mass[0] += b->centre_of_mass[0];
  a->centre_of_mass[1] += b->centre_of_mass[1];
  a->centre_of_mass[2] += b->centre_of_mass[2];
  a->gas_H_mass += b->gas_H_mass;
  a->gas_H2_mass += b->gas_H2_mass;
  a->gas_HI_mass += b->gas_HI_mass;
  a->gas_He_mass += b->gas_He_mass;
  a->E_mag += b->E_mag;
  a->divB_error += b->divB_error;
  a->H_cross += b->H_cross;
  a->H_mag += b->H_mag;
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
  const struct space_index_data *data = (struct space_index_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = s->e;
  const int with_ext_grav = (e->policy & engine_policy_external_gravity);
  const int with_self_grav = (e->policy & engine_policy_self_gravity);
  const double time = e->time;
  const struct part *const parts = (struct part *)map_data;
  const struct xpart *const xparts = s->xparts + (ptrdiff_t)(parts - s->parts);
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

    /* Ignore non-existing particles */
    if (p->time_bin == time_bin_inhibited ||
        p->time_bin == time_bin_not_created)
      continue;

    /* Get position and velocity */
    double x[3];
    float v[3];
    convert_part_pos(e, p, xp, x);
    convert_part_vel(e, p, xp, v);

    const float m = hydro_get_mass(p);
    const float entropy = hydro_get_drifted_physical_entropy(p, cosmo);
    const float u_inter = hydro_get_drifted_physical_internal_energy(p, cosmo);

    /* Collect mass */
    stats.gas_mass += m;

    /* Collect metal mass */
    stats.gas_Z_mass += chemistry_get_total_metal_mass_for_stats(p);

#if defined(CHEMISTRY_EAGLE)
#if defined(COOLING_EAGLE) || defined(COOLING_PS2020)

    const struct unit_system *us = e->internal_units;
    const struct hydro_props *hydro_props = e->hydro_properties;
    const struct entropy_floor_properties *floor_props = e->entropy_floor;
    const struct cooling_function_data *cooling = e->cooling_func;

    /* Collect H and He species */
    const float H_mass_frac =
        p->chemistry_data.metal_mass_fraction[chemistry_element_H];
    const float He_mass_frac =
        p->chemistry_data.metal_mass_fraction[chemistry_element_He];
    const float H_mass = m * H_mass_frac;
    const float He_mass = m * He_mass_frac;
    const float HI_frac = cooling_get_particle_subgrid_HI_fraction(
        us, phys_const, cosmo, hydro_props, floor_props, cooling, p, xp);
    const float H2_frac = cooling_get_particle_subgrid_H2_fraction(
        us, phys_const, cosmo, hydro_props, floor_props, cooling, p, xp);

    const float HI_mass = H_mass * HI_frac;
    const float H2_mass = H_mass * H2_frac * 2.;

    stats.gas_H_mass += H_mass;
    stats.gas_HI_mass += HI_mass;
    stats.gas_H2_mass += H2_mass;
    stats.gas_He_mass += He_mass;

#endif
#endif

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

    /* Collect magnetic energy */
    stats.E_mag += mhd_get_magnetic_energy(p, xp);

    /* Collect helicity */
    stats.H_mag += mhd_get_magnetic_helicity(p, xp);
    stats.H_cross += mhd_get_cross_helicity(p, xp);

    /* Collect div B error */
    stats.divB_error += mhd_get_divB_error(p, xp);
  }

  /* Now write back to memory */
  if (lock_lock(&global_stats->lock) == 0) stats_add(global_stats, &stats);
  if (lock_unlock(&global_stats->lock) != 0) error("Failed to unlock stats.");
}

/**
 * @brief The #threadpool mapper function used to collect statistics for #spart.
 *
 * @param map_data Pointer to the particles.
 * @param nr_sparts The number of particles in this chunk
 * @param extra_data The #statistics aggregator.
 */
void stats_collect_spart_mapper(void *map_data, int nr_sparts,
                                void *extra_data) {

  /* Unpack the data */
  const struct space_index_data *data = (struct space_index_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = s->e;
  const int with_ext_grav = (e->policy & engine_policy_external_gravity);
  const int with_self_grav = (e->policy & engine_policy_self_gravity);
  const double time = e->time;
  const struct spart *const sparts = (struct spart *)map_data;
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
  for (int k = 0; k < nr_sparts; k++) {

    /* Get the particle */
    const struct spart *sp = &sparts[k];
    const struct gpart *gp = sp->gpart;

    /* Ignore non-existing particles */
    if (sp->time_bin == time_bin_inhibited ||
        sp->time_bin == time_bin_not_created)
      continue;

    /* Get position and velocity */
    double x[3];
    float v[3];
    convert_spart_pos(e, sp, x);
    convert_spart_vel(e, sp, v);

    const float m = sp->mass;

    /* Collect mass */
    stats.star_mass += m;

    /* Collect metal mass */
    stats.star_Z_mass += chemistry_get_star_total_metal_mass_for_stats(sp);

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
    if (gp != NULL && with_self_grav)
      stats.E_pot_self += 0.5f * m * gravity_get_physical_potential(gp, cosmo);
    if (gp != NULL && with_ext_grav)
      stats.E_pot_ext += m * external_gravity_get_potential_energy(
                                 time, potential, phys_const, gp);
  }

  /* Now write back to memory */
  if (lock_lock(&global_stats->lock) == 0) stats_add(global_stats, &stats);
  if (lock_unlock(&global_stats->lock) != 0) error("Failed to unlock stats.");
}

/**
 * @brief The #threadpool mapper function used to collect statistics for #sink.
 *
 * @param map_data Pointer to the particles.
 * @param nr_sinks The number of particles in this chunk
 * @param extra_data The #statistics aggregator.
 */
void stats_collect_sink_mapper(void *map_data, int nr_sinks, void *extra_data) {

  /* Unpack the data */
  const struct space_index_data *data = (struct space_index_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = s->e;
  const int with_ext_grav = (e->policy & engine_policy_external_gravity);
  const int with_self_grav = (e->policy & engine_policy_self_gravity);
  const double time = e->time;
  const struct sink *const sinks = (struct sink *)map_data;
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
  for (int k = 0; k < nr_sinks; k++) {

    /* Get the particle */
    const struct sink *sp = &sinks[k];
    const struct gpart *gp = sp->gpart;

    /* Ignore non-existing particles */
    if (sp->time_bin == time_bin_inhibited ||
        sp->time_bin == time_bin_not_created)
      continue;

    /* Get position and velocity */
    double x[3];
    float v[3];
    convert_sink_pos(e, sp, x);
    convert_sink_vel(e, sp, v);

    const float m = sp->mass;

    /* Collect mass */
    stats.star_mass += m;

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
    if (gp != NULL && with_self_grav)
      stats.E_pot_self += 0.5f * m * gravity_get_physical_potential(gp, cosmo);
    if (gp != NULL && with_ext_grav)
      stats.E_pot_ext += m * external_gravity_get_potential_energy(
                                 time, potential, phys_const, gp);
  }

  /* Now write back to memory */
  if (lock_lock(&global_stats->lock) == 0) stats_add(global_stats, &stats);
  if (lock_unlock(&global_stats->lock) != 0) error("Failed to unlock stats.");
}

/**
 * @brief The #threadpool mapper function used to collect statistics for #bpart.
 *
 * @param map_data Pointer to the particles.
 * @param nr_bparts The number of particles in this chunk
 * @param extra_data The #statistics aggregator.
 */
void stats_collect_bpart_mapper(void *map_data, int nr_bparts,
                                void *extra_data) {

  /* Unpack the data */
  const struct space_index_data *data = (struct space_index_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = s->e;
  const int with_ext_grav = (e->policy & engine_policy_external_gravity);
  const int with_self_grav = (e->policy & engine_policy_self_gravity);
  const double time = e->time;
  const struct bpart *const bparts = (struct bpart *)map_data;
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
  for (int k = 0; k < nr_bparts; k++) {

    /* Get the particle */
    const struct bpart *bp = &bparts[k];
    const struct gpart *gp = bp->gpart;

    /* Ignore non-existing particles */
    if (bp->time_bin == time_bin_inhibited ||
        bp->time_bin == time_bin_not_created)
      continue;

    /* Get position and velocity */
    double x[3];
    float v[3];
    convert_bpart_pos(e, bp, x);
    convert_bpart_vel(e, bp, v);

    const float m = bp->mass;

    /* Collect mass */
    stats.bh_mass += m;

    /* Collect subgrid mass */
    stats.bh_subgrid_mass += black_holes_get_subgrid_mass(bp);

    /* Collect metal mass */
    stats.bh_Z_mass += chemistry_get_bh_total_metal_mass_for_stats(bp);

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
    if (gp != NULL && with_self_grav)
      stats.E_pot_self += 0.5f * m * gravity_get_physical_potential(gp, cosmo);
    if (gp != NULL && with_ext_grav)
      stats.E_pot_ext += m * external_gravity_get_potential_energy(
                                 time, potential, phys_const, gp);

    /* Collect accretion data. */
    stats.bh_accretion_rate += black_holes_get_accretion_rate(bp);
    stats.bh_accreted_mass += black_holes_get_accreted_mass(bp);

    /* Collect bolometric luminosity and jet powers. */
    stats.bh_bolometric_luminosity +=
        black_holes_get_bolometric_luminosity(bp, phys_const);
    stats.bh_jet_power += black_holes_get_jet_power(bp, phys_const);
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
  const struct space_index_data *data = (struct space_index_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = s->e;
  const int with_ext_grav = (e->policy & engine_policy_external_gravity);
  const int with_self_grav = (e->policy & engine_policy_self_gravity);
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

    /* Ignore the hydro particles as they are already computed and skip
     * neutrinos */
    if (gp->type != swift_type_dark_matter &&
        gp->type != swift_type_dark_matter_background)
      continue;

    /* Ignore non-existing particles */
    if (gp->time_bin == time_bin_inhibited ||
        gp->time_bin == time_bin_not_created)
      continue;

    /* Get position and velocity */
    double x[3];
    float v[3];
    convert_gpart_pos(e, gp, x);
    convert_gpart_vel(e, gp, v);

    const float m = gravity_get_mass(gp);

    /* Collect mass */
    stats.dm_mass += m;

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
  struct space_index_data extra_data;
  extra_data.s = s;
  extra_data.stats = stats;

  /* Run parallel collection of statistics for parts */
  if (s->nr_parts > 0)
    threadpool_map(&s->e->threadpool, stats_collect_part_mapper, s->parts,
                   s->nr_parts, sizeof(struct part), threadpool_auto_chunk_size,
                   &extra_data);

  /* Run parallel collection of statistics for sparts */
  if (s->nr_sparts > 0)
    threadpool_map(&s->e->threadpool, stats_collect_spart_mapper, s->sparts,
                   s->nr_sparts, sizeof(struct spart),
                   threadpool_auto_chunk_size, &extra_data);

  /* Run parallel collection of statistics for sparts */
  if (s->nr_sinks > 0)
    threadpool_map(&s->e->threadpool, stats_collect_sink_mapper, s->sinks,
                   s->nr_sinks, sizeof(struct sink), threadpool_auto_chunk_size,
                   &extra_data);

  /* Run parallel collection of statistics for sparts */
  if (s->nr_bparts > 0)
    threadpool_map(&s->e->threadpool, stats_collect_bpart_mapper, s->bparts,
                   s->nr_bparts, sizeof(struct bpart),
                   threadpool_auto_chunk_size, &extra_data);

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

  stats->total_mass = stats->gas_mass + stats->dm_mass + stats->sink_mass +
                      stats->star_mass + stats->bh_mass;

  if (stats->total_mass > 0.) {
    stats->centre_of_mass[0] /= stats->total_mass;
    stats->centre_of_mass[1] /= stats->total_mass;
    stats->centre_of_mass[2] /= stats->total_mass;
  }
}

void stats_write_file_header(FILE *file, const struct unit_system *restrict us,
                             const struct phys_const *phys_const) {

  fprintf(file, "# Global statistics file\n");
  fprintf(file, "######################################################\n");
  fprintf(file, "# The quantities are all given in internal units!\n");
  fprintf(file, "#\n");
  fprintf(file, "# (0)  Simulation step\n");
  fprintf(file, "#      Unit = dimensionless\n");
  fprintf(file,
          "# (1)  Time since Big Bang (cosmological run), Time since start of "
          "the simulation (non-cosmological run).\n");
  fprintf(file, "#      Unit = %e s\n", us->UnitTime_in_cgs);
  fprintf(file, "#      Unit = %e yr \n", 1.f / phys_const->const_year);
  fprintf(file, "#      Unit = %e Myr \n", 1.f / phys_const->const_year / 1e6);
  fprintf(file, "# (2)  Scale factor\n");
  fprintf(file, "#      Unit = dimensionless\n");
  fprintf(file, "# (3)  Redshift\n");
  fprintf(file, "#      Unit = dimensionless\n");
  fprintf(file, "# (4)  Total mass in the simulation. \n");
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file,
          "# (5)  Total gas mass in the simulation (Particle type %d). \n",
          swift_type_gas);
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file,
          "# (6)  Total dark matter mass in the simulation (Particle type %d & "
          "%d). \n",
          swift_type_dark_matter, swift_type_dark_matter_background);
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file,
          "# (7)  Total sink mass in the simulation (Particle type %d). \n",
          swift_type_sink);
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file,
          "# (8)  Total stellar mass in the simulation (Particle type %d). \n",
          swift_type_stars);
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(
      file,
      "# (9)  Total black hole mass in the simulation (Particle type %d). \n",
      swift_type_black_hole);
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file, "# (10) Total metal mass in the gas phase. \n");
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file, "# (11) Total metal mass locked in stars. \n");
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file, "# (12) Total metal mass locked in black holes. \n");
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file, "# (13) Total kinetic energy (physical). \n");
  fprintf(file, "#      Unit = %e erg\n",
          units_cgs_conversion_factor(us, UNIT_CONV_ENERGY));
  fprintf(file,
          "# (14) Total internal (thermal) energy of the gas (physical). \n");
  fprintf(file, "#      Unit = %e erg\n",
          units_cgs_conversion_factor(us, UNIT_CONV_ENERGY));
  fprintf(file, "# (15) Total potential energy (physical). \n");
  fprintf(file, "#      Unit = %e erg\n",
          units_cgs_conversion_factor(us, UNIT_CONV_ENERGY));
  fprintf(file, "# (16) Total radiated energy of the gas (physical). \n");
  fprintf(file, "#      Unit = %e erg\n",
          units_cgs_conversion_factor(us, UNIT_CONV_ENERGY));
  fprintf(file, "# (17) Total gas entropy (physical). \n");
  fprintf(file, "#      Unit = %e gram**(%.3f) * cm**(%.3f) * s**(%.3f)\n",
          units_cgs_conversion_factor(us, UNIT_CONV_ENTROPY), 2.f - hydro_gamma,
          3.f * hydro_gamma - 1.f, -2.f);
  fprintf(
      file,
      "# (18) Comoving centre of mass of the simulation (x coordinate). \n");
  fprintf(file, "#      Unit = %e cm\n", us->UnitLength_in_cgs);
  fprintf(file, "#      Unit = %e pc\n", 1. / phys_const->const_parsec);
  fprintf(file, "#      Unit = %e Mpc\n", 1. / phys_const->const_parsec / 1e6);
  fprintf(
      file,
      "# (19) Comoving centre of mass of the simulation (y coordinate). \n");
  fprintf(file, "#      Unit = %e cm\n", us->UnitLength_in_cgs);
  fprintf(file, "#      Unit = %e pc\n", 1. / phys_const->const_parsec);
  fprintf(file, "#      Unit = %e Mpc\n", 1. / phys_const->const_parsec / 1e6);
  fprintf(
      file,
      "# (20) Comoving centre of mass of the simulation (z coordinate). \n");
  fprintf(file, "#      Unit = %e cm\n", us->UnitLength_in_cgs);
  fprintf(file, "#      Unit = %e pc\n", 1. / phys_const->const_parsec);
  fprintf(file, "#      Unit = %e Mpc\n", 1. / phys_const->const_parsec / 1e6);
  fprintf(file,
          "# (21) Comoving momentum of the simulation (x coordinate). \n");
  fprintf(file, "#      Unit = %e gram * cm * s**-1\n",
          units_cgs_conversion_factor(us, UNIT_CONV_MOMENTUM));
  fprintf(file,
          "# (22) Comoving momentum of the simulation (y coordinate). \n");
  fprintf(file, "#      Unit = %e gram * cm * s**-1\n",
          units_cgs_conversion_factor(us, UNIT_CONV_MOMENTUM));
  fprintf(file,
          "# (23) Comoving momentum of the simulation (z coordinate). \n");
  fprintf(file, "#      Unit = %e gram * cm * s**-1\n",
          units_cgs_conversion_factor(us, UNIT_CONV_MOMENTUM));
  fprintf(
      file,
      "# (24) Comoving angular momentum of the simulation (x coordinate). \n");
  fprintf(file, "#      Unit = %e gram * cm**2 * s**-1\n",
          units_cgs_conversion_factor(us, UNIT_CONV_ANGULAR_MOMENTUM));
  fprintf(
      file,
      "# (25) Comoving angular momentum of the simulation (y coordinate). \n");
  fprintf(file, "#      Unit = %e gram * cm**2 * s**-1\n",
          units_cgs_conversion_factor(us, UNIT_CONV_ANGULAR_MOMENTUM));
  fprintf(
      file,
      "# (26) Comoving angular momentum of the simulation (z coordinate). \n");
  fprintf(file, "#      Unit = %e gram * cm**2 * s**-1\n",
          units_cgs_conversion_factor(us, UNIT_CONV_ANGULAR_MOMENTUM));
  fprintf(file,
          "# (27) Sum of instantaneous accretion rate of all black holes in "
          "the simulation. \n");
  fprintf(file, "#      Unit = %e gram * s**-1\n",
          units_cgs_conversion_factor(us, UNIT_CONV_MASS_PER_UNIT_TIME));
  fprintf(file, "#      Unit = %e Msun/yr\n",
          phys_const->const_year / phys_const->const_solar_mass);
  fprintf(file,
          "# (28)  Total mass accreted by black holes in the simulation (not "
          "including mass accreted by progenitors that have merged in the "
          "current BHs). \n");
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file, "# (29)  Total black hole subgrid mass in the simulation. \n");
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file,
          "# (30) Total Hydrogen (all species) mass in the gas phase of the "
          "simulation. \n");
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file,
          "# (31) Total Molecular Hydrogen mass in the gas phase of the "
          "simulation. \n");
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file,
          "# (32) Total Atomic Hydrogen mass in the gas phase of the "
          "simulation. \n");
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file,
          "# (33) Total Helium (all species) mass in the gas phase of the "
          "simulation. \n");
  fprintf(file, "#      Unit = %e gram\n", us->UnitMass_in_cgs);
  fprintf(file, "#      Unit = %e Msun\n", 1. / phys_const->const_solar_mass);
  fprintf(file,
          "# (34) Total Magnetic Energy B2/(2*mu0) in the"
          "simulation. \n");
  fprintf(file, "#      Unit = %e erg\n",
          units_cgs_conversion_factor(us, UNIT_CONV_ENERGY));
  fprintf(file,
          "# (35) Total DivB error in the"
          "simulation. \n");
  fprintf(file, "#      Unit = dimensionless\n");
  fprintf(file,
          "# (36) Total Cross Helicity :: sum(V.B) in the"
          "simulation. \n");
  fprintf(file, "#      Unit = %e gram * cm * s**-3 * A**-1 \n",
          units_cgs_conversion_factor(us, UNIT_CONV_MAGNETIC_CROSS_HELICITY));
  fprintf(file,
          "# (37) Total Magnetic Helicity :: sum(A.B) in the"
          "simulation. \n");
  fprintf(file, "#      Unit = %e gram**2 * cm * s**-4 * A**-2\n",
          1. / units_cgs_conversion_factor(us, UNIT_CONV_MAGNETIC_HELICITY));
  fprintf(file, "# (38) Total bolometric luminosity of the BHs. \n");
  fprintf(file, "#      Unit = %e erg * s**-1\n",
          units_cgs_conversion_factor(us, UNIT_CONV_POWER));
  fprintf(file, "# (39) Total jet power of the BHs. \n");
  fprintf(file, "#      Unit = %e erg * s**-1\n",
          units_cgs_conversion_factor(us, UNIT_CONV_POWER));

  fprintf(file, "#\n");
  fprintf(
      file,
      "#%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s "
      "%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s "
      "%14s %14s %14s %14s %14s %14s %14s  %14s  %14s  %14s %14s  %14s \n",
      "(0)", "(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)", "(9)",
      "(10)", "(11)", "(12)", "(13)", "(14)", "(15)", "(16)", "(17)", "(18)",
      "(19)", "(20)", "(21)", "(22)", "(23)", "(24)", "(25)", "(26)", "(27)",
      "(28)", "(29)", "(30)", "(31)", "(32)", "(33)", "(34)", "(35)", "(36)",
      "(37)", "(38)", "(39)");
  fprintf(
      file,
      "#%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s "
      "%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s "
      "%14s %14s %14s %14s %14s %14s %14s  %14s  %14s  %14s %14s %14s \n",
      "Step", "Time", "a", "z", "Total mass", "Gas mass", "DM mass",
      "Sink mass", "Star mass", "BH mass", "Gas Z mass", "Star Z mass",
      "BH Z mass", "Kin. Energy", "Int. Energy", "Pot. energy", "Rad. energy",
      "Gas Entropy", "CoM x", "CoM y", "CoM z", "Mom. x", "Mom. y", "Mom. z",
      "Ang. mom. x", "Ang. mom. y", "Ang. mom. z", "BH acc. rate",
      "BH acc. mass", "BH sub. mass", "Gas H mass", "Gas H2 mass",
      "Gas HI mass", "Gas He mass", "Mag. Energy", "DivB err", "Cr. Helicity",
      "Mag. Helicity", "BH bol. lum.", "BH jet power");

  fflush(file);
}

/**
 * @brief Prints the content of a #statistics aggregator to a file
 *
 * @param file File to write to.
 * @param stats The #statistics object to write to the file
 * @param time The current physical time.
 * @param a The current scale-factor.
 * @param z The current redshift.
 * @param step The current time-step.
 */
void stats_write_to_file(FILE *file, const struct statistics *stats,
                         const double time, const double a, const double z,
                         const int step) {

  /* Compute the total potential */
  double E_pot = stats->E_pot_ext + stats->E_pot_self;

  /* Write to the file */
  fprintf(
      file,
      " %14d %14e %14.7f %14.7f %14e %14e %14e %14e %14e %14e %14e %14e %14e "
      "%14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e "
      "%14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e %14e\n",
      step, time, a, z, stats->total_mass, stats->gas_mass, stats->dm_mass,
      stats->sink_mass, stats->star_mass, stats->bh_mass, stats->gas_Z_mass,
      stats->star_Z_mass, stats->bh_Z_mass, stats->E_kin, stats->E_int, E_pot,
      stats->E_rad, stats->entropy, stats->centre_of_mass[0],
      stats->centre_of_mass[1], stats->centre_of_mass[2], stats->mom[0],
      stats->mom[1], stats->mom[2], stats->ang_mom[0], stats->ang_mom[1],
      stats->ang_mom[2], stats->bh_accretion_rate, stats->bh_accreted_mass,
      stats->bh_subgrid_mass, stats->gas_H_mass, stats->gas_H2_mass,
      stats->gas_HI_mass, stats->gas_He_mass, stats->E_mag, stats->divB_error,
      stats->H_cross, stats->H_mag, stats->bh_bolometric_luminosity,
      stats->bh_jet_power);

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
