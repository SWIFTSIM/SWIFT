/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_SINK_H
#define SWIFT_GEAR_SINK_H

#include <float.h>

/* Put pragma if gsl around here */
#ifdef HAVE_LIBGSL
#include <gsl/gsl_cdf.h>
#endif

/* Local includes */
#include "active.h"
#include "chemistry.h"
#include "cooling.h"
#include "feedback.h"
#include "minmax.h"
#include "random.h"
#include "sink_part.h"
#include "sink_properties.h"

/**
 * @brief Computes the time-step of a given sink particle.
 *
 * @param sp Pointer to the sink-particle data.
 */
__attribute__((always_inline)) INLINE static float sink_compute_timestep(
    const struct sink* const sp) {

  return FLT_MAX;
}

/**
 * @brief Update the target mass of the sink particle.
 *
 * @param sink the sink particle.
 * @param sink_props the sink properties to use.
 * @param phys_const the physical constants in internal units.
 * @param e The #engine
 * @param star_counter The star loop counter.
 */
INLINE static void sink_update_target_mass(struct sink* sink,
                                           const struct sink_props* sink_props,
                                           const struct engine* e,
                                           int star_counter) {

  float random_number = random_unit_interval_part_ID_and_index(
      sink->id, star_counter, e->ti_current, random_number_sink_formation);

  const struct feedback_props* feedback_props = e->feedback_props;

  /* Pick the correct table. (if only one table, threshold is < 0) */
  const float metal =
      chemistry_get_sink_total_iron_mass_fraction_for_feedback(sink);
  const float threshold = feedback_props->metallicity_max_first_stars;

  /* If metal < threshold, then the sink generates first star particles. */
  const int is_first_star = metal < threshold;
  const struct stellar_model* model;
  double minimal_discrete_mass_Msun;

  /* Take the correct values if your are a first star or not. */
  if (!is_first_star) /* (metal >= threshold)*/ {
    model = &feedback_props->stellar_model;
    minimal_discrete_mass_Msun = sink_props->minimal_discrete_mass_Msun;
  } else {
    model = &feedback_props->stellar_model_first_stars;
    minimal_discrete_mass_Msun =
        sink_props->minimal_discrete_mass_first_stars_Msun;
  }

  const struct initial_mass_function* imf = &model->imf;

  if (random_number < imf->sink_Pc) {
    /* We are dealing with the continous part of the IMF. */
    sink->target_mass_Msun = imf->stellar_particle_mass_Msun;
    sink->target_type = star_population_continuous_IMF;
  } else {
    /* We are dealing with the discrete part of the IMF. */
    random_number = random_unit_interval_part_ID_and_index(
        sink->id, star_counter + 1, e->ti_current,
        random_number_sink_formation);
    double m = initial_mass_function_sample_power_law(
        minimal_discrete_mass_Msun, imf->mass_max, imf->exp[imf->n_parts - 1],
        random_number);
    sink->target_mass_Msun = m;
    sink->target_type = single_star;
  }
}

/**
 * @brief Initialises the sink-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The #sink particle to act upon.
 * @param sink_props The properties of the sink particles scheme.
 * @param e The #engine
 */
__attribute__((always_inline)) INLINE static void sink_first_init_sink(
    struct sink* sp, const struct sink_props* sink_props,
    const struct engine* e) {

  sp->r_cut = sink_props->cut_off_radius;
  sp->time_bin = 0;

  sp->number_of_gas_swallows = 0;
  sp->number_of_direct_gas_swallows = 0;
  sp->number_of_sink_swallows = 0;
  sp->number_of_direct_sink_swallows = 0;
  sp->swallowed_angular_momentum[0] = 0.f;
  sp->swallowed_angular_momentum[1] = 0.f;
  sp->swallowed_angular_momentum[2] = 0.f;
  sp->n_stars = 0;

  sp->has_IMF_changed_from_popIII_to_popII = 0;

  sink_mark_sink_as_not_swallowed(&sp->merger_data);

  /* Bug fix: Setup the target mass for sink formation after reading the
     ICs. Otherwise sink->target_mass_Msun = 0.0 and a sink present in the IC
     spawn a star of mass 0.0... */
  sink_update_target_mass(sp, sink_props, e, 0);

  /* Initialize to the mass of the sink */
  sp->mass_tot_before_star_spawning = sp->mass;
}

/**
 * @brief Initialisation of particle data before the hydro density loop.
 * Note: during initalisation (space_init)
 *
 * @param p The #part to act upon.
 * @param sink_props The properties of the sink particles scheme.
 */
__attribute__((always_inline)) INLINE static void sink_init_part(
    struct part* restrict p, const struct sink_props* sink_props) {

  struct sink_part_data* cpd = &p->sink_data;

  if (sink_props->disable_sink_formation) {
    cpd->can_form_sink = 0;
  } else {
    cpd->can_form_sink = 1;
  }
  cpd->E_kin_neighbours = 0.f;
  cpd->E_int_neighbours = 0.f;
  cpd->E_rad_neighbours = 0.f;
  cpd->E_pot_self_neighbours = 0.f;
  cpd->E_pot_ext_neighbours = 0.f;
  cpd->E_mag_neighbours = 0.f;
  cpd->E_rot_neighbours[0] = 0.f;
  cpd->E_rot_neighbours[1] = 0.f;
  cpd->E_rot_neighbours[2] = 0.f;
  cpd->potential = 0.f;
  cpd->E_mec_bound = 0.f; /* Gravitationally bound particles will have
                             E_mec_bound < 0. This is checked before comparing
                             any other value with this one. So no need to put
                             it to the max of float. */
  cpd->is_overlapping_sink = 0;
}

/**
 * @brief Initialisation of sink particle data before sink loops.
 * Note: during initalisation (space_init_sinks)
 *
 * @param sp The #sink particle to act upon.
 */
__attribute__((always_inline)) INLINE static void sink_init_sink(
    struct sink* sp) {

  /* Reset to the mass of the sink */
  sp->mass_tot_before_star_spawning = sp->mass;

#ifdef DEBUG_INTERACTIONS_SINKS
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_SINKS; ++i)
    sp->ids_ngbs_accretion[i] = -1;
  sp->num_ngb_accretion = 0;

  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_SINKS; ++i)
    sp->ids_ngbs_merger[i] = -1;
  sp->num_ngb_merger = 0;

  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_SINKS; ++i)
    sp->ids_ngbs_formation[i] = -1;
  sp->num_ngb_formation = 0;
#endif
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param sp The #sink.
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void sink_predict_extra(
    struct sink* restrict sp, float dt_drift) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param sp The #sink particle.
 */
__attribute__((always_inline)) INLINE static void sink_reset_predicted_values(
    struct sink* restrict sp) {}

/**
 * @brief Kick the additional variables
 *
 * @param sp The #sink particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void sink_kick_extra(
    struct sink* sp, float dt) {}

/**
 * @brief Compute the rotational energy of the neighbouring gas particles.
 *
 * Note: This function must be used after having computed the rotational energy
 * per components, i.e. after sink_prepare_part_sink_formation().
 *
 * @param p The gas particle.
 *
 */
INLINE static double sink_compute_neighbour_rotation_energy_magnitude(
    const struct part* restrict p) {
  double E_rot_x = p->sink_data.E_rot_neighbours[0];
  double E_rot_y = p->sink_data.E_rot_neighbours[1];
  double E_rot_z = p->sink_data.E_rot_neighbours[2];
  double E_rot =
      sqrtf(E_rot_x * E_rot_x + E_rot_y * E_rot_y + E_rot_z * E_rot_z);
  return E_rot;
}

/**
 * @brief Retrieve the physical velocity divergence from the gas particle.
 *
 * @param p The gas particles.
 *
 */
INLINE static float sink_get_physical_div_v_from_part(
    const struct part* restrict p) {

  float div_v = 0.0;

  /* The implementation of div_v depends on the Hydro scheme. Furthermore, some
     add a Hubble flow term, some do not. We need to take care of this */
#ifdef SPHENIX_SPH
  /* SPHENIX is already including the Hubble flow. */
  div_v = hydro_get_div_v(p);
#elif GADGET2_SPH
  div_v = p->density.div_v;

  /* Add the missing term */
  div_v += hydro_dimension * cosmo->H;
#elif MINIMAL_SPH
  div_v = hydro_get_div_v(p);

  /* Add the missing term */
  div_v += hydro_dimension * cosmo->H;
#elif GASOLINE_SPH
  /* Copy the velocity divergence */
  div_v = (1. / 3.) * (p->viscosity.velocity_gradient[0][0] +
                       p->viscosity.velocity_gradient[1][1] +
                       p->viscosity.velocity_gradient[2][2]);
#elif HOPKINS_PU_SPH
  div_v = p->density.div_v;
#else
#error \
    "This scheme is not implemented. Note that Different scheme apply the Hubble flow in different places. Be careful about it."
#endif
  return div_v;
}

/**
 * @brief Calculate if the gas has the potential of becoming
 * a sink.
 *
 * Return 0 if no sink formation should occur.
 * Note: called in runner_do_sink_formation
 *
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sink_props the sink properties to use.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 * @param entropy_floor The entropy_floor properties.
 *
 */
INLINE static int sink_is_forming(
    const struct part* restrict p, const struct xpart* restrict xp,
    const struct sink_props* sink_props, const struct phys_const* phys_const,
    const struct cosmology* cosmo,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    const struct entropy_floor_properties* restrict entropy_floor) {

  /* the particle is not elligible */
  if (!p->sink_data.can_form_sink) return 0;

  const struct sink_part_data* sink_data = &p->sink_data;

  const float temperature_max = sink_props->maximal_temperature;
  const float temperature = cooling_get_temperature(phys_const, hydro_props, us,
                                                    cosmo, cooling, p, xp);

  const float density_threshold = sink_props->density_threshold;
  const float density = hydro_get_physical_density(p, cosmo);

  const float div_v = sink_get_physical_div_v_from_part(p);

  const float h = p->h;
  const float sink_cut_off_radius = sink_props->cut_off_radius;

  double E_grav = sink_data->E_pot_self_neighbours;
  double E_rot_neighbours = sink_compute_neighbour_rotation_energy_magnitude(p);
  double E_tot = sink_data->E_kin_neighbours + sink_data->E_int_neighbours +
                 E_grav + sink_data->E_mag_neighbours;

  /* Density and temperature criterion */
  if (density <= density_threshold || temperature >= temperature_max) {
    return 0;
  }

  /* Contracting gas criterion */
  if ((sink_props->sink_formation_contracting_gas_criterion) && (div_v > 0)) {
    return 0;
  }

  /* Smoothing length criterion */
  if ((sink_props->sink_formation_smoothing_length_criterion) &&
      (kernel_gamma * h >= sink_cut_off_radius)) {
    return 0;
  }

  /* Active neighbours criterion */
  /* This is checked on the fly in runner_do_sink_formation(). The part is
     flagged to not form sink through p->sink_data.can_form_sink */

  /* Jeans instability criterion */
  if ((sink_props->sink_formation_jeans_instability_criterion) &&
      (sink_data->E_int_neighbours >= 0.5f * fabs(E_grav))) {
    return 0;
  }

  if ((sink_props->sink_formation_jeans_instability_criterion) &&
      (sink_data->E_int_neighbours + E_rot_neighbours >= fabs(E_grav))) {
    return 0;
  }

  /* Bound state criterion */
  if ((sink_props->sink_formation_bound_state_criterion) && (E_tot >= 0)) {
    return 0;
  }

  /* Minimum of the potential criterion */
  /* Done in density loop. The gas is then flagged through
     sink_data.can_form_sink to not form sink. The check is done at the
     beginning. */

  /* Overlapping existing sinks criterion */
  if (sink_props->sink_formation_overlapping_sink_criterion &&
      sink_data->is_overlapping_sink) {
    return 0;
  }

  message("Gas particle %lld can form a sink !", p->id);
  return 1;
}

/**
 * @brief Decides whether a particle should be converted into a
 * sink or not.
 *
 * No SF should occur, so return 0.
 * Note: called in runner_do_sink_formation
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param sink_props The properties of the sink model.
 * @param e The #engine (for random numbers).
 * @param dt_sink The time-step of this particle
 * @return 1 if a conversion should be done, 0 otherwise.
 */
INLINE static int sink_should_convert_to_sink(
    const struct part* p, const struct xpart* xp,
    const struct sink_props* sink_props, const struct engine* e,
    const double dt_sink) {

  /* We do not use a stockastic approach.
   * Once elligible (sink_is_forming), the gas particle form a sink */

  return 1;
}

/**
 * @brief Copies the properties of the gas particle over to the
 * sink particle.
 *
 * @param p The gas particles.
 * @param xp The additional properties of the gas particles.
 * @param sink the new created #sink particle.
 * @param e The #engine.
 * @param sink_props The sink properties to use.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 * @param phys_const The physical constants in internal units.
 * @param hydro_props The hydro properties to use.
 * @param us The internal unit system.
 * @param cooling The cooling function to use.
 */
INLINE static void sink_copy_properties(
    const struct part* p, const struct xpart* xp, struct sink* sink,
    const struct engine* e, const struct sink_props* sink_props,
    const struct cosmology* cosmo, const int with_cosmology,
    const struct phys_const* phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling) {

  /* First initialisation */
  sink_init_sink(sink);

  /* Set a smoothing length */
  sink->r_cut = sink_props->cut_off_radius;

  /* Flag it as not swallowed */
  sink_mark_sink_as_not_swallowed(&sink->merger_data);

  /* Additional initialisation */
  sink->number_of_gas_swallows = 0;
  sink->number_of_direct_gas_swallows = 0;
  sink->number_of_sink_swallows = 0;
  sink->number_of_direct_sink_swallows = 0;
  sink->swallowed_angular_momentum[0] = 0.f;
  sink->swallowed_angular_momentum[1] = 0.f;
  sink->swallowed_angular_momentum[2] = 0.f;
  sink->n_stars = 0;
  sink->has_IMF_changed_from_popIII_to_popII = 0;

  /* setup the target mass for sink star formation */
  sink_update_target_mass(sink, sink_props, e, 0);

  /* Copy the chemistry properties */
  chemistry_copy_sink_properties(p, xp, sink);

  /* Note, we do not need to update sp->mass_tot_before_star_spawning because
     it is performed within the 'sink_init_sink()' function. */
}

/**
 * @brief Update the properties of a sink particles by swallowing
 * a gas particle.
 *
 * @param sp The #sink to update.
 * @param p The #part that is swallowed.
 * @param xp The #xpart that is swallowed.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sink_swallow_part(
    struct sink* sp, const struct part* p, const struct xpart* xp,
    const struct cosmology* cosmo) {

  /* Get the current dynamical masses */
  const float gas_mass = hydro_get_mass(p);
  const float sink_mass = sp->mass;

  /* store the mass of the sink part i */
  const float msp_old = sp->mass;

  /* Increase the dynamical mass of the sink. */
  sp->mass += gas_mass;
  sp->gpart->mass += gas_mass;

  /* Comoving and physical distance between the particles */
  const float dx[3] = {sp->x[0] - p->x[0], sp->x[1] - p->x[1],
                       sp->x[2] - p->x[2]};
  const float dx_physical[3] = {dx[0] * cosmo->a, dx[1] * cosmo->a,
                                dx[2] * cosmo->a};

  /* Relative velocity between the sink and the part */
  const float dv[3] = {sp->v[0] - p->v[0], sp->v[1] - p->v[1],
                       sp->v[2] - p->v[2]};

  const float a = cosmo->a;
  const float H = cosmo->H;
  const float a2H = a * a * H;

  /* Calculate the velocity with the Hubble flow */
  const float v_plus_H_flow[3] = {a2H * dx[0] + dv[0], a2H * dx[1] + dv[1],
                                  a2H * dx[2] + dv[2]};

  /* Compute the physical relative velocity between the particles */
  const float dv_physical[3] = {v_plus_H_flow[0] * cosmo->a_inv,
                                v_plus_H_flow[1] * cosmo->a_inv,
                                v_plus_H_flow[2] * cosmo->a_inv};

  /* Collect the swallowed angular momentum */
  sp->swallowed_angular_momentum[0] +=
      gas_mass *
      (dx_physical[1] * dv_physical[2] - dx_physical[2] * dv_physical[1]);
  sp->swallowed_angular_momentum[1] +=
      gas_mass *
      (dx_physical[2] * dv_physical[0] - dx_physical[0] * dv_physical[2]);
  sp->swallowed_angular_momentum[2] +=
      gas_mass *
      (dx_physical[0] * dv_physical[1] - dx_physical[1] * dv_physical[0]);

  /* Update the sink momentum */
  const float sink_mom[3] = {sink_mass * sp->v[0] + gas_mass * p->v[0],
                             sink_mass * sp->v[1] + gas_mass * p->v[1],
                             sink_mass * sp->v[2] + gas_mass * p->v[2]};

  sp->v[0] = sink_mom[0] / sp->mass;
  sp->v[1] = sink_mom[1] / sp->mass;
  sp->v[2] = sink_mom[2] / sp->mass;
  sp->gpart->v_full[0] = sp->v[0];
  sp->gpart->v_full[1] = sp->v[1];
  sp->gpart->v_full[2] = sp->v[2];

  /* Update the sink metal masses fraction */
  chemistry_add_part_to_sink(sp, p, msp_old);

  /* This sink swallowed a gas particle */
  sp->number_of_gas_swallows++;
  sp->number_of_direct_gas_swallows++;

  /* Update the total mass before star spawning */
  sp->mass_tot_before_star_spawning = sp->mass;

#ifdef SWIFT_DEBUG_CHECKS
  const float dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  message(
      "sink %lld swallow gas particle %lld. "
      "(Mass = %e, "
      "Delta_v = [%f, %f, %f] U_V, "
      "Delta_x = [%f, %f, %f] U_L, "
      "Delta_v_rad = %f)",
      sp->id, p->id, sp->mass, -dv[0], -dv[1], -dv[2], -dx[0], -dx[1], -dx[2],
      (dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2]) / dr);
#endif
}

/**
 * @brief Update the properties of a sink particles by swallowing
 * a sink particle.
 *
 * @param spi The #sink to update.
 * @param spj The #sink that is swallowed.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sink_swallow_sink(
    struct sink* spi, const struct sink* spj, const struct cosmology* cosmo) {

  /* Get the current dynamical masses */
  const float spi_dyn_mass = spi->mass;
  const float spj_dyn_mass = spj->mass;

  /* store the mass of the sink part i */
  const float mi_old = spi->mass;

  /* Increase the masses of the sink. */
  spi->mass += spj->mass;
  spi->gpart->mass += spj->mass;

  /* Collect the swallowed angular momentum */
  spi->swallowed_angular_momentum[0] += spj->swallowed_angular_momentum[0];
  spi->swallowed_angular_momentum[1] += spj->swallowed_angular_momentum[1];
  spi->swallowed_angular_momentum[2] += spj->swallowed_angular_momentum[2];

  /* Update the sink momentum */
  const float sink_mom[3] = {
      spi_dyn_mass * spi->v[0] + spj_dyn_mass * spj->v[0],
      spi_dyn_mass * spi->v[1] + spj_dyn_mass * spj->v[1],
      spi_dyn_mass * spi->v[2] + spj_dyn_mass * spj->v[2]};

  spi->v[0] = sink_mom[0] / spi->mass;
  spi->v[1] = sink_mom[1] / spi->mass;
  spi->v[2] = sink_mom[2] / spi->mass;
  spi->gpart->v_full[0] = spi->v[0];
  spi->gpart->v_full[1] = spi->v[1];
  spi->gpart->v_full[2] = spi->v[2];

  /* Update the sink metal masses fraction */
  chemistry_add_sink_to_sink(spi, spj, mi_old);

  /* This sink swallowed a sink particle */
  spi->number_of_sink_swallows++;
  spi->number_of_direct_sink_swallows++;

  /* Add all other swallowed particles swallowed by the swallowed sink */
  spi->number_of_sink_swallows += spj->number_of_sink_swallows;
  spi->number_of_gas_swallows += spj->number_of_gas_swallows;

  /* Add the stars spawned by the swallowed sink */
  spi->n_stars += spj->n_stars;

  /* Update the total mass before star spawning */
  spi->mass_tot_before_star_spawning = spi->mass;

  message("sink %lld swallow sink particle %lld. New mass: %e.", spi->id,
          spj->id, spi->mass);
}

/**
 * @brief Should the sink spawn a star particle?
 *
 * @param sink the sink particle.
 * @param e The #engine
 * @param sink_props The sink properties to use.
 * @param cosmo The cosmological parameters and properties.
 * @param with_cosmology If we run with cosmology.
 * @param phys_const The physical constants in internal units.
 * @param us The internal unit system.
 */
INLINE static int sink_spawn_star(struct sink* sink, const struct engine* e,
                                  const struct sink_props* sink_props,
                                  const struct cosmology* cosmo,
                                  const int with_cosmology,
                                  const struct phys_const* phys_const,
                                  const struct unit_system* restrict us) {

  if (sink->mass > sink->target_mass_Msun * phys_const->const_solar_mass)
    return 1;
  else
    return 0;
}

/**
 * @brief Give the #spart a new position.
 *
 * In GEAR: Positions are set by randomly sampling coordinates in an homogeneous
 * sphere centered on the #sink with radius  the sink's r_cut.
 *
 * @param e The #engine.
 * @param si The #sink generating a star.
 * @param sp The #spart generated.
 */
INLINE static void sink_star_formation_give_new_position(const struct engine* e,
                                                         struct sink* si,
                                                         struct spart* sp) {
#ifdef SWIFT_DEBUG_CHECKS
  if (si->x[0] != sp->x[0] || si->x[1] != sp->x[1] || si->x[2] != sp->x[2]) {
    error(
        "Moving particles that are not at the same location."
        " (%g, %g, %g) - (%g, %g, %g)",
        si->x[0], si->x[1], si->x[2], sp->x[0], sp->x[1], sp->x[2]);
  }
#endif

  /* Put the star randomly within the accretion radius of the sink */
  const double phi =
      2 * M_PI *
      random_unit_interval(sp->id, e->ti_current, (enum random_number_type)3);
  const double r = si->r_cut * random_unit_interval(sp->id, e->ti_current,
                                                    (enum random_number_type)4);
  const double cos_theta =
      1.0 - 2.0 * random_unit_interval(sp->id, e->ti_current,
                                       (enum random_number_type)5);
  const double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

  double new_pos[3] = {r * sin_theta * cos(phi), r * sin_theta * sin(phi),
                       r * cos_theta};

  /* Assign this new position to the star and its gpart */
  sp->x[0] += new_pos[0];
  sp->x[1] += new_pos[1];
  sp->x[2] += new_pos[2];
  sp->gpart->x[0] = sp->x[0];
  sp->gpart->x[1] = sp->x[1];
  sp->gpart->x[2] = sp->x[2];
}

/**
 * @brief Give a velocity to the #spart.
 *
 * In GEAR: Currently, a gaussian centered on 0 is used. The standard deviation
 * is computed based on the local gravitational dynamics of the system.
 *
 * @param e The #engine.
 * @param si The #sink generating a star.
 * @param sp The new #spart.
 * @param sink_props The sink properties to use.
 */
INLINE static void sink_star_formation_give_new_velocity(
    const struct engine* e, struct sink* si, struct spart* sp,
    const struct sink_props* sink_props) {

#ifdef HAVE_LIBGSL
  /* Those intermediate variables are the values that will be given to the star
     and subtracted from the sink. */
  double v_given[3] = {0.0, 0.0, 0.0};
  const double G_newton = e->physical_constants->const_newton_G;
  const double sigma_2 =
      G_newton * si->mass_tot_before_star_spawning / si->r_cut;
  const double sigma = sink_props->star_spawning_sigma_factor * sqrt(sigma_2);

  for (int i = 0; i < 3; ++i) {

    /* Draw a random value in unform interval (0, 1] */
    const double random_number = random_unit_interval_part_ID_and_index(
        sp->id, i, e->ti_current, (enum random_number_type)1);

    /* Sample a gaussian with mu=0 and sigma=sigma */
    double v_i_random = gsl_cdf_gaussian_Pinv(random_number, sigma);
    v_given[i] = v_i_random;
  }

  /* Update the star velocity. Do not forget to update the gpart velocity */
  sp->v[0] = si->v[0] + v_given[0];
  sp->v[1] = si->v[1] + v_given[1];
  sp->v[2] = si->v[2] + v_given[2];
  sp->gpart->v_full[0] = sp->v[0];
  sp->gpart->v_full[1] = sp->v[1];
  sp->gpart->v_full[2] = sp->v[2];
  message(
      "New star velocity: v = (%lf %lf %lf). Sink velocity: v = (%lf %lf %lf). "
      "Sigma = %lf",
      sp->v[0], sp->v[1], sp->v[2], si->v[0], si->v[1], si->v[2], sigma);
#else
  error("Code not compiled with GSL. Can't compute Star new velocity.");
#endif
}

/**
 * @brief Copy the properties of the sink particle towards the new star. Also,
 * give the stars some properties such as position and velocity.
 *
 * This function also needs to update the sink particle.
 *
 * @param sink The #sink particle.
 * @param sp The star particle.
 * @param e The #engine
 * @param sink_props The sink properties to use.
 * @param cosmo The cosmological parameters and properties.
 * @param with_cosmology If we run with cosmology.
 * @param phys_const The physical constants in internal units.
 * @param us The internal unit system.
 */
INLINE static void sink_copy_properties_to_star(
    struct sink* sink, struct spart* sp, const struct engine* e,
    const struct sink_props* sink_props, const struct cosmology* cosmo,
    const int with_cosmology, const struct phys_const* phys_const,
    const struct unit_system* restrict us) {

  /* Give the stars a new position */
  sink_star_formation_give_new_position(e, sink, sp);

  /* Set the mass (do not forget the sink's gpart friend!) */
  sp->mass = sink->target_mass_Msun * phys_const->const_solar_mass;
  sp->gpart->mass = sp->mass;

  /* Give a new velocity to the stars */
  sink_star_formation_give_new_velocity(e, sink, sp, sink_props);

  /* Sph smoothing length */
  sp->h = sink->r_cut;

  /* Feedback related initialisation */
  /* ------------------------------- */

  /* Initialize the feedback */
  feedback_init_after_star_formation(sp, e->feedback_props, sink->target_type);

  /* Star formation related initalisation */
  /* ------------------------------------ */

  /* Note: The sink module need to be compiled with GEAR SF as we store data
     in the SF struct. However, we do not need to run with --star-formation */

  /* Store either the birth_scale_factor or birth_time depending  */
  if (with_cosmology) {
    sp->birth_scale_factor = cosmo->a;
  } else {
    sp->birth_time = e->time;
  }

  /* Copy the chemistry properties */
  chemistry_copy_sink_properties_to_star(sink, sp);

  /* Copy the progenitor id */
  sp->sf_data.progenitor_id = sink->id;
}

/**
 * @brief Update the #sink particle properties before spawning a star.
 *
 * In GEAR, we check if the sink had an IMF change from pop III to pop II
 * during the last gas/sink accretion loops. If so, we draw a new target mass
 * with the correct IMF so that stars have the right metallicities.
 *
 * @param sink The #sink particle.
 * @param e The #engine
 * @param sink_props The sink properties to use.
 * @param phys_const The physical constants in internal units.
 */
INLINE static void sink_update_sink_properties_before_star_formation(
    struct sink* sink, const struct engine* e,
    const struct sink_props* sink_props, const struct phys_const* phys_const) {

  /* Has the sink accumulated enough metallicity so that the target mass
     should be updated before spawning stars?
     Between the last update of the target_mass_Msun, the sink may have accreted
     gas with metallicities that that are higher than those of population III
     stars. However, the target mass was set with the pop
     III IMF. */

  const struct feedback_props* feedback_props = e->feedback_props;

  /* Pick the correct table. (if only one table, threshold is < 0) */
  const float metal =
      chemistry_get_sink_total_iron_mass_fraction_for_feedback(sink);
  const float threshold = feedback_props->metallicity_max_first_stars;

  /* If metal < threshold, then the sink generate first star particles. */
  const int is_first_star = metal < threshold;

  /* If the sink has not changed its IMF yet
     (has_IMF_changed_from_popIII_to_popII = 0)
     but is eligible to (sink metal > threshold), get a target_mass_Msun of the
     pop II stars. */
  if (!(sink->has_IMF_changed_from_popIII_to_popII) && !is_first_star) {
    sink_update_target_mass(sink, sink_props, e, 0);

    /* Flag the sink to have made the transition of IMF. This ensures that next
    time we do not update the target_mass_Msun because metal > threshold
    (otherwise we would update it without needing to) */
    sink->has_IMF_changed_from_popIII_to_popII = 1;
    message("IMF transition : Sink %lld will now spawn Pop II stars.",
            sink->id);
  }
}

/**
 * @brief Update the #sink particle properties right after spawning a star.
 *
 * In GEAR: Important properties that are updated are the sink mass and the
 * sink->target_mass_Msun to draw the next star mass.
 *
 * @param sink The #sink particle that spawed stars.
 * @param sp The #spart particle spawned.
 * @param e The #engine
 * @param sink_props the sink properties to use.
 * @param phys_const the physical constants in internal units.
 * @param star_counter The star loop counter.
 */
INLINE static void sink_update_sink_properties_during_star_formation(
    struct sink* sink, const struct spart* sp, const struct engine* e,
    const struct sink_props* sink_props, const struct phys_const* phys_const,
    int star_counter) {

  /* count the number of stars spawned by this particle */
  sink->n_stars++;

  /* Update the mass */
  sink->mass =
      sink->mass - sink->target_mass_Msun * phys_const->const_solar_mass;

  /* Bug fix: Do not forget to update the sink gpart's mass. */
  sink->gpart->mass = sink->mass;

  /* This message must be put carefully after giving the star its mass,
     updated the sink mass and before changing the target_type */
  message(
      "%010lld spawn a star (%010lld) with mass %8.2f Msol type=%d  "
      "star_counter=%03d. Sink remaining mass: %e Msol.",
      sink->id, sp->id, sp->mass / phys_const->const_solar_mass,
      sink->target_type, star_counter,
      sink->mass / phys_const->const_solar_mass);

  /* Sample the IMF to the get next target mass */
  sink_update_target_mass(sink, sink_props, e, star_counter);
}

/**
 * @brief Update the #sink particle properties after star formation.
 *
 * In GEAR, this is unused.
 *
 * @param sink The #sink particle.
 * @param e The #engine
 * @param sink_props The sink properties to use.
 * @param phys_const The physical constants in internal units.
 */
INLINE static void sink_update_sink_properties_after_star_formation(
    struct sink* sink, const struct engine* e,
    const struct sink_props* sink_props, const struct phys_const* phys_const) {}

/**
 * @brief Store the gravitational potential of a particle by copying it from
 * its #gpart friend.
 *
 * @param p_data The sink data of a gas particle.
 * @param gp The part's #gpart.
 */
__attribute__((always_inline)) INLINE static void sink_store_potential_in_part(
    struct sink_part_data* p_data, const struct gpart* gp) {
  p_data->potential = gp->potential;
}

/**
 * @brief Compute all quantities required for the formation of a sink such as
 * kinetic energy, potential energy, etc. This function works on the
 * neighbouring gas particles.
 *
 * @param e The #engine.
 * @param p The #part for which we compute the quantities.
 * @param xp The #xpart data of the particle #p.
 * @param pi A neighbouring #part of #p.
 * @param xpi The #xpart data of the particle #pi.
 * @param cosmo The cosmological parameters and properties.
 * @param sink_props The sink properties to use.
 */
INLINE static void sink_prepare_part_sink_formation_gas_criteria(
    struct engine* e, struct part* restrict p, struct xpart* restrict xp,
    struct part* restrict pi, struct xpart* restrict xpi,
    const struct cosmology* cosmo, const struct sink_props* sink_props) {

  /* If for some reason the particle has been flagged to not form sink,
     do not continue and save some computationnal ressources. */
  if (!p->sink_data.can_form_sink) {
    return;
  }

  const int with_self_grav = (e->policy & engine_policy_self_gravity);

  /* Physical accretion radius of part p */
  const float r_acc_p = sink_props->cut_off_radius * cosmo->a;

  /* Comoving distance of particl p */
  const float px[3] = {(float)(p->x[0]), (float)(p->x[1]), (float)(p->x[2])};

  /* No need to check if the particle has been flagged to form a sink or
     not. This is done in runner_prepare_part_sink_formation(). */

  /* Compute the pairwise physical distance */
  const float pix[3] = {(float)(pi->x[0]), (float)(pi->x[1]),
                        (float)(pi->x[2])};

  const float dx[3] = {px[0] - pix[0], px[1] - pix[1], px[2] - pix[2]};
  const float dx_physical[3] = {dx[0] * cosmo->a, dx[1] * cosmo->a,
                                dx[2] * cosmo->a};
  const float r2_physical = dx_physical[0] * dx_physical[0] +
                            dx_physical[1] * dx_physical[1] +
                            dx_physical[2] * dx_physical[2];

  /* Checks that this part is a neighbour */
  if ((r2_physical > r_acc_p * r_acc_p) || (r2_physical == 0.0)) {
    return;
  }

  /* Do not form sinks if some neighbours are not active */
  if (!part_is_active(pi, e)) {
    p->sink_data.can_form_sink = 0;
    return;
  }

  const float mi = hydro_get_mass(p);
  const float u_inter_i = hydro_get_drifted_physical_internal_energy(p, cosmo);

  /* Compute the relative comoving velocity between p and pi */
  const float dv[3] = {pi->v[0] - p->v[0], pi->v[1] - p->v[1],
                       pi->v[2] - p->v[2]};

  const float a = cosmo->a;
  const float H = cosmo->H;
  const float a2H = a * a * H;

  /* Calculate the velocity with the Hubble flow */
  const float v_plus_H_flow[3] = {a2H * dx[0] + dv[0], a2H * dx[1] + dv[1],
                                  a2H * dx[2] + dv[2]};

  /* Compute the physical relative velocity between the particles */
  const float dv_physical[3] = {v_plus_H_flow[0] * cosmo->a_inv,
                                v_plus_H_flow[1] * cosmo->a_inv,
                                v_plus_H_flow[2] * cosmo->a_inv};

  const float dv_physical_squared = dv_physical[0] * dv_physical[0] +
                                    dv_physical[1] * dv_physical[1] +
                                    dv_physical[2] * dv_physical[2];

  /* Compute specific physical angular momentum between pk and pi */
  const float specific_angular_momentum[3] = {
      dx_physical[1] * dv_physical[2] - dx_physical[2] * dv_physical[1],
      dx_physical[2] * dv_physical[0] - dx_physical[0] * dv_physical[2],
      dx_physical[0] * dv_physical[1] - dx_physical[1] * dv_physical[0]};

  /* Updates the energies */
  p->sink_data.E_kin_neighbours += 0.5f * mi * dv_physical_squared;
  p->sink_data.E_int_neighbours += mi * u_inter_i;
  p->sink_data.E_rad_neighbours += cooling_get_radiated_energy(xpi);

  /* Notice that we skip the potential of the current particle here
     instead of subtracting it later */
  if ((with_self_grav) && (pi != p))
    p->sink_data.E_pot_self_neighbours +=
        0.5 * mi * pi->sink_data.potential * cosmo->a_inv;

  /* No external potential for now */
  /* if (gpi != NULL && with_ext_grav)	 */
  /* p->sink_data.E_pot_ext_neighbours +=  mi *
   * external_gravity_get_potential_energy( */
  /* time, potential, phys_const, gpi); */

  /* Need to include mhd header */
  /* p->sink_data.E_mag_neighbours += mhd_get_magnetic_energy(p, xpi); */

  /* Compute rotation energies per component */
  p->sink_data.E_rot_neighbours[0] +=
      0.5 * mi * specific_angular_momentum[0] * specific_angular_momentum[0] /
      sqrtf(dx_physical[1] * dx_physical[1] + dx_physical[2] * dx_physical[2]);
  p->sink_data.E_rot_neighbours[1] +=
      0.5 * mi * specific_angular_momentum[1] * specific_angular_momentum[1] /
      sqrtf(dx_physical[0] * dx_physical[0] + dx_physical[2] * dx_physical[2]);
  p->sink_data.E_rot_neighbours[2] +=
      0.5 * mi * specific_angular_momentum[2] * specific_angular_momentum[2] /
      sqrtf(dx_physical[0] * dx_physical[0] + dx_physical[1] * dx_physical[1]);

  /* Shall we reset the values of the energies for the next timestep? No, it is
     done in cell_drift.c and space_init.c, for active particles. The
     potential is set in runner_others.c->runner_do_end_grav_force() */
}

/**
 * @brief Compute all quantities required for the formation of a sink. This
 * function works on the neighbouring sink particles.
 *
 * @param e The #engine.
 * @param p The #part for which we compute the quantities.
 * @param xp The #xpart data of the particle #p.
 * @param si A neighbouring #sink of #p.
 * @param cosmo The cosmological parameters and properties.
 * @param sink_props The sink properties to use.
 */
INLINE static void sink_prepare_part_sink_formation_sink_criteria(
    struct engine* e, struct part* restrict p, struct xpart* restrict xp,
    struct sink* restrict si, const struct cosmology* cosmo,
    const struct sink_props* sink_props) {
  /* Do not continue if the gas cannot form sink for any reason */
  if (!p->sink_data.can_form_sink) {
    return;
  }

  /* Physical accretion radius of part p */
  const float r_acc_p = sink_props->cut_off_radius * cosmo->a;

  /* Physical accretion radius of sink si */
  const float r_acc_si = si->r_cut * cosmo->a;

  /* Comoving distance of particl p */
  const float px[3] = {(float)(p->x[0]), (float)(p->x[1]), (float)(p->x[2])};

  /* Compute the pairwise physical distance */
  const float six[3] = {(float)(si->x[0]), (float)(si->x[1]),
                        (float)(si->x[2])};

  const float dx[3] = {(px[0] - six[0]) * cosmo->a, (px[1] - six[1]) * cosmo->a,
                       (px[2] - six[2]) * cosmo->a};
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  /* If forming a sink from this particle will create a sink overlapping an
     existing sink's accretion radius, do not form a sink. This criterion can
     be disabled. */
  if (r2 < (r_acc_si + r_acc_p) * (r_acc_si + r_acc_p)) {
    p->sink_data.is_overlapping_sink = 1;
  }
}

#endif /* SWIFT_GEAR_SINK_H */
