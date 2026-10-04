/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Include header */
#include "feedback.h"

/* Local includes */
#include "../GEAR/radiation.h"
#include "../GEAR/radiation_iact.h"
#include "../GEAR/radiation_propagation_iact.h"
#include "../GEAR/stellar_evolution.h"
#include "chemistry.h"
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "minmax.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

#include <float.h>
#include <math.h>
#include <strings.h>

/**
 * @brief Update the properties of the particle due to a supernovae.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 */
void feedback_update_part(struct part *p, struct xpart *xp,
                          const struct engine *e) {

  /* Did the particle receive an event? delta_mass is tested on its own:
     ejecta can arrive with no hit flag set, and the mass must still be
     applied. */
  /* TODO: Remove the ionization part from here and move it to cooling */
  if (!xp->feedback_data.hit_by_SN && !xp->feedback_data.hit_by_winds &&
      !xp->feedback_data.hit_by_radiation &&
      xp->feedback_data.delta_mass == 0.f &&
      !radiation_is_part_tagged_as_ionized(p, xp))
    return;

  const struct cosmology *cosmo = e->cosmology;
  const struct pressure_floor_props *pressure_floor = e->pressure_floor_props;

  /* Turn off the cooling only when SN are involved */
  if (xp->feedback_data.hit_by_SN) {
    cooling_set_part_time_cooling_off(p, xp, e->time);
  }

  /* Update mass */
  const float old_mass = hydro_get_mass(p);
  const float new_mass = old_mass + xp->feedback_data.delta_mass;

  if (xp->feedback_data.delta_mass < 0.) {
    error("Delta mass smaller than 0");
  }

  hydro_set_mass(p, new_mass);

  xp->feedback_data.delta_mass = 0;

  /* Update the density */
  p->rho *= new_mass / old_mass;

  /* Update internal energy (the final mass shares the energy of all events) */
  const float new_mass_inv = 1.0f / new_mass;
  const float u = hydro_get_physical_internal_energy(p, xp, cosmo) * old_mass *
                  new_mass_inv;
  const float u_new = u + xp->feedback_data.delta_E_th * new_mass_inv;

  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo, pressure_floor, u_new);

  xp->feedback_data.delta_E_th = 0.0f;

  /* Update the velocities */
  for (int i = 0; i < 3; i++) {
    const float dv = xp->feedback_data.delta_p[i] * new_mass_inv;

    xp->v_full[i] += dv;
    p->v[i] += dv;

    xp->feedback_data.delta_p[i] = 0;
  }

  /*----------------------------------------*/
  /* Update the radiation fields */
  feedback_update_part_radiation(p, xp, e, old_mass);

  /* Update the wind fields */
  xp->feedback_data.hit_by_SN = 0;
  xp->feedback_data.hit_by_winds = 0;
}

/**
 * @brief Finishes the #part density calculation: caches the LW/PE
 * propagation speed, see #radiation_end_density_propagation.
 *
 * @param p The particle to act upon
 * @param xp The extra particle to act upon
 * @param e The #engine.
 */
void feedback_end_density(struct part *p, struct xpart *xp,
                          const struct engine *e) {
  radiation_end_density_propagation(p, e);
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0
 * neighbours, see #radiation_part_has_no_neighbours.
 *
 * @param p The particle to act upon.
 * @param xp The extra particle to act upon.
 * @param e The #engine.
 */
void feedback_part_has_no_neighbours(struct part *p, struct xpart *xp,
                                     const struct engine *e) {
  radiation_part_has_no_neighbours(p, e);
}

/**
 * @brief Finishes the #part gradient calculation, see
 * #radiation_end_gradient_propagation.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void feedback_end_gradient(struct part *p, const struct engine *e) {
  radiation_end_gradient_propagation(p, e);
}

/**
 * @brief Finishes the #part force calculation, see
 * #radiation_end_force_propagation.
 *
 * @param p The particle to act upon.
 * @param e The #engine.
 */
void feedback_end_force(struct part *p, const struct engine *e) {
  radiation_end_force_propagation(p, e);
}

/**
 * @brief Radiation timestep bound of a particle, see
 * #radiation_isrf_part_timestep.
 *
 * Stops the run if the bound is below TimeIntegration:dt_min.
 *
 * @param p The particle to consider.
 * @param e The #engine.
 * @return The radiation timestep bound (before the cosmology factor), or
 *     FLT_MAX if none applies.
 */
float feedback_compute_part_timestep(const struct part *restrict p,
                                     const struct engine *e) {
  const float dt_isrf = radiation_isrf_part_timestep(p, e);
  /* Compared to dt_min after the cosmology factor, like the other
   * candidates in get_part_timestep(). */
  const float dt_isrf_scaled = dt_isrf * e->cosmology->time_step_factor;
  if (dt_isrf_scaled < e->dt_min)
    error(
        "part (id=%lld) wants an ISRF radiation time-step (%e, %e after "
        "the cosmology factor) below TimeIntegration:dt_min (%e): "
        "GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c=%g forces dt_rad = "
        "C_hyp*h/(f*c) below dt_min for this particle's h. Lower the "
        "fraction (dt_rad grows as 1/f), or lower dt_min.",
        p->id, dt_isrf, dt_isrf_scaled, e->dt_min,
        e->feedback_props->ISRF_c_hyp_fixed_fraction_of_c);
  return dt_isrf;
}

/**
 * @brief Reset the feedback fields of a gas particle once per step, before
 * the density loop's h-iterations.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param e The #engine.
 */
void feedback_reset_part(struct part *p, struct xpart *xp,
                         const struct engine *e) {
  radiation_snapshot_part_propagation(p, e);
  radiation_reset_part_ISRF_illumination_tag(p, e);
  /* Must stay after the tag reset: an expired tag can zero `u`. */
  radiation_cache_m1_closure_part(p);
}

/**
 * @brief Re-initialise the feedback fields of a gas particle at the start of
 * each density h-iteration.
 *
 * @param p The particle.
 * @param e The #engine.
 */
void feedback_init_part(struct part *p, const struct engine *e) {
  radiation_init_part_propagation(p);
}

/**
 * @brief First-init of a #part's feedback state, see
 * #radiation_first_init_part.
 *
 * @param p The #part to initialise.
 */
void feedback_first_init_part(struct part *restrict p) {
  radiation_first_init_part(p);
}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
int feedback_is_active(const struct spart *sp, const struct engine *e) {

  /* the particle is inactive if its birth_scale_factor or birth_time is
   * negative */

  /* If the spart is dead, don't do anything */
  if (sp->birth_scale_factor < 0.0 || sp->birth_time < 0.0) return 0;

  return sp->feedback_data.will_do_feedback;
}

/**
 * @brief Is this star particle done evolving, i.e. finished with its
 * feedback-relevant lifetime?
 *
 * @param sp The #spart to query.
 * @return sp->feedback_data.is_dead.
 */
int feedback_is_star_dead(const struct spart *sp) {

  return sp->feedback_data.is_dead;
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * Note: In GEAR, this function must not reset the data as the are computed
 * at the end of the tasks by the stellar_evolution functions.
 *
 * @param sp The particle to act upon
 */
void feedback_init_spart(struct spart *sp) {

  sp->feedback_data.enrichment_weight = 0.f;
  sp->feedback_data.num_ngbs = 0;

  /* mass_HII_region is not reset here: the HII search only reruns on a
     rebuild step. It is reset in feedback_will_do_feedback(). */

  sp->feedback_data.grad_rho_star[0] = 0.0;
  sp->feedback_data.grad_rho_star[1] = 0.0;
  sp->feedback_data.grad_rho_star[2] = 0.0;

  sp->feedback_data.Z_star = 0.0;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 *
 * This is called in the stars ghost.
 */
void feedback_reset_feedback(struct spart *sp,
                             const struct feedback_props *feedback_props) {}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
void feedback_prepare_spart(struct spart *sp,
                            const struct feedback_props *feedback_props) {}

/**
 * @brief Prepare a #spart for the feedback task.
 *
 * This is called in the stars ghost task.
 *
 * In here, we only need to add the missing coefficients.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 * @param time The physical time in internal units.
 * @param ti_begin The integer time at the beginning of the step.
 * @param with_cosmology Are we running with cosmology on?
 */
void feedback_prepare_feedback(struct spart *restrict sp,
                               const struct feedback_props *feedback_props,
                               const struct cosmology *cosmo,
                               const struct unit_system *us,
                               const struct phys_const *phys_const,
                               const double star_age_beg_step, const double dt,
                               const double time, const integertime_t ti_begin,
                               const int with_cosmology) {
  /* Add missing h factor */
  const float hi_inv = 1.f / sp->h;
  const float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
  sp->feedback_data.enrichment_weight *= hi_inv_dim;

  /* Do radiation feedback */
  feedback_prepare_radiation_feedback(sp, feedback_props, cosmo, us, phys_const,
                                      star_age_beg_step, dt, time, ti_begin,
                                      with_cosmology);
}
