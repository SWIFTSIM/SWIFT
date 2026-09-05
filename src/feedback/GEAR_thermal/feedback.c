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
#include "../GEAR/stellar_evolution.h"
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

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

  /* Did the particle receive an event? */
  /* TODO: Remove the ionization part from here and move it to cooling */
  if (!xp->feedback_data.hit_by_SN && !xp->feedback_data.hit_by_winds &&
      !xp->feedback_data.hit_by_radiation &&
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

  /* Update internal energy */
  const float u =
      hydro_get_physical_internal_energy(p, xp, cosmo) * old_mass / new_mass;
  const float u_new = u + xp->feedback_data.delta_u;

  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo, pressure_floor, u_new);

  xp->feedback_data.delta_u = 0.;

  /* Update the velocities */
  for (int i = 0; i < 3; i++) {
    const float dv = xp->feedback_data.delta_p[i] / new_mass;

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
 * @brief Finishes the #part density calculation.
 *
 * Nothing to do here.
 *
 * @param p The particle to act upon
 * @param xp The extra particle to act upon
 */
__attribute__((always_inline)) INLINE void feedback_end_density(
    struct part *p, struct xpart *xp) {}

/**
 * @brief Reset the gas particle-carried fields related to feedback at the
 * start of a step.
 *
 * Nothing to do here: u_FUV/u_LW are consumed and zeroed on the cooling
 * side instead (#cooling_expire_LW_FUV_dose_subgrid, called right after
 * cooling_copy_to_grackle has read them for this step's solve). Drift
 * (hence this function) always runs before this step's own cooling call,
 * so a reset placed here would fire before the value it is supposed to
 * follow has ever been read, wiping out an injection landed since the
 * last cooling call before Grackle ever saw it. is_illuminated_LW_FUV is
 * likewise never cleared here: Phase 1 has no decay model to end an
 * illumination episode, so the first-touch flag is set once
 * (radiation_iact_nonsym_feedback_apply) and never reset.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
void feedback_reset_part(struct part *p, struct xpart *xp) {}

/**
 * @brief First-init of a #part's feedback-model state (S3.3/F3-style
 * negative sentinel, mirroring #feedback_part_data.neutral_H_frac's own
 * treatment).
 *
 * u_FUV/u_LW are set to -1.f (never a physical specific energy) rather
 * than left at struct part's zero-init default, 0.0f: 0.0f is itself a
 * legitimate "genuinely unilluminated" value, so it cannot distinguish that
 * from "not yet consumed by this particle's own first cooling call".
 * cooling_expire_LW_FUV_dose_subgrid clears the sentinel to a legitimate
 * value the first time cooling_new_energy() runs for this particle; every
 * reader in between (#radiation_get_part_isrf_habing and friends) clamps
 * a negative value to 0 rather than assume that has already happened.
 *
 * @param p The #part to initialise.
 */
void feedback_first_init_part(struct part *restrict p) {
  p->feedback_data.u_FUV = -1.f;
  p->feedback_data.u_LW = -1.f;
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

  /* mass_HII_region is NOT reset here: this function runs every active
     step for every star, but the HII search only actually reruns on a
     rebuild step (feedback_will_do_feedback's need_HII_region_rebuild).
     Resetting unconditionally here wiped the accumulated mass long before
     any snapshot dump could read it. It is reset instead in
     feedback_will_do_feedback(), exactly when a fresh rebuild is about to
     recompute it. */

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
