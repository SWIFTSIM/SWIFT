/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_SPHM1RT_H
#define SWIFT_RT_SPHM1RT_H

#include "rt_properties.h"
#include "rt_struct.h"

#include <float.h>

/**
 * @file src/rt/SPHM1RT/rt.h
 * @brief Main header file for no radiative transfer scheme.
 * SPHM1RT method described in Chan+21: 2102.08404
 */

/**
 * @brief Returns the comoving radiation energy per mass of a particle
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param p Pointer to the particle data.
 * @param urad The comoving radiation energy per mass.
 *
 */
__attribute__((always_inline)) INLINE static void
radiation_get_comoving_urad_multifrequency(const struct part* restrict p,
                                           float urad[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    urad[g] = p->rt_data.conserved[g].urad;
  }
}

/**
 * @brief Returns the physical radiation energy per mass of a particle
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param p Pointer to the particle data.
 * @param cosmo Cosmology data structure.
 * @param urad The physical radiation energy.
 *
 */
__attribute__((always_inline)) INLINE static void
radiation_get_physical_urad_multifrequency(const struct part* restrict p,
                                           const struct cosmology* cosmo,
                                           float urad[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    urad[g] = p->rt_data.conserved[g].urad;
  }
}

/**
 * @brief Sets the comoving radiation energy per mass of a particle
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param p The particle of interest.
 * @param urad The comoving radiation energy per mass
 *
 */
__attribute__((always_inline)) INLINE static void
radiation_set_comoving_urad_multifrequency(struct part* p,
                                           const float urad[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].urad = urad[g];
  }
}

/**
 * @brief Sets the physical radiation energy per mass of a particle
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param urad The physical radiation energy per mass
 */
__attribute__((always_inline)) INLINE static void
radiation_set_physical_urad_multifrequency(struct part* p,
                                           const struct cosmology* cosmo,
                                           const float urad[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].urad = urad[g];
  }
}

/**
 * @brief Returns the comoving radiation flux per gas density of a particle
 * (note that the comoving and physical flux per gas density are the same in our
 * convention)
 *
 * @param p Pointer to the particle data.
 * @param fradtemp The comoving radiation flux per gas density
 */
__attribute__((always_inline)) INLINE static void
radiation_get_comoving_frad_multifrequency(const struct part* restrict p,
                                           float fradtemp[RT_NGROUPS][3]) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    fradtemp[g][0] = p->rt_data.conserved[g].frad[0];
    fradtemp[g][1] = p->rt_data.conserved[g].frad[1];
    fradtemp[g][2] = p->rt_data.conserved[g].frad[2];
  }
}

/**
 * @brief Returns the physical radiation flux per gas density of a particle
 * (note that the comoving and physical flux per gas density are the same in our
 * convention)
 *
 * @param p Pointer to the particle data.
 * @param cosmo Cosmology data structure
 * @param fradtemp The comoving radiation flux per gas density
 */
__attribute__((always_inline)) INLINE static void
radiation_get_physical_frad_multifrequency(const struct part* restrict p,
                                           const struct cosmology* cosmo,
                                           float fradtemp[RT_NGROUPS][3]) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    fradtemp[g][0] = p->rt_data.conserved[g].frad[0];
    fradtemp[g][1] = p->rt_data.conserved[g].frad[1];
    fradtemp[g][2] = p->rt_data.conserved[g].frad[2];
  }
}

/**
 * @brief Sets the comoving radiation flux per density of a particle
 * (note that the comoving and physical flux per density are the same in our
 * convention)
 *
 * @param p The particle of interest.
 * @param frad The comoving radiation flux
 */
__attribute__((always_inline)) INLINE static void
radiation_set_comoving_frad_multifrequency(struct part* p,
                                           const float frad[RT_NGROUPS][3]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].frad[0] = frad[g][0];
    p->rt_data.conserved[g].frad[1] = frad[g][1];
    p->rt_data.conserved[g].frad[2] = frad[g][2];
  }
}

/**
 * @brief Sets the physical radiation flux of a particle
 * (note that the comoving and physical flux are the same in our convention)
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param frad The comoving radiation flux
 */
__attribute__((always_inline)) INLINE static void
radiation_set_physical_radiation_flux_multifrequency(
    struct part* p, const struct cosmology* cosmo, float frad[RT_NGROUPS][3]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].frad[0] = frad[g][0];
    p->rt_data.conserved[g].frad[1] = frad[g][1];
    p->rt_data.conserved[g].frad[2] = frad[g][2];
  }
}

/**
 * @brief Initialisation of the RT density loop related particle data.
 * Note: during initalisation (space_init), rt_reset_part and rt_init_part
 * are both called individually.
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_init_part(
    struct part* restrict p) {}

/**
 * @brief Reset of the RT hydro particle data not related to the density.
 * Note: during initalisation (space_init), rt_reset_part and rt_init_part
 * are both called individually.
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_reset_part(
    struct part* restrict p) {

  struct rt_part_data* rpd = &p->rt_data;

  for (int g = 0; g < RT_NGROUPS; g++) {
    rpd->dconserved_dt[g].urad = 0.0f;
    rpd->dconserved_dt[g].frad[0] = 0.0f;
    rpd->dconserved_dt[g].frad[1] = 0.0f;
    rpd->dconserved_dt[g].frad[2] = 0.0f;
  }

  for (int g = 0; g < RT_NGROUPS; g++) {
    rpd->viscosity[g].divf = 0.0f;
    rpd->diffusion[g].graduradc[0] = 0.0f;
    rpd->diffusion[g].graduradc[1] = 0.0f;
    rpd->diffusion[g].graduradc[2] = 0.0f;
  }

  float fradmag, flux_max, correct;

  /* To avoid radiation reaching other dimension and violating conservation */
  for (int g = 0; g < RT_NGROUPS; g++) {
    if (hydro_dimension < 1.001f) {
      rpd->conserved[g].frad[1] = 0.0f;
    }
    if (hydro_dimension < 2.001f) {
      rpd->conserved[g].frad[2] = 0.0f;
    }
  }

  for (int g = 0; g < RT_NGROUPS; g++) {
    /* TK: avoid the radiation flux to violate causality. Impose a limit: F<Ec
     */

    if ((rpd->conserved[g].frad[0] == 0.f) &&
        (rpd->conserved[g].frad[1] == 0.f) &&
        (rpd->conserved[g].frad[2] == 0.f)) {
      fradmag = 0.f;
    } else {
      fradmag = sqrtf(rpd->conserved[g].frad[0] * rpd->conserved[g].frad[0] +
                      rpd->conserved[g].frad[1] * rpd->conserved[g].frad[1] +
                      rpd->conserved[g].frad[2] * rpd->conserved[g].frad[2]);
    }

    flux_max = rpd->conserved[g].urad * rpd->params.cred;
    if (fradmag > 0.f) {
      if (fradmag > flux_max) {
        correct = flux_max / fradmag;
        rpd->conserved[g].frad[0] = rpd->conserved[g].frad[0] * correct;
        rpd->conserved[g].frad[1] = rpd->conserved[g].frad[1] * correct;
        rpd->conserved[g].frad[2] = rpd->conserved[g].frad[2] * correct;
      }
    } else {
      rpd->conserved[g].frad[0] = 0.0f;
      rpd->conserved[g].frad[1] = 0.0f;
      rpd->conserved[g].frad[2] = 0.0f;
    }
  }
}

/**
 * @brief First initialisation of the RT hydro particle data.
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_first_init_part(
    struct part* restrict p, const struct rt_props* restrict rt_props) {

  struct rt_part_data* rpd = &p->rt_data;

  for (int g = 0; g < RT_NGROUPS; g++) {
    rpd->viscosity[g].alpha = 1.0f;
    rpd->diffusion[g].alpha = 1.0f;
    rpd->params.chi[g] = rt_props->initialchi[g];
    rpd->viscosity[g].divf_previous_step = 0.0f;
  }

  /* We can get parameters for diffusion (force loop) */

  rpd->params.cred = rt_props->cred;

  rpd->force.f = 1.0f;

  rpd->dt = 1.0f;

  rt_init_part(p);
  rt_reset_part(p);
}

/**
 * @brief Initialises particle quantities that can't be set
 * otherwise before the zeroth step is finished. E.g. because
 * they require the particle density and time step to be known.
 *
 * @param p particle to work on
 * @param rt_props RT properties struct
 */
__attribute__((always_inline)) INLINE static void
rt_init_part_after_zeroth_step(struct part* restrict p,
                               const struct rt_props* rt_props) {}

/**
 * @brief Initialisation of the RT density loop related star particle data.
 * Note: during initalisation (space_init), rt_reset_spart and rt_init_spart
 * are both called individually.
 * @param sp star particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_init_spart(
    struct spart* restrict sp) {}

/**
 * @brief Reset of the RT star particle data not related to the density.
 * Note: during initalisation (space_init), rt_reset_spart and rt_init_spart
 * are both called individually.
 * @param sp star particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_reset_spart(
    struct spart* restrict sp) {}

/**
 * @brief First initialisation of the RT star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_spart(
    struct spart* restrict sp) {}

/**
 * @brief Initialises particle quantities that can't be set
 * otherwise before the zeroth step is finished. E.g. because
 * they require the star density and time step to be known.
 * @param sp star particle to work on
 * @param time current system time
 * @param star_age age of the star *at the end of the step*
 * @param dt star time step
 * @param rt_props RT properties struct
 * @param phys_const physical constants struct
 * @param internal_units struct holding internal units
 */
__attribute__((always_inline)) INLINE static void
rt_init_star_after_zeroth_step(struct spart* restrict sp, double time,
                               double star_age, double dt,
                               const struct rt_props* rt_props,
                               const struct phys_const* phys_const,
                               const struct unit_system* internal_units) {}

/**
 * @brief Split the RT data of a particle into n pieces
 *
 * @param p The #part.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void rt_split_part(struct part* p,
                                                                double n) {
  error("RT can't run with split particles for now.");
}

/**
 * @brief Exception handle a hydro part not having any neighbours in ghost task
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void rt_part_has_no_neighbours(
    struct part* p) {
  message("WARNING: found particle without neighbours");
};

/**
 * @brief Exception handle a star part not having any neighbours in ghost task
 *
 * @param sp The #spart.
 */
__attribute__((always_inline)) INLINE static void rt_spart_has_no_neighbours(
    struct spart* sp){};

/**
 * @brief Do checks/conversions on particles on startup.
 *
 * @param p The particle to work on
 * @param rtp The RT properties struct
 * @param phys_const physical constants struct
 * @param us unit_system struct
 * @param cosmo cosmology struct
 */
__attribute__((always_inline)) INLINE static void rt_convert_quantities(
    struct part* restrict p, const struct rt_props* rt_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo) {

  struct rt_part_data* rpd = &p->rt_data;
  /* Note that in the input, we read radiation energy and flux
   * then we convert these quantities to radiation energy per mass and flux per
   * mass
   */
  for (int g = 0; g < RT_NGROUPS; g++) {
    rpd->conserved[g].urad = rpd->conserved[g].urad / p->mass;
    rpd->conserved[g].frad[0] = rpd->conserved[g].frad[0] / p->mass;
    rpd->conserved[g].frad[1] = rpd->conserved[g].frad[1] / p->mass;
    rpd->conserved[g].frad[2] = rpd->conserved[g].frad[2] / p->mass;
  }
}

/**
 * @brief Computes the next radiative transfer time step size
 * of a given particle (during timestep tasks)
 *
 * @param p particle to work on
 * @param rt_props the RT properties struct
 * @param cosmo the cosmology
 */
__attribute__((always_inline)) INLINE static float rt_compute_timestep(
    const struct part* restrict p, const struct rt_props* restrict rt_props,
    const struct cosmology* restrict cosmo) {
  float dt = p->h * cosmo->a / (p->rt_data.params.cred + FLT_MIN) *
             rt_props->CFL_condition;

  return dt;
}

/**
 * @brief Computes the next radiative transfer time step size
 * of a given star particle (during timestep tasks).
 *
 * @param sp spart to work on
 * @param rt_props the RT properties struct
 * @param cosmo the cosmology
 */
__attribute__((always_inline)) INLINE static float rt_compute_spart_timestep(
    const struct spart* restrict sp, const struct rt_props* restrict rt_props,
    const struct cosmology* restrict cosmo) {

  return FLT_MAX;
}

/**
 * @brief Compute the time-step length for an RT step of a particle from given
 * integer times ti_beg and ti_end. This time-step length is then used to
 * compute the actual time integration of the transport/force step and the
 * thermochemistry. This is not used to determine the time-step length during
 * the time-step tasks.
 *
 * @param ti_beg Start of the time-step (on the integer time-line).
 * @param ti_end End of the time-step (on the integer time-line).
 * @param time_base Minimal time-step size on the time-line.
 * @param with_cosmology Are we running with cosmology integration?
 * @param cosmo The #cosmology object.
 *
 * @return The time-step size for the rt integration. (internal units).
 */
__attribute__((always_inline)) INLINE static double rt_part_dt(
    const integertime_t ti_beg, const integertime_t ti_end,
    const double time_base, const int with_cosmology,
    const struct cosmology* cosmo) {
  if (with_cosmology) {
    error("SPHM1RT with cosmology not implemented yet! :(");
    return 0.f;
  } else {
    return (ti_end - ti_beg) * time_base;
  }
}

/**
 * @brief This function finalises the injection step.
 *
 * @param p particle to work on
 * @param props struct #rt_props that contains global RT properties
 */
__attribute__((always_inline)) INLINE static void rt_finalise_injection(
    struct part* restrict p, struct rt_props* props) {}

/**
 * @brief Compute the photon emission rates for this stellar particle
 *        This function is called every time the spart is being reset
 *        (during start-up and during stars ghost if spart is active)
 *        and assumes that the photon emission rate is an intrinsic
 *        stellar property, i.e. doesn't depend on the environment.
 *
 * @param sp star particle to work on
 * @param time current system time
 * @param star_age age of the star *at the end of the step*
 * @param dt star time step
 * @param rt_props RT properties struct
 * @param phys_const physical constants struct
 * @param internal_units struct holding internal units
 */
__attribute__((always_inline)) INLINE static void
rt_compute_stellar_emission_rate(struct spart* restrict sp, double time,
                                 double star_age, double dt,
                                 const struct rt_props* rt_props,
                                 const struct phys_const* phys_const,
                                 const struct unit_system* internal_units) {}

/**
 * @brief finishes up the gradient computation
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_end_gradient(
    struct part* restrict p) {
  struct rt_part_data* rpd = &p->rt_data;
  /* artificial diffusion for shock capturing */
  const float vsig_diss = rpd->params.cred;
  /* similar to Cullen & Dehnen 2010 switch */
  float divf, divf_previous_step, urad, viscosity_alpha, diffusion_alpha;
  float divf_dt, shockest, alphaflim, alpha_f_diss, alpha_f_diss_loc;
  float alpha_diss_loc, alpha_diss;
  if (rpd->dt == 0) return;

  for (int g = 0; g < RT_NGROUPS; g++) {
    divf = rpd->viscosity[g].divf;
    divf_previous_step = rpd->viscosity[g].divf_previous_step;
    urad = rpd->conserved[g].urad;
    viscosity_alpha = rpd->viscosity[g].alpha;
    diffusion_alpha = rpd->diffusion[g].alpha;
    divf_dt = (divf - divf_previous_step) / (rpd->dt);

    if (urad == 0.f) {
      shockest = FLT_MAX;
    } else {
      shockest = -p->h * p->h / (vsig_diss) / (vsig_diss)*divf_dt * 200.f;
      shockest /= urad;
    }
    alphaflim = max(shockest, 0.0f); /* should be positive or 0 */
    alpha_f_diss = viscosity_alpha;
    alpha_f_diss_loc = 0.0f;

    /* f diffusion only operates in compression */
    if (divf < 0.0f) {
      /* limit the diffusivity to Courant time step */
      alpha_f_diss_loc = min(alphaflim, 1.0f);
    }

    if (alpha_f_diss_loc > alpha_f_diss) {
      /* Reset the value of alpha to the appropriate value */
      alpha_f_diss = alpha_f_diss_loc;
    } else {
      /* Integrate the alpha forward in time to decay back to alpha = alpha_loc
       */
      alpha_f_diss = alpha_f_diss_loc +
                     (alpha_f_diss - alpha_f_diss_loc) *
                         expf(-rpd->dt * vsig_diss *
                              (1.f / p->h + rpd->params.chi[g] * p->rho));
    }

    /* alpha inspired by Price 2010: it should vanish where radiation energy
     * difference is small */
    alpha_diss_loc = 1.0f;
    alpha_diss = diffusion_alpha;
    if (alpha_diss_loc > alpha_diss) {
      /* Reset the value of alpha to the appropriate value */
      alpha_diss = alpha_diss_loc;
    } else {
      /* Integrate the alpha forward in time to decay back to alpha = alpha_loc
       */
      alpha_diss = alpha_diss_loc +
                   (alpha_diss - alpha_diss_loc) *
                       expf(-rpd->dt * vsig_diss *
                            (0.01f / p->h + rpd->params.chi[g] * p->rho));
    }

    /* Cap the dissipation to avoid instabilities */
    alpha_diss = min(alpha_diss, 1.0f);
    alpha_diss = max(alpha_diss, 0.0f);

    alpha_f_diss = min(alpha_f_diss, 1.0f);
    alpha_f_diss = max(alpha_f_diss, 0.0f);

    rpd->diffusion[g].alpha = alpha_diss;
    rpd->viscosity[g].alpha = alpha_f_diss;
  }
}

/**
 * @brief finishes up the transport step
 *
 * @param p particle to work on
 * @param dt the current time step of the particle
 */
__attribute__((always_inline)) INLINE static void rt_finalise_transport(
    struct part* restrict p, const double dt) {
  struct rt_part_data* rpd = &p->rt_data;

  for (int g = 0; g < RT_NGROUPS; g++) {
    rpd->conserved[g].urad += rpd->dconserved_dt[g].urad * dt;
    rpd->conserved[g].frad[0] += rpd->dconserved_dt[g].frad[0] * dt;
    rpd->conserved[g].frad[1] += rpd->dconserved_dt[g].frad[1] * dt;
    rpd->conserved[g].frad[2] += rpd->dconserved_dt[g].frad[2] * dt;
  }

  /* add frad source term implicitly */
  float dfrac;

  for (int g = 0; g < RT_NGROUPS; g++) {
    dfrac = -rpd->params.chi[g] * p->rho * rpd->params.cred;
    rpd->conserved[g].frad[0] *= expf(dfrac * dt);
    rpd->conserved[g].frad[1] *= expf(dfrac * dt);
    rpd->conserved[g].frad[2] *= expf(dfrac * dt);

    /* update urad */
    /* limiter to avoid negative urad */
    /* negative urad will make the dissipation (diffusion) unstable) */
    if (rpd->conserved[g].urad < 0.0f) {
      rpd->conserved[g].urad = 0.0f;
      rpd->conserved[g].frad[0] = 0.0f;
      rpd->conserved[g].frad[1] = 0.0f;
      rpd->conserved[g].frad[2] = 0.0f;
    }

    /* save next time step */
    rpd->viscosity[g].divf_previous_step = rpd->viscosity[g].divf;
  }

  rpd->dt = dt;

  /* To avoid radiation reaching other dimension and violating conservation */
  for (int g = 0; g < RT_NGROUPS; g++) {
    if (hydro_dimension < 1.001f) {
      rpd->conserved[g].frad[1] = 0.0f;
    }
    if (hydro_dimension < 2.001f) {
      rpd->conserved[g].frad[2] = 0.0f;
    }
  }
}

/**
 * @brief Do the thermochemistry on a particle.
 *
 * @param p Particle to work on.
 * @param xp Pointer to the particle' extended data.
 * @param rt_props RT properties struct
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE static void rt_tchem(
    struct part* restrict p, struct xpart* restrict xp,
    struct rt_props* rt_props, const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us, const double dt) {}

/**
 * @brief Extra operations done during the kick.
 *
 * @param p Particle to act upon.
 * @param dt_therm Thermal energy time-step @f$\frac{dt}{a^2}@f$.
 * @param dt_grav Gravity time-step @f$\frac{dt}{a}@f$.
 * @param dt_hydro Hydro acceleration time-step
 * @f$\frac{dt}{a^{3(\gamma{}-1)}}@f$.
 * @param dt_kick_corr Gravity correction time-step @f$adt@f$.
 * @param cosmo Cosmology.
 * @param hydro_props Additional hydro properties.
 */
__attribute__((always_inline)) INLINE static void rt_kick_extra(
    struct part* p, float dt_therm, float dt_grav, float dt_hydro,
    float dt_kick_corr, const struct cosmology* cosmo,
    const struct hydro_props* hydro_props) {}

/**
 * @brief Prepare a particle for the !HYDRO! force calculation.
 * E.g. for the meshless schemes, we need to take into account the
 * mass fluxes of the ionizing species between particles.
 * NOTE: don't call this during rt_init_part or rt_reset_part,
 * follow the hydro_prepare_force logic.
 *
 * @param p particle to work on
 **/
__attribute__((always_inline)) INLINE static void rt_prepare_force(
    struct part* p) {

  struct rt_part_data* rpd = &p->rt_data;

  /* Some smoothing length multiples. */
  const float rho = hydro_get_comoving_density(p);
  const float rho_inv = 1.0f / rho; /* 1 / rho */

  /* Compute the "grad h" term */
  float rho_dh = p->density.rho_dh;

  const float omega_inv =
      1.f / (1.f + hydro_dimension_inv * p->h * rho_dh * rho_inv);

  /* Update variables. */
  rpd->force.f = omega_inv;
}

/**
 * @brief Clean the allocated memory inside the RT properties struct.
 *
 * @param props the #rt_props.
 * @param restart did we restart?
 */
__attribute__((always_inline)) INLINE static void rt_clean(
    struct rt_props* props, int restart) {}

#endif /* SWIFT_RT_SPHM1RT_H */
