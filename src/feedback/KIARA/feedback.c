/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2022 Doug Rennehan (douglas.rennehan@gmail.com)
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

/* This file's header */
#include "feedback.h"

/* Local includes. */
#include "hydro_properties.h"
#include "inline.h"
#include "random.h"
#include "timers.h"
#include "timestep_sync_part.h"


/**
 * @brief Determine if the star forming gas should kick a particle,
 *        and then kick it.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param fb_props The feedback properties.
 * @param ti_current The current timestep.
 */
double feedback_wind_probability(struct part* p, struct xpart* xp, 
                                 const struct engine* e, 
                                 const struct cosmology* cosmo,
                                 const struct feedback_props* fb_props, 
                                 const integertime_t ti_current, 
                                 const double dt_part,
                                 double *rand_for_sf_wind,
                                 double *wind_mass) {

  /* First thing we will do is generate a random number */
  *rand_for_sf_wind = random_unit_interval(p->id, ti_current,
                                           random_number_stellar_feedback_1);

  /* This is done in the RUNNER files. Therefore, we have access
   * to the gpart. */
  double galaxy_stellar_mass = p->gpart->fof_data.group_stellar_mass;
  if (galaxy_stellar_mass <= 0.) return 0.;

  const double stellar_mass_this_step = p->sf_data.SFR * dt_part;
  if (stellar_mass_this_step <= 0.) return 0.;

  /* If M* is non-zero, make sure it is at least resolved in the
   * following calculations.
   */
  if (galaxy_stellar_mass < fb_props->minimum_galaxy_stellar_mass) {
    galaxy_stellar_mass = fb_props->minimum_galaxy_stellar_mass;
  }

  /* When early wind suppression is enabled, we alter the minimum
   * stellar mass to be safe.
   */
  if (fb_props->early_wind_suppression_enabled) {
    const double early_minimum_stellar_mass =
        fb_props->early_stellar_mass_norm *
        exp(
          -1. *
          (
            (cosmo->a / fb_props->early_wind_suppression_scale_factor) *
            (cosmo->a / fb_props->early_wind_suppression_scale_factor)
          )
        );
    if (cosmo->a < fb_props->early_wind_suppression_scale_factor) {
      galaxy_stellar_mass = early_minimum_stellar_mass;
    }
  }

  *wind_mass = 
      fb_props->FIRE_eta_normalization * stellar_mass_this_step;
  if (galaxy_stellar_mass < fb_props->FIRE_eta_break) {
    (*wind_mass) *= pow(
      galaxy_stellar_mass / fb_props->FIRE_eta_break, 
      fb_props->FIRE_eta_lower_slope /*-0.317*/
    );
  } else {
    (*wind_mass) *= pow(
      galaxy_stellar_mass / fb_props->FIRE_eta_break, 
      fb_props->FIRE_eta_upper_slope /*-0.761*/
    );
  }

  /* Suppress stellar feedback in the early universe when galaxies are
   * too small. Star formation can destroy unresolved galaxies, so
   * we must suppress the stellar feedback.
   */
  if (fb_props->early_wind_suppression_enabled) {
    if (cosmo->a < fb_props->early_wind_suppression_scale_factor) {
      (*wind_mass) *= pow(cosmo->a / fb_props->early_wind_suppression_scale_factor, 
                       fb_props->early_wind_suppression_slope);
    }
  }

  const float probability_to_kick = 1. - exp(-(*wind_mass) / hydro_get_mass(p));
  return probability_to_kick;

}

void feedback_kick_and_decouple_part(struct part* p, struct xpart* xp, 
                                     const struct engine* e, 
                                     const struct cosmology* cosmo,
                                     const struct feedback_props* fb_props, 
                                     const integertime_t ti_current,
                                     const int with_cosmology,
                                     const double dt_part,
                                     const double wind_mass) {

  const double galaxy_stellar_mass = 
      p->gpart->fof_data.group_stellar_mass;
  const double galaxy_stellar_mass_Msun =
      galaxy_stellar_mass * fb_props->mass_to_solar_mass;
  /* This is done in the RUNNER files. Therefore, we have
   * access to the gpart */
  const double galaxy_gas_stellar_mass_Msun = 
      p->gpart->fof_data.group_mass * fb_props->mass_to_solar_mass;
  if (galaxy_gas_stellar_mass_Msun <= 0. || galaxy_stellar_mass <= 0.) return;

  /* Physical circular velocity km/s */
  const double v_circ_km_s = 
      pow(galaxy_gas_stellar_mass_Msun / 102.329, 0.26178) *
      pow(cosmo->H / cosmo->H0, 1. / 3.);
  const double rand_for_scatter = random_unit_interval(p->id, ti_current,
                                      random_number_stellar_feedback_2);

  /* The wind velocity in internal units */
  double wind_velocity =
      fb_props->FIRE_velocity_normalization *
      pow(v_circ_km_s / 200., fb_props->FIRE_velocity_slope) *
      (
        1. - fb_props->kick_velocity_scatter + 
        2. * fb_props->kick_velocity_scatter * rand_for_scatter
      ) *
      v_circ_km_s *
      fb_props->kms_to_internal;

  /* Now we have wind_velocity in internal units, determine how much should go to heating */
  const double u_wind = 0.5 * wind_velocity * wind_velocity;
  
  /* Metal mass fraction (Z) of the gas particle */
  const double Z = p->chemistry_data.metal_mass_fraction_total;

  /* Supernova energy in internal units */
  double u_SN = ((1.e51 * (0.0102778 / fb_props->solar_mass_in_g) * 
                    (p->sf_data.SFR * dt_part / wind_mass)) /
                    (fb_props->kms_to_cms * fb_props->kms_to_cms)) *
		                (fb_props->kms_to_internal * fb_props->kms_to_internal);
  if (Z > 1.e-9) {
    u_SN *= pow(10., -0.0029 * pow(log10(Z) + 9., 2.5) + 0.417694);
  } else {
    u_SN *= 2.61634;
  }

  /* Limit the kinetic energy in the winds to the available SN energy */
  if (u_wind > u_SN) wind_velocity *= sqrt(u_SN / u_wind);

  /* 0.2511886 = pow(10., -0.6) */
  float pandya_slope = 0.f;
  if (galaxy_stellar_mass_Msun > 3.16e10) {
    pandya_slope = -2.1f;
  } else {
    pandya_slope = -0.1f;
  }

  const double f_warm = 0.2511886 * pow(galaxy_stellar_mass_Msun / 3.16e10, pandya_slope);
  const double hot_wind_fraction = max(0., 0.9 - f_warm); /* additional 10% removed for cold phase */
  const double rand_for_hot = random_unit_interval(p->id, ti_current,
                                                   random_number_stellar_feedback_3);
  const double rand_for_spread = random_unit_interval(p->id, ti_current,
                                                      random_number_stellar_feedback);

  /* We want these for logging purposes */
  const double u_init = hydro_get_physical_internal_energy(p, xp, cosmo);
  double u_new = 0.;
  if (u_SN > u_wind && rand_for_hot < hot_wind_fraction) {
    u_new = u_init + (u_SN - u_wind) * (0.5 + rand_for_spread);
  } else {
    u_new = fb_props->cold_wind_internal_energy;
  }

  if (u_new / u_init > 1000) {
    warning("Wind heating too large! T0=%g Tnew=%g fw=%g hwf=%g TSN=%g Tw=%g vw=%g ms=%g mwind=%g", 
            u_init / fb_props->temp_to_u_factor, 
            u_new / fb_props->temp_to_u_factor, 
            f_warm, 
            hot_wind_fraction, 
            u_SN / fb_props->temp_to_u_factor, 
            u_wind / fb_props->temp_to_u_factor, 
            wind_velocity, 
            p->sf_data.SFR * dt_part, 
            wind_mass);

    u_new = u_init * 1000;
  }

  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo, u_new);

  const double dir[3] = {
    p->gpart->a_grav[1] * p->gpart->v_full[2] - 
        p->gpart->a_grav[2] * p->gpart->v_full[1],
    p->gpart->a_grav[2] * p->gpart->v_full[0] - 
        p->gpart->a_grav[0] * p->gpart->v_full[2],
    p->gpart->a_grav[0] * p->gpart->v_full[1] - 
        p->gpart->a_grav[1] * p->gpart->v_full[0]
  };
  const double norm = sqrt(
    dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]
  );
  /* No norm, no wind */
  if (norm <= 0.) return;
  const double prefactor = cosmo->a * wind_velocity / norm;

  p->v_full[0] += prefactor * dir[0];
  p->v_full[1] += prefactor * dir[1];
  p->v_full[2] += prefactor * dir[2];

  /* Update the signal velocity of the particle based on the velocity kick. */
  hydro_set_v_sig_based_on_velocity_kick(p, cosmo, wind_velocity);

  /* Impose maximal viscosity */
  hydro_diffusive_feedback_reset(p);

  /* Synchronize the particle on the timeline */
  timestep_sync_part(p);

  /* Decouple the particles from the hydrodynamics */
  p->feedback_data.decoupling_delay_time = 
      fb_props->wind_decouple_time_factor * 
      cosmology_get_time_since_big_bang(cosmo, cosmo->a);

  p->feedback_data.number_of_times_decoupled += 1;

#ifdef WITH_FOF_GALAXIES
  /* Wind particles are never grouppable. This is done in the
   * RUNNER files. Therefore, we have access to the gpart. */
  p->gpart->fof_data.is_grouppable = 0;
#endif

  /* Wind cannot be star forming */
  if (p->sf_data.SFR > 0.f) {

    /* Record the current time as an indicator of when this particle was last
       star-forming. */
    if (with_cosmology) {
      p->sf_data.SFR = -e->cosmology->a;
    } else {
      p->sf_data.SFR = -e->time;
    }

  }

  /**
   * z pid dt M* Mb vkick vkx vky vkz h x y z vx vy vz T rho v_sig decoupletime 
   * Ndecouple
   */
  const float length_convert = cosmo->a * fb_props->length_to_kpc;
  const float velocity_convert = cosmo->a_inv / fb_props->kms_to_internal;
  const float rho_convert = cosmo->a3_inv * fb_props->rho_to_n_cgs;
  const float u_convert = 
      cosmo->a_factor_internal_energy / fb_props->temp_to_u_factor;
  printf("WIND_LOG %.3f %lld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %g\n",
          cosmo->z,
          p->id, 
          dt_part * fb_props->time_to_Myr,
          galaxy_stellar_mass * fb_props->mass_to_solar_mass,
          galaxy_gas_stellar_mass_Msun,
          wind_velocity / fb_props->kms_to_internal,
          prefactor * dir[0] * velocity_convert,
          prefactor * dir[1] * velocity_convert,
          prefactor * dir[2] * velocity_convert,
          p->h * cosmo->a * fb_props->length_to_kpc,
          p->x[0] * length_convert, 
          p->x[1] * length_convert, 
          p->x[2] * length_convert,
          p->v_full[0] * velocity_convert, 
          p->v_full[1] * velocity_convert, 
          p->v_full[2] * velocity_convert,
          p->u * u_convert, 
          p->rho * rho_convert, 
          p->viscosity.v_sig * velocity_convert,
          p->feedback_data.decoupling_delay_time * fb_props->time_to_Myr, 
          p->feedback_data.number_of_times_decoupled,
          u_new / u_init);
}

/**
 * @brief Run the Chem5 module that interpolates the yield tables and returns
 *        the ejected mass, metals, and unprocessed materials.
 *
 * @param sp The #spart to consider.
 * @param age The stellar age in code units.
 * @param fb_props The feedback properties.
 * @param dt The current timestep.
 * @param ejecta_energy The total ejected energy in code units.
 * @param ejecta_mass The total ejected mass in code units.
 * @param ejecta_unprocessed The unprocessed mass in code units.
 * @param ejecta_metal_mass The metal masses for each element in chem5_element_count in code units.
 */
void feedback_get_ejecta_from_star_particle(const struct spart* sp,
                                            double age,
                                            const struct feedback_props* fb_props,
                                            double dt,
                                            float *ejecta_energy,
                                            float *ejecta_mass,
                                            float *ejecta_unprocessed,
                                            float ejecta_metal_mass[chem5_element_count]) {
  int j, k, j1, j2, l, l1, l2, ll1, ll2, lll1, lll2;
  double SW_R, SNII_R, SNII_U, SNII_E, SNII_Z[chem5_element_count];
  double SNII_ENE, SNIa_R, SNIa_E, SNIa_Z[chem5_element_count];
  double SNn, SWn;
  double SNIIa, SNIIb, z, lz;

  /* Convert to yr for code below */
  age *= fb_props->time_to_yr;
  dt *= fb_props->time_to_yr;

  *ejecta_energy = 0.f;
  *ejecta_mass = 0.f;
  *ejecta_unprocessed = 0.f;
  for (k = 0; k < chem5_element_count; k++) ejecta_metal_mass[k] = 0.f;

  /* @TODO What does "fb" mean? fb stage? */
  int fb = 0;
  int fb_first = 0;

  if (sp->mass_init == sp->mass) fb_first = 1;

  z = sp->chemistry_data.metal_mass_fraction_total;
  feedback_set_turnover_mass(z);

  /* [Fe/H] */
  float feh = -10.f;
  if (z < 1.e-10f) {
    lz = -10.f;
  } 
  else {
    lz = log10f(z);
    feh = sp->chemistry_data.metal_mass_fraction[chemistry_element_Fe] / sp->chemistry_data.metal_mass_fraction[chemistry_element_H];
    if (feh > 0.f) feh = log10f((feh / fb_props->Fe_mf) * fb_props->H_mf);
  }

  float tm1 = feedback_get_turnover_mass(fb_props, age);

  if (tm1 >= fb_props->M_u3) return;

  float ltm = log10f(tm1 * fb_props->solar_mass_to_mass);

  for (j = 1; j < NM; j++) {
    j1 = j - 1;
    j2 = j;
    if (SNLM[j] < ltm) break;
  }

  /* @TODO D. Rennehan: Pretty sure this is NEVER true? */
  if (z <= fb_props->zmax3) {
    SNII_U = feedback_linear_interpolation(SNLM[j1], SN2E[1][0][j1], SNLM[j2], SN2E[1][0][j2], ltm);
    SNII_E = feedback_linear_interpolation(SNLM[j1], SN2E[2][0][j1], SNLM[j2], SN2E[2][0][j2], ltm);
    SNII_ENE = feedback_linear_interpolation(SNLM[j1], SN2E[0][0][j1], SNLM[j2], SN2E[0][0][j2], ltm);

    for (k = 0; k < chem5_element_count; k++) {
      SNII_Z[k] = feedback_linear_interpolation(SNLM[j1], SN2E[k + 3][0][j1], SNLM[j2], SN2E[k+3][0][j2], ltm);
    }

    SNII_R = feedback_linear_interpolation(SNLM[j1], SN2R[0][j1], SNLM[j2], SN2R[0][j2], ltm);
    SW_R = feedback_linear_interpolation(SNLM[j1], SWR[0][j1], SNLM[j2], SWR[0][j2], ltm);
    SNIa_R = 0.0;
  } 
  else {
    for (l = 2; l < NZSN; l++) {
      l1 = l - 1;
      l2 = l;
      if (SNLZ[l] > lz) break;
    }
    for (l = 1; l < NZSN1R; l++) {
      ll1 = l - 1;
      ll2 = l;
      if (SNLZ1R[l] > feh) break;
    }
    for (l = 1; l < NZSN1Y; l++) {
      lll1 = l - 1;
      lll2 = l;
      if (SNLZ1Y[l] > lz) break;
    }

    SNIIa = feedback_linear_interpolation(SNLM[j1], SN2E[1][l1][j1], SNLM[j2], SN2E[1][l1][j2], ltm);
    SNIIb = feedback_linear_interpolation(SNLM[j1], SN2E[1][l2][j1], SNLM[j2], SN2E[1][l2][j2], ltm);
    SNII_U = feedback_linear_interpolation(SNLZ[l1], SNIIa, SNLZ[l2], SNIIb, lz);
    SNIIa = feedback_linear_interpolation(SNLM[j1], SN2E[2][l1][j1], SNLM[j2], SN2E[2][l1][j2], ltm);
    SNIIb = feedback_linear_interpolation(SNLM[j1], SN2E[2][l2][j1], SNLM[j2], SN2E[2][l2][j2], ltm);
    SNII_E = feedback_linear_interpolation(SNLZ[l1], SNIIa, SNLZ[l2], SNIIb,lz);

    if (l2 == NZSN - 1) {
      SNII_ENE = feedback_linear_interpolation(SNLM[j1], SN2E[0][l2][j1], SNLM[j2], SN2E[0][l2][j2], ltm);
    }
    else {
      SNIIa = feedback_linear_interpolation(SNLZ[l1], SN2E[0][l1][j1], SNLZ[l2], SN2E[0][l2][j1], lz);
      SNIIb = feedback_linear_interpolation(SNLZ[l1], SN2E[0][l1][j2], SNLZ[l2], SN2E[0][l2][j2], lz);
      SNII_ENE = feedback_linear_interpolation(SNLM[j1], SNIIa,SNLM[j2], SNIIb, ltm);
    }

    for (k = 0; k < chem5_element_count; k++) {
      SNIIa = feedback_linear_interpolation(SNLM[j1], SN2E[k+3][l1][j1], SNLM[j2], SN2E[k+3][l1][j2], ltm);
      SNIIb = feedback_linear_interpolation(SNLM[j1], SN2E[k+3][l2][j1], SNLM[j2], SN2E[k+3][l2][j2], ltm);
      SNII_Z[k] = feedback_linear_interpolation(SNLZ[l1], SNIIa, SNLZ[l2], SNIIb, lz);
    }

    SNIIa = feedback_linear_interpolation(SNLM[j1], SN2R[l1][j1], SNLM[j2], SN2R[l1][j2], ltm);
    SNIIb = feedback_linear_interpolation(SNLM[j1], SN2R[l2][j1], SNLM[j2], SN2R[l2][j2], ltm);
    SNII_R = feedback_linear_interpolation(SNLZ[l1], SNIIa, SNLZ[l2], SNIIb, lz);
    SNIIa = feedback_linear_interpolation(SNLM[j1], SWR[l1][j1], SNLM[j2], SWR[l1][j2], ltm);
    SNIIb = feedback_linear_interpolation(SNLM[j1], SWR[l2][j1], SNLM[j2], SWR[l2][j2], ltm);
    SW_R = feedback_linear_interpolation(SNLZ[l1], SNIIa, SNLZ[l2], SNIIb, lz);

    if (feh < SNLZ1R[0]) {
      SNIa_R = 0.0;
    }
    else {
      SNIa_E = feedback_linear_interpolation(SNLZ1Y[lll1], SN1E[2][lll1], SNLZ1Y[lll2], SN1E[2][lll2], lz);

      for (k = 0; k < chem5_element_count; k++) SNIa_Z[k] = feedback_linear_interpolation(SNLZ1Y[lll1], SN1E[k+3][lll1], SNLZ1Y[lll2], SN1E[k+3][lll2], lz);

      SNIIa = feedback_linear_interpolation(SNLM[j1], SN1R[ll1][j1], SNLM[j2], SN1R[ll1][j2], ltm);
      SNIIb = feedback_linear_interpolation(SNLM[j1], SN1R[ll2][j1], SNLM[j2], SN1R[ll2][j2], ltm);
      SNIa_R = feedback_linear_interpolation(SNLZ1R[ll1], SNIIa, SNLZ1R[ll2], SNIIb, feh);
    }
  }

  float tm2 = 1.e-10f;

  if (age > dt) {
    tm2 = feedback_get_turnover_mass(fb_props, age - dt);

    ltm = log10(tm2 * fb_props->solar_mass_to_mass);
    for (j = 1 ; j < NM; j++) {
      j1 = j - 1;
      j2 = j;
      if (SNLM[j] < ltm) break;
    }

    if (z <= fb_props->zmax3) {
      SNII_U -= feedback_linear_interpolation(SNLM[j1], SN2E[1][0][j1], SNLM[j2], SN2E[1][0][j2], ltm);
      SNII_E -= feedback_linear_interpolation(SNLM[j1], SN2E[2][0][j1], SNLM[j2], SN2E[2][0][j2], ltm);
      SNII_ENE -= feedback_linear_interpolation(SNLM[j1], SN2E[0][0][j1], SNLM[j2], SN2E[0][0][j2], ltm);

      for (k = 0; k < chem5_element_count; k++) SNII_Z[k] -= feedback_linear_interpolation(SNLM[j1], SN2E[k+3][0][j1], SNLM[j2], SN2E[k+3][0][j2], ltm);

      SNII_R -= feedback_linear_interpolation(SNLM[j1], SN2R[0][j1], SNLM[j2], SN2R[0][j2], ltm);
      SW_R -= feedback_linear_interpolation(SNLM[j1], SWR[0][j1], SNLM[j2], SWR[0][j2], ltm);
    } 
    else {
      SNIIa = feedback_linear_interpolation(SNLM[j1], SN2E[1][l1][j1], SNLM[j2], SN2E[1][l1][j2], ltm);
      SNIIb = feedback_linear_interpolation(SNLM[j1], SN2E[1][l2][j1], SNLM[j2], SN2E[1][l2][j2], ltm);
      SNII_U -= feedback_linear_interpolation(SNLZ[l1], SNIIa, SNLZ[l2], SNIIb, lz);
      SNIIa = feedback_linear_interpolation(SNLM[j1], SN2E[2][l1][j1], SNLM[j2], SN2E[2][l1][j2], ltm);
      SNIIb = feedback_linear_interpolation(SNLM[j1], SN2E[2][l2][j1], SNLM[j2], SN2E[2][l2][j2], ltm);
      SNII_E -= feedback_linear_interpolation(SNLZ[l1], SNIIa,SNLZ[l2], SNIIb, lz);

      if (l2 == NZSN - 1) {
        SNII_ENE -= feedback_linear_interpolation(SNLM[j1], SN2E[0][l2][j1], SNLM[j2], SN2E[0][l2][j2], ltm);
      }
      else {
        SNIIa = feedback_linear_interpolation(SNLZ[l1], SN2E[0][l1][j1], SNLZ[l2], SN2E[0][l2][j1], lz);
        SNIIb = feedback_linear_interpolation(SNLZ[l1], SN2E[0][l1][j2], SNLZ[l2], SN2E[0][l2][j2], lz);
        SNII_ENE -= feedback_linear_interpolation(SNLM[j1], SNIIa,SNLM[j2], SNIIb, ltm);
      }
    
      for (k = 0; k < chem5_element_count; k++) {
        SNIIa = feedback_linear_interpolation(SNLM[j1], SN2E[k+3][l1][j1], SNLM[j2], SN2E[k+3][l1][j2], ltm);
        SNIIb = feedback_linear_interpolation(SNLM[j1], SN2E[k+3][l2][j1], SNLM[j2], SN2E[k+3][l2][j2], ltm);
        SNII_Z[k] -= feedback_linear_interpolation(SNLZ[l1], SNIIa, SNLZ[l2], SNIIb, lz);
      }

      SNIIa = feedback_linear_interpolation(SNLM[j1], SN2R[l1][j1], SNLM[j2], SN2R[l1][j2], ltm);
      SNIIb = feedback_linear_interpolation(SNLM[j1], SN2R[l2][j1], SNLM[j2], SN2R[l2][j2], ltm);
      SNII_R -= feedback_linear_interpolation(SNLZ[l1], SNIIa, SNLZ[l2], SNIIb, lz);
      SNIIa = feedback_linear_interpolation(SNLM[j1], SWR[l1][j1], SNLM[j2], SWR[l1][j2], ltm);
      SNIIb = feedback_linear_interpolation(SNLM[j1], SWR[l2][j1], SNLM[j2], SWR[l2][j2], ltm);
      SW_R -= feedback_linear_interpolation(SNLZ[l1], SNIIa, SNLZ[l2], SNIIb, lz);

      if (feh < SNLZ1R[0]) {
        SNIa_R = 0.0;
      }
      else {
        SNIIa = feedback_linear_interpolation(SNLM[j1], SN1R[ll1][j1], SNLM[j2], SN1R[ll1][j2], ltm);
        SNIIb = feedback_linear_interpolation(SNLM[j1], SN1R[ll2][j1], SNLM[j2], SN1R[ll2][j2], ltm);
        SNIa_R  -= feedback_linear_interpolation(SNLZ1R[ll1], SNIIa, SNLZ1R[ll2], SNIIb, feh);
      }
    }
  }

  *ejecta_unprocessed = max(0.f, sp->mass_init * SNII_U);
  *ejecta_mass = max(0.f, sp->mass_init * SNII_E);

  if (tm1 > fb_props->M_u2) {
    fb = 3;
    if (z == 0.f) {
      *ejecta_energy = 0.f;
    }
    else {
      SWn = sp->mass_init * SW_R;
      if (fb_props->with_HN_energy_from_chem5) {
        *ejecta_energy = SWn * fb_props->E_sw * powf(z / fb_props->Z_mf, 0.8f);
      } 
      else {
        *ejecta_energy += 0.f;
      }
    }
  } else {
    if (tm2 > fb_props->M_l2 || fb_first == 1) {
      fb = 2;

      SWn = sp->mass_init * SW_R;
      SNn = sp->mass_init * SNII_R;
      if (fb_props->with_SNII_energy_from_chem5) {
        *ejecta_energy = SWn * fb_props->E_sw;
        *ejecta_energy += sp->mass_init * SNII_ENE;
      }
      else {
        *ejecta_energy += 0.f;
      }
    }

    for (k = 0; k < chem5_element_count; k++) ejecta_metal_mass[k] = sp->mass_init * SNII_Z[k];

    if (tm1 <= fb_props->M_l2) {
      if (feh < SNLZ1R[0]) {
        SNIa_R = 0.f;
      }
      else {
        fb = 1;
        SNn = sp->mass_init * SNIa_R;
        if (fb_props->with_SNIa_energy_from_chem5) {
          *ejecta_energy += SNn * fb_props->E_sn1;
        }
        else {
          *ejecta_energy += 0.f;
        }

        *ejecta_mass += SNn * SNIa_E;
        for (k = 0; k < chem5_element_count; k++) {
          ejecta_metal_mass[k] += SNn * SNIa_Z[k];
        }
      }
    }
  }
}

float feedback_life_time(const float m, const float z) {
  int j, j1, j2, l, l1, l2;
  float lm, lz, ta, tb, t;

  lm = log10f(m);
  j1 = 0;
  j2 = 1;
  for (j = 1; j < NMLF; j++) {
    j1 = j - 1;
    j2 = j;
    if (LFLM[j] > lm) break;
  }

  if (z == 0.) {
    lz = -990.f;
  }
  else {
    lz = log10f(z);
  }

  if (lz <= LFLZ[0]) {
    t = feedback_linear_interpolation(LFLM[j1], LFLT[0][j1], LFLM[j2], LFLT[0][j2], lm);
  }
  else {
    l1 = 0;
    l2 = 1;
    for (l = 1; l < NZLF; l++) {
      l1 = l - 1;
      l2 = l;
      if (LFLZ[l] > lz) break;
    }
    ta = feedback_linear_interpolation(LFLZ[l1], LFLT[l1][j1], LFLZ[l2], LFLT[l2][j1], lz);
    tb = feedback_linear_interpolation(LFLZ[l1], LFLT[l1][j2], LFLZ[l2], LFLT[l2][j2], lz);
    t = feedback_linear_interpolation(LFLM[j1], ta, LFLM[j2], tb, lm);
  }

  return powf(10.f, t);
}

float feedback_imf(const struct feedback_props* fb_props, const float m) {
  if (fb_props->imf == 0) { /* Kroupa */
    if (m >= 0.5) {
      return powf(m, -fb_props->ximf) * 0.5f;
    }
    else if (m >= 0.08f) {
      return powf(m, -0.3f);
    }
    else {
      return powf(m, 0.7f) / 0.08f;
    }
  }
  else {
    return pow(m, -fb_props->ximf);
  }
}

void feedback_set_turnover_mass(const float z) {
  float lz;
  int j, l, l1, l2;

  if (z == 0.) {
    lz = -4.1f;
  }
  else {
    lz = log10f(z);
  }

  if (lz < -4.1f) lz = -4.1f;

  l1 = 0;
  l2 = 1;
  for (l = 1; l < NZLF; l++) {
    l1 = l - 1;
    l2 = l;
    if (LFLZ[l] > lz) break;
  }

  for (j = 0; j < NMLF; j++) LFLT2[j] = feedback_linear_interpolation(LFLZ[l1], LFLT[l1][j], LFLZ[l2], LFLT[l2][j], lz);
  return;
}

float feedback_get_turnover_mass(const struct feedback_props* fb_props, const float t) {
  if (t == 0.0) return fb_props->M_u3;

  float result, m, lt;
  int j, j1, j2;

  lt = log10f(t);
  j1 = 0;
  j2 = 1;
  for (j = 1; j < NMLF; j++) {
    j1 = j - 1;
    j2 = j;
    if (LFLT2[j] < lt) break;
  }

  m = feedback_linear_interpolation(LFLT2[j1], LFLM[j1], LFLT2[j2], LFLM[j2], lt);
  result = powf(10.f, m);

  if (result < fb_props->M_l) return fb_props->M_l;
  if (result > fb_props->M_u3) return fb_props->M_u3;

  return result;
}

void feedback_prepare_interpolation_tables(const struct feedback_props* fb_props) {
  int i, j, k, j1, j2, l;
  double sni[NXSNall][NZSN1Y], sn[NXSNall][NZSN][NMSN - 2], hn[NXSNall][NZSN][4];
  double sniilm[NMSN], snii[chem5_NXSN][NZSN][NMSN], snii2[chem5_NXSN][NZSN][NM];
  double SN1wd[NZSN1R][NM], SN1ms[NZSN1R][NM], SN1rg[NZSN1R][NM];
  double m[NM], imf[NZSN][NM];
  float dlm, norm, norm3;
  double m_l;
  FILE *fp;
  char buf[1000];
  double a1, a2, a3, a4, a5, a6, a7;
  double a8, a9, a10, a11, a12, a13, a14, a15, a16;
  double a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27;
  double a28, a29, a30;
  double effSNz[NZSN], temp, tempz;
  const double lfz[NZLF] = {.0001, .0002, .0005, .0010, .0020, .0040, .0080, .0200, .0500};
  double temp_ms, temp_rg;

  /* Massloss (Kobayashi et al. 2000)
   * sw[2][24]: progenitor mass, He core mass = NSorWD mass 
   */
  const double sw[2][NMSN] = {
      {40., 30., 25., 20., 18., 15., 13., 10., 9., 8.5, 8., 7.5, 7., 6.5, 6.0, 5.5, 5.0, 
       4.5, 4., 3.5, 3., 2.5, 2.25, 2., 1.9, 1.75, 1.5, 1.25, 1.0, 0.9, 0.7, 0.05},
      {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0.473, 0.459, 0.05}
  };
  const double sniie[9] = {30, 20, 10, 10, 1, 1, 1, 1, 1};
  const double sniiz[NZSN] = {0., 0., .001, .004, .008, .02, .05};
  const double sniz[NZSN1Y] = {0., .002, .01, .02, .04, .06, .10};
  const double effHNz[NZSN] = {0.5, 0.5, 0.5, 0.4, 0.232036142, 0.01, 0.01};

  const double feh_ia[NZSN1R] = {-1.1, -1.0, -0.69896996, 0., 0.39794001}; /* [Fe/H] */
  const double M_l1rg[NZSN1R] = {0.9, 0.9, 0.9, 0.9, 0.8};
  const double M_l1ms[NZSN1R] = {1.8, 1.8, 1.8, 1.8, 1.8};
  const double M_u1rg[NZSN1R] = {0.9, 1.5, 2.0, 3.0, 3.5};
  const double M_u1ms[NZSN1R] = {1.8, 2.6, 4.0, 5.5, 6.0};
  const double M_l1wd[NZSN1R] = {2.4, 2.5, 2.80103004, 3.5, 3.89794001};
  const double M_u1wd[NZSN1R] = {6.75, 6.75, 7.05, 7.95, 7.95};

  for (l = 0; l < NZSN1R; l++) SNLZ1R[l] = feh_ia[l];

  message("\n\nset nucleosynthesis yields for %i elements...", chem5_element_count);
  message("\nZ-dependent HN efficiency !!! ");
  message("effHN = %f %f", effHNz[0], effHNz[NZSN - 1]);
  message("\nZ-dependent SNIa model !!! ");
  message("b=(%.3f %.3f) [Fe/H] > %f", fb_props->b_rg, fb_props->b_ms, SNLZ1R[0]);
  message("\nZ-dependent SAGB!!!\n");

  sprintf(buf, "SN2SAGBYIELD.DAT");
  if ((fp = fopen(buf, "r")) == NULL) {
    fprintf(stderr, "Can not open File %s\n", buf);
    exit(-1);
  }

  for (j = 1; j < NZSN; j++) {
    fgets(buf, 1000, fp); /* metallicity */
    fgets(buf, 1000, fp); /* mass */
    fgets(buf, 1000, fp); /* energy */
    for (k = 0; k < NXSNall; k++) { /* k=0: masscut */
      fgets(buf, 1000, fp);
      sscanf(buf, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n",
            &a1, &a2, &a3, &a4, &a5, &a6, &a7,
            &a8, &a9, &a10, &a11, &a12, &a13, &a14, &a15, &a16,
            &a17, &a18, &a19, &a20, &a21, &a22, &a23, &a24, &a25, &a26, &a27, &a28, &a29, &a30);
      sn[k][j][0] = a30;
      sn[k][j][1] = a29;
      sn[k][j][2] = a28;
      sn[k][j][3] = a27;
      sn[k][j][4] = a26;
      sn[k][j][5] = a25;
      sn[k][j][6] = a24;
      sn[k][j][7] = a23;
      sn[k][j][8] = a22;
      sn[k][j][9] = a21;
      sn[k][j][10] = a20;
      sn[k][j][11] = a19;
      sn[k][j][12] = a18;
      sn[k][j][13] = a17;
      sn[k][j][14] = a16;
      sn[k][j][15] = a15;
      sn[k][j][16] = a14;
      sn[k][j][17] = a13;
      sn[k][j][18] = a12;
      sn[k][j][19] = a11;
      sn[k][j][20] = a10;
      sn[k][j][21] = a9;
      sn[k][j][22] = a8;
      sn[k][j][23] = a7;
      sn[k][j][24] = a6;
      sn[k][j][25] = a5;
      sn[k][j][26] = a4;
      sn[k][j][27] = a3;
      sn[k][j][28] = a2;
      sn[k][j][29] = a1;
    }
  }

  for (j = 1; j < NZSN; j++) {
    fgets(buf, 200, fp);
    fgets(buf, 200, fp);
    fgets(buf, 200, fp);
    for (k = 0; k < NXSNall; k++) {
        fgets(buf, 200, fp);
        sscanf(buf, "%lf%lf%lf%lf\n", &a1, &a2, &a3, &a4);
        hn[k][j][0] = a4;
        hn[k][j][1] = a3;
        hn[k][j][2] = a2;
        hn[k][j][3] = a1;
    }
  }
  
  fclose(fp);

  for (i = 0; i < NMSN-2; i++) {
    for (k = 0; k < NXSNall; k++) sn[k][0][i] = sn[k][1][i]; /* pop3 */
  }

  for (i = 0; i < 4; i++) {
    for (k = 0; k < NXSNall; k++) hn[k][0][i] = hn[k][1][i]; /* pop3 */
  }

  for (j = 0; j < NZSN; j++) effSNz[j] = 1.0 - effHNz[j];

  for (i = 0; i < 4; i++) {
    for (j = 0; j < NZSN; j++) {
      for (k = 0; k < NXSNall; k++) sn[k][j][i] = effSNz[j] * sn[k][j][i] + effHNz[j] * hn[k][j][i];
    }
  }

  for (i = 0; i < NMSN - 2; i++) {
    sniilm[i] = log10(sw[0][i]);
    for (j = 0; j < NZSN; j++) {
      temp = tempz = 0.0;
      for (k = 0; k < NXSNall; k++) {
        sn[k][j][i] /= sw[0][i];
        if (k > 0) temp += sn[k][j][i];
        if (k > 9) tempz += sn[k][j][i];
      }

      snii[2][j][i] = 1.0 - sn[0][j][i]; /* ejected mass */
      snii[1][j][i] = snii[2][j][i] - temp; /* unprocessed mass */
      if (snii[1][j][i] < 0.0) {
        snii[2][j][i] = temp;
        snii[1][j][i] = 0.0;
      }
      snii[3][j][i] = tempz; /* Z */
      snii[4][j][i] = sn[1][j][i] + sn[2][j][i]; /* H */
      snii[5][j][i] = sn[3][j][i] + sn[4][j][i]; /* He */
      snii[6][j][i] = sn[5][j][i] + sn[6][j][i]; /* Li */
      snii[7][j][i] = sn[7][j][i];             /* Be */
      snii[8][j][i] = sn[8][j][i] + sn[9][j][i]; /* B */
      snii[9][j][i] = sn[10][j][i] + sn[11][j][i]; /* C */
      snii[10][j][i] = sn[12][j][i] + sn[13][j][i]; /* N */
      snii[11][j][i] = sn[14][j][i] + sn[15][j][i] 
                        + sn[16][j][i]; /* O */
      snii[12][j][i] = sn[17][j][i];              /* F */
      snii[13][j][i] = sn[18][j][i] + sn[19][j][i] 
                        + sn[20][j][i]; /* Ne */
      snii[14][j][i] = sn[21][j][i];              /* Na */
      snii[15][j][i] = sn[22][j][i] + sn[23][j][i] 
                        + sn[24][j][i]; /* Mg */
      snii[16][j][i] = sn[25][j][i];              /* Al */
      snii[17][j][i] = sn[26][j][i] + sn[27][j][i] 
                        + sn[28][j][i]; /* Si */
      snii[18][j][i] = sn[29][j][i];              /* P */
      snii[19][j][i] = sn[30][j][i] + sn[31][j][i] 
                      + sn[32][j][i] + sn[33][j][i]; /* S */
      snii[20][j][i] = sn[34][j][i] + sn[35][j][i]; /* Cl */
      snii[21][j][i] = sn[36][j][i] + sn[37][j][i] 
                      + sn[38][j][i]; /* Ar */
      snii[22][j][i] = sn[39][j][i] + sn[40][j][i] 
                      + sn[41][j][i]; /* K */
      snii[23][j][i] = sn[42][j][i] + sn[43][j][i] 
                      + sn[44][j][i] + sn[45][j][i] 
                      + sn[46][j][i] + sn[47][j][i]; /* Ca */
      snii[24][j][i] = sn[48][j][i]; /* Sc */
      snii[25][j][i] = sn[49][j][i] + sn[50][j][i] 
                      + sn[51][j][i] + sn[52][j][i] 
                      + sn[53][j][i]; /* Ti */
      snii[26][j][i] = sn[54][j][i] + sn[55][j][i]; /* V */
      snii[27][j][i] = sn[56][j][i] + sn[57][j][i] 
                      + sn[58][j][i] + sn[59][j][i]; /* Cr */
      snii[28][j][i] = sn[60][j][i]; /* Mn */
      snii[29][j][i] = sn[61][j][i] + sn[62][j][i] 
                      + sn[63][j][i] + sn[64][j][i]; /* Fe */
      snii[30][j][i] = sn[65][j][i]; /* Co */
      snii[31][j][i] = sn[66][j][i] + sn[67][j][i] 
                      + sn[68][j][i] + sn[69][j][i] 
                      + sn[70][j][i]; /* Ni */
      snii[32][j][i] = sn[71][j][i] + sn[72][j][i]; /* Cu */
      snii[33][j][i] = sn[73][j][i] + sn[74][j][i] 
                      + sn[75][j][i] + sn[76][j][i] 
                      + sn[77][j][i]; /* Zn */
      snii[34][j][i] = sn[78][j][i] + sn[79][j][i]; /* Ga */
      snii[35][j][i] = sn[80][j][i] + sn[81][j][i] 
                      + sn[82][j][i] 
                      + sn[83][j][i]; /* Ge */
      snii[36][j][i] = 0.0;
    }
  }

  for (i = 9; i < NMSN-2; i++) {
    for (j = 0; j < NZSN; j++) {
      for (k = 29; k < 36; k++) snii[k][j][i] = 0.0; /* Fe in AGB 09.10.12 */
    }
  }

  for (i = NMSN - 2; i < NMSN; i++) /* 1 < Msun */{
    sniilm[i] = log10(sw[0][i]);
    for (j = 0; j < NZSN; j++) {
      snii[1][j][i] = 1.0 - sw[1][i]/sw[0][i];
      snii[2][j][i] = snii[1][j][i];
      for (k = 3; k < chem5_NXSN; k++) snii[k][j][i] = 0.0;
    }
  }

  for (j = 0; j < NZSN; j++) {
    for (i = 0; i < 9; i++) snii[0][j][i] = effSNz[j] * 1.0 + effHNz[j] * sniie[i]; /* Energy, 10-40 Msun */
    for (i = 9; i < NMSN; i++) snii[0][j][i] = 0.0;
  }

  SNLZ[0] = -999.; /* z=0 */
  SNLZ[1] = -10.; /* z=0 */
  for (j = 2; j < NZSN; j++) SNLZ[j] = log10(sniiz[j]);

  SNLZ1Y[0] = -999.; /* z=0 */
  for (j = 1; j < NZSN1Y; j++) SNLZ1Y[j] = log10(sniz[j]);

  /* SNIa yield (Kobayashi et al. 2020a, DDT) */
  sprintf(buf, "SN1YIELD_Z.DAT");
  if ((fp = fopen(buf, "r")) == NULL) {
      fprintf(stderr, "Can not open File %s\n", buf);
      exit(-1);
  }

  fgets(buf, 1000, fp); /* metallicity */
  for (k = 0; k < NXSNall; k++) { /* k=0: ejected mass */
    fgets(buf, 1000, fp);
    sscanf(buf, "%lf%lf%lf%lf%lf%lf%lf\n", &a1, &a2, &a3, &a4, &a5, &a6, &a7);
    sni[k][0] = a1;
    sni[k][1] = a2;
    sni[k][2] = a3;
    sni[k][3] = a4;
    sni[k][4] = a5;
    sni[k][5] = a6;
    sni[k][6] = a7;
  }
  fclose(fp);

  for (j = 0; j < NZSN1Y; j++) {
    temp = tempz = 0.0;
    for (k = 0; k < NXSNall; k++) {
      if (k > 0) temp += sni[k][j];
      if (k > 9) tempz += sni[k][j];
    }

    SN1E[0][j] = 1.3; /* energy */
    SN1E[1][j] = 0.0; /* unprocessed */
    SN1E[2][j] = temp; /* ejected */
    SN1E[3][j] = tempz; /* Z */
    SN1E[4][j] = sni[1][j] + sni[2][j]; /* H */
    SN1E[5][j] = sni[3][j] + sni[4][j]; /* He */
    SN1E[6][j] = sni[5][j] + sni[6][j]; /* Li */
    SN1E[7][j] = sni[7][j];             /* Be */
    SN1E[8][j] = sni[8][j] + sni[9][j]; /* B */
    SN1E[9][j] = sni[10][j] + sni[11][j]; /* C */
    SN1E[10][j] = sni[12][j] + sni[13][j]; /* N */
    SN1E[11][j] = sni[14][j] + sni[15][j] 
                  + sni[16][j]; /* O */
    SN1E[12][j] = sni[17][j];              /* F */
    SN1E[13][j] = sni[18][j] + sni[19][j] 
                  + sni[20][j]; /* Ne */
    SN1E[14][j] = sni[21][j];              /* Na */
    SN1E[15][j] = sni[22][j] + sni[23][j] 
                  + sni[24][j]; /* Mg */
    SN1E[16][j] = sni[25][j];              /* Al */
    SN1E[17][j] = sni[26][j] + sni[27][j] 
                  + sni[28][j]; /* Si */
    SN1E[18][j] = sni[29][j];              /* P */
    SN1E[19][j] = sni[30][j] + sni[31][j] 
                  + sni[32][j] + sni[33][j]; /* S */
    SN1E[20][j] = sni[34][j] + sni[35][j]; /* Cl */
    SN1E[21][j] = sni[36][j] + sni[37][j] 
                  + sni[38][j]; /* Ar */
    SN1E[22][j] = sni[39][j] + sni[40][j] 
                  + sni[41][j]; /* K */
    SN1E[23][j] = sni[42][j] + sni[43][j] 
                  + sni[44][j] + sni[45][j] 
                  + sni[46][j] + sni[47][j]; /* Ca */
    SN1E[24][j] = sni[48][j]; /* Sc */
    SN1E[25][j] = sni[49][j] + sni[50][j] 
                + sni[51][j] + sni[52][j] 
                + sni[53][j]; /* Ti */
    SN1E[26][j] = sni[54][j] + sni[55][j]; /* V */
    SN1E[27][j] = sni[56][j] + sni[57][j] 
                + sni[58][j] + sni[59][j]; /* Cr */
    SN1E[28][j] = sni[60][j]; /* Mn */
    SN1E[29][j] = sni[61][j] + sni[62][j] 
                + sni[63][j] + sni[64][j]; /* Fe */
    SN1E[30][j] = sni[65][j]; /* Co */
    SN1E[31][j] = sni[66][j] + sni[67][j] 
                + sni[68][j] + sni[69][j] 
                + sni[70][j]; /* Ni */
    SN1E[32][j] = sni[71][j] + sni[72][j]; /* Cu */
    SN1E[33][j] = sni[73][j] + sni[74][j] 
                + sni[75][j] + sni[76][j] 
                + sni[77][j]; /* Zn */
    SN1E[34][j] = sni[78][j] + sni[79][j]; /* Ga */
    SN1E[35][j] = sni[80][j] + sni[81][j] 
                + sni[82][j] + sni[83][j]; /* Ge */
    SN1E[36][j] = SN1E[29][j];
    for (k = 1; k < chem5_NXSN; k++) SN1E[k][j] *= fb_props->solar_mass_to_mass;
  }

  /* lifetime (Kobayashi et al. 2000) */
  sprintf(buf, "LIFETIME.DAT");
  if ((fp = fopen(buf, "r")) == NULL) {
      fprintf(stderr, "Can not open File %s\n", buf);
      exit(-1);
  }
  fgets(buf, 1000, fp);
  for (j = 0; j < NZLF; j++) {
    LFLZ[j] = log10(lfz[j]);
    fgets(buf, 1000, fp);
    fgets(buf, 1000, fp);
    fgets(buf, 1000, fp);
    for (i = 0; i < NMLF; i++) {
      fgets(buf, 1000, fp);
      sscanf(buf, "%lf%lf\n", &a1, &a2);
      LFLM[i] = log10(a1); /* sm */
      LFLT[j][i] = log10(a2); /* yr */
    }
  }
  fclose(fp);

  message(
    "total: %.2f %.1f  %.2e %.2e\n",
    fb_props->M_l,
    fb_props->M_u,
    feedback_life_time(fb_props->M_l, 0.02f),
    feedback_life_time(fb_props->M_u, 0.02f)
  );
  message(
    "SN2:   %.2f %.1f  %.2e %.2e  x=%.2f\n",
    fb_props->M_l2,
    fb_props->M_u2,
    feedback_life_time(fb_props->M_l2, 0.02f),
    feedback_life_time(fb_props->M_u2, 0.02f),
    fb_props->ximf
  );
  if (fb_props->zmax3 >= 0.0f) {
    message(
      "Pop3:  %.2f %.1f  %.2e %.2e  x=%.2f\n", 
      fb_props->M_l3, 
      fb_props->M_u3, 
      feedback_life_time(fb_props->M_l3, 0.02f), 
      feedback_life_time(fb_props->M_u3, 0.02f), 
      fb_props->ximf3
    );
  }

  /* IMF is normalized to 1 solar mass */
  if (fb_props->imf == 0) { /* Kroupa */
    if (fb_props->ximf == 1.) {
      norm = log10f(fb_props->M_u / 0.5f) * 0.5f 
              + (powf(0.5f, 0.7f) - powf(0.08f, 0.7f)) / 0.7f 
              + (powf(0.08f, 1.7f) - powf(fb_props->M_l, 1.7f)) / 1.7f / 0.08f;
    }
    else {
      norm = (powf(fb_props->M_u, 1.f - fb_props->ximf) - powf(0.5f, 1.f - fb_props->ximf)) 
              / (1.f - fb_props->ximf) * 0.5f + (powf(0.5f, 0.7f) - powf(0.08f, 0.7f)) 
              / 0.7f + (powf(0.08f, 1.7f) - powf(fb_props->M_l, 1.7f)) / 1.7f / 0.08f;
    }

    norm = 1.f / norm;
  }
  else { /* Chabrier, anything else */
    if (fb_props->ximf == 1.) {
      norm = 1.f / log(fb_props->M_u / fb_props->M_l);
    }
    else {
      norm = (1.f - fb_props->ximf) 
              / (pow(fb_props->M_u, (1.f - fb_props->ximf)) 
                  - pow(fb_props->M_l, (1.f - fb_props->ximf)));
    }
  }

  if (fb_props->ximf3 == 1.) {
    norm3 = 1.f / logf(fb_props->M_u3 / fb_props->M_l3);
  }
  else {
    norm3 = (1.f - fb_props->ximf3) / (powf(fb_props->M_u3, (1.f - fb_props->ximf3)) - powf(fb_props->M_l3, (1.f - fb_props->ximf3)));
  }

  /* IMF integration */
  dlm = (log10f(fb_props->M_u3) - log10f(fb_props->M_l)) / NM;
  for (i = 0; i < NM; i++) {
    SNLM[i] = log10f(fb_props->M_u3) - dlm * i;
    m[i] = powf(10.f, SNLM[i]);

    if (m[i] >= fb_props->M_l3) {
      imf[0][i] = powf(m[i], -fb_props->ximf3) * norm3;
    }
    else {
      imf[0][i] = 0.;
    }

    if (m[i] <= fb_props->M_u) {
      imf[1][i] = feedback_imf(fb_props, m[i]) * norm;
    }
    else {
      imf[1][i] = 0.;
    }

    for (l = 2; l < NZSN; l++) {
      imf[l][i] = imf[1][i];
    }

    j1 = 0;
    j2 = 1;
    for (j = 1; j < NMSN; j++) {
      j1 = j-1;
      j2 = j;
      if (sniilm[j] < SNLM[i]) break;
    }

    for (l = 0; l < NZSN; l++) {
      if (m[i] < fb_props->M_u2) { /* avoid extrapolation */
        for (k = 0; k < chem5_NXSN; k++) {
          snii2[k][l][i] = feedback_linear_interpolation(sniilm[j1], snii[k][l][j1], sniilm[j2], snii[k][l][j2], SNLM[i]);
          if (m[i] > fb_props->M_l2 && snii2[k][l][i] < 0.) snii2[k][l][i] = 0.0;
        }
      } else {
        snii2[0][l][i] = 0.0;
        snii2[1][l][i] = feedback_linear_interpolation(sniilm[j1], snii[1][l][j1], sniilm[j2], snii[1][l][j2], SNLM[i]);
        if (snii2[1][l][i] < 0.) snii2[1][l][i]=0.0;
        snii2[2][l][i] = snii2[1][l][i];

        for (k = 3; k < chem5_NXSN; k++) snii2[k][l][i] = 0.0;
      }
    }
  }

  for (i = 0; i < NM; i++) {
    for (l = 0; l < NZSN1R; l++) {
      SN1wd[l][i] = 0.;
      SN1ms[l][i] = 0.;
      SN1rg[l][i] = 0.;
      SN1R[l][i] = 0.;
    }

    for (l = 0; l < NZSN; l++) {
      SN2R[l][i] = 0.;
      SWR[l][i] = 0.;
      for (k = 0; k < chem5_NXSN; k++) SN2E[k][l][i] = 0.;
    }
  }
  for (i = 1; i < NM; i++) {
    for (l = 0; l < NZSN; l++) {
      if (l == 0) {
        m_l = max(fb_props->M_l2, fb_props->M_l3);
      }
      else {
        m_l = fb_props->M_l2;
      }

      if (m[i] > m_l) {
        SWR[l][i] = SWR[l][i-1] + sqrt(imf[l][i] * imf[l][i-1]) * dlm * log(10.);
      }
      else {
        SWR[l][i] = SWR[l][i-1];
      }

      if (l == 0) {
        m_l = fb_props->M_l3;
      }
      else {
        m_l = 0.0;
      }

      if (m[i] > m_l) {
        for (k = 1; k < 3; k++) {
          SN2E[k][l][i] = SN2E[k][l][i-1] + (snii2[k][l][i] + snii2[k][l][i-1]) / 2. 
                          * sqrt(m[i] * m[i-1] * imf[l][i] * imf[l][i-1]) * dlm * log(10.);
        }
      } else {
        for (k = 1; k < 3; k++) {
          SN2E[k][l][i] = SN2E[k][l][i-1];
        }
      }

      if (m[i] > m_l && m[i] < fb_props->M_u2) {
        for (k = 3; k < chem5_NXSN; k++) {
          SN2E[k][l][i] = SN2E[k][l][i-1] + (snii2[k][l][i] + snii2[k][l][i-1]) / 2. 
                          * sqrt(m[i] * m[i-1] * imf[l][i] * imf[l][i-1]) * dlm * log(10.);
        }
      } else {
        for (k = 3; k < chem5_NXSN; k++) {
          SN2E[k][l][i] = SN2E[k][l][i-1];
        }
      }

      if (l == 0) {
        m_l = max(fb_props->M_l2, fb_props->M_l3);
      }
      else {
        m_l = fb_props->M_l2;
      }

      if (m[i] > m_l && m[i] < fb_props->M_u2) {
        SN2R[l][i] = SN2R[l][i-1] + sqrt(imf[l][i] * imf[l][i-1]) * dlm * log(10.);
        SN2E[0][l][i] = SN2E[0][l][i-1] + (snii2[0][l][i] + snii2[0][l][i-1]) / 2. 
                        * sqrt(imf[l][i] * imf[l][i-1]) * dlm * log(10.);
      } else {
        SN2R[l][i] = SN2R[l][i-1];
        SN2E[0][l][i] = SN2E[0][l][i-1];
      }
    }
    for (l = 1; l < NZSN1R; l++) {
      if (m[i] > M_l1wd[l] && m[i] < M_u1wd[l]) {
        SN1wd[l][i] = SN1wd[l][i-1] + sqrt(imf[l][i] * imf[l][i-1]) * dlm * log(10.);
      }
      else {
        SN1wd[l][i] = SN1wd[l][i-1];
      }

      if (m[i] > M_l1ms[l] && m[i] < M_u1ms[l]) {
        SN1ms[l][i] = SN1ms[l][i-1] + pow(sqrt(m[i] * m[i-1]), -0.35) * dlm * log(10.);
      }
      else {
        SN1ms[l][i] = SN1ms[l][i-1];
      }

      if (m[i] > M_l1rg[l] && m[i] < M_u1rg[l]) {
        SN1rg[l][i] = SN1rg[l][i-1] + pow(sqrt(m[i] * m[i-1]), -0.35) * dlm * log(10.);
      }
      else {
        SN1rg[l][i] = SN1rg[l][i-1];
      }
    }
  }
  temp_ms = SN1ms[2][NM-1]; /* normalized at Z=0.004 */
  temp_rg = SN1rg[2][NM-1]; /* normalized at Z=0.004 */

  for (i = 0; i < NM; i++) {
    for (l = 1; l < NZSN1R; l++) {
      SN1ms[l][i] *= fb_props->b_ms / temp_ms;
      SN1rg[l][i] *= fb_props->b_rg / temp_rg;
      SN1R[l][i] = SN1wd[l][i] * (SN1ms[l][i] + SN1rg[l][i]);
      SN1R[l][i] /= fb_props->solar_mass_to_mass;
    }

    for (l = 1; l < NZSN1Y; l++) SN1E[0][1] *= (1.e51 / fb_props->energy_to_cgs);
    /* convert gadget unit */
    SNLM[i] += log10(fb_props->solar_mass_to_mass);

    for (l = 0; l < NZSN; l++) {
      SN2R[l][i] /= fb_props->solar_mass_to_mass;
      SWR[l][i] /= fb_props->solar_mass_to_mass;
      SN2E[0][l][i] *= (1.e51 / fb_props->energy_to_cgs / fb_props->solar_mass_to_mass);
    }
  }

  message("Done Chem5 setup.");
}

double feedback_linear_interpolation(const double x1, 
                                     const double y1, 
                                     const double x2, 
                                     const double y2, 
                                     const double x) {
  return (y2 - y1) / (x2 - x1) * (x - x1) + y1;
}

/**
 * @brief calculates stellar mass in spart that died over the timestep, calls
 * functions to calculate feedback due to SNIa, SNII and AGB
 *
 * @param feedback_props feedback_props data structure
 * @param phys_const The physical constants in internal units.
 * @param cosmo The cosmological model.
 * @param sp spart that we're evolving
 * @param us unit_system data structure
 * @param age age of spart at beginning of step
 * @param dt length of current timestep
 * @param ti_begin The current integer time (for random number hashing).
 */
void compute_stellar_evolution(const struct feedback_props* feedback_props,
                               const struct phys_const* phys_const,
                               const struct cosmology* cosmo, struct spart* sp,
                               const struct unit_system* us, const double age,
                               const double dt, const integertime_t ti_begin) {

  TIMER_TIC;

#ifdef SIMBA_DEBUG_CHECKS
  if (sp->feedback_data.to_collect.ngb_rho <= 0) {
    warning("Star %lld with mass %g has no neighbors!",
            sp->id, sp->mass);
    return;
  }
#endif

#ifdef SWIFT_DEBUG_CHECKS
  if (age < 0.f) error("Negative age for a star.");
  if (sp->feedback_data.to_collect.ngb_rho <= 0)
    error("Star %lld with mass %g has no neighbors!",
            sp->id, sp->mass);
  if (sp->count_since_last_enrichment != 0 && engine_current_step > 0)
    error("Computing feedback on a star that should not");
#endif

  /* Properties collected in the stellar density loop. */
  const float ngb_gas_mass = sp->feedback_data.to_collect.ngb_mass;

  /* Check if there are neighbours, otherwise exit */
  if (ngb_gas_mass == 0.f || sp->density.wcount * pow_dimension(sp->h) < 1e-4) {
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  /* Update the enrichment weights */
  const float enrichment_weight_inv =
      sp->feedback_data.to_collect.enrichment_weight_inv;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_collect.enrichment_weight_inv < 0.)
    error("Negative inverse weight!");
#endif

  /* Now we start filling the data structure for information to apply to the
   * particles. Do _NOT_ read from the to_collect substructure any more. */

  /* Zero all the output fields */
  feedback_reset_feedback(sp, feedback_props);

  /* Update the weights used for distribution */
  const float enrichment_weight =
      (enrichment_weight_inv != 0.f) ? 1.f / enrichment_weight_inv : 0.f;
  sp->feedback_data.to_distribute.enrichment_weight = enrichment_weight;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.enrichment_weight < 0.)
    error("Negative weight!");
#endif

  int elem;
  float ejecta_energy = 0.f;
  float ejecta_mass = 0.f;
  float ejecta_unprocessed = 0.f;
  float ejecta_metal_mass[chemistry_element_count];
  for (elem = 0; elem < chemistry_element_count; elem++) ejecta_metal_mass[elem] = 0.f;

  feedback_get_ejecta_from_star_particle(sp, age, feedback_props, dt, 
                                         &ejecta_energy, 
                                         &ejecta_mass, &ejecta_unprocessed, 
                                         ejecta_metal_mass);

#ifdef SIMBA_DEBUG_CHECKS
  if (ejecta_mass > sp->mass) {
    error("Star particle %lld with mass %g is trying to give away more mass (Mejecta=%g) than it has!",
          sp->id, sp->mass, ejecta_mass);
  }
#endif

  /* D. Rennehan: Do some magic that I still don't understand 
   * https://www.youtube.com/watch?v=cY2xBNWrBZ4
   */
  float dum = 0.f;
  int flag_negative = 0;
  for (elem = 0; elem < chemistry_element_count; elem++) {
    dum = ejecta_unprocessed * sp->chemistry_data.metal_mass_fraction[elem];
    ejecta_metal_mass[feedback_props->element_index_conversions[elem]] += dum;
    if (ejecta_metal_mass[feedback_props->element_index_conversions[elem]] < 0.f) {
      ejecta_metal_mass[feedback_props->element_index_conversions[elem]] = 0.f;
      flag_negative = 1;
      /* Do not break here, we need the zeroed elements where negative */
    }
  }
  
  /* If ANY element ended up negative we recompute everything */
  if (flag_negative) {
    ejecta_mass = 0.f;
    ejecta_metal_mass[chem5_element_Z] = 0.f;
    for (elem = chem5_element_H; elem < chem5_element_Zn; elem++) ejecta_mass += ejecta_metal_mass[elem];
    for (elem = chem5_element_C; elem < chem5_element_Zn; elem++) ejecta_metal_mass[chem5_element_Z] += ejecta_metal_mass[elem];
  }

  sp->feedback_data.to_distribute.total_metal_mass = ejecta_mass;
  for (elem = 0; elem < chemistry_element_count; elem++) {
    sp->feedback_data.to_distribute.metal_mass[elem] = ejecta_metal_mass[feedback_props->element_index_conversions[elem]];

    /* Only count real metals in total_metal_mass */
    if (elem != chemistry_element_H && elem != chemistry_element_He) {
      sp->feedback_data.to_distribute.total_metal_mass += ejecta_metal_mass[feedback_props->element_index_conversions[elem]];
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.mass != 0.f)
    error("Injected mass will be lost");
#endif

  /* Compute the total mass to distribute */
  sp->feedback_data.to_distribute.mass = ejecta_mass;

  /* Compute energy change due to kinetic energy of ejectas */
  sp->feedback_data.to_distribute.energy += ejecta_energy;

  /* Compute energy change due to kinetic energy of the star */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass * 0.5f *
      (sp->v[0] * sp->v[0] + sp->v[1] * sp->v[1] + sp->v[2] * sp->v[2]) *
      cosmo->a2_inv;

  TIMER_TOC(timer_do_star_evol);
}

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
void feedback_props_init(struct feedback_props* fp,
                         const struct phys_const* phys_const,
                         const struct unit_system* us,
                         struct swift_params* params,
                         const struct hydro_props* hydro_props,
                         const struct cosmology* cosmo) {

  /* Common conversions ------------------------------------------------- */

  /* Calculate internal mass to solar mass conversion factor */
  const double Msun_cgs = phys_const->const_solar_mass *
                          units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  fp->mass_to_solar_mass = unit_mass_cgs / Msun_cgs;
  fp->solar_mass_in_g = Msun_cgs;
  fp->solar_mass_to_mass = 1. / fp->mass_to_solar_mass;

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  fp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);

  /* Calculate conversion factor from rho to n_H
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  fp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  fp->kms_to_internal = 1.0e5f / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  fp->kms_to_cms = 1.e5;

  fp->time_to_Myr = units_cgs_conversion_factor(us, UNIT_CONV_TIME) /
      (1.e6f * 365.25f * 24.f * 60.f * 60.f);

  /* Convert to Myr first, then multiply by a factor of 1e6 yr / 1 Myr */
  fp->time_to_yr = fp->time_to_Myr * 1.e6f;

  fp->length_to_kpc = 
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH) / 3.08567758e21f;

  fp->energy_to_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Constant Chem5 parameters ----------------------------------------------- */

  /* Solar values */
  fp->H_mf = 7.35224e-1f;
  fp->He_mf = 2.50274e-1f;
  fp->Z_mf = 0.0144404378f;
  fp->O_mf = 0.675327e-2f + 0.272594e-5f + 0.152311e-4f;
  fp->Fe_mf = 0.733849e-4f + 0.119465e-2f + 0.280824e-4f + 0.380282e-5f;

  /* supernova energy in foe (D. Rennehan: What does foe mean??) */
  fp->E_sw = 0.2 * (1.e51 / fp->energy_to_cgs);
  fp->E_sn1 = 1.3 * (1.e51 / fp->energy_to_cgs);

  fp->imf = parser_get_param_int(params, "KIARAFeedback:imf");
  
  /* Kroupa IMF || Chabrier IMF */
  if (fp->imf == 0 || fp->imf == 1) {
    fp->ximf = 1.3f;
    fp->M_u = 120.f;
    fp->M_l = 0.01f;
  }
  else { /* No idea what this is? */
    fp->ximf = 1.35f;
    fp->M_u = 120.f;
    fp->M_l = 0.07f;
  }

  fp->ximf3 = 1.35f;
  fp->M_u3 = 120.f; /* >= M_u */
  fp->M_l3 = 20.f; /* >= M_l */
  fp->zmax3 = -999.f;
  fp->M_u2 = 50.f;
  fp->M_l2 = 8.f;
  fp->b_rg = 0.02f; /* binary parameter for SNIa */
  fp->b_ms = 0.04f; /* binary parameter for SNIa */

  /* Chem5 indices to Swift */
  fp->element_index_conversions[chemistry_element_H] = chem5_element_H;
  fp->element_index_conversions[chemistry_element_He] = chem5_element_He;
  fp->element_index_conversions[chemistry_element_C] = chem5_element_C;
  fp->element_index_conversions[chemistry_element_N] = chem5_element_N;
  fp->element_index_conversions[chemistry_element_O] = chem5_element_O;
  fp->element_index_conversions[chemistry_element_Ne] = chem5_element_Ne;
  fp->element_index_conversions[chemistry_element_Mg] = chem5_element_Mg;
  fp->element_index_conversions[chemistry_element_Si] = chem5_element_Si;
  fp->element_index_conversions[chemistry_element_Fe] = chem5_element_Fe;

  /* Main operation modes ------------------------------------------------- */

  fp->with_HN_energy_from_chem5 =
      parser_get_param_int(params, "KIARAFeedback:use_HN_energy_from_chem5");

  fp->with_SNII_energy_from_chem5 =
      parser_get_param_int(params, "KIARAFeedback:use_SNII_energy_from_chem5");

  fp->with_SNIa_energy_from_chem5 =
      parser_get_param_int(params, "KIARAFeedback:use_SNIa_energy_from_chem5");

  /* Properties of Simba kinetic winds -------------------------------------- */

  fp->FIRE_velocity_normalization =
      parser_get_param_double(params, "KIARAFeedback:FIRE_velocity_normalization");
  fp->FIRE_velocity_slope =
      parser_get_param_double(params, "KIARAFeedback:FIRE_velocity_slope");
  fp->FIRE_eta_normalization =
      parser_get_param_double(params, "KIARAFeedback:FIRE_eta_normalization");
  fp->FIRE_eta_break =
      parser_get_param_double(params, "KIARAFeedback:FIRE_eta_break_Msun");
  fp->FIRE_eta_break *= fp->solar_mass_to_mass;
  fp->FIRE_eta_lower_slope =
      parser_get_param_double(params, "KIARAFeedback:FIRE_eta_lower_slope");
  fp->FIRE_eta_upper_slope =
      parser_get_param_double(params, "KIARAFeedback:FIRE_eta_upper_slope");

  fp->early_stellar_mass_norm =
      parser_get_param_double(params, "KIARAFeedback:early_stellar_mass_norm_Msun");
  fp->early_stellar_mass_norm *= fp->solar_mass_to_mass;
  fp->early_wind_suppression_enabled = 
      parser_get_param_int(params, "KIARAFeedback:early_wind_suppression_enabled");
  fp->early_wind_suppression_scale_factor =
      parser_get_param_double(params, 
          "KIARAFeedback:early_wind_suppression_scale_factor");
  fp->early_wind_suppression_slope =
      parser_get_param_double(params, "KIARAFeedback:early_wind_suppression_slope");

  fp->minimum_galaxy_stellar_mass =
      parser_get_param_double(params, "KIARAFeedback:minimum_galaxy_stellar_mass_Msun");
  fp->minimum_galaxy_stellar_mass *= fp->solar_mass_to_mass;

  fp->kick_velocity_scatter =
      parser_get_param_double(params, "KIARAFeedback:kick_velocity_scatter");

  fp->wind_decouple_time_factor = parser_get_param_double(
      params, "KIARAFeedback:wind_decouple_time_factor");

  fp->cold_wind_internal_energy = parser_get_opt_param_double(
      params, "KIARAFeedback:cold_wind_temperature_K", 1.e4);
  if (fp->cold_wind_internal_energy <= 0.) {
    error("KIARAFeedback:cold_wind_temperature_K must be strictly positive.");
  }
  /* Convert Kelvin to internal energy and internal units */
  fp->cold_wind_internal_energy *= fp->temp_to_u_factor / 
                                   units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Initialise the yield/mass tables ------------------------------------------- */

  feedback_prepare_interpolation_tables(fp);


  /* Output some information to the people -------------------------------------- */
  
  if (engine_rank == 0) {
    message("Feedback model is KIARA");
    message("Feedback FIRE velocity normalization: %g", 
            fp->FIRE_velocity_normalization);
    message("Feedback FIRE velocity slope: %g", fp->FIRE_velocity_slope);
    message("Feedback velocity scatter: %g", fp->kick_velocity_scatter);
    message("Feedback FIRE eta normalization: %g", fp->FIRE_eta_normalization);
    message("Feedback FIRE eta break: %g", fp->FIRE_eta_break);
    message("Feedback FIRE eta upper slope: %g", fp->FIRE_eta_upper_slope);
    message("Feedback FIRE eta lower slope: %g", fp->FIRE_eta_lower_slope);

    message("Feedback early suppression enabled: %d", 
            fp->early_wind_suppression_enabled);
    
    if (fp->early_wind_suppression_enabled) {
      message("Feedback early stellar mass norm: %g Msun", 
              fp->early_stellar_mass_norm);
      message("Feedback early suppression scale factor: %g", 
              fp->early_wind_suppression_scale_factor);
      message("Feedback early suppression slope: %g", fp->early_wind_suppression_slope);
    }

    message("Feedback use Chem5 SNII energy: %d", fp->with_SNII_energy_from_chem5);
    message("Feedback use Chem5 SNIa energy: %d", fp->with_SNIa_energy_from_chem5);
  }
}

/**
 * @brief Clean-up the memory allocated for the feedback routines
 *
 * We simply free all the arrays.
 *
 * @param fp the feedback data structure.
 */
void feedback_clean(struct feedback_props* fp) { }

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream) {
  /* To make sure everything is restored correctly, we zero all the pointers to
     tables. If they are not restored correctly, we would crash after restart on
     the first call to the feedback routines. Helps debugging. */
  struct feedback_props feedback_copy = *feedback;

  restart_write_blocks((void*)&feedback_copy, sizeof(struct feedback_props), 1,
                       stream, "feedback", "feedback function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Read the structure from the stream and restore the feedback tables by
 * re-reading them.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream) {
  restart_read_blocks((void*)feedback, sizeof(struct feedback_props), 1, stream,
                      NULL, "feedback function");
}
