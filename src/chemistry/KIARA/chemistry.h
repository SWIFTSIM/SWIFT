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
#ifndef SWIFT_CHEMISTRY_KIARA_H
#define SWIFT_CHEMISTRY_KIARA_H

/**
 * @file src/chemistry/KIARA/chemistry.h
 * @brief Empty infrastructure for the cases without chemistry function
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <signal.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "timestep_sync_part.h"
#include "units.h"

/**
 * @brief Initializes summed particle quantities for the firehose wind model
 *
 * This is called from chemistry_init_part.
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void
firehose_init_ambient_quantities(struct part *restrict p,
                                 const struct chemistry_global_data *cd) {

  struct chemistry_part_data *cpd = &p->chemistry_data;

  cpd->w_ambient = 0.f;
  cpd->rho_ambient = 0.f;
  cpd->u_ambient = 0.f;
}

__attribute__((always_inline)) INLINE static void logger_windprops_printprops(
    struct part *pi, const struct cosmology *cosmo,
    const struct chemistry_global_data *cd) {

  /* Ignore COUPLED particles */
  if (!pi->decoupled) return;

#ifdef CHEMISTRY_OUTPUT_FIREHOSE_LOG
  /* Print wind properties */
  const float length_convert = cosmo->a * cd->length_to_kpc;
  const float rho_convert = cosmo->a3_inv * cd->rho_to_n_cgs;
  const float u_convert =
      cosmo->a_factor_internal_energy / cd->temp_to_u_factor;

  message(
      "FIREHOSE: z=%.3f id=%lld Mgal=%g h=%g T=%g rho=%g Rs=%g Z=%g tdel=%g "
      "Ndec=%d rhoamb=%g Tamb=%g tcmix=%g\n",
      cosmo->z, pi->id,
      (pi->galaxy_data.gas_mass + pi->galaxy_data.stellar_mass) *
          cd->mass_to_solar_mass,
      pi->h * cosmo->a * cd->length_to_kpc,
      hydro_get_drifted_comoving_internal_energy(pi) * u_convert,
      pi->rho * rho_convert, pi->chemistry_data.radius_stream * length_convert,
      pi->chemistry_data.metal_mass_fraction_total,
      pi->chemistry_data.decoupling_delay_time * cd->time_to_Myr,
      pi->chemistry_data.number_of_times_decoupled,
      pi->chemistry_data.rho_ambient * cd->rho_to_n_cgs * cosmo->a3_inv,
      pi->chemistry_data.u_ambient * cosmo->a_factor_internal_energy /
          cd->temp_to_u_factor,
      pi->cooling_data.mixing_layer_cool_time);
#endif

  return;
}

/**
 * @brief Finishes up ambient quantity calculation for the firehose wind model
 *
 * This is called from chemistry_end_density
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void
firehose_end_ambient_quantities(struct part *restrict p,
                                struct xpart *restrict xp,
                                const struct chemistry_global_data *cd,
                                const struct cosmology *cosmo) {

  const float u_floor = cd->firehose_u_floor / cosmo->a_factor_internal_energy;
  const float rho_max =
      cd->firehose_ambient_rho_max * cosmo->a * cosmo->a * cosmo->a;

  /* No ambient properties for non-wind particles */
  if (p->decoupled) {

    /* Some smoothing length multiples. */
    const float h = p->h;
    const float h_inv = 1.0f / h;                 /* 1/h */
    const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

#ifdef FIREHOSE_DEBUG_CHECKS
    message(
        "FIREHOSE_prelim: id=%lld rhoamb=%g wamb=%g uamb=%g"
        "h=%g h_inv=%g h_inv_dim=%g",
        p->id, p->chemistry_data.rho_ambient, p->chemistry_data.w_ambient,
        p->chemistry_data.u_ambient, h, h_inv, h_inv_dim);
#endif

    /* h_inv_dim can sometimes lead to rho_ambient -> 0 after normalization */
    p->chemistry_data.rho_ambient *= h_inv_dim;

    if (p->chemistry_data.rho_ambient > 0.f) {
      p->chemistry_data.u_ambient *= h_inv_dim / p->chemistry_data.rho_ambient;
    } else {
      p->chemistry_data.rho_ambient = hydro_get_comoving_density(p);
      p->chemistry_data.u_ambient = u_floor;
    }

#ifdef FIREHOSE_DEBUG_CHECKS
    message("FIREHOSE_lim: id=%lld rhoamb=%g wamb=%g uamb=%g ufloor=%g\n",
            p->id, p->chemistry_data.rho_ambient, p->chemistry_data.w_ambient,
            p->chemistry_data.u_ambient,
            cd->firehose_u_floor / cd->temp_to_u_factor);
#endif
  } else {
    /* Set them to reasonable values for non-wind, just in case */
    p->chemistry_data.rho_ambient = hydro_get_comoving_density(p);
    p->chemistry_data.u_ambient = hydro_get_drifted_comoving_internal_energy(p);
  }

  /* Limit ambient density to the user settings */
  p->chemistry_data.rho_ambient = min(p->chemistry_data.rho_ambient, rho_max);
  p->chemistry_data.u_ambient = max(p->chemistry_data.u_ambient, u_floor);
#ifdef FIREHOSE_DEBUG_CHECKS
  if (p->decoupled) {
    message(
        "FIREHOSE_AMB: z=%g id=%lld nH=%g nHamb=%g u=%g uamb=%g T=%g "
        "Tamb=%g tcool=%g",
        cosmo->z, p->id, p->rho * cd->rho_to_n_cgs * cosmo->a3_inv,
        p->chemistry_data.rho_ambient * cd->rho_to_n_cgs * cosmo->a3_inv,
        hydro_get_drifted_comoving_internal_energy(p),
        p->chemistry_data.u_ambient,
        hydro_get_drifted_comoving_internal_energy(p) *
            cosmo->a_factor_internal_energy / cd->temp_to_u_factor,
        p->chemistry_data.u_ambient * cosmo->a_factor_internal_energy /
            cd->temp_to_u_factor,
        p->cooling_data.mixing_layer_cool_time);
  }
#endif

#ifdef FIREHOSE_DEBUG_CHECKS
  logger_windprops_printprops(p, cosmo, cd);
#endif
}

/**
 * @brief Return a string containing the name of a given #chemistry_element.
 */
__attribute__((always_inline)) INLINE static const char *
chemistry_get_element_name(enum chemistry_element elem) {

  static const char *chemistry_element_names[chemistry_element_count] = {
      "Hydrogen",  "Helium",  "Carbon", "Nitrogen", "Oxygen", "Neon",
      "Magnesium", "Silicon", "Sulfur", "Calcium",  "Iron"};

  return chemistry_element_names[elem];
}

/**
 * @brief Prepares a particle for the smooth metal calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various smooth metallicity tasks
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void chemistry_init_part(
    struct part *restrict p, const struct chemistry_global_data *cd) {

  struct chemistry_part_data *cpd = &p->chemistry_data;

  /* Reset the shear tensor */
  for (int i = 0; i < 3; i++) {
    cpd->shear_tensor[i][0] = 0.f;
    cpd->shear_tensor[i][1] = 0.f;
    cpd->shear_tensor[i][2] = 0.f;

    /* Accumulated velocity from the firehose model */
    cpd->dv[i] = 0.f;
  }

  /* Reset the diffusion. */
  cpd->diffusion_coefficient = 0.f;

  /* Reset the changes for the accumulated properties */
  cpd->dZ_dt_total = 0.f;
  cpd->du = 0.;
  cpd->dm = 0.f;
  cpd->dm_dust = 0.f;
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    cpd->dZ_dt[elem] = 0.f;
    cpd->dm_Z[elem] = 0.f;
    cpd->dm_dust_Z[elem] = 0.f;
  }

#if COOLING_GRACKLE_MODE >= 2
  cpd->local_sfr_density = 0.f;
#endif

  if (cd->use_firehose_wind_model) {
    firehose_init_ambient_quantities(p, cd);
  }
}

/**
 * @brief Finishes the smooth metal calculation.
 *
 * Multiplies the metallicity and number of neighbours by the
 * appropiate constants and add the self-contribution term.
 *
 * This function requires the #hydro_end_density to have been called.
 *
 * @param p The particle to act upon.
 * @param cd #chemistry_global_data containing chemistry informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_density(
    struct part *restrict p, struct xpart *restrict xp,
    const struct chemistry_global_data *cd, const struct cosmology *cosmo) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  struct chemistry_part_data *cpd = &p->chemistry_data;

  /* If diffusion is on, finish up shear tensor & particle's diffusion coeff */
  if (cd->diffusion_flag == 1 && cd->C_Smagorinsky > 0.f) {
    const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */
    const float rho = hydro_get_comoving_density(p);

    /* convert the shear factor into physical */
    const float factor_shear = h_inv_dim_plus_one * cosmo->a2_inv / rho;
    for (int k = 0; k < 3; k++) {
      cpd->shear_tensor[k][0] *= factor_shear;
      cpd->shear_tensor[k][1] *= factor_shear;
      cpd->shear_tensor[k][2] *= factor_shear;
    }

    /* Compute the trace over 3 and add the hubble flow. */
    float trace_3 = 0.f;
    for (int i = 0; i < 3; i++) {
      cpd->shear_tensor[i][i] += cosmo->H;
      trace_3 += cpd->shear_tensor[i][i];
    }
    trace_3 /= 3.f;

    float shear_tensor[3][3] = {{0.f}};
    for (int i = 0; i < 3; i++) {
      /* Make the tensor symmetric. */
      float avg = 0.5f * (cpd->shear_tensor[i][0] + cpd->shear_tensor[0][i]);
      shear_tensor[i][0] = avg;
      shear_tensor[0][i] = avg;

      avg = 0.5f * (cpd->shear_tensor[i][1] + cpd->shear_tensor[1][i]);
      shear_tensor[i][1] = avg;
      shear_tensor[1][i] = avg;

      avg = 0.5f * (cpd->shear_tensor[i][2] + cpd->shear_tensor[2][i]);
      shear_tensor[i][2] = avg;
      shear_tensor[2][i] = avg;

      /* Remove the trace. */
      shear_tensor[i][i] -= trace_3;
    }

    /* Compute the norm. */
    float velocity_gradient_norm = 0.f;
    for (int i = 0; i < 3; i++) {
      velocity_gradient_norm += shear_tensor[i][0] * shear_tensor[i][0];
      velocity_gradient_norm += shear_tensor[i][1] * shear_tensor[i][1];
      velocity_gradient_norm += shear_tensor[i][2] * shear_tensor[i][2];

      /* Copy the final values into the particle quantity */
      cpd->shear_tensor[i][0] = shear_tensor[i][0];
      cpd->shear_tensor[i][1] = shear_tensor[i][1];
      cpd->shear_tensor[i][2] = shear_tensor[i][2];
    }

    velocity_gradient_norm = sqrtf(velocity_gradient_norm);

    /* Never set D for wind, or ISM particles */
    if (!(p->decoupled) && !(p->cooling_data.subgrid_temp > 0.f)) {

      /* Rennehan: Limit to maximum resolvable velocity scale */
      const float v_phys =
          sqrtf(p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2]) *
          cosmo->a_inv;
      const float h_phys = cosmo->a * p->h * kernel_gamma;
      const float vel_norm_phys_max = 0.5f * v_phys / h_phys;
      if (velocity_gradient_norm > vel_norm_phys_max) {
        velocity_gradient_norm = vel_norm_phys_max;
      }

      /* Compute the diffusion coefficient in physical coordinates.
       * The norm is already in physical coordinates.
       * kernel_gamma is necessary (see Rennehan 2021)
       */
      const float rho_phys = hydro_get_physical_density(p, cosmo);
      const float smag_length_scale = cd->C_Smagorinsky * h_phys;
      const float D_phys = rho_phys * smag_length_scale * smag_length_scale *
                           velocity_gradient_norm;

      cpd->diffusion_coefficient = D_phys;
    }
  } /* end Smagorinsky diffusion */

#if COOLING_GRACKLE_MODE >= 2
  /* Finish SFR density calculation */
  cpd->local_sfr_density *= h_inv_dim;
#endif

  if (cd->use_firehose_wind_model) {
    firehose_end_ambient_quantities(p, xp, cd, cosmo);
  }
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_has_no_neighbours(struct part *restrict p,
                                 struct xpart *restrict xp,
                                 const struct chemistry_global_data *cd,
                                 const struct cosmology *cosmo) {

  /* Just make all the smoothed fields default to the un-smoothed values */
  struct chemistry_part_data *cpd = &p->chemistry_data;
  /* Reset the shear tensor */
  for (int i = 0; i < 3; i++) {
    cpd->shear_tensor[i][0] = 0.f;
    cpd->shear_tensor[i][1] = 0.f;
    cpd->shear_tensor[i][2] = 0.f;
  }

  /* Reset the diffusion. */
  cpd->diffusion_coefficient = 0.f;

  /* Reset the change in metallicity */
  cpd->dZ_dt_total = 0.f;
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    cpd->dZ_dt[elem] = 0.f;
  }

#if COOLING_GRACKLE_MODE >= 2
  /* If there is no nearby SF, set to zero */
  cpd->local_sfr_density = 0.f;
#endif
}

/**
 * @brief Sets the chemistry properties of the (x-)particles to a valid start
 * state.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The current cosmological model.
 * @param data The global chemistry information.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct chemistry_global_data *data, struct part *restrict p,
    struct xpart *restrict xp) {

  /* Initialize mass fractions for total metals and each metal individually */
  if (data->initial_metal_mass_fraction_total != -1) {
    p->chemistry_data.metal_mass_fraction_total =
        data->initial_metal_mass_fraction_total;

    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      p->chemistry_data.metal_mass_fraction[elem] =
          data->initial_metal_mass_fraction[elem];
    }
  }
  chemistry_init_part(p, data);

  if (data->use_firehose_wind_model) {
    firehose_init_ambient_quantities(p, data);
  }
}

/**
 * @brief Sets the chemistry properties of the sparticles to a valid start
 * state.
 *
 * @param data The global chemistry information.
 * @param sp Pointer to the sparticle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_spart(
    const struct chemistry_global_data *data, struct spart *restrict sp) {

  /* Initialize mass fractions for total metals and each metal individually */
  if (data->initial_metal_mass_fraction_total != -1) {
    sp->chemistry_data.metal_mass_fraction_total =
        data->initial_metal_mass_fraction_total;

    for (int elem = 0; elem < chemistry_element_count; ++elem)
      sp->chemistry_data.metal_mass_fraction[elem] =
          data->initial_metal_mass_fraction[elem];
  }
}

/**
 * @brief Sets the chemistry properties of the sink particles to a valid start
 * state.
 *
 * @param data The global chemistry information.
 * @param sink Pointer to the sink particle data.
 * Required by space_first_init.c
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_sink(
    const struct chemistry_global_data *data, struct sink *restrict sink) {}

/**
 * @brief Initialises the chemistry properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param data The properties to initialise.
 */
static INLINE void chemistry_init_backend(struct swift_params *parameter_file,
                                          const struct unit_system *us,
                                          const struct phys_const *phys_const,
                                          struct chemistry_global_data *data) {

  /* Set some useful unit conversions */
  const double Msun_cgs = phys_const->const_solar_mass *
                          units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  data->mass_to_solar_mass = unit_mass_cgs / Msun_cgs;
  data->temp_to_u_factor =
      phys_const->const_boltzmann_k /
      (hydro_gamma_minus_one * phys_const->const_proton_mass *
       units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE));
  data->T_to_internal =
      1. / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  const double X_H = 0.752;
  data->rho_to_n_cgs =
      (X_H / phys_const->const_proton_mass) *
      units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);
  data->kms_to_internal =
      1.0e5 / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);
  data->time_to_Myr = units_cgs_conversion_factor(us, UNIT_CONV_TIME) /
                      (1.e6 * 365.25 * 24. * 60. * 60.);
  data->length_to_kpc =
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH) / 3.08567758e21;

  /* Is metal diffusion turned on? */
  data->diffusion_flag =
      parser_get_param_int(parameter_file, "KIARAChemistry:diffusion_on");

  /* Read the diffusion coefficient */
  data->C_Smagorinsky = parser_get_opt_param_float(
      parameter_file, "KIARAChemistry:diffusion_coefficient", 0.23f);

  /* Time-step restriction to <= 0.15*rho*h^2 / D from
     Parshikov & Medin 2002 equation 41 */
  data->diffusion_beta = parser_get_opt_param_float(
      parameter_file, "KIARAChemistry:diffusion_beta", 0.1f);
  if (data->diffusion_beta < 0.f || data->diffusion_beta > 0.1f) {
    error("diffusion_beta must be >= 0 and <= 0.1");
  }

  data->time_step_min = parser_get_opt_param_float(
      parameter_file, "KIARAChemistry:minimum_timestep_Myr", 0.1f);
  data->time_step_min /= data->time_to_Myr;
  if (data->time_step_min < 0.f) {
    error("time_step_min must be > 0");
  }

  data->max_fractional_Z_transfer = parser_get_opt_param_float(
      parameter_file, "KIARAChemistry:max_fractional_Z_transfer", 0.25f);
  if (data->max_fractional_Z_transfer < 0.f ||
      data->max_fractional_Z_transfer > 1.f) {
    error("diffusion_beta must be >= 0 and <= 1");
  }

  /* Are we using the firehose wind model? */
  data->use_firehose_wind_model = parser_get_opt_param_int(
      parameter_file, "KIARAChemistry:use_firehose_wind_model", 0);

  if (data->use_firehose_wind_model) {
    /* Firehose model parameters */
    data->firehose_ambient_rho_max = parser_get_opt_param_float(
        parameter_file, "KIARAChemistry:firehose_nh_ambient_max_cgs", 0.1f);
    data->firehose_ambient_rho_max /= data->rho_to_n_cgs;

    data->firehose_u_floor = parser_get_opt_param_float(
        parameter_file, "KIARAChemistry:firehose_temp_floor", 1.e4f);
    data->firehose_u_floor *= data->temp_to_u_factor * data->T_to_internal;

    /* Firehose recoupling parameters */
    data->firehose_recoupling_mach = parser_get_opt_param_float(
        parameter_file, "KIARAChemistry:firehose_recoupling_mach", 0.5f);

    data->firehose_recoupling_u_factor = parser_get_opt_param_float(
        parameter_file, "KIARAChemistry:firehose_recoupling_u_factor", 0.5f);

    data->firehose_recoupling_fmix = parser_get_opt_param_float(
        parameter_file, "KIARAChemistry:firehose_recoupling_fmix", 0.9f);

    data->firehose_max_velocity = parser_get_opt_param_float(
        parameter_file, "KIARAChemistry:firehose_max_velocity", 2000.f);
    data->firehose_max_velocity *= data->kms_to_internal;

    data->firehose_max_fmix_per_step = parser_get_opt_param_float(
        parameter_file, "KIARAChemistry:firehose_max_fmix_per_step", 0.1f);
  }

  /* Read the total metallicity */
  data->initial_metal_mass_fraction_total = parser_get_opt_param_float(
      parameter_file, "KIARAChemistry:init_abundance_metal", -1);

  if (data->initial_metal_mass_fraction_total != -1) {
    /* Read the individual mass fractions */
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      char buffer[50];
      sprintf(buffer, "KIARAChemistry:init_abundance_%s",
              chemistry_get_element_name((enum chemistry_element)elem));

      data->initial_metal_mass_fraction[elem] =
          parser_get_param_float(parameter_file, buffer);
    }

    /* Let's check that things make sense (broadly) */

    /* H + He + Z should be ~1 */
    float total_frac = data->initial_metal_mass_fraction[chemistry_element_H] +
                       data->initial_metal_mass_fraction[chemistry_element_He] +
                       data->initial_metal_mass_fraction_total;

    if (total_frac < 0.98 || total_frac > 1.02)
      error("The abundances provided seem odd! H + He + Z = %f =/= 1.",
            total_frac);

    /* Sum of metal elements should be <= Z */
    total_frac = 0.f;
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      if (elem != chemistry_element_H && elem != chemistry_element_He) {
        total_frac += data->initial_metal_mass_fraction[elem];
      }
    }

    if (!data->diffusion_flag) {
      if (total_frac > 1.02 * data->initial_metal_mass_fraction_total) {
        error(
            "The abundances provided seem odd! \\sum metal elements (%f) > Z "
            "(%f)",
            total_frac, data->initial_metal_mass_fraction_total);
      }
    } else {
      /* If diffusion is on, need a metallicity floor so reset Z total */
      if (total_frac > 1.02 * data->initial_metal_mass_fraction_total) {
        warning("Resetting total Z to the sum of all available metals.");
        data->initial_metal_mass_fraction_total = total_frac;

        /* H + He + Z should be ~1 */
        float total_frac_check =
            data->initial_metal_mass_fraction[chemistry_element_H] +
            data->initial_metal_mass_fraction[chemistry_element_He] +
            data->initial_metal_mass_fraction_total;

        if (total_frac_check < 0.98 || total_frac_check > 1.02) {
          error(
              "After resetting, the abundances provided seem odd! "
              "H + He + Z = %f =/= 1.",
              total_frac_check);
        }
      }
    }

    /* Sum of all elements should be <= 1 */
    total_frac = 0.f;
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      total_frac += data->initial_metal_mass_fraction[elem];
    }

    if (total_frac > 1.02) {
      error("The abundances provided seem odd! \\sum elements (%f) > 1",
            total_frac);
    }
  }
}

/**
 * @brief Prints the properties of the chemistry model to stdout.
 *
 * @brief The #chemistry_global_data containing information about the current
 * model.
 */
static INLINE void chemistry_print_backend(
    const struct chemistry_global_data *data) {

  if (data->use_firehose_wind_model) {
    if (data->diffusion_flag) {
      message(
          "Chemistry model is 'KIARA' tracking %d elements with the "
          "firehose wind model and metal diffusion.",
          chemistry_element_count);
    } else {
      message(
          "Chemistry model is 'KIARA' tracking %d elements with "
          "the firehose wind model.",
          chemistry_element_count);
    }
  } else {
    if (data->diffusion_flag) {
      message(
          "Chemistry model is 'KIARA' tracking %d elements with "
          " metal diffusion on.",
          chemistry_element_count);
    } else {
      message("Chemistry model is 'KIARA' tracking %d elements.",
              chemistry_element_count);
    }
  }
}

/**
 * @brief Check recoupling criterion for firehose stream particle .
 * Returns negative value if it should recouple.
 * Actual recoupling is done in feedback.h.
 *
 * @param pi Wind particle (not updated).
 * @param Mach Stream Mach number vs ambient
 * @param r_stream Current radius of stream
 * @param cd #chemistry_global_data containing chemistry information.
 *
 */
__attribute__((always_inline)) INLINE static float
firehose_recoupling_criterion(struct part *p, const float Mach,
                              const float r_stream,
                              const struct chemistry_global_data *cd) {

  if (!cd->use_firehose_wind_model) return 0.f;

  float rs = r_stream;
  const double u = hydro_get_drifted_comoving_internal_energy(p);
  const double u_max = max(u, p->chemistry_data.u_ambient);
  const double u_diff = fabs(u - p->chemistry_data.u_ambient) / u_max;
  if (Mach < cd->firehose_recoupling_mach &&
      u_diff < cd->firehose_recoupling_u_factor)
    rs = -1.f;

  const float exchanged_mass_frac =
      p->chemistry_data.exchanged_mass / hydro_get_mass(p);

  if (exchanged_mass_frac > cd->firehose_recoupling_fmix) rs = -1.f;
  if (r_stream == 0.f) rs = -1.f;

  return rs;
}

/**
 * @brief Finishes the gradient calculation.
 *
 * Nothing to do here.
 *
 * @param p The particle to act upon.
 * @param cd The global properties of the chemistry scheme.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_gradient(
    struct part *p, const struct chemistry_global_data *cd) {}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * Nothing to do here.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The current cosmological model.
 * @param dt_alpha The time-step used to evolve non-cosmological quantities such
 *                 as the artificial viscosity.
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 * @param cd The global properties of the chemistry scheme.
 */
__attribute__((always_inline)) INLINE static void chemistry_prepare_force(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const float dt_alpha, const float dt_therm,
    const struct chemistry_global_data *cd) {}

/**
 * @brief Updates to the chemistry data after the hydro force loop.
 *
 * Finish off the diffusion by actually exchanging the metals
 *
 * @param p The particle to act upon.
 * @param cosmo The current cosmological model.
 * @param with_cosmology Are we running with the cosmology?
 * @param time Current time of the simulation.
 * @param dt Time step (in physical units).
 */
__attribute__((always_inline)) INLINE static void chemistry_end_force(
    struct part *restrict p, struct xpart *xp, const struct cosmology *cosmo,
    const int with_cosmology, const double time, const double dt,
    const struct chemistry_global_data *cd) {

  if (dt == 0.) return;

  struct chemistry_part_data *ch = &p->chemistry_data;

  const float h_inv = 1.f / p->h;
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  /* Missing factors in iact. */
  const float factor = h_inv_dim * h_inv;

  if (cd->use_firehose_wind_model && ch->dm > 0.f) {
    struct cooling_part_data *co = &p->cooling_data;

    const float m = hydro_get_mass(p);
    const float v =
          sqrtf(p->v[0] * p->v[0] + p->v[1] * p->v[1] +
                p->v[2] * p->v[2]);
    const float dv = sqrtf(ch->dv[0] * ch->dv[0] + ch->dv[1] * ch->dv[1] +
                           ch->dv[2] * ch->dv[2]);
    float dv_phys = dv * cosmo->a_inv;

    /* Use this to limit energy change in v^2 and u */
    float alpha = 1.f;

    if (dv >= FIREHOSE_EPSILON_TOLERANCE * v) {
      const float v_new[3] = {p->v[0] + ch->dv[0],
                              p->v[1] + ch->dv[1],
                              p->v[2] + ch->dv[2]};
      const float v_new_norm = sqrtf(v_new[0] * v_new[0] + v_new[1] * v_new[1] +
                                     v_new[2] * v_new[2]);

      /* Apply a kinetic energy limiter */
      const double v2 = v * v;
      double dv2 = dv * dv;
      const double v_new2 = v_new_norm * v_new_norm;
      const double KE_ratio = (v > 0.) ? v_new2 / v2 : 1.;
      const int KE_low_flag = (KE_ratio < FIREHOSE_COOLLIM);
      const int KE_high_flag = (KE_ratio > FIREHOSE_HEATLIM);
      const int KE_out_of_bounds = KE_low_flag || KE_high_flag;

      if (KE_out_of_bounds && dv2 > 0.) {
        /* Solve the same scaling equation, just with a different target */
        const float target_KE_factor =
            (KE_low_flag) ? FIREHOSE_COOLLIM : FIREHOSE_HEATLIM;

        const float v_dot_dv = p->v[0] * ch->dv[0] +
                               p->v[1] * ch->dv[1] +
                               p->v[2] * ch->dv[2];

        /* How to scale all components equally? Solve quadratic:
         * v^2 + 2 * alpha * v * dv + alpha^2 * dv^2 = target_KE_factor * v^2
         *
         * Or equivalently:
         *  A * alpha^2 + B * alpha + C = 0
         *
         * where A = 1
         *       B = 2 * (v * dv) / (dv^2))
         *       C = (v / dv)^2 * (1 - target_KE_factor)
         */
        const float B = 2.f * v_dot_dv / dv2;
        const float C = (v2 / dv2) * (1.f - target_KE_factor);
        const float discriminant = B * B - 4.f * C;
        /* For logging */
        const double u_drift = hydro_get_drifted_comoving_internal_energy(p);

        if (discriminant >= 0.) {
          const float alpha1 = (-B - sqrtf(discriminant)) / 2.f;
          const float alpha2 = (-B + sqrtf(discriminant)) / 2.f;

          /* Minimize alpha1 and alpha2 between (0, 1) */
          if (alpha1 > 0.f && alpha1 < 1.f) alpha = alpha1;
          if (alpha2 < alpha && alpha2 > 0.f) alpha = alpha2;

          /* If there is predicted to be no change, just cancel the operation */
          if (alpha == 1.f) alpha = 0.f;

          ch->dv[0] *= alpha;
          ch->dv[1] *= alpha;
          ch->dv[2] *= alpha;

          message(
              "FIREHOSE_KE_LIMIT p=%lld alpha=%.4g KE_ratio=%.4g v=%.4g "
              "dv=%g m=%g dm=%g u=%g du=%g "
              "dv[0]=%g dv[1]=%g dv[2]=%g "
              "v[0]=%g v[1]=%g v[2]=%g",
              p->id, alpha, KE_ratio, v, dv, m, ch->dm, u_drift, ch->du,
              ch->dv[0], ch->dv[1], ch->dv[2], p->v[0], p->v[1], p->v[2]);
        } else {
          ch->dv[0] = 0.f;
          ch->dv[1] = 0.f;
          ch->dv[2] = 0.f;

          message(
              "FIREHOSE_KE_LIMIT p=%lld alpha=INVALID KE_ratio=%.4g v=%.4g "
              "dv=%g m=%g dm=%g u=%g du=%g "
              "dv[0]=%g dv[1]=%g dv[2]=%g "
              "v[0]=%g v[1]=%g v[2]=%g",
              p->id, KE_ratio, v, dv, m, ch->dm, u_drift, ch->du, 
              ch->dv[0], ch->dv[1], ch->dv[2], p->v[0], p->v[1], p->v[2]);
        }

        /* Recompute the new updated limited values to set v_sig */
        dv2 = ch->dv[0] * ch->dv[0] + ch->dv[1] * ch->dv[1] +
              ch->dv[2] * ch->dv[2];
        dv_phys = sqrtf(dv2) * cosmo->a_inv;
      }
    } else {
      /* Cancel everything if the kick is so small it doesn't matter */
      dv_phys = 0.f;
    }

    /* Make sure there were also no problems with the KE of the particle.
     * Skip all exchanges if there was! */
    if (dv_phys > 0.f) {

      hydro_set_v_sig_based_on_velocity_kick(p, cosmo, dv_phys);

      p->v[0] += ch->dv[0];
      p->v[1] += ch->dv[1];
      p->v[2] += ch->dv[2];

      /* Grab the comoving internal energy at last kick */
      const double u = hydro_get_drifted_comoving_internal_energy(p);

      /* Reset du based on previously calculated alpha limiter */
      ch->du *= alpha * alpha;

      double u_new = u + ch->du;
      const double u_floor =
          cd->firehose_u_floor / cosmo->a_factor_internal_energy;
      if (u_new < u_floor) {
        u_new = u_floor;
        ch->du = u_new - u;
      }

      /* Ignore small changes to the internal energy */
      const double u_eps = fabs(ch->du) / u;

      if (u_eps > FIREHOSE_EPSILON_TOLERANCE) {
#ifdef FIREHOSE_DEBUG_CHECKS
        if (!isfinite(u) || !isfinite(ch->du)) {
          message("FIREHOSE_BAD p=%lld u=%g du=%g dv_phys=%g m=%g dm=%g", p->id,
                  u, ch->du, dv_phys, m, ch->dm);
        }
#endif

        const double energy_frac = (u > 0.) ? u_new / u : 1.;
        if (energy_frac > FIREHOSE_HEATLIM) u_new = FIREHOSE_HEATLIM * u;
        if (energy_frac < FIREHOSE_COOLLIM) u_new = FIREHOSE_COOLLIM * u;

        /* If it's in subgrid ISM mode, use additional heat to
         * lower ISM cold fraction */
        const int firehose_add_heat_to_ISM =
            (p->cooling_data.subgrid_temp > 0.f &&
             p->cooling_data.subgrid_fcold > 0.f && ch->du > 0.);

        if (firehose_add_heat_to_ISM) {

          /* 0.8125 is mu for a fully neutral gas with XH=0.75;
           * approximate but good enough */
          const double T_conv =
              cd->temp_to_u_factor / cosmo->a_factor_internal_energy;
          const double u_cold = 0.8125 * p->cooling_data.subgrid_temp * T_conv;

          const double delta_u = u - u_cold;
          double f_evap = 0.;

          if (delta_u > FIREHOSE_EPSILON_TOLERANCE * u) {
            f_evap = ch->du / delta_u;
            f_evap = min(f_evap, 1.0);
          } else {
            f_evap = 1.0;
          }

          /* Clip values in case of overflow */
          if (f_evap > 0.) {

            p->cooling_data.subgrid_fcold *= 1. - f_evap;

            /* Make sure any extra heat goes into the particle */
            const double u_remaining = ch->du - f_evap * delta_u;
            u_new = u + max(u_remaining, 0.);

            if (p->cooling_data.subgrid_fcold <= 0.f) {
              p->cooling_data.subgrid_temp = 0.f;
              p->cooling_data.subgrid_dens =
                  hydro_get_physical_density(p, cosmo);
              p->cooling_data.subgrid_fcold = 0.f;
            }
          }
        }

        double u_phys = u_new * cosmo->a_factor_internal_energy;

        hydro_set_physical_internal_energy(p, xp, cosmo, u_phys);
        hydro_set_drifted_physical_internal_energy(p, cosmo, NULL, u_phys);
      } else {
        ch->du = 0.;
      }

      /* Check dust change */
      float dust_eps = 0.f;

      /* Check dust change */
      if (co->dust_mass > 0.f) {
        dust_eps = fabs(ch->dm_dust) / co->dust_mass;
      }

      float new_dust_mass = co->dust_mass;
      if (dust_eps >= FIREHOSE_EPSILON_TOLERANCE) new_dust_mass += ch->dm_dust;

      ch->metal_mass_fraction_total = 0.f;
      for (int elem = 0; elem < chemistry_element_count; ++elem) {
        const float old_mass_Z = ch->metal_mass_fraction[elem] * m;
        if (old_mass_Z > 0.f) {
          const float Z_eps = fabs(ch->dm_Z[elem]) / old_mass_Z;

          if (Z_eps >= FIREHOSE_EPSILON_TOLERANCE) {
            ch->metal_mass_fraction[elem] = (old_mass_Z + ch->dm_Z[elem]) / m;
          }
        }

        /* Recompute Z */
        if (elem != chemistry_element_H && elem != chemistry_element_He) {
          ch->metal_mass_fraction_total += ch->metal_mass_fraction[elem];
        }

        if (dust_eps >= FIREHOSE_EPSILON_TOLERANCE) {
          const float old_dust_mass_Z =
              co->dust_mass_fraction[elem] * co->dust_mass;
          co->dust_mass_fraction[elem] =
              (old_dust_mass_Z + ch->dm_dust_Z[elem]) / new_dust_mass;
        }
      }

      /* Set the new dust mass from the exchange */
      co->dust_mass = (new_dust_mass > 0.f) ? new_dust_mass : 0.f;
      if (co->dust_mass <= 0.f) {
        for (int elem = 0; elem < chemistry_element_count; ++elem) {
          co->dust_mass_fraction[elem] = 0.f;
        }

        co->dust_mass = 0.f;
      }

      /* Make sure that X + Y + Z = 1 */
      const float Y_He = ch->metal_mass_fraction[chemistry_element_He];
      ch->metal_mass_fraction[chemistry_element_H] =
          1.f - Y_He - ch->metal_mass_fraction_total;

      /* Make sure H fraction does not go out of bounds */
      if (ch->metal_mass_fraction[chemistry_element_H] > 1.f ||
          ch->metal_mass_fraction[chemistry_element_H] < 0.f) {
        for (int i = chemistry_element_H; i < chemistry_element_count; i++) {
          warning("\telem[%d] is %g", i, ch->metal_mass_fraction[i]);
        }

        error(
            "Hydrogen fraction exeeds unity or is negative for"
            " particle id=%lld due to firehose exchange",
            p->id);
      }

      /* Update stream radius for stream particle */
      if (p->decoupled) {
        const float stream_growth_factor = 1.f + ch->dm / hydro_get_mass(p);
        ch->radius_stream *= sqrtf(stream_growth_factor);

        const double c_s =
            sqrt(ch->u_ambient * hydro_gamma * hydro_gamma_minus_one);
        const float Mach = dv_phys / (cosmo->a_factor_sound_speed * c_s);
        ch->radius_stream =
            firehose_recoupling_criterion(p, Mach, ch->radius_stream, cd);
      }
    }
  }

  /* Are we a decoupled wind? Skip diffusion. */
  if (p->decoupled) return;

  /* Check if we are hypersonic*/
  /* Reset dZ_dt and return? */
  bool reset_time_derivatives = false;

  /* Add diffused metals to particle */
  const float dZ_tot = ch->dZ_dt_total * dt * factor;
  const float new_metal_mass_fraction_total =
      ch->metal_mass_fraction_total + dZ_tot;
  if (ch->metal_mass_fraction_total > 0.f) {
    const float abs_fractional_change =
        fabs(dZ_tot) / ch->metal_mass_fraction_total;
    /* Check if dZ is bigger than 1/4 of the Z */
    if (abs_fractional_change > cd->max_fractional_Z_transfer) {
      reset_time_derivatives = true;
    }
  }

  /* Handle edge case where diffusion leads to negative metallicity */
  if (new_metal_mass_fraction_total < 0.f) {
    warning(
        "Metal diffusion led to negative metallicity!\n"
        "\tpid=%lld\n\tdt=%g\n\tZ=%g\n\tdZ_dt=%g\n"
        "\tdZtot=%g\n\tZnewtot=%g\n\tfactor=%g",
        p->id, dt, ch->metal_mass_fraction_total, ch->dZ_dt_total, dZ_tot,
        new_metal_mass_fraction_total, factor);
    reset_time_derivatives = true;
  }

  /* Handle edge case where diffusion leads to super-unity metallicity */
  if (new_metal_mass_fraction_total > 1.f) {
    warning(
        "Metal diffusion led to metal fractions above unity!\n"
        "pid=%lld\n\tdt=%g\n\tZ=%g\n\tdZ_dt=%g\n"
        "\tdZtot=%g\n\tZnewtot=%g\n\tfactor=%g",
        p->id, dt, ch->metal_mass_fraction_total, ch->dZ_dt_total, dZ_tot,
        new_metal_mass_fraction_total, factor);
    reset_time_derivatives = true;
  }

  /* Add individual element contributions from diffusion */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const float dZ = ch->dZ_dt[elem] * dt * factor;
    const float new_metal_fraction_elem = ch->metal_mass_fraction[elem] + dZ;

    if (ch->metal_mass_fraction[elem] > 0.f) {
      const float abs_fractional_change =
          fabs(dZ) / ch->metal_mass_fraction[elem];
      if (abs_fractional_change > cd->max_fractional_Z_transfer) {
        reset_time_derivatives = true;
      }
    }

    /* Make sure that the metallicity is 0 <= x <= 1 */
    if (new_metal_fraction_elem < 0.f) {
      warning(
          "Z[elem] < 0! pid=%lld, dt=%g, elem=%d, Z=%g, dZ_dt=%g, dZ=%g, "
          "dZtot=%g Ztot=%g Zdust=%g.",
          p->id, dt, elem, ch->metal_mass_fraction[elem], ch->dZ_dt[elem], dZ,
          dZ_tot, ch->metal_mass_fraction_total,
          p->cooling_data.dust_mass_fraction[elem]);
      reset_time_derivatives = true;
    }

    if (new_metal_fraction_elem > 1.f) {
      warning(
          "Z[elem] > 1! pid=%lld, dt=%g, elem=%d, Z=%g, dZ_dt=%g, "
          "dZ=%g, dZtot=%g Ztot=%g.",
          p->id, dt, elem, ch->metal_mass_fraction[elem], ch->dZ_dt[elem], dZ,
          dZ_tot, ch->metal_mass_fraction_total);
      reset_time_derivatives = true;
    }
  }

  /* Found weird dZ_dt values so we should reset everything and exit */
  if (reset_time_derivatives) {
    ch->dZ_dt_total = 0.f;
    for (int elem = 0; elem < chemistry_element_count; elem++) {
      ch->dZ_dt[elem] = 0.f;
    }
    return;
  } else {
#if COOLING_GRACKLE_MODE >= 2
    if (ch->metal_mass_fraction_total > 0.f) {
      /* Add diffused dust to particle, in proportion to added metals */
      p->cooling_data.dust_mass *= 1.f + dZ_tot / ch->metal_mass_fraction_total;
    }
#endif

    /* Reset the total metallicity Z */
    ch->metal_mass_fraction_total = new_metal_mass_fraction_total;

    /* Add individual element contributions from diffusion */
    for (int elem = 0; elem < chemistry_element_count; elem++) {
      const float dZ = ch->dZ_dt[elem] * dt * factor;
      const float new_metal_fraction_elem = ch->metal_mass_fraction[elem] + dZ;

#if COOLING_GRACKLE_MODE >= 2
      /* Add diffused dust to particle, in proportion to added metals */
      if (ch->metal_mass_fraction[elem] > 0.f) {
        p->cooling_data.dust_mass_fraction[elem] *=
            1.f + dZ / ch->metal_mass_fraction[elem];
      }
#endif

      /* Treating Z like a passive scalar */
      ch->metal_mass_fraction[elem] = new_metal_fraction_elem;
    }
  }

  /* Make sure that X + Y + Z = 1 */
  const float Y_He = ch->metal_mass_fraction[chemistry_element_He];
  ch->metal_mass_fraction[chemistry_element_H] =
      1.f - Y_He - ch->metal_mass_fraction_total;

  /* Make sure H fraction does not go out of bounds */
  if (ch->metal_mass_fraction[chemistry_element_H] > 1.f ||
      ch->metal_mass_fraction[chemistry_element_H] < 0.f) {
    for (int i = chemistry_element_H; i < chemistry_element_count; i++) {
      warning("\telem[%d] is %g", i, ch->metal_mass_fraction[i]);
    }

    error(
        "Hydrogen fraction exeeds unity or is negative for"
        " particle id=%lld due to metal diffusion",
        p->id);
  }
}

/**
 * @brief Computes the chemistry-related time-step constraint.
 *
 * Only constraint in KIARA is the diffusion time-step.
 *
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param us The internal system of units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cd The global properties of the chemistry scheme.
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float chemistry_timestep(
    const struct phys_const *restrict phys_const,
    const struct cosmology *restrict cosmo,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props,
    const struct chemistry_global_data *cd, const struct part *restrict p) {

  float dt_chem = FLT_MAX;
  if (cd->diffusion_flag) {
    if (p->chemistry_data.diffusion_coefficient > 0.f) {
      const struct chemistry_part_data *ch = &p->chemistry_data;

      /* Parshikov & Medin 2002 equation 41 */
      const float h_phys = p->h * cosmo->a * kernel_gamma;
      const float D_phys = ch->diffusion_coefficient;
      const float rho_phys = hydro_get_physical_density(p, cosmo);
      dt_chem = cd->diffusion_beta * rho_phys * h_phys * h_phys / D_phys;
      if (dt_chem < cd->time_step_min) {
        message(
            "Warning! dZ_dt timestep low: id=%lld (%g Myr) is below "
            "time_step_min (%g Myr).",
            p->id, dt_chem * cd->time_to_Myr,
            cd->time_step_min * cd->time_to_Myr);
      }

      dt_chem = max(dt_chem, cd->time_step_min);
    }
  }

  if (cd->use_firehose_wind_model) {
    /* About-to-recouple winds need the hydro time-step. */
    if (p->decoupled == 2) {
      const float CFL_condition = hydro_props->CFL_condition;
      const float h = kernel_gamma * cosmo->a * p->h;
      const float v_sig = 2.f * hydro_get_physical_soundspeed(p, cosmo);
      const float dt_cfl = 2.f * CFL_condition * h / v_sig;

      /* The actual minimum time-step is handled in the runner file. */
      dt_chem = min(dt_chem, dt_cfl);
    }
  }

  return dt_chem;
}

/**
 * @brief Initialise the chemistry properties of a black hole with
 * the chemistry properties of the gas it is born from.
 *
 * Black holes don't store fractions so we need to use element masses.
 *
 * @param bp_data The black hole data to initialise.
 * @param p_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_bpart_from_part(
    struct chemistry_bpart_data *bp_data,
    const struct chemistry_part_data *p_data, const double gas_mass) {

  bp_data->metal_mass_total = p_data->metal_mass_fraction_total * gas_mass;
  for (int i = 0; i < chemistry_element_count; ++i) {
    bp_data->metal_mass[i] = p_data->metal_mass_fraction[i] * gas_mass;
  }

  bp_data->formation_metallicity = p_data->metal_mass_fraction_total;
}

/**
 * @brief Add the chemistry data of a gas particle to a black hole.
 *
 * Black holes don't store fractions so we need to add element masses.
 *
 * @param bp_data The black hole data to add to.
 * @param p_data The gas data to use.
 * @param gas_mass The mass of the gas particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_part_to_bpart(
    struct chemistry_bpart_data *bp_data,
    const struct chemistry_part_data *p_data, const double gas_mass) {

  bp_data->metal_mass_total += p_data->metal_mass_fraction_total * gas_mass;
  for (int i = 0; i < chemistry_element_count; ++i) {
    bp_data->metal_mass[i] += p_data->metal_mass_fraction[i] * gas_mass;
  }
}

/**
 * @brief Transfer chemistry data of a gas particle to a black hole.
 *
 * In contrast to `chemistry_add_part_to_bpart`, only a fraction of the
 * masses stored in the gas particle are transferred here. Absolute masses
 * of the gas particle are adjusted as well.
 * Black holes don't store fractions so we need to add element masses.
 *
 * We expect the nibble_mass to be the gas particle mass multiplied by the
 * nibble_fraction.
 *
 * @param bp_data The black hole data to add to.
 * @param p_data The gas data to use.
 * @param nibble_mass The mass to be removed from the gas particle.
 * @param nibble_fraction The fraction of the (original) mass of the gas
 *        particle that is removed.
 */
__attribute__((always_inline)) INLINE static void
chemistry_transfer_part_to_bpart(struct chemistry_bpart_data *bp_data,
                                 struct chemistry_part_data *p_data,
                                 const double nibble_mass,
                                 const double nibble_fraction) {

  bp_data->metal_mass_total += p_data->metal_mass_fraction_total * nibble_mass;
  for (int i = 0; i < chemistry_element_count; ++i)
    bp_data->metal_mass[i] += p_data->metal_mass_fraction[i] * nibble_mass;
}

/**
 * @brief Add the chemistry data of a black hole to another one.
 *
 * @param bp_data The black hole data to add to.
 * @param swallowed_data The black hole data to use.
 */
__attribute__((always_inline)) INLINE static void chemistry_add_bpart_to_bpart(
    struct chemistry_bpart_data *bp_data,
    const struct chemistry_bpart_data *swallowed_data) {

  bp_data->metal_mass_total += swallowed_data->metal_mass_total;
  for (int i = 0; i < chemistry_element_count; ++i) {
    bp_data->metal_mass[i] += swallowed_data->metal_mass[i];
  }
}

/**
 * @brief Split the metal content of a particle into n pieces
 *
 * We only need to split the fields that are not fractions.
 *
 * @param p The #part.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void chemistry_split_part(
    struct part *p, const double n) {}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * We return the un-smoothed quantity here as the star will smooth
 * over its gas neighbours.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_feedback(
    const struct part *restrict p) {

  return p->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * We return the un-smoothed quantity here as the star will smooth
 * over its gas neighbours.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const *
chemistry_get_metal_mass_fraction_for_feedback(const struct part *restrict p) {

  return p->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * star particle to be used in feedback/enrichment related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_fraction_for_feedback(
    const struct spart *restrict sp) {

  return sp->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * star particle to be used in feedback/enrichment related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const *
chemistry_get_star_metal_mass_fraction_for_feedback(
    const struct spart *restrict sp) {

  return sp->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in cooling related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_cooling(
    const struct part *restrict p) {

  return p->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in cooling related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const *
chemistry_get_metal_mass_fraction_for_cooling(const struct part *restrict p) {

  return p->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in star formation related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_star_formation(
    const struct part *restrict p) {

  return p->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in star formation related routines.
 *
 * KIARA uses smooth abundances for everything.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const *
chemistry_get_metal_mass_fraction_for_star_formation(
    const struct part *restrict p) {

  return p->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metal mass of the
 * gas particle to be used in the stats related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_for_stats(const struct part *restrict p) {

  return p->chemistry_data.metal_mass_fraction_total * hydro_get_mass(p);
}

/**
 * @brief Returns the total metal mass of the
 * star particle to be used in the stats related routines.
 *
 * @param sp Pointer to the star particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_for_stats(const struct spart *restrict sp) {

  return sp->chemistry_data.metal_mass_fraction_total * sp->mass;
}

/**
 * @brief Returns the total metal mass of the
 * black hole particle to be used in the stats related routines.
 *
 * @param bp Pointer to the BH particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_bh_total_metal_mass_for_stats(const struct bpart *restrict bp) {

  return bp->chemistry_data.metal_mass_total;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * star particle to be used in the luminosity calculations.
 *
 * @param sp Pointer to the star particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_fraction_for_luminosity(
    const struct spart *restrict sp) {

  return sp->chemistry_data.metal_mass_fraction_total;
}

/**
 * @brief Extra chemistry operations to be done during the drift.
 *
 * @param p Particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param dt_drift The drift time-step for positions.
 * @param dt_therm The drift time-step for thermal quantities.
 * @param cosmo The current cosmological model.
 * @param chem_data The global properties of the chemistry scheme.
 */
__attribute__((always_inline)) INLINE static void chemistry_predict_extra(
    struct part *p, struct xpart *xp, float dt_drift, float dt_therm,
    const struct cosmology *cosmo,
    const struct chemistry_global_data *chem_data) {}

#endif /* SWIFT_CHEMISTRY_KIARA_H */
