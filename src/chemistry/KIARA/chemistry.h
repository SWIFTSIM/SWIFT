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
#include "units.h"

/**
 * @brief Return a string containing the name of a given #chemistry_element.
 */
__attribute__((always_inline)) INLINE static const char*
chemistry_get_element_name(enum chemistry_element elem) {

  static const char* chemistry_element_names[chemistry_element_count] = {
      "Hydrogen", "Helium",    "Carbon",  "Nitrogen", "Oxygen",
      "Neon",     "Magnesium", "Silicon", "Sulfur", "Calcium", "Iron"};

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
    struct part* restrict p, const struct chemistry_global_data* cd) {

  struct chemistry_part_data* cpd = &p->chemistry_data;

  /* Reset the shear tensor */
  for (int i = 0; i < 3; i++) {
    cpd->shear_tensor[i][0] = 0.f;
    cpd->shear_tensor[i][1] = 0.f;
    cpd->shear_tensor[i][2] = 0.f;
  }

  /* Reset the diffusion. */
  cpd->diffusion_coefficient = 0.f;

#if COOLING_GRACKLE_MODE >= 2
  cpd->local_sfr_density = 0.f;
#endif
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
    struct part* restrict p, const struct chemistry_global_data* cd,
    const struct cosmology* cosmo) {


  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h; /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */

  struct chemistry_part_data* cpd = &p->chemistry_data;

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

    float shear_tensor[3][3];
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
    }
    velocity_gradient_norm = sqrtf(velocity_gradient_norm);

    /* Compute the diffusion coefficient in physical coordinates.
     * The norm is already in physical coordinates.
     * kernel_gamma is necessary (see Rennehan 2021)
     */
    const float rho_phys = hydro_get_physical_density(p, cosmo);
    const float h_phys = cosmo->a * p->h * kernel_gamma;
    const float smag_length_scale = cd->C_Smagorinsky * h_phys;
    cpd->diffusion_coefficient = 2.f * rho_phys * smag_length_scale 
                                * smag_length_scale * velocity_gradient_norm;
  }
#if COOLING_GRACKLE_MODE >= 2
  /* Add self contribution to SFR */
  cpd->local_sfr_density += max(0.f, p->sf_data.SFR);
  const float vol_factor = 0.238732 * h_inv_dim; /* 1./(4/3 pi) */
  /* Convert to physical density from comoving */
  cpd->local_sfr_density *= vol_factor * cosmo->a3_inv;
#endif
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
chemistry_part_has_no_neighbours(struct part* restrict p,
                                 struct xpart* restrict xp,
                                 const struct chemistry_global_data* cd,
                                 const struct cosmology* cosmo) {

  /* Just make all the smoothed fields default to the un-smoothed values */
  struct chemistry_part_data* cpd = &p->chemistry_data;
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
  for (int elem = 0; elem < chemistry_element_count; ++elem) cpd->dZ_dt[elem] = 0.f;

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
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct chemistry_global_data* data, struct part* restrict p,
    struct xpart* restrict xp) {

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
}

/**
 * @brief Sets the chemistry properties of the sparticles to a valid start
 * state.
 *
 * @param data The global chemistry information.
 * @param sp Pointer to the sparticle data.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_spart(
    const struct chemistry_global_data* data, struct spart* restrict sp) {

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
 * @brief Initialises the chemistry properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param data The properties to initialise.
 */
static INLINE void chemistry_init_backend(struct swift_params* parameter_file,
                                          const struct unit_system* us,
                                          const struct phys_const* phys_const,
                                          struct chemistry_global_data* data) {

  /* Is metal diffusion turned on? */
  data->diffusion_flag = parser_get_param_int(parameter_file,
                                              "KIARAChemistry:diffusion_on");

  /* Read the diffusion coefficient */
  data->C_Smagorinsky = parser_get_opt_param_float(parameter_file,
                                                   "KIARAChemistry:diffusion_coefficient",
                                                   0.23);

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

    if (total_frac > 1.02 * data->initial_metal_mass_fraction_total)
      error(
          "The abundances provided seem odd! \\sum metal elements (%f) > Z "
          "(%f)",
          total_frac, data->initial_metal_mass_fraction_total);

    /* Sum of all elements should be <= 1 */
    total_frac = 0.f;
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      total_frac += data->initial_metal_mass_fraction[elem];
    }

    if (total_frac > 1.02)
      error("The abundances provided seem odd! \\sum elements (%f) > 1",
            total_frac);
  }
}

/**
 * @brief Prints the properties of the chemistry model to stdout.
 *
 * @brief The #chemistry_global_data containing information about the current
 * model.
 */
static INLINE void chemistry_print_backend(
    const struct chemistry_global_data* data) {

  message("Chemistry model is 'KIARA' tracking %d elements.",
          chemistry_element_count);
}

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
    struct part* restrict p, const struct cosmology* cosmo,
    const int with_cosmology, const double time, const double dt) {

  if (dt == 0.) return;

  struct chemistry_part_data* ch = &p->chemistry_data;

  /* Nothing to do if no change in metallicity */
  if (ch->dZ_dt_total == 0.f) return;

  const float h_inv = cosmo->a / p->h;
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  /* Missing factors in iact. */
  const float factor = h_inv_dim * h_inv;
  //const float mass = hydro_get_mass(p);

  /* Add diffused metals to particle */
  const float dZ_tot = ch->dZ_dt_total * dt * factor;
#if COOLING_GRACKLE_MODE >= 2
  /* Add diffused dust to particle, in proportion to added metals */
  p->cooling_data.dust_mass *= 1. + dZ_tot / ch->metal_mass_fraction_total;
#endif
  ch->metal_mass_fraction_total += dZ_tot;

  /* Handle edge case where diffusion leads to <0 metallicity */
  if (ch->metal_mass_fraction_total <= 0.f) {
    ch->metal_mass_fraction_total = 0.f;
    p->cooling_data.dust_mass = 0.f;
    for (int elem = 0; elem < chemistry_element_count; elem++) {
      ch->metal_mass_fraction[elem] = 0.f;
      p->cooling_data.dust_mass_fraction[elem] = 0.f;
    }
    return;
  }

  /* Add individual element contributions from diffusion */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const float dZ = ch->dZ_dt[elem] * dt * factor;
#if COOLING_GRACKLE_MODE >= 2
  /* Add diffused dust to particle, in proportion to added metals */
    p->cooling_data.dust_mass_fraction[elem] *= 1. + dZ / ch->metal_mass_fraction[elem];
#endif
    /* Treating Z like a passive scalar */
    ch->metal_mass_fraction[elem] += dZ;

    /* Make sure that the metallicity is 0 <= x <= 1 */
    if (ch->metal_mass_fraction[elem] < 0.f ) {
      warning("Z<0! pid=%lld, dt=%g, elem=%d, Z=%g, dZ_dt=%g, dZ=%g, dZtot=%g Ztot=%g Zdust=%g.", 
            p->id, dt, elem, ch->metal_mass_fraction[elem], ch->dZ_dt[elem], dZ,
	    dZ_tot, ch->metal_mass_fraction_total, p->cooling_data.dust_mass_fraction[elem]);
      ch->metal_mass_fraction[elem] = 0.f;
      p->cooling_data.dust_mass_fraction[elem] = 0.f;
    }
    if (ch->metal_mass_fraction[elem] > 1.f) {
      error("met frac>1! pid=%lld, dt=%g, elem=%d, Z=%g, dZ_dt=%g, dZ=%g, dZtot=%g Ztot=%g.", 
            p->id, dt, elem, ch->metal_mass_fraction[elem], ch->dZ_dt[elem], dZ,
	    dZ_tot, ch->metal_mass_fraction_total);
    }
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
    const struct phys_const* restrict phys_const,
    const struct cosmology* restrict cosmo,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_props,
    const struct chemistry_global_data* cd, const struct part* restrict p) {

  if (cd->diffusion_flag) {
    const struct chemistry_part_data* ch = &p->chemistry_data;
    float max_dZ_dt = FLT_MIN;
    for (int elem = 0; elem < chemistry_element_count; elem++) {
      const float abs_dZ_dt_elem = fabs(ch->dZ_dt[elem]);
      if (abs_dZ_dt_elem > max_dZ_dt) {
        max_dZ_dt = abs_dZ_dt_elem;
      }
    }

    if (max_dZ_dt > FLT_MIN) {
      return 1.f / max_dZ_dt;
    }
  }

  return FLT_MAX;
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
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_part_data* p_data, const double gas_mass) {

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
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_part_data* p_data, const double gas_mass) {

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
chemistry_transfer_part_to_bpart(struct chemistry_bpart_data* bp_data,
                                 struct chemistry_part_data* p_data,
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
    struct chemistry_bpart_data* bp_data,
    const struct chemistry_bpart_data* swallowed_data) {

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
    struct part* p, const double n) { }

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
    const struct part* restrict p) {

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
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_feedback(const struct part* restrict p) {

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
    const struct spart* restrict sp) {

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
__attribute__((always_inline)) INLINE static float const*
chemistry_get_star_metal_mass_fraction_for_feedback(
    const struct spart* restrict sp) {

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
    const struct part* restrict p) {

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
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_cooling(const struct part* restrict p) {

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
    const struct part* restrict p) {

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
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_star_formation(
    const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metal mass of the
 * gas particle to be used in the stats related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_for_stats(const struct part* restrict p) {

  return p->chemistry_data.metal_mass_fraction_total * hydro_get_mass(p);
}

/**
 * @brief Returns the total metal mass of the
 * star particle to be used in the stats related routines.
 *
 * @param sp Pointer to the star particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_for_stats(const struct spart* restrict sp) {

  return sp->chemistry_data.metal_mass_fraction_total * sp->mass;
}

/**
 * @brief Returns the total metal mass of the
 * black hole particle to be used in the stats related routines.
 *
 * @param bp Pointer to the BH particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_bh_total_metal_mass_for_stats(const struct bpart* restrict bp) {

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
    const struct spart* restrict sp) {

  return sp->chemistry_data.metal_mass_fraction_total;
}

#endif /* SWIFT_CHEMISTRY_KIARA_H */
