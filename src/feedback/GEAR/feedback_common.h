/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_FEEDBACK_GEAR_COMMON_H
#define SWIFT_FEEDBACK_GEAR_COMMON_H

/* We need to explicitly point to the src/ file to ensure the correct file is
   included for each feedback */
#include "../../feedback_properties.h"
#include "cooling.h"
#include "hydro_properties.h"
#include "part.h"
#include "radiation_selection.h"
#include "units.h"

/**
 * @file src/feedback/GEAR/feedback_common.h
 * @brief Header file with common functions for GEAR and GEAR-mechanical
 * feedback modules.
 */

void feedback_compute_spart_timestep(
    const struct spart *const sp, const struct feedback_props *feedback_props,
    const struct phys_const *phys_const, const struct unit_system *us,
    const int with_cosmology, const struct cosmology *cosmo,
    const integertime_t ti_current, const double time, const double time_base,
    const timebin_t old_time_bin, float *dt_event_side,
    float *dt_evolution_ssp);

void feedback_will_do_feedback(
    struct spart *sp, const struct feedback_props *feedback_props,
    const int with_cosmology, const struct cosmology *cosmo, const double time,
    const struct unit_system *us, const struct phys_const *phys_const,
    const integertime_t ti_current, const double time_base,
    const timebin_t old_time_bin);

void compute_time(const struct spart *sp, const int with_cosmology,
                  const struct cosmology *cosmo, double *star_age_beg_of_step,
                  double *dt_enrichment, integertime_t *ti_begin_star,
                  const integertime_t ti_current, const double time_base,
                  const double time, const timebin_t old_time_bin);

double compute_star_age_end_of_step(const struct spart *sp,
                                    const int with_cosmology,
                                    const struct cosmology *cosmo,
                                    const double time);

double feedback_get_enrichment_timestep(const struct spart *sp,
                                        const int with_cosmology,
                                        const struct cosmology *cosmo,
                                        const double time,
                                        const double dt_star);

#ifdef GEAR_SUBGRID_RADIATION_HII
void feedback_will_do_HII_ionization(
    struct spart *sp, const struct feedback_props *feedback_props,
    const double star_age_beg_step, const double star_age_end_step);
int feedback_is_HII_ionization_active(const struct spart *sp,
                                      const struct engine *e);
double feedback_get_star_ionization_rate(const struct spart *sp, int pixel);
double feedback_get_star_ionization_budget(const struct spart *sp, int pixel);
double feedback_get_star_ionization_budget_max(const struct spart *sp);
double feedback_get_star_ionization_budget_total(const struct spart *sp);
int feedback_get_star_HII_pixel_count(const struct spart *sp);
double feedback_get_star_HII_last_rebuild(const struct spart *sp);
double feedback_get_star_HII_nominal_interval(
    const struct feedback_props *feedback_props, const double dt_enrichment);
void feedback_open_star_ionizing_photon_budget(struct spart *sp,
                                               double dt_back);
void feedback_resync_star_ionizing_photon_rate_cache(struct spart *sp);
void feedback_set_star_HII_last_rebuild(struct spart *sp,
                                        double star_age_beg_step);
double feedback_get_star_HII_last_attempt(const struct spart *sp);
void feedback_set_star_HII_last_attempt(struct spart *sp,
                                        double star_age_beg_step);

float feedback_get_star_HII_mass(const struct spart *sp);
#else
/* Without the part: the star carries no HII state. */
__attribute__((always_inline)) INLINE static void
feedback_will_do_HII_ionization(struct spart *sp,
                                const struct feedback_props *feedback_props,
                                const double star_age_beg_step,
                                const double star_age_end_step) {}
__attribute__((always_inline)) INLINE static double
feedback_get_star_HII_last_rebuild(const struct spart *sp) {
  return 0.;
}
__attribute__((always_inline)) INLINE static float feedback_get_star_HII_mass(
    const struct spart *sp) {
  return 0.f;
}
#endif /* GEAR_SUBGRID_RADIATION_HII */
double feedback_get_star_L_PE(const struct spart *sp);
double feedback_get_star_L_LW(const struct spart *sp);
float feedback_get_star_teff(const struct spart *sp);

void feedback_init_after_star_formation(
    struct spart *sp, const struct feedback_props *feedback_props,
    enum stellar_type star_type);

void feedback_first_init_spart(struct spart *sp,
                               const struct feedback_props *feedback_props);

float feedback_get_comoving_gas_density_at_star(const struct spart *sp);

/*! ISRF layout marker written ahead of #feedback_props: 1 reduced flux, 2
 * pending fields at a = 0 only, 3 phi-weighted pending at any a. */
#define FEEDBACK_RESTART_ISRF_PART_LAYOUT 3

/*! The restart layout marker adds this factor times
    #radiation_selection_absent_mask(), so that a build with all the parts
    writes #FEEDBACK_RESTART_ISRF_PART_LAYOUT unchanged. */
#define FEEDBACK_RESTART_SELECTION_FACTOR 16

/**
 * @brief The particle layout marker this build writes to its restart files.
 */
__attribute__((always_inline)) INLINE static int
feedback_restart_particle_layout(void) {
  return FEEDBACK_RESTART_ISRF_PART_LAYOUT +
         FEEDBACK_RESTART_SELECTION_FACTOR * radiation_selection_absent_mask();
}

void feedback_struct_dump(const struct feedback_props *feedback, FILE *stream);
void feedback_struct_restore(struct feedback_props *feedback, FILE *stream,
                             const struct unit_system *us,
                             const struct phys_const *phys_const);
void feedback_clean(struct feedback_props *feedback);

#endif /* SWIFT_FEEDBACK_GEAR_COMMON_H */
