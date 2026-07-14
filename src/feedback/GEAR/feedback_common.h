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

/* We need to explicitely point to the src/ file to ensure the correct file is
   included for each feedback */
#include "../../feedback_properties.h"
#include "cooling.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

/**
 * @file src/feebback/GEAR/feedback_common.h
 * @brief Header file with common functions for GEAR and GEAR-mechanical
 * feedback modules.
 */

float feedback_compute_spart_timestep(
    const struct spart *const sp, const struct feedback_props *feedback_props,
    const struct phys_const *phys_const, const struct unit_system *us,
    const int with_cosmology, const struct cosmology *cosmo,
    const integertime_t ti_current, const double time, const double time_base);

void feedback_will_do_feedback(
    struct spart *sp, const struct feedback_props *feedback_props,
    const int with_cosmology, const struct cosmology *cosmo, const double time,
    const struct unit_system *us, const struct phys_const *phys_const,
    const integertime_t ti_current, const double time_base);

void compute_time(const struct spart *sp, const int with_cosmology,
                  const struct cosmology *cosmo, double *star_age_beg_of_step,
                  double *dt_enrichment, integertime_t *ti_begin_star,
                  const integertime_t ti_current, const double time_base,
                  const double time);

double compute_star_age_end_of_step(const struct spart *sp,
                                    const int with_cosmology,
                                    const struct cosmology *cosmo,
                                    const double time);

double feedback_get_enrichment_timestep(const struct spart *sp,
                                        const int with_cosmology,
                                        const struct cosmology *cosmo,
                                        const double time,
                                        const double dt_star);

int feedback_is_HII_ionization_active(const struct spart *sp,
                                      const struct engine *e);
double feedback_get_star_ionization_rate(const struct spart *sp, int pixel);
double feedback_get_star_ionization_rate_max(const struct spart *sp);
char feedback_part_can_be_ionized(const struct part *p, const struct xpart *xp,
                                  const struct engine *e);
void feedback_iact_HII_ionization(
    struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xpj, float r2, int pixel,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const integertime_t ti_begin);

void feedback_init_after_star_formation(
    struct spart *sp, const struct feedback_props *feedback_props,
    enum stellar_type star_type);

void feedback_first_init_spart(struct spart *sp,
                               const struct feedback_props *feedback_props);

void feedback_struct_dump(const struct feedback_props *feedback, FILE *stream);
void feedback_struct_restore(struct feedback_props *feedback, FILE *stream);
void feedback_clean(struct feedback_props *feedback);

#endif /* SWIFT_FEEDBACK_GEAR_COMMON_H */
