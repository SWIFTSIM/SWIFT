/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_FEEDBACK_GEAR_MECHANICAL_H
#define SWIFT_FEEDBACK_GEAR_MECHANICAL_H

#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "stellar_evolution.h"
#include "units.h"

#include <strings.h>

void feedback_update_part(struct part* p, struct xpart* xp,
                          const struct engine* e);

void feedback_reset_part(struct part* p, struct xpart* xp);

void feedback_will_do_feedback(
    struct spart* sp, const struct feedback_props* feedback_props,
    const int with_cosmology, const struct cosmology* cosmo, const double time,
    const struct unit_system* us, const struct phys_const* phys_const,
    const integertime_t ti_current, const double time_base);

int feedback_is_active(const struct spart* sp, const struct engine* e);

double feedback_get_enrichment_timestep(const struct spart* sp,
                                        const int with_cosmology,
                                        const struct cosmology* cosmo,
                                        const double time,
                                        const double dt_star);
void feedback_init_spart(struct spart* sp);

void feedback_init_after_star_formation(
    struct spart* sp, const struct feedback_props* feedback_props);
void feedback_reset_feedback(struct spart* sp,
                             const struct feedback_props* feedback_props);
void feedback_first_init_spart(struct spart* sp,
                               const struct feedback_props* feedback_props);
void feedback_prepare_spart(struct spart* sp,
                            const struct feedback_props* feedback_props);
void feedback_prepare_feedback(struct spart* restrict sp,
                               const struct feedback_props* feedback_props,
                               const struct cosmology* cosmo,
                               const struct unit_system* us,
                               const struct phys_const* phys_const,
                               const double star_age_beg_step, const double dt,
                               const double time, const integertime_t ti_begin,
                               const int with_cosmology);
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream);
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream);
void feedback_clean(struct feedback_props* feedback);

/**
 * @brief Writes the current model of feedback to the file
 *
 * @param feedback The #feedback_props.
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void feedback_write_flavour(struct feedback_props* feedback,
                                          hid_t h_grp) {

  io_write_attribute_s(h_grp, "Feedback Model", "GEAR-mechanical");
};

double feedback_compute_scalar_weight(const float r2, const float* dx,
                                      const float hi, const float hj,
                                      const struct spart* restrict si,
                                      const struct part* restrict pj,
                                      double* dx_ij_plus, double* dx_ij_minus);

void feedback_compute_vector_weight_non_normalized(
    const float r2, const float* dx, const float hi, const float hj,
    const struct spart* restrict si, const struct part* restrict pj,
    double* f_plus_i, double* f_minus_i, double* w_j);

void feedback_compute_vector_weight_normalized(const float r2, const float *dx,
					       const float hi, const float hj,
					       const struct spart *restrict si,
					       const struct part *restrict pj,
					       double* w_j_bar);

double feedback_get_SN_terminal_momentum(const struct spart* restrict sp,
					 const struct part* restrict p,
					 const struct xpart* restrict xp,
					 const struct phys_const* phys_const,
					 const struct unit_system* us);

#endif /* SWIFT_FEEDBACK_GEAR_MECHANICAL_H */
