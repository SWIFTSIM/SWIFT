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
#ifndef SWIFT_FEEDBACK_GEAR_H
#define SWIFT_FEEDBACK_GEAR_H

#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "stellar_evolution.h"
#include "units.h"

#include <strings.h>

void feedback_update_part(struct part* restrict p, struct xpart* restrict xp,
                          const struct engine* restrict e);

int feedback_will_do_feedback(const struct spart* sp,
                              const struct feedback_props* feedback_props,
                              const int with_cosmology,
                              const struct cosmology* cosmo, const double time);

int feedback_is_active(const struct spart* sp, const double time,
                       const struct cosmology* cosmo, const int with_cosmology);
double feedback_get_enrichment_timestep(const struct spart* sp,
                                        const int with_cosmology,
                                        const struct cosmology* cosmo,
                                        const double time,
                                        const double dt_star);
void feedback_init_spart(struct spart* sp);

void feedback_reset_feedback(struct spart* sp,
                             const struct feedback_props* feedback_props);
void feedback_first_init_spart(struct spart* sp,
                               const struct feedback_props* feedback_props);
void feedback_prepare_spart(struct spart* sp,
                            const struct feedback_props* feedback_props);
void feedback_evolve_spart(struct spart* restrict sp,
                           const struct feedback_props* feedback_props,
                           const struct cosmology* cosmo,
                           const struct unit_system* us,
                           const struct phys_const* phys_const,
                           const double star_age_beg_step, const double dt,
                           const double time, const integertime_t ti_begin,
                           const int with_cosmology);
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream);
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream,
                             struct engine* e);
void feedback_clean(struct feedback_props* feedback);

/**
 * @brief Writes the current model of feedback to the file
 *
 * @param feedback The #feedback_props.
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void feedback_write_flavour(struct feedback_props* feedback,
                                          hid_t h_grp) {

  io_write_attribute_s(h_grp, "Feedback Model", "GEAR");
};

#endif /* SWIFT_FEEDBACK_GEAR_H */
