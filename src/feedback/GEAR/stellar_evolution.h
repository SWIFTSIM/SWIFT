/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_STELLAR_EVOLUTION_GEAR_H
#define SWIFT_STELLAR_EVOLUTION_GEAR_H

#include "hdf5_functions.h"
#include "initial_mass_function.h"
#include "lifetime.h"
#include "random.h"
#include "stellar_evolution_struct.h"
#include "supernovae_ia.h"
#include "supernovae_ii.h"

#include <math.h>
#include <stddef.h>

void stellar_model_print(const struct stellar_model* sm);
int stellar_evolution_compute_integer_number_supernovae(
    struct spart* restrict sp, float number_supernovae_f,
    const integertime_t ti_begin, enum random_number_type random_type);

void stellar_evolution_compute_continuous_feedback_properties(
    struct spart* restrict sp, const struct stellar_model* sm,
    const struct phys_const* phys_const, const float log_m_beg_step,
    const float log_m_end_step, const float m_beg_step, const float m_end_step,
    const float m_init, const float number_snia_f, const float number_snii_f);
void stellar_evolution_compute_discrete_feedback_properties(
    struct spart* restrict sp, const struct stellar_model* sm,
    const struct phys_const* phys_const, const float log_m_beg_step,
    const float m_end_step, const float m_init, const int number_snia,
    const int number_snii);

void stellar_evolution_evolve_spart(
    struct spart* restrict sp, const struct stellar_model* sm,
    const struct cosmology* cosmo, const struct unit_system* us,
    const struct phys_const* phys_const, const integertime_t ti_begin,
    const double star_age_beg_step, const double dt);

const char* stellar_evolution_get_element_name(const struct stellar_model* sm,
                                               int i);
int stellar_evolution_get_element_index(const struct stellar_model* sm,
                                        const char* element_name);

void stellar_evolution_read_elements(struct stellar_model* sm,
                                     struct swift_params* params);
void stellar_evolution_props_init(struct stellar_model* sm,
                                  const struct phys_const* phys_const,
                                  const struct unit_system* us,
                                  struct swift_params* params,
                                  const struct cosmology* cosmo);

void stellar_evolution_dump(const struct stellar_model* sm, FILE* stream);
void stellar_evolution_restore(struct stellar_model* sm, FILE* stream);

void stellar_evolution_clean(struct stellar_model* sm);

#endif  // SWIFT_STELLAR_EVOLUTION_GEAR_H
