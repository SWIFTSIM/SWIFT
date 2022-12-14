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
#ifndef SWIFT_SUPERNOVAE_II_GEAR_H
#define SWIFT_SUPERNOVAE_II_GEAR_H

#include "hdf5_functions.h"
#include "interpolation.h"
#include "stellar_evolution_struct.h"

/**
 * @brief Print the supernovae II model.
 *
 * @param snii The #supernovae_ii.
 */
void supernovae_ii_print(const struct supernovae_ii *snii);
int supernovae_ii_can_explode(const struct supernovae_ii *snii, float m_low,
                              float m_high);
float supernovae_ii_get_number_per_unit_mass(const struct supernovae_ii *snii,
                                             float m1, float m2);
void supernovae_ii_get_yields_from_integral(const struct supernovae_ii *snii,
                                            float log_m1, float log_m2,
                                            float *yields);
void supernovae_ii_get_yields_from_raw(const struct supernovae_ii *snii,
                                       float log_m, float *yields);
float supernovae_ii_get_ejected_mass_fraction_non_processed_from_integral(
    const struct supernovae_ii *snii, float log_m1, float log_m2);
float supernovae_ii_get_ejected_mass_fraction_non_processed_from_raw(
    const struct supernovae_ii *snii, float log_m);

float supernovae_ii_get_ejected_mass_fraction_processed_from_integral(
    const struct supernovae_ii *snii, float log_m1, float log_m2);
float supernovae_ii_get_ejected_mass_fraction_processed_from_raw(
    const struct supernovae_ii *snii, float log_m);
float supernovae_ii_get_energy_from_progenitor_mass(
    const struct supernovae_ii *snii, float mass);
void supernovae_ii_read_yields_array(
    struct supernovae_ii *snii, struct interpolation_1d *interp_raw,
    struct interpolation_1d *interp_int, const struct stellar_model *sm,
    hid_t group_id, const char *hdf5_dataset_name, hsize_t *previous_count,
    int interpolation_size);
void supernovae_ii_read_yields(struct supernovae_ii *snii,
                               struct swift_params *params,
                               const struct stellar_model *sm,
                               const int restart);
void supernovae_ii_read_from_params(struct supernovae_ii *snii,
                                    struct swift_params *params);

void supernovae_ii_read_from_tables(struct supernovae_ii *snii,
                                    struct swift_params *params,
                                    const char *filename);
void supernovae_ii_init(struct supernovae_ii *snii, struct swift_params *params,
                        const struct stellar_model *sm,
                        const struct unit_system *us);
void supernovae_ii_dump(const struct supernovae_ii *snii, FILE *stream,
                        const struct stellar_model *sm);
void supernovae_ii_restore(struct supernovae_ii *snii, FILE *stream,
                           const struct stellar_model *sm);
void supernovae_ii_clean(struct supernovae_ii *snii);
#endif  // SWIFT_SUPERNOVAE_II_GEAR_H
