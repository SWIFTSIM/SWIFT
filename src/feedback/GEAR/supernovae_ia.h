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
#ifndef SWIFT_SUPERNOVAE_IA_GEAR_H
#define SWIFT_SUPERNOVAE_IA_GEAR_H

#include "hdf5_functions.h"
#include "stellar_evolution.h"
#include "stellar_evolution_struct.h"

void supernovae_ia_print(const struct supernovae_ia *snia);
int supernovae_ia_can_explode(const struct supernovae_ia *snia, float m_low,
                              float m_high);
const float *supernovae_ia_get_yields(const struct supernovae_ia *snia);
float supernovae_ia_get_ejected_mass_processed(
    const struct supernovae_ia *snia);
float supernovae_ia_get_companion_fraction(const struct supernovae_ia *snia,
                                           float m1, float m2,
                                           int companion_type);
float supernovae_ia_get_number_per_unit_mass(const struct supernovae_ia *snia,
                                             float m1, float m2);

void supernovae_ia_read_yields(struct supernovae_ia *snia,
                               struct swift_params *params,
                               const struct stellar_model *sm,
                               const char *filename);

void supernovae_ia_init_companion(struct supernovae_ia *snia);

void supernovae_ia_read_from_tables(struct supernovae_ia *snia,
                                    struct swift_params *params,
                                    const char *filename);

void supernovae_ia_init(struct supernovae_ia *snia,
                        const struct phys_const *phys_const,
                        const struct unit_system *us,
                        struct swift_params *params,
                        const struct stellar_model *sm);

void supernovae_ia_dump(const struct supernovae_ia *snia, FILE *stream,
                        const struct stellar_model *sm);

void supernovae_ia_restore(struct supernovae_ia *snia, FILE *stream,
                           const struct stellar_model *sm);
void supernovae_ia_clean(struct supernovae_ia *snia);

#endif  // SWIFT_SUPERNOVAE_IA_GEAR_H
