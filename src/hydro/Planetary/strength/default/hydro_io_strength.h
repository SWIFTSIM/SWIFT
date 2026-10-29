/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

#ifndef SWIFT_PLANETARY_STRENGTH_IO_DEFAULT_H
#define SWIFT_PLANETARY_STRENGTH_IO_DEFAULT_H

/**
 * @file Planetary/strength/default/hydro_io_strength.h
 * @brief Planetary implementation of SPH with default material strength
 */

#include "adiabatic_index.h"
#include "hydro.h"
#include "hydro_parameters.h"
#include "io_properties.h"
#include "kernel_hydro.h"

/**
 * @brief Specifies which particle strength fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void hydro_read_particles_strength(struct part* parts,
                                                 struct io_props* list,
                                                 int* num_fields) {

  list[7] = io_make_input_field("Density", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_DENSITY, parts, rho);

#ifdef STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG
  list[*num_fields] = io_make_input_field("NumFlaws", INT, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, strength_data.number_of_flaws);
  list[*num_fields+1] = io_make_input_field("FirstFlawIndices", LONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, strength_data.first_flaw_index);

  *num_fields += 2;
#endif
#ifdef STRENGTH_OBJECT_IDS
  list[*num_fields] = io_make_input_field("ObjectIDs", INT, 1, COMPULSORY,
                                          UNIT_CONV_NO_UNITS, parts,
                                          strength_data.object_id);
  *num_fields += 1;
#endif
}

#ifdef STRENGTH_YIELD_STRESS
/**
 * @brief Computes the yield stress of a particle for the snapshot output.
 *
 * @param e The #engine.
 * @param p The particle.
 * @param xp The extended particle.
 * @param ret The yield stress.
 */
INLINE static void convert_yield_stress(const struct engine* e,
                                        const struct part* p,
                                        const struct xpart* xp, float* ret) {

  const struct mat_params *mp = material_params(p->mat_id);
  struct strength_state s;
  hydro_strength_state_init_drift(&s, p, mp);

  ret[0] = strength_state_get_yield_stress(&s, mp);
}
#endif

/**
 * @brief Specifies which particle strength fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static void hydro_write_particles_strength(const struct part* parts,
                                         const struct xpart* xparts,
                                         struct io_props* list,
                                         int* num_fields) {

  list[*num_fields] = io_make_output_field(
      "DeviatoricStressTensors", FLOAT, 6, UNIT_CONV_PRESSURE, -3.f * hydro_gamma,
      parts, strength_data.deviatoric_stress_tensor,
      "Deviatoric stress tensors of the particles (xx, yy, zz, xy, xz, yz)");
  *num_fields += 1;

#ifdef STRENGTH_YIELD_STRESS
  list[*num_fields] = io_make_output_field_convert_part(
      "YieldStresses", FLOAT, 1, UNIT_CONV_PRESSURE, -3.f * hydro_gamma, parts,
      xparts, convert_yield_stress, "Yield stresses of the particles");
  *num_fields += 1;
#endif
#ifdef STRENGTH_DAMAGE
  list[*num_fields] = io_make_output_field(
      "Damage", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts, strength_data.damage.total, "Damage of the particles");
  *num_fields += 1;
#endif
/* The damage from each source, only needed if both contribute to the total. */
#if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG) && defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
  list[*num_fields] = io_make_output_field(
      "TensileDamage", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      strength_data.damage.tensile, "Damage of the particles due to tension");
  list[*num_fields + 1] = io_make_output_field(
      "ShearDamage", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      strength_data.damage.shear, "Damage of the particles due to shear");
  *num_fields += 2;
#endif
#ifdef STRENGTH_OBJECT_IDS
  list[*num_fields] = io_make_output_field(
      "ObjectIDs", INT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      strength_data.object_id, "IDs of the objects that the particles belong to");
  *num_fields += 1;
#endif
}

#endif /* SWIFT_PLANETARY_STRENGTH_IO_DEFAULT_H */
