/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "potentials.h"

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param potential The external potential properties to initialize
 */
void potential_init(const struct swift_params* parameter_file,
                    const struct phys_const* phys_const,
                    const struct UnitSystem* us,
                    struct external_potential* potential) {

#ifdef EXTERNAL_POTENTIAL_POINTMASS

  potential->point_mass.x =
      parser_get_param_double(parameter_file, "PointMass:position_x");
  potential->point_mass.y =
      parser_get_param_double(parameter_file, "PointMass:position_y");
  potential->point_mass.z =
      parser_get_param_double(parameter_file, "PointMass:position_z");
  potential->point_mass.mass =
      parser_get_param_double(parameter_file, "PointMass:mass");
  potential->point_mass.timestep_mult =
      parser_get_param_float(parameter_file, "PointMass:timestep_mult");

#endif /* EXTERNAL_POTENTIAL_POINTMASS */

#ifdef EXTERNAL_POTENTIAL_ISOTHERMALPOTENTIAL

  potential->isothermal_potential.x =
      parser_get_param_double(parameter_file, "IsothermalPotential:position_x");
  potential->isothermal_potential.y =
      parser_get_param_double(parameter_file, "IsothermalPotential:position_y");
  potential->isothermal_potential.z =
      parser_get_param_double(parameter_file, "IsothermalPotential:position_z");
  potential->isothermal_potential.vrot =
      parser_get_param_double(parameter_file, "IsothermalPotential:vrot");
  potential->isothermal_potential.timestep_mult = parser_get_param_float(
      parameter_file, "IsothermalPotential:timestep_mult");

#endif /* EXTERNAL_POTENTIAL_ISOTHERMALPOTENTIAL */
#ifdef EXTERNAL_POTENTIAL_DISK_PATCH
  potential->disk_patch_potential.surface_density = parser_get_param_double(
      parameter_file, "Disk-PatchPotential:surface_density");
  potential->disk_patch_potential.scale_height = parser_get_param_double(
      parameter_file, "Disk-PatchPotential:scale_height");
  potential->disk_patch_potential.z_disk =
      parser_get_param_double(parameter_file, "Disk-PatchPotential:z_disk");
  potential->disk_patch_potential.timestep_mult = parser_get_param_double(
      parameter_file, "Disk-PatchPotential:timestep_mult");
  potential->disk_patch_potential.dynamical_time =
      sqrt(potential->disk_patch_potential.scale_height /
           (phys_const->const_newton_G *
            potential->disk_patch_potential.surface_density));
#endif /* EXTERNAL_POTENTIAL_DISK_PATCH */
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
void potential_print(const struct external_potential* potential) {

#ifdef EXTERNAL_POTENTIAL_POINTMASS

  message(
      "Point mass properties are (x,y,z) = (%e, %e, %e), M = %e timestep "
      "multiplier = %e",
      potential->point_mass.x, potential->point_mass.y, potential->point_mass.z,
      potential->point_mass.mass, potential->point_mass.timestep_mult);

#endif /* EXTERNAL_POTENTIAL_POINTMASS */

#ifdef EXTERNAL_POTENTIAL_ISOTHERMALPOTENTIAL

  message(
      "Isothermal potential properties are (x,y,z) = (%e, %e, %e), vrot = %e "
      "timestep multiplier= %e",
      potential->isothermal_potential.x, potential->isothermal_potential.y,
      potential->isothermal_potential.z, potential->isothermal_potential.vrot,
      potential->isothermal_potential.timestep_mult);

#endif /* EXTERNAL_POTENTIAL_ISOTHERMALPOTENTIAL */
#ifdef EXTERNAL_POTENTIAL_DISK_PATCH
  message(
      "Disk-patch potential properties are surface_density = %e disk height= "
      "%e scale height= %e timestep multiplier= %e",
      potential->disk_patch_potential.surface_density,
      potential->disk_patch_potential.z_disk,
      potential->disk_patch_potential.scale_height,
      potential->disk_patch_potential.timestep_mult);
#ifdef EXTERNAL_POTENTIAL_DISK_PATCH_ICS
  message(
      "Disk-patch potential: imposing growth of gravity over time, and adding "
      "viscous force to gravity");
#endif
#endif /* EXTERNAL_POTENTIAL_DISK_PATCH */
}
