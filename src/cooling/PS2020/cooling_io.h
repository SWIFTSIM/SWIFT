/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_PS2020_IO_H
#define SWIFT_COOLING_PS2020_IO_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "cooling.h"
#include "engine.h"
#include "io_properties.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of cooling to the file.
 *
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param cooling The #cooling_function_data
 */
__attribute__((always_inline)) INLINE static void cooling_write_flavour(
    hid_t h_grp, hid_t h_grp_columns,
    const struct cooling_function_data* cooling) {

  io_write_attribute_s(h_grp, "Cooling Model", "PS2020");

  const int number_of_species = 3;
  const int species_name_length = 4;
  char species_names[number_of_species][species_name_length];
  sprintf(species_names[0], "%s", "HI");
  sprintf(species_names[1], "%s", "HII");
  sprintf(species_names[2], "%s", "H2");

  /* Add the species names to the named columns */
  hsize_t dims[1] = {number_of_species};
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, species_name_length);
  hid_t space = H5Screate_simple(1, dims, NULL);
  hid_t dset = H5Dcreate(h_grp_columns, "SpeciesFractions", type, space,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, species_names[0]);
  H5Dclose(dset);

  H5Tclose(type);
  H5Sclose(space);
}
#endif

INLINE static void convert_part_T(const struct engine* e, const struct part* p,
                                  const struct xpart* xp, float* ret) {

  ret[0] = cooling_get_temperature(e->physical_constants, e->hydro_properties,
                                   e->internal_units, e->cosmology,
                                   e->cooling_func, p, xp);
}

INLINE static void convert_part_sub_T(const struct engine* e,
                                      const struct part* p,
                                      const struct xpart* xp, float* ret) {

  ret[0] = cooling_get_particle_subgrid_temperature(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);
}

INLINE static void convert_part_sub_rho(const struct engine* e,
                                        const struct part* p,
                                        const struct xpart* xp, float* ret) {

  ret[0] = cooling_get_particle_subgrid_density(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);
}

INLINE static void convert_part_sub_species_frac(const struct engine* e,
                                                 const struct part* p,
                                                 const struct xpart* xp,
                                                 float* ret) {

  ret[0] = cooling_get_particle_subgrid_HI_fraction(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);

  ret[1] = cooling_get_particle_subgrid_HII_fraction(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);

  ret[2] = cooling_get_particle_subgrid_H2_fraction(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);

  /* normalize the sum of the hydrogen fractions to 1 */
  const float sum = ret[0] + ret[1] + 2. * ret[2];
  ret[0] /= sum;
  ret[1] /= sum;
  ret[2] /= sum;
}

INLINE static void convert_part_HI_mass(const struct engine* e,
                                        const struct part* p,
                                        const struct xpart* xp, float* ret) {

  const float X_H =
      chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];

  const float HI_frac = cooling_get_particle_subgrid_HI_fraction(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);

  *ret = hydro_get_mass(p) * X_H * HI_frac;
}

INLINE static void convert_part_H2_mass(const struct engine* e,
                                        const struct part* p,
                                        const struct xpart* xp, float* ret) {

  const float X_H =
      chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];

  const float H2_frac = cooling_get_particle_subgrid_H2_fraction(
      e->internal_units, e->physical_constants, e->cosmology,
      e->hydro_properties, e->entropy_floor, e->cooling_func, p, xp);

  *ret = hydro_get_mass(p) * X_H * H2_frac * 2.f;
}

INLINE static void convert_part_e_density(const struct engine* e,
                                          const struct part* p,
                                          const struct xpart* xp, double* ret) {

  *ret = cooling_get_electron_density(e->physical_constants,
                                      e->hydro_properties, e->internal_units,
                                      e->cosmology, e->cooling_func, p, xp);
}

INLINE static void convert_part_y_compton(const struct engine* e,
                                          const struct part* p,
                                          const struct xpart* xp, double* ret) {

  *ret = cooling_get_ycompton(e->physical_constants, e->hydro_properties,
                              e->internal_units, e->cosmology, e->cooling_func,
                              p, xp);
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int cooling_write_particles(
    const struct part* parts, const struct xpart* xparts,
    struct io_props* list) {

  list[0] = io_make_output_field_convert_part(
      "Temperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, parts, xparts,
      convert_part_T, "Temperatures of the gas particles");

  list[1] = io_make_output_field_convert_part(
      "SubgridTemperatures", FLOAT, 1, UNIT_CONV_TEMPERATURE, 0.f, parts,
      xparts, convert_part_sub_T,
      "The subgrid temperatures if the particles are within deltaT of the "
      "entropy floor the subgrid temperature is calculated assuming a "
      "pressure equilibrium on the entropy floor, if the particles are "
      "above deltaT of the entropy floor the subgrid temperature is "
      "identical to the SPH temperature.");

  list[2] = io_make_output_field_convert_part(
      "SubgridPhysicalDensities", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, parts,
      xparts, convert_part_sub_rho,
      "The subgrid physical density if the particles are within deltaT of the "
      "entropy floor the subgrid density is calculated assuming a pressure "
      "equilibrium on the entropy floor, if the particles are above deltaT "
      "of the entropy floor the subgrid density is identical to the "
      "physical SPH density.");

  list[3] = io_make_output_field_convert_part(
      "SpeciesFractions", FLOAT, 3, UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
      convert_part_sub_species_frac,
      "Fractions of neutral, ionized and molecular hydrogen: [nHI/nH, nHII/nH, "
      "nH2/nH], assuming equilibrium "
      "tables. If the particles are within deltaT of the entropy floor the "
      "fractions are calculated using the subgrid quantities, i.e. assuming a "
      "pressure equilibrium on the entropy floor. If the particles are "
      "above deltaT of the entropy floor, the normal hydro quantities are "
      "used.");

  list[4] = io_make_output_field_convert_part(
      "AtomicHydrogenMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
      convert_part_HI_mass,
      "Atomic hydrogen masses containted in the particles. This quantity is "
      "obtained from the cooling tables and, if the particle is on the entropy "
      "floor, by extrapolating to the equilibrium curve assuming constant "
      "pressure.");

  list[5] = io_make_output_field_convert_part(
      "MolecularHydrogenMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts, xparts,
      convert_part_H2_mass,
      "Molecular hydrogen masses containted in the particles. This quantity is "
      "obtained from the cooling tables and, if the particle is on the entropy "
      "floor, by extrapolating to the equilibrium curve assuming constant "
      "pressure.");

  list[6] = io_make_output_field_convert_part(
      "ElectronNumberDensities", DOUBLE, 1, UNIT_CONV_NUMBER_DENSITY, 0.f,
      parts, xparts, convert_part_e_density,
      "Electron number densities in the physical frame computed based on the "
      "cooling tables. This is 0 for star-forming particles.");

  list[7] = io_make_output_field_convert_part(
      "ComptonYParameters", DOUBLE, 1, UNIT_CONV_AREA, 0.f, parts, xparts,
      convert_part_y_compton,
      "Compton y parameters in the physical frame computed based on the "
      "cooling tables. This is 0 for star-forming particles.");

  return 8;
}

#endif /* SWIFT_COOLING_PS2020_IO_H */
