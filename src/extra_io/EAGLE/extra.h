/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EXTRA_EAGLE_H
#define SWIFT_EXTRA_EAGLE_H

#include "chemistry.h"
#include "cooling.h"
#include "cooling/PS2020/cooling_tables.h"
#include "engine.h"
#include "star_formation.h"

#define xray_table_date_string 20240406

#define xray_emission_N_temperature 46
#define xray_emission_N_density 71
#define xray_emission_N_helium 10
#define xray_emission_N_element 10
#define xray_emission_N_redshift 45

/**
 * @brief X-ray bands available in the interpolation tables
 */
enum xray_band_types {
  xray_band_types_erosita_low_intrinsic_photons,   /*< eROSITA 0.2 - 2.3 keV */
  xray_band_types_erosita_high_intrinsic_photons,  /*< eROSITA 2.3 - 8.0 keV */
  xray_band_types_ROSAT_intrinsic_photons,         /*< ROSAT 0.5 - 2.0 keV */
  xray_band_types_erosita_low_intrinsic_energies,  /*< eROSITA 0.2 - 2.3 keV */
  xray_band_types_erosita_high_intrinsic_energies, /*< eROSITA 2.3 - 8.0 keV */
  xray_band_types_ROSAT_intrinsic_energies,        /*< ROSAT 0.5 - 2.0 keV */
  xray_band_types_count
};

/**
 * @brief The general properties required for the extra i/o fields.
 */
struct extra_io_properties {

  struct xray_properties {

    /* Element masses for the chemistry elements (cgs) */
    float *element_mass;

    /* Temperature bins from xray table (cgs) */
    float *Temperatures;

    /* Minimum and maximum temperature the table exists for */
    float Temperature_min;
    float Temperature_max;

    /* Density bins from xray table (physical cgs) */
    float *Densities;

    /* Minimum and maximum density the table exists for */
    float Density_min;
    float Density_max;

    /* Helium fraction bins from xray table */
    float *He_bins;

    /* Redshift bins from xray table */
    float *Redshifts;

    /* Maximum redshift the table exists for */
    float Redshift_max;

    /* Solar metallicites from xray table */
    float *Solar_metallicity;

    /* Log of solar metallicites from xray table */
    float *Log10_solar_metallicity;

    /* Integrated photon emissivity in the erosita-low band (0.2-2.3 keV)
     * (physical) */
    float *emissivity_erosita_low_intrinsic_photons;

    /* Integrated photon emissivity in the erosita-high band (2.3-8.0 keV)
     * (physical) */
    float *emissivity_erosita_high_intrinsic_photons;

    /* Integrated photon emissivity in the ROSAT band (0.5-2.0 keV) (physical)
     */
    float *emissivity_ROSAT_intrinsic_photons;

    /* Integrated emissivity in the erosita-low band (0.2-2.3 keV)
     * (physical) */
    float *emissivity_erosita_low_intrinsic_energies;

    /* Integrated emissivity in the erosita-high band (2.3-8.0 keV)
     * (physical) */
    float *emissivity_erosita_high_intrinsic_energies;

    /* Integrated emissivity in the ROSAT band (0.5-2.0 keV) (physical)
     */
    float *emissivity_ROSAT_intrinsic_energies;

    /* Path to the xray table */
    char xray_table_path[500];

    /* Photon emissivity unit conversion factor */
    double xray_photon_emissivity_unit_conversion;

    /* Energy emissivity unit conversion factor */
    double xray_energy_emissivity_unit_conversion;
  } xray_data;
};

/**
 * @brief Reads in xray table header. Consists of tables
 * of values for temperature, hydrogen number density, helium fraction,
 * solar metallicity, redshifts, and element masses.
 *
 * @param xrays Xray data structure
 * @param fname Xray table path
 */
INLINE static void read_xray_header(struct xray_properties *xrays,
                                    const char *fname) {

  hid_t tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (tempfile_id < 0) error("unable to open file %s\n", fname);

  /* Check whether the correct table version is being used */
  int datestring;

  hid_t dataset = H5Dopen(tempfile_id, "Date_String", H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, &datestring);
  if (status < 0) error("error reading the date string");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
  if (datestring != xray_table_date_string)
    error(
        "The table and code version do not match, please use table version %i",
        xray_table_date_string);

  /* Read temperature bins */
  if (posix_memalign((void **)&xrays->Temperatures, SWIFT_STRUCT_ALIGNMENT,
                     xray_emission_N_temperature * sizeof(float)) != 0)
    error("Failed to allocate temperatures array\n");

  dataset = H5Dopen(tempfile_id, "Bins/Temperature_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   xrays->Temperatures);
  if (status < 0) error("error reading temperatures");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* Read density bins */
  if (posix_memalign((void **)&xrays->Densities, SWIFT_STRUCT_ALIGNMENT,
                     xray_emission_N_density * sizeof(float)) != 0)
    error("Failed to allocate densities array\n");

  dataset = H5Dopen(tempfile_id, "Bins/Density_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   xrays->Densities);
  if (status < 0) error("error reading densities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* Read Helium bins */
  if (posix_memalign((void **)&xrays->He_bins, SWIFT_STRUCT_ALIGNMENT,
                     xray_emission_N_helium * sizeof(float)) != 0)
    error("Failed to allocate He_bins array\n");

  dataset = H5Dopen(tempfile_id, "Bins/He_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   xrays->He_bins);
  if (status < 0) error("error reading Helium massfractions");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* Read solar metallicity */
  if (posix_memalign((void **)&xrays->Log10_solar_metallicity,
                     SWIFT_STRUCT_ALIGNMENT,
                     chemistry_element_count * sizeof(float)) != 0)
    error("Failed to allocate Solar_metallicity array\n");

  dataset = H5Dopen(tempfile_id, "Bins/Solar_metallicities", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   xrays->Log10_solar_metallicity);
  if (status < 0) error("error reading solar metalicities");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* Get Solar metallicities from log solar metallicities */
  if (posix_memalign((void **)&xrays->Solar_metallicity, SWIFT_STRUCT_ALIGNMENT,
                     chemistry_element_count * sizeof(float)) != 0)
    error("Failed to allocate Solar_metallicity array\n");

  for (int i = 0; i < chemistry_element_count; ++i)
    xrays->Solar_metallicity[i] = exp10f(xrays->Log10_solar_metallicity[i]);

  /* Read redshift bins */
  if (posix_memalign((void **)&xrays->Redshifts, SWIFT_STRUCT_ALIGNMENT,
                     xray_emission_N_redshift * sizeof(float)) != 0)
    error("Failed to allocate Redshifts array\n");

  dataset = H5Dopen(tempfile_id, "Bins/Redshift_bins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   xrays->Redshifts);
  if (status < 0) error("error reading redshift bins");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* Read element mass */
  if (posix_memalign((void **)&xrays->element_mass, SWIFT_STRUCT_ALIGNMENT,
                     10 * sizeof(float)) != 0)
    error("Failed to allocate element_mass array\n");

  dataset = H5Dopen(tempfile_id, "Bins/Element_masses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   xrays->element_mass);
  if (status < 0) error("error reading element masses");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");
}

/**
 * @brief Reads in xray table. Consists of tables
 * of values for xray emissivities in different bands.
 * We read the erosita-low, erosita-high and ROSAT band
 * in their intrinsic forms (no multiplication with response function)
 *
 * @param xrays Xray data structure
 * @param fname Xray table path
 */
INLINE static void read_xray_table(struct xray_properties *xrays,
                                   const char *fname) {

  /* Open File */
  hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  // erosita-low intrinsic photons
  if (swift_memalign("xrays_table_erosita_low_photons",
                     (void **)&xrays->emissivity_erosita_low_intrinsic_photons,
                     SWIFT_STRUCT_ALIGNMENT,
                     xray_emission_N_redshift * xray_emission_N_helium *
                         xray_emission_N_element * xray_emission_N_temperature *
                         xray_emission_N_density * sizeof(float)) != 0)
    error(
        "Failed to allocate xray emissivity_erosita_low_intrinsic_photons "
        "array\n");

  /* Read full table */
  hid_t dataset =
      H5Dopen(file_id, "erosita-low/photons_intrinsic", H5P_DEFAULT);
  herr_t status =
      H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
              xrays->emissivity_erosita_low_intrinsic_photons);
  if (status < 0) error("error reading X-Ray table\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* erosita-high intrinsic photons */
  if (swift_memalign("xrays_table_erosita_high_photons",
                     (void **)&xrays->emissivity_erosita_high_intrinsic_photons,
                     SWIFT_STRUCT_ALIGNMENT,
                     xray_emission_N_redshift * xray_emission_N_helium *
                         xray_emission_N_element * xray_emission_N_temperature *
                         xray_emission_N_density * sizeof(float)) != 0)
    error(
        "Failed to allocate xray emissivity_erosita_high_intrinsic_photons "
        "array\n");

  /* Read full table */
  dataset = H5Dopen(file_id, "erosita-high/photons_intrinsic", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   xrays->emissivity_erosita_high_intrinsic_photons);
  if (status < 0) error("error reading X-Ray table\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* ROSAT intrinsic photons */
  if (swift_memalign("xray_table_ROSAT_photons",
                     (void **)&xrays->emissivity_ROSAT_intrinsic_photons,
                     SWIFT_STRUCT_ALIGNMENT,
                     xray_emission_N_redshift * xray_emission_N_helium *
                         xray_emission_N_element * xray_emission_N_temperature *
                         xray_emission_N_density * sizeof(float)) != 0)
    error("Failed to allocate xray emissivity_ROSAT_intrinsic_photons array\n");

  /* Read full table */
  dataset = H5Dopen(file_id, "ROSAT/photons_intrinsic", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   xrays->emissivity_ROSAT_intrinsic_photons);
  if (status < 0) error("error reading X-Ray table\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  // erosita-low intrinsic energies
  if (swift_memalign("xrays_table_erosita_low_energies",
                     (void **)&xrays->emissivity_erosita_low_intrinsic_energies,
                     SWIFT_STRUCT_ALIGNMENT,
                     xray_emission_N_redshift * xray_emission_N_helium *
                         xray_emission_N_element * xray_emission_N_temperature *
                         xray_emission_N_density * sizeof(float)) != 0)
    error(
        "Failed to allocate xray emissivity_erosita_low_intrinsic_energies "
        "array\n");

  /* Read full table */
  dataset = H5Dopen(file_id, "erosita-low/energies_intrinsic", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   xrays->emissivity_erosita_low_intrinsic_energies);
  if (status < 0) error("error reading X-Ray table\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* erosita-high intrinsic energies */
  if (swift_memalign(
          "xrays_table_erosita_high_energies",
          (void **)&xrays->emissivity_erosita_high_intrinsic_energies,
          SWIFT_STRUCT_ALIGNMENT,
          xray_emission_N_redshift * xray_emission_N_helium *
              xray_emission_N_element * xray_emission_N_temperature *
              xray_emission_N_density * sizeof(float)) != 0)
    error(
        "Failed to allocate xray emissivity_erosita_high_intrinsic_energies "
        "array\n");

  /* Read full table */
  dataset = H5Dopen(file_id, "erosita-high/energies_intrinsic", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   xrays->emissivity_erosita_high_intrinsic_energies);
  if (status < 0) error("error reading X-Ray table\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* ROSAT intrinsic energies */
  if (swift_memalign("xray_table_ROSAT_energies",
                     (void **)&xrays->emissivity_ROSAT_intrinsic_energies,
                     SWIFT_STRUCT_ALIGNMENT,
                     xray_emission_N_redshift * xray_emission_N_helium *
                         xray_emission_N_element * xray_emission_N_temperature *
                         xray_emission_N_density * sizeof(float)) != 0)
    error(
        "Failed to allocate xray emissivity_ROSAT_intrinsic_energies array\n");

  /* Read full table */
  dataset = H5Dopen(file_id, "ROSAT/energies_intrinsic", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   xrays->emissivity_ROSAT_intrinsic_energies);
  if (status < 0) error("error reading X-Ray table\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing dataset");

  /* Close file */
  status = H5Fclose(file_id);
  if (status < 0) error("error closing file");
}

/**
 * @brief Find the 1d index for a table dimension
 *
 * @param table 1d array with binned values
 * @param size dimensionality
 * @param x value for which we aim to find the index
 * @param i (return) index
 * @param dx (return) offset from index bin
 */
INLINE static void get_index_1d(const float *restrict table, const int size,
                                const float x, int *i, float *restrict dx) {

  const float epsilon = 1e-4f;

  const float delta = (size - 1) / (table[size - 1] - table[0]);

  if (x < table[0] + epsilon) {
    /* We are below the first element */
    *i = 0;
    *dx = 0.f;
  } else if (x < table[size - 1] - epsilon) {
    /* Normal case */
    *i = (x - table[0]) * delta;
    *dx = (x - table[*i]) * delta;
  } else {
    /* We are after the last element */
    *i = size - 2;
    *dx = 1.f;
  }
}

/**
 * @brief Find the 1d index for a table dimension with irregularly spaced bins
 *
 * @param table 1d array with monotonically increasing binned values
 * @param size dimensionality
 * @param x value for which we aim to find the index
 * @param i (return) index
 * @param dx (return) offset from index bin
 */
INLINE static void get_index_1d_irregular(const float *restrict table,
                                          const int size, const float x, int *i,
                                          float *restrict dx) {
  const float epsilon = 1e-6f;

  if (x < table[0] + epsilon) {

    *i = 0;
    *dx = 0.f;

  } else if (x < table[size - 1] - epsilon) {

    int min_idx = -1;

    /* Do this the hard way: Search the table
     * for the largest index i in table[] such
     * that table[i] < x */
    for (int idx = 0; idx < size; idx++) {

      if (x - table[idx] <= 0.f) {

        /* Found the first entry that is larger than x, go back by 1. */
        min_idx = idx - 1;
        break;
      }
    }

    *i = min_idx;
    *dx = (x - table[min_idx]) / (table[min_idx + 1] - table[min_idx]);

  } else {

    *i = size - 2;
    *dx = 1.f;
  }
}

/**
 * @brief Returns the 1d index of element with 5d indices x,y,z,w
 * from a flattened 5d array in row major order
 *
 * @param x, y, z, v, w Indices of element of interest
 * @param Nx, Ny, Nz, Nv, Nw Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_5d(
    const int x, const int y, const int z, const int w, const int v,
    const int Nx, const int Ny, const int Nz, const int Nw, const int Nv) {

#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
  assert(z < Nz);
  assert(w < Nw);
  assert(v < Nv);
#endif

  return x * Ny * Nz * Nw * Nv + y * Nz * Nw * Nv + z * Nw * Nv + w * Nv + v;
}

/**
 * @brief 4d interpolation of the Xray table
 *
 * @param emissivity xray table
 * @param element number table index for missing element
 * @param n_H_index Index along the Hydrogen number density dimension
 * @param He_index Index along the Helium abundance dimension
 * @param T_index Index along temperature dimension
 * @param z_index Index along redshift dimension
 * @param d_nH Offset between Hydrogen density and table[n_H_index]
 * @param d_He Offset between Helium abundance and table[He_index]
 * @param d_T Offset between temperture and table[T_index]
 * @param d_z Offset between redshift and table[red_index]
 *
 * @return The log10 of the emssisivity
 */
INLINE static float interpolate_xray(const float *emissivity,
                                     const int element_number,
                                     const int nH_index, const int He_index,
                                     const int T_index, const int z_index,
                                     const float d_nH, const float d_He,
                                     const float d_T, const float d_z) {
  const float t_nH = 1.f - d_nH;
  const float t_He = 1.f - d_He;
  const float t_T = 1.f - d_T;
  const float t_z = 1.f - d_z;

  float result = 0.f;

  result += t_nH * t_He * t_T * t_z *
            emissivity[row_major_index_5d(
                z_index + 0, He_index + 0, element_number, T_index + 0,
                nH_index + 0, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += t_nH * t_He * d_T * t_z *
            emissivity[row_major_index_5d(
                z_index + 0, He_index + 0, element_number, T_index + 1,
                nH_index + 0, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += t_nH * d_He * t_T * t_z *
            emissivity[row_major_index_5d(
                z_index + 0, He_index + 1, element_number, T_index + 0,
                nH_index + 0, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += d_nH * t_He * t_T * t_z *
            emissivity[row_major_index_5d(
                z_index + 0, He_index + 0, element_number, T_index + 0,
                nH_index + 1, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += t_nH * d_He * d_T * t_z *
            emissivity[row_major_index_5d(
                z_index + 0, He_index + 1, element_number, T_index + 1,
                nH_index + 0, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += d_nH * t_He * d_T * t_z *
            emissivity[row_major_index_5d(
                z_index + 0, He_index + 0, element_number, T_index + 1,
                nH_index + 1, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += d_nH * d_He * t_T * t_z *
            emissivity[row_major_index_5d(
                z_index + 0, He_index + 1, element_number, T_index + 0,
                nH_index + 1, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += d_nH * d_He * d_T * t_z *
            emissivity[row_major_index_5d(
                z_index + 0, He_index + 1, element_number, T_index + 1,
                nH_index + 1, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += t_nH * t_He * t_T * d_z *
            emissivity[row_major_index_5d(
                z_index + 1, He_index + 0, element_number, T_index + 0,
                nH_index + 0, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += t_nH * t_He * d_T * d_z *
            emissivity[row_major_index_5d(
                z_index + 1, He_index + 0, element_number, T_index + 1,
                nH_index + 0, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += t_nH * d_He * t_T * d_z *
            emissivity[row_major_index_5d(
                z_index + 1, He_index + 1, element_number, T_index + 0,
                nH_index + 0, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += d_nH * t_He * t_T * d_z *
            emissivity[row_major_index_5d(
                z_index + 1, He_index + 0, element_number, T_index + 0,
                nH_index + 1, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += t_nH * d_He * d_T * d_z *
            emissivity[row_major_index_5d(
                z_index + 1, He_index + 1, element_number, T_index + 1,
                nH_index + 0, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += d_nH * t_He * d_T * d_z *
            emissivity[row_major_index_5d(
                z_index + 1, He_index + 0, element_number, T_index + 1,
                nH_index + 1, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += d_nH * d_He * t_T * d_z *
            emissivity[row_major_index_5d(
                z_index + 1, He_index + 1, element_number, T_index + 0,
                nH_index + 1, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  result += d_nH * d_He * d_T * d_z *
            emissivity[row_major_index_5d(
                z_index + 1, He_index + 1, element_number, T_index + 1,
                nH_index + 1, xray_emission_N_redshift, xray_emission_N_helium,
                xray_emission_N_element, xray_emission_N_temperature,
                xray_emission_N_density)];

  return result;
}

/**
 * @brief Find index of the Xray table to interpolate between and compute
 * emissivities
 *
 * @param xrays Xray data structure
 * @param He_fraction Helium fraction
 * @param log_nH_cgs physical number density in CGS
 * @param log_T temperature
 * @param redshift redshift
 * @param solar_ratio abundance ratio relative to solar
 * @param band xray band to use
 *
 * @return The X-ray emmisivity in the corresponding band in CGS units.
 */
INLINE static float do_xray_interpolation(
    const struct xray_properties *xrays, const float log10_He_fraction,
    const float log_nH_cgs, const float log_T, const float redshift,
    const float solar_ratio[colibre_cooling_N_elementtypes],
    const enum xray_band_types band) {

  /* Get indices in the interpolation table along the He, nH, T
   * and z dimensions */
  int He_index, log_nH_cgs_index, log_T_index, z_index;
  float d_He, d_log_nH_cgs, d_log_T, d_z;
  get_index_1d_irregular(xrays->He_bins, xray_emission_N_helium,
                         log10_He_fraction, &He_index, &d_He);

  get_index_1d(xrays->Densities, xray_emission_N_density, log_nH_cgs,
               &log_nH_cgs_index, &d_log_nH_cgs);

  get_index_1d(xrays->Temperatures, xray_emission_N_temperature, log_T,
               &log_T_index, &d_log_T);

  get_index_1d(xrays->Redshifts, xray_emission_N_redshift, redshift, &z_index,
               &d_z);

  /* Select the table corresponding to this band */
  float *table;
  switch (band) {
    case xray_band_types_erosita_low_intrinsic_photons:
      table = xrays->emissivity_erosita_low_intrinsic_photons;
      break;
    case xray_band_types_erosita_high_intrinsic_photons:
      table = xrays->emissivity_erosita_high_intrinsic_photons;
      break;
    case xray_band_types_ROSAT_intrinsic_photons:
      table = xrays->emissivity_ROSAT_intrinsic_photons;
      break;
    case xray_band_types_erosita_low_intrinsic_energies:
      table = xrays->emissivity_erosita_low_intrinsic_energies;
      break;
    case xray_band_types_erosita_high_intrinsic_energies:
      table = xrays->emissivity_erosita_high_intrinsic_energies;
      break;
    case xray_band_types_ROSAT_intrinsic_energies:
      table = xrays->emissivity_ROSAT_intrinsic_energies;
      break;
    default:
      error("Band doesn't exist");
  }

  /* The total flux is computed in two steps:
   * - First, we compute the contribution excluding all elements.
   * - Next, we loop over all the metals we trace and compute the flux
   *   in a case one metal is added
   * - Finally, for each metal, we add a portion of the difference
   *   between the (one metal) and (no metals) case based on the
   *   the abundance of each metal
   *
   * The interpolation table structure is as follows:
   * First we have the individual element contributions in the same order
   * as the PS2020 cooling model. As we only include metals, the first
   * entry is Carbon, second Nitrogen, etc.
   * The contribution of no metals is the last entry.
   *
   * The table size is hence:
   * colibre_cooling_N_elementtypes - 3 (H + He + OA) + 1 (no metals)
   * which is equal to xray_emission_N_element
   */

  /* Perform the interpolation of the no metal case
   * Note: That is stored as the last entry of the table */
  const float log10_x_ray_no_metals_cgs = interpolate_xray(
      table, xray_emission_N_element - 1, log_nH_cgs_index, He_index,
      log_T_index, z_index, d_log_nH_cgs, d_He, d_log_T, d_z);

  const float x_ray_no_metals_cgs = exp10f(log10_x_ray_no_metals_cgs);
  float x_ray_cgs = x_ray_no_metals_cgs;

  /* Loop over the *individual metals* used in the PS2020 cooling */
  for (int elem = element_C; elem <= element_Fe; elem++) {

    /* Note: we deduct 2 since the interpolation tables do not include H and He
     * and start straight with the metals */
    const float log10_x_ray_elem_cgs = interpolate_xray(
        table, elem - 2, log_nH_cgs_index, He_index, log_T_index, z_index,
        d_log_nH_cgs, d_He, d_log_T, d_z);

    const float x_ray_elem_cgs = exp10f(log10_x_ray_elem_cgs);

    /* Add the difference multiplied by the  abundance to solar-abundance
     * ratio */
    x_ray_cgs += x_ray_elem_cgs * solar_ratio[elem];
  }

  /* Convert from cm^3 to cm^-3 (i.e. multiply by nH^2) */
  x_ray_cgs *= exp10f(2.f * log_nH_cgs);

  return x_ray_cgs;
}

/**
 * @brief Compute the emmisivity of a particle in a given X-ray band.
 *
 * Particles that are star-forming or have a (rho,T) pair outside
 * the table range return a flux of 0.
 *
 * @param p The #part.
 * @param xp The corresponding #xpart.
 * @param e The #engine.
 * @param band xray band to use
 *
 * @return The emissivity in internal units.
 */
INLINE static double extra_io_get_xray_fluxes(const struct part *p,
                                              const struct xpart *xp,
                                              const struct engine *e,
                                              const enum xray_band_types band) {

  /* Get gas particle temperature */
  const float T = cooling_get_temperature(
      e->physical_constants, e->hydro_properties, e->internal_units,
      e->cosmology, e->cooling_func, p, xp);
  const float log10_T = log10f(T);

  /* Get gas particle Hydrogen number density in cgs */
  const float rho_phys = hydro_get_physical_density(p, e->cosmology);
  const float XH =
      chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];
  const float nH = rho_phys * XH / e->physical_constants->const_proton_mass *
                   e->cooling_func->number_density_to_cgs;
  const float log10_nH_cgs = log10f(nH);

  /* If the particle is not in the table range or star-forming, we return a flux
   * of 0 */
  if ((log10_T < e->io_extra_props->xray_data.Temperature_min ||
       log10_T > e->io_extra_props->xray_data.Temperature_max) ||
      (log10_nH_cgs < e->io_extra_props->xray_data.Density_min ||
       log10_nH_cgs > e->io_extra_props->xray_data.Density_max) ||
      e->cosmology->z > e->io_extra_props->xray_data.Redshift_max ||
      star_formation_get_SFR(p, xp) > 0.)
    return 0.;

  /* Get gas particle element mass fractions */
  const float *const mass_fractions =
      chemistry_get_metal_mass_fraction_for_cooling(p);

  /* Convert to abundances. For now, ignore Ca and S that are not tracked */
  float abundances[chemistry_element_count];
  for (int el = 0; el < chemistry_element_count; el++) {
    abundances[el] = (mass_fractions[el] / mass_fractions[0]) *
                     (e->io_extra_props->xray_data.element_mass[0] /
                      e->io_extra_props->xray_data.element_mass[el]);
  }

  /* We now need to convert the array we received from the chemistry
   * module (likely EAGLE) into the PS2020-cooling format.
   * This means adding un-tracked elements and changing their order */

  /* Finally onvert to abundances relative to solar */
  float abundance_ratio[colibre_cooling_N_elementtypes];
  for (int el = 0; el < colibre_cooling_N_elementtypes; el++) {

    /* Treat all regular elements */
    if (el <= element_Si) {

      abundance_ratio[el] =
          abundances[el] / e->io_extra_props->xray_data.Solar_metallicity[el];

      /* Special case for the two elements not traced in the chemistry */
    } else if (el == element_S || el == element_Ca) {

      /* S and Ca are fixed to have the same abundance ratio as Si */
      abundance_ratio[el] = abundance_ratio[element_Si];

      /* Final special case: Iron. */
    } else if (el == element_Fe) {

      /* We need to fish it out of the chemistry where it was at a different
       * location in the array */
      abundance_ratio[el] =
          abundances[chemistry_element_Fe] /
          e->io_extra_props->xray_data.Solar_metallicity[chemistry_element_Fe];
    } else {

      /* Any other element is not used in the Xray interpolation */
      abundance_ratio[el] = 0.f;
    }
  }

  /* Extract the (log of) Helium abundance */
  const float log10_He_fraction = log10f(abundances[chemistry_element_He]);

  /* Compute the X-ray emission in the given band */
  const double xray_em_cgs = do_xray_interpolation(
      &e->io_extra_props->xray_data, log10_He_fraction, log10_nH_cgs, log10_T,
      e->cosmology->z, abundance_ratio, band);

  /* Convert back to internal units */
  double xray_em;
  switch (band) {
    case xray_band_types_erosita_low_intrinsic_photons:
    case xray_band_types_erosita_high_intrinsic_photons:
    case xray_band_types_ROSAT_intrinsic_photons:
      xray_em =
          xray_em_cgs /
          e->io_extra_props->xray_data.xray_photon_emissivity_unit_conversion;
      break;
    case xray_band_types_erosita_low_intrinsic_energies:
    case xray_band_types_erosita_high_intrinsic_energies:
    case xray_band_types_ROSAT_intrinsic_energies:
      xray_em =
          xray_em_cgs /
          e->io_extra_props->xray_data.xray_energy_emissivity_unit_conversion;
      break;
    default:
      error("Band doesn't exist");
  }

  /* Now compute the luminosity from the emissivity
   *  To do so, we multiply by the particle volume
   *  luminosity = emissivity * (mass / density)
   */
  const double xray_lum = xray_em * (hydro_get_mass(p) / rho_phys);

  return xray_lum;
}

/**
 * @brief Initialises properties stored for the extra i/o fields
 *
 * @param parameter_file The parsed parameter file
 * @param us Internal system of units data structure
 * @param phys_const #phys_const data structure
 * @param cosmo The cosmology model
 * @param props #extra_io_properties struct to initialize
 */
INLINE static void extra_io_init(struct swift_params *parameter_file,
                                 const struct unit_system *us,
                                 const struct phys_const *phys_const,
                                 const struct cosmology *cosmo,
                                 struct extra_io_properties *props) {

  parser_get_param_string(parameter_file, "XrayEmissivity:xray_table_path",
                          props->xray_data.xray_table_path);

  read_xray_header(&props->xray_data, props->xray_data.xray_table_path);
  read_xray_table(&props->xray_data, props->xray_data.xray_table_path);

  /* Find the minimum and maximum density and temperature and the maximum
     redshift Print this information to the screen*/
  props->xray_data.Density_min = props->xray_data.Densities[0];
  props->xray_data.Density_max =
      props->xray_data.Densities[xray_emission_N_density - 1];

  props->xray_data.Temperature_min = props->xray_data.Temperatures[0];
  props->xray_data.Temperature_max =
      props->xray_data.Temperatures[xray_emission_N_temperature - 1];

  props->xray_data.Redshift_max =
      props->xray_data.Redshifts[xray_emission_N_redshift - 1];

  message(
      "X-ray broad band interpolation for particles between densities of "
      "nH=%f-%f cm-3"
      "temperature of logT=%f-%f K"
      "and redshift less than z<%f",
      props->xray_data.Density_min, props->xray_data.Density_max,
      props->xray_data.Temperature_min, props->xray_data.Temperature_max,
      props->xray_data.Redshift_max);

  /* Compute unit conversions only once and use them throughout */
  props->xray_data.xray_photon_emissivity_unit_conversion =
      units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY_PER_TIME);
  props->xray_data.xray_energy_emissivity_unit_conversion =
      units_cgs_conversion_factor(us, UNIT_CONV_POWER_DENSITY);
}

/**
 * @brief Free the memory allocated for the extra i/o fields.
 *
 * @param props #extra_io_properties struct to clean
 */
INLINE static void extra_io_clean(struct extra_io_properties *props) {

  free(props->xray_data.Temperatures);
  free(props->xray_data.Densities);
  free(props->xray_data.He_bins);
  free(props->xray_data.Solar_metallicity);
  free(props->xray_data.Log10_solar_metallicity);
  free(props->xray_data.Redshifts);
  free(props->xray_data.element_mass);

  swift_free("xrays_table_erosita_low_photons",
             props->xray_data.emissivity_erosita_low_intrinsic_photons);
  swift_free("xrays_table_erosita_high_photons",
             props->xray_data.emissivity_erosita_high_intrinsic_photons);
  swift_free("xray_table_ROSAT_photons",
             props->xray_data.emissivity_ROSAT_intrinsic_photons);
  swift_free("xrays_table_erosita_low_energies",
             props->xray_data.emissivity_erosita_low_intrinsic_energies);
  swift_free("xrays_table_erosita_high_energies",
             props->xray_data.emissivity_erosita_high_intrinsic_energies);
  swift_free("xray_table_ROSAT_energies",
             props->xray_data.emissivity_ROSAT_intrinsic_energies);
}

/**
 * @brief Write a extra i/o struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
INLINE static void extra_io_struct_dump(const struct extra_io_properties *props,
                                        FILE *stream) {

  struct extra_io_properties props_copy = *props;

  props_copy.xray_data.Temperatures = NULL;
  props_copy.xray_data.Densities = NULL;
  props_copy.xray_data.He_bins = NULL;
  props_copy.xray_data.Solar_metallicity = NULL;
  props_copy.xray_data.Log10_solar_metallicity = NULL;
  props_copy.xray_data.Redshifts = NULL;
  props_copy.xray_data.element_mass = NULL;

  props_copy.xray_data.emissivity_erosita_low_intrinsic_photons = NULL;
  props_copy.xray_data.emissivity_erosita_high_intrinsic_photons = NULL;
  props_copy.xray_data.emissivity_ROSAT_intrinsic_photons = NULL;
  props_copy.xray_data.emissivity_erosita_low_intrinsic_energies = NULL;
  props_copy.xray_data.emissivity_erosita_high_intrinsic_energies = NULL;
  props_copy.xray_data.emissivity_ROSAT_intrinsic_energies = NULL;

  restart_write_blocks((void *)&props_copy, sizeof(struct extra_io_properties),
                       1, stream, "extra_io", "extra i/o properties");
}

/**
 * @brief Restore a extra_io_properties struct from the given FILE as a
 * stream of bytes.
 *
 * Read the structure from the stream and restore the extra i/o tables by
 * re-reading them.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
INLINE static void extra_io_struct_restore(struct extra_io_properties *props,
                                           FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct extra_io_properties), 1,
                      stream, NULL, "extra i/o properties");

  read_xray_header(&props->xray_data, props->xray_data.xray_table_path);
  read_xray_table(&props->xray_data, props->xray_data.xray_table_path);
}

#endif /* SWIFT_EXTRA_EAGLE_H */
