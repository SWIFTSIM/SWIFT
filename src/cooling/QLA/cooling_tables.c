/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/**
 * @file src/cooling/QLA/cooling_tables.c
 * @brief Functions to read QLA tables
 */

/* Config parameters. */
#include <config.h>

/* This file's header */
#include "cooling_tables.h"

/* Standard includes */
#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "cooling_properties.h"
#include "error.h"
#include "exp10.h"
#include "interpolate.h"

/**
 * @brief Reads in PS2020 cooling table header. Consists of tables
 * of values for temperature, hydrogen number density, metallicity,
 * abundance ratios, and elements used to index the cooling tables.
 *
 * @param cooling Cooling data structure
 */
void read_cooling_header(struct cooling_function_data *cooling) {

#ifdef HAVE_HDF5

  hid_t dataset;
  herr_t status;

  /* read sizes of array dimensions */
  hid_t tempfile_id =
      H5Fopen(cooling->cooling_table_path, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (tempfile_id < 0)
    error("unable to open file %s\n", cooling->cooling_table_path);

  /* allocate arrays of bins */
  if (posix_memalign((void **)&cooling->Temp, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_temperature * sizeof(float)) != 0)
    error("Failed to allocate temperature table\n");

  if (posix_memalign((void **)&cooling->Redshifts, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * sizeof(float)) != 0)
    error("Failed to allocate redshift table\n");

  if (posix_memalign((void **)&cooling->nH, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_density * sizeof(float)) != 0)
    error("Failed to allocate density table\n");

  if (posix_memalign((void **)&cooling->Metallicity, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_metallicity * sizeof(float)) != 0)
    error("Failed to allocate metallicity table\n");

  if (posix_memalign((void **)&cooling->Therm, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_internalenergy * sizeof(float)) != 0)
    error("Failed to allocate internal energy table\n");

  if (posix_memalign((void **)&cooling->LogAbundances, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_metallicity * qla_cooling_N_elementtypes *
                         sizeof(float)) != 0)
    error("Failed to allocate abundance array\n");

  if (posix_memalign((void **)&cooling->Abundances, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_metallicity * qla_cooling_N_elementtypes *
                         sizeof(float)) != 0)
    error("Failed to allocate abundance array\n");

  if (posix_memalign((void **)&cooling->Abundances_inv, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_metallicity * qla_cooling_N_elementtypes *
                         sizeof(float)) != 0)
    error("Failed to allocate abundance array\n");

  if (posix_memalign((void **)&cooling->atomicmass, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_elementtypes * sizeof(float)) != 0)
    error("Failed to allocate atomic masses array\n");

  if (posix_memalign((void **)&cooling->atomicmass_inv, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_elementtypes * sizeof(float)) != 0)
    error("Failed to allocate inverse atomic masses array\n");

  if (posix_memalign((void **)&cooling->Zsol, SWIFT_STRUCT_ALIGNMENT,
                     1 * sizeof(float)) != 0)
    error("Failed to allocate solar metallicity array\n");

  if (posix_memalign((void **)&cooling->Zsol_inv, SWIFT_STRUCT_ALIGNMENT,
                     1 * sizeof(float)) != 0)
    error("Failed to allocate inverse solar metallicity array\n");

  if (posix_memalign((void **)&cooling->LogMassFractions,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_metallicity * qla_cooling_N_elementtypes *
                         sizeof(float)) != 0)
    error("Failed to allocate log mass fraction array\n");

  if (posix_memalign((void **)&cooling->MassFractions, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_metallicity * qla_cooling_N_elementtypes *
                         sizeof(float)) != 0)
    error("Failed to allocate mass fraction array\n");

  /* read in bins and misc information */
  dataset = H5Dopen(tempfile_id, "/TableBins/TemperatureBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Temp);
  if (status < 0) error("error reading temperature bins\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TableBins/RedshiftBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Redshifts);
  if (status < 0) error("error reading redshift bins\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TableBins/DensityBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->nH);
  if (status < 0) error("error reading density bins\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TableBins/MetallicityBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Metallicity);
  if (status < 0) error("error reading metallicity bins\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TableBins/InternalEnergyBins", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Therm);
  if (status < 0) error("error reading internal energy bins\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TotalAbundances", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->LogAbundances);
  if (status < 0) error("error reading total abundances\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/TotalMassFractions", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->LogMassFractions);
  if (status < 0) error("error reading total mass fractions\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/ElementMasses", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->atomicmass);
  if (status < 0) error("error reading element masses\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  dataset = H5Dopen(tempfile_id, "/SolarMetallicity", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->Zsol);
  if (status < 0) error("error reading solar metallicity \n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Close the file */
  H5Fclose(tempfile_id);

  cooling->Zsol_inv[0] = 1.f / cooling->Zsol[0];

  /* find the metallicity bin that refers to solar metallicity */
  const float tol = 1.e-3;
  for (int i = 0; i < qla_cooling_N_metallicity; i++) {
    if (fabsf(cooling->Metallicity[i]) < tol) {
      cooling->indxZsol = i;
    }
  }

#if defined(__ICC)
#pragma novector
#endif
  for (int i = 0; i < qla_cooling_N_elementtypes; i++) {
    cooling->atomicmass_inv[i] = 1.f / cooling->atomicmass[i];
  }

  /* set some additional useful abundance arrays */
  for (int i = 0; i < qla_cooling_N_metallicity; i++) {

#if defined(__ICC)
#pragma novector
#endif
    for (int j = 0; j < qla_cooling_N_elementtypes; j++) {
      const int indx1d = row_major_index_2d(i, j, qla_cooling_N_metallicity,
                                            qla_cooling_N_elementtypes);
      cooling->Abundances[indx1d] = exp10f(cooling->LogAbundances[indx1d]);
      cooling->Abundances_inv[indx1d] = 1.f / cooling->Abundances[indx1d];
      cooling->MassFractions[indx1d] =
          exp10f(cooling->LogMassFractions[indx1d]);
    }
  }

#else
  error("Need HDF5 to read cooling tables");
#endif
}

/**
 * @brief Allocate space for cooling tables and read them
 *
 * @param cooling #cooling_function_data structure
 */
void read_cooling_tables(struct cooling_function_data *restrict cooling) {

#ifdef HAVE_HDF5
  hid_t dataset;
  herr_t status;

  /* open hdf5 file */
  hid_t tempfile_id =
      H5Fopen(cooling->cooling_table_path, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (tempfile_id < 0)
    error("unable to open file %s\n", cooling->cooling_table_path);

  /* Allocate and read arrays to store cooling tables. */

  /* Mean particle mass (temperature) */
  if (swift_memalign("cooling_table.Tmu", (void **)&cooling->table.Tmu,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_temperature *
                         qla_cooling_N_metallicity * qla_cooling_N_density *
                         sizeof(float)) != 0)
    error("Failed to allocate Tmu array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/MeanParticleMass", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.Tmu);
  if (status < 0) error("error reading Tmu\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing mean particle mass dataset");

  /* Mean particle mass (internal energy) */
  if (swift_memalign("cooling_table.Umu", (void **)&cooling->table.Umu,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_internalenergy *
                         qla_cooling_N_metallicity * qla_cooling_N_density *
                         sizeof(float)) != 0)
    error("Failed to allocate Umu array\n");

  dataset = H5Dopen(tempfile_id, "/Udep/MeanParticleMass", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.Umu);
  if (status < 0) error("error reading Umu\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing mean particle mass dataset");

  /* Cooling (temperature) */
  if (swift_memalign("cooling_table.Tcooling",
                     (void **)&cooling->table.Tcooling, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_temperature *
                         qla_cooling_N_metallicity * qla_cooling_N_density *
                         qla_cooling_N_cooltypes * sizeof(float)) != 0)
    error("Failed to allocate Tcooling array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/Cooling", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.Tcooling);
  if (status < 0) error("error reading Tcooling\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Cooling (internal energy) */
  if (swift_memalign("cooling_table.Ucooling",
                     (void **)&cooling->table.Ucooling, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_internalenergy *
                         qla_cooling_N_metallicity * qla_cooling_N_density *
                         qla_cooling_N_cooltypes * sizeof(float)) != 0)
    error("Failed to allocate Ucooling array\n");

  dataset = H5Dopen(tempfile_id, "/Udep/Cooling", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.Ucooling);
  if (status < 0) error("error reading Ucooling\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Heating (temperature) */
  if (swift_memalign("cooling_table.Theating",
                     (void **)&cooling->table.Theating, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_temperature *
                         qla_cooling_N_metallicity * qla_cooling_N_density *
                         qla_cooling_N_heattypes * sizeof(float)) != 0)
    error("Failed to allocate Theating array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/Heating", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.Theating);
  if (status < 0) error("error reading Theating\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Heating (internal energy) */
  if (swift_memalign("cooling_table.Uheating",
                     (void **)&cooling->table.Uheating, SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_internalenergy *
                         qla_cooling_N_metallicity * qla_cooling_N_density *
                         qla_cooling_N_heattypes * sizeof(float)) != 0)
    error("Failed to allocate Uheating array\n");

  dataset = H5Dopen(tempfile_id, "/Udep/Heating", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.Uheating);
  if (status < 0) error("error reading Uheating\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Electron fraction (temperature) */
  if (swift_memalign("cooling_table.Tefrac",
                     (void **)&cooling->table.Telectron_fraction,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_temperature *
                         qla_cooling_N_metallicity * qla_cooling_N_density *
                         qla_cooling_N_electrontypes * sizeof(float)) != 0)
    error("Failed to allocate Telectron_fraction array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/ElectronFractionsVol", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.Telectron_fraction);
  if (status < 0) error("error reading electron_fraction (temperature)\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Electron fraction (internal energy) */
  if (swift_memalign("cooling_table.Uefrac",
                     (void **)&cooling->table.Uelectron_fraction,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_internalenergy *
                         qla_cooling_N_metallicity * qla_cooling_N_density *
                         qla_cooling_N_electrontypes * sizeof(float)) != 0)
    error("Failed to allocate Uelectron_fraction array\n");

  dataset = H5Dopen(tempfile_id, "/Udep/ElectronFractionsVol", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.Uelectron_fraction);
  if (status < 0) error("error reading electron_fraction (internal energy)\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Internal energy from temperature */
  if (swift_memalign("cooling_table.UfromT", (void **)&cooling->table.U_from_T,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_temperature *
                         qla_cooling_N_metallicity * qla_cooling_N_density *
                         sizeof(float)) != 0)
    error("Failed to allocate U_from_T array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/U_from_T", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.U_from_T);
  if (status < 0) error("error reading U_from_T array\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Temperature from interal energy */
  if (swift_memalign("cooling_table.TfromU", (void **)&cooling->table.T_from_U,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_internalenergy *
                         qla_cooling_N_metallicity * qla_cooling_N_density *
                         sizeof(float)) != 0)
    error("Failed to allocate T_from_U array\n");

  dataset = H5Dopen(tempfile_id, "/Udep/T_from_U", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.T_from_U);
  if (status < 0) error("error reading T_from_U array\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing cooling dataset");

  /* Thermal equilibrium temperature */
  if (swift_memalign("cooling_table.Teq", (void **)&cooling->table.logTeq,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_metallicity *
                         qla_cooling_N_density * sizeof(float)) != 0)
    error("Failed to allocate logTeq array\n");

  dataset = H5Dopen(tempfile_id, "/ThermEq/Temperature", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.logTeq);
  if (status < 0) error("error reading Teq array\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing logTeq dataset");

  /* Mean particle mass at thermal equilibrium temperature */
  if (swift_memalign("cooling_table.mueq",
                     (void **)&cooling->table.meanpartmass_Teq,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_metallicity *
                         qla_cooling_N_density * sizeof(float)) != 0)
    error("Failed to allocate mu array\n");

  dataset = H5Dopen(tempfile_id, "/ThermEq/MeanParticleMass", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.meanpartmass_Teq);
  if (status < 0) error("error reading mu array\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing mu dataset");

  /* Hydrogen fractions at thermal equilibirum temperature */
  if (swift_memalign("cooling_table.Hfracs",
                     (void **)&cooling->table.logHfracs_Teq,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_metallicity *
                         qla_cooling_N_density * 3 * sizeof(float)) != 0)
    error("Failed to allocate hydrogen fractions array\n");

  dataset = H5Dopen(tempfile_id, "/ThermEq/HydrogenFractionsVol", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.logHfracs_Teq);
  if (status < 0) error("error reading hydrogen fractions array\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing hydrogen fractions dataset");

  /* All hydrogen fractions */
  if (swift_memalign("cooling_table.Hfracs",
                     (void **)&cooling->table.logHfracs_all,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_temperature *
                         qla_cooling_N_metallicity * qla_cooling_N_density * 3 *
                         sizeof(float)) != 0)
    error("Failed to allocate big hydrogen fractions array\n");

  dataset = H5Dopen(tempfile_id, "/Tdep/HydrogenFractionsVol", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   cooling->table.logHfracs_all);
  if (status < 0) error("error reading big hydrogen fractions array\n");
  status = H5Dclose(dataset);
  if (status < 0) error("error closing big hydrogen fractions dataset");

  /* Close the file */
  H5Fclose(tempfile_id);

  /* Pressure at thermal equilibrium temperature */
  if (swift_memalign("cooling_table.Peq", (void **)&cooling->table.logPeq,
                     SWIFT_STRUCT_ALIGNMENT,
                     qla_cooling_N_redshifts * qla_cooling_N_metallicity *
                         qla_cooling_N_density * sizeof(float)) != 0)
    error("Failed to allocate logPeq array\n");

  const float log10_kB_cgs = cooling->log10_kB_cgs;

  /* Compute the pressures at thermal eq. */
  for (int ired = 0; ired < qla_cooling_N_redshifts; ired++) {
    for (int imet = 0; imet < qla_cooling_N_metallicity; imet++) {

      const int index_XH = row_major_index_2d(
          imet, 0, qla_cooling_N_metallicity, qla_cooling_N_elementtypes);

      const float log10_XH = cooling->LogMassFractions[index_XH];

      for (int iden = 0; iden < qla_cooling_N_density; iden++) {

        const int index_Peq = row_major_index_3d(
            ired, imet, iden, qla_cooling_N_redshifts,
            qla_cooling_N_metallicity, qla_cooling_N_density);

        cooling->table.logPeq[index_Peq] =
            cooling->nH[iden] + cooling->table.logTeq[index_Peq] - log10_XH -
            log10(cooling->table.meanpartmass_Teq[index_Peq]) + log10_kB_cgs;
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  message("Done reading in general cooling table");
#endif

#else
  error("Need HDF5 to read cooling tables");
#endif
}
