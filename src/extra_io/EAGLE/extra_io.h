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
#ifndef SWIFT_EXTRA_IO_EAGLE_H
#define SWIFT_EXTRA_IO_EAGLE_H

#include "extra.h"
#include "io_properties.h"

INLINE static void convert_part_Xray_photons(const struct engine *e,
                                             const struct part *p,
                                             const struct xpart *xp,
                                             double *ret) {

  ret[0] = extra_io_get_xray_fluxes(
      p, xp, e, xray_band_types_erosita_low_intrinsic_photons);
  ret[1] = extra_io_get_xray_fluxes(
      p, xp, e, xray_band_types_erosita_high_intrinsic_photons);
  ret[2] = extra_io_get_xray_fluxes(p, xp, e,
                                    xray_band_types_ROSAT_intrinsic_photons);
}

INLINE static void convert_part_Xray_energies(const struct engine *e,
                                              const struct part *p,
                                              const struct xpart *xp,
                                              double *ret) {

  ret[0] = extra_io_get_xray_fluxes(
      p, xp, e, xray_band_types_erosita_low_intrinsic_energies);
  ret[1] = extra_io_get_xray_fluxes(
      p, xp, e, xray_band_types_erosita_high_intrinsic_energies);
  ret[2] = extra_io_get_xray_fluxes(p, xp, e,
                                    xray_band_types_ROSAT_intrinsic_energies);
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extra particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology?
 *
 * @return Returns the number of fields to write.
 */
INLINE static int extra_io_write_particles(const struct part *parts,
                                           const struct xpart *xparts,
                                           struct io_props *list,
                                           const int with_cosmology) {

  list[0] = io_make_output_field_convert_part(
      "XrayPhotonLuminosities", DOUBLE, 3, UNIT_CONV_PHOTONS_PER_TIME, 0.f,
      parts, xparts, convert_part_Xray_photons,
      "Intrinsic X-ray photon luminosities in various bands. This is 0 for "
      "star-forming particles.");

  list[1] = io_make_output_field_convert_part(
      "XrayLuminosities", DOUBLE, 3, UNIT_CONV_POWER, 0.f, parts, xparts,
      convert_part_Xray_energies,
      "Intrinsic X-ray luminosities in various bands. This is 0 for "
      "star-forming particles.");

  return 2;
}

/**
 * @brief Specifies which star particle fields to write to a dataset
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int extra_io_write_sparticles(const struct spart *sparts,
                                            struct io_props *list,
                                            const int with_cosmology) {

  return 0;
}

/**
 * @brief Specifies which black hole particle fields to write to a dataset
 *
 * @param bparts The BH particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int extra_io_write_bparticles(const struct bpart *bparts,
                                            struct io_props *list,
                                            const int with_cosmology) {
  return 0;
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of extra-io to the file
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 */
INLINE static void extra_io_write_flavour(hid_t h_grp, hid_t h_grp_columns) {

  /* Write the extra-io model */
  io_write_attribute_s(h_grp, "Extra-io", "EAGLE");

  /* Create an array of xray band names */
  static const char xrayband_names[xray_band_types_count / 2][32] = {
      "erosita_low", "erosita_high", "ROSAT"};

  /* Add to the named columns. We do it twice. Once for
   * the energies and once for the photons.
   * The columns use the same names for both arrays. */
  hsize_t dims[1] = {xray_band_types_count / 2};
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, 32);
  hid_t space = H5Screate_simple(1, dims, NULL);
  hid_t dset = H5Dcreate(h_grp_columns, "XrayLuminosities", type, space,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, xrayband_names[0]);
  H5Dclose(dset);
  dset = H5Dcreate(h_grp_columns, "XrayPhotonLuminosities", type, space,
                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, xrayband_names[0]);
  H5Dclose(dset);

  H5Tclose(type);
  H5Sclose(space);
}
#endif

/*
  Extra lightcone map types
*/
/*
   Healpix map of intrinsic erosita-low photons band
*/
double lightcone_map_xray_erosita_low_intrinsic_photons_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);

/*
   Healpix map of intrinsic erosita-low energy band
*/
double lightcone_map_xray_erosita_low_intrinsic_energy_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);

/*
   Healpix map of intrinsic erosita-high photons band
*/
double lightcone_map_xray_erosita_high_intrinsic_photons_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);

/*
   Healpix map of intrinsic erosita-high energy band
*/
double lightcone_map_xray_erosita_high_intrinsic_energy_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);

/*
   Healpix map of intrinsic ROSAT photons band
*/
double lightcone_map_xray_rosat_intrinsic_photons_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);

/*
   Healpix map of intrinsic ROSAT energy band
*/
double lightcone_map_xray_rosat_intrinsic_energy_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);
/*
   Healpix map of compton y
*/
int lightcone_map_compton_y_type_contributes(int ptype);

double lightcone_map_compton_y_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);
/*
   Healpix map of doppler b
*/
int lightcone_map_doppler_b_type_contributes(int ptype);

double lightcone_map_doppler_b_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);
/*
   Healpix map of dispersion meassure
*/
int lightcone_map_dispersion_meassure_type_contributes(int ptype);

double lightcone_map_dispersion_meassure_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]);

static const struct lightcone_map_type extra_lightcone_map_types[] = {
    {
        /* .name = */ "XrayErositaLowIntrinsicPhotons",
        /* .update_map = */
        lightcone_map_xray_erosita_low_intrinsic_photons_get_value,
        /* .ptype_contributes = */ lightcone_map_gas_only,
        /* .baseline_func = */ NULL,
        /* .units = */ UNIT_CONV_PHOTON_FLUX_PER_UNIT_SURFACE,
        /* .smoothing = */ map_smoothed,
        /* .compression = */ compression_write_lossless,
        /* .buffer_scale_factor = */ 1.0e-62, /* to keep in range of a float */
    },
    {
        /* .name = */ "XrayErositaLowIntrinsicEnergies",
        /* .update_map = */
        lightcone_map_xray_erosita_low_intrinsic_energy_get_value,
        /* .ptype_contributes = */ lightcone_map_gas_only,
        /* .baseline_func = */ NULL,
        /* .units = */ UNIT_CONV_ENERGY_FLUX_PER_UNIT_SURFACE,
        /* .smoothing = */ map_smoothed,
        /* .compression = */ compression_write_lossless,
        /* .buffer_scale_factor = */ 1.0,
    },
    {
        /* .name = */ "XrayErositaHighIntrinsicPhotons",
        /* .update_map = */
        lightcone_map_xray_erosita_high_intrinsic_photons_get_value,
        /* .ptype_contributes = */ lightcone_map_gas_only,
        /* .baseline_func = */ NULL,
        /* .units = */ UNIT_CONV_PHOTON_FLUX_PER_UNIT_SURFACE,
        /* .smoothing = */ map_smoothed,
        /* .compression = */ compression_write_lossless,
        /* .buffer_scale_factor = */ 1.0e-62, /* to keep in range of a float */
    },
    {
        /* .name = */ "XrayErositaHighIntrinsicEnergies",
        /* .update_map = */
        lightcone_map_xray_erosita_high_intrinsic_energy_get_value,
        /* .ptype_contributes = */ lightcone_map_gas_only,
        /* .baseline_func = */ NULL,
        /* .units = */ UNIT_CONV_ENERGY_FLUX_PER_UNIT_SURFACE,
        /* .smoothing = */ map_smoothed,
        /* .compression = */ compression_write_lossless,
        /* .buffer_scale_factor = */ 1.0,
    },
    {
        /* .name = */ "XrayROSATIntrinsicPhotons",
        /* .update_map = */
        lightcone_map_xray_rosat_intrinsic_photons_get_value,
        /* .ptype_contributes = */ lightcone_map_gas_only,
        /* .baseline_func = */ NULL,
        /* .units = */ UNIT_CONV_PHOTON_FLUX_PER_UNIT_SURFACE,
        /* .smoothing = */ map_smoothed,
        /* .compression = */ compression_write_lossless,
        /* .buffer_scale_factor = */ 1.0e-62, /* to keep in range of a float */
    },
    {
        /* .name = */ "XrayROSATIntrinsicEnergies",
        /* .update_map = */ lightcone_map_xray_rosat_intrinsic_energy_get_value,
        /* .ptype_contributes = */ lightcone_map_gas_only,
        /* .baseline_func = */ NULL,
        /* .units = */ UNIT_CONV_ENERGY_FLUX_PER_UNIT_SURFACE,
        /* .smoothing = */ map_smoothed,
        /* .compression = */ compression_write_lossless,
        /* .buffer_scale_factor = */ 1.0,
    },
    {
        /* .name = */ "ComptonY",
        /* .update_map = */ lightcone_map_compton_y_get_value,
        /* .ptype_contributes = */ lightcone_map_gas_only,
        /* .baseline_func = */ NULL,
        /* .units = */ UNIT_CONV_NO_UNITS,
        /* .smoothing = */ map_smoothed,
        /* .compression = */ compression_write_lossless,
        /* .buffer_scale_factor = */ 1.0,
    },
    {
        /* .name = */ "DopplerB",
        /* .update_map = */ lightcone_map_doppler_b_get_value,
        /* .ptype_contributes = */ lightcone_map_gas_only,
        /* .baseline_func = */ NULL,
        /* .units = */ UNIT_CONV_NO_UNITS,
        /* .smoothing = */ map_smoothed,
        /* .compression = */ compression_write_lossless,
        /* .buffer_scale_factor = */ 1.0,
    },
    {
        /* .name = */ "DM",
        /* .update_map = */ lightcone_map_dispersion_meassure_get_value,
        /* .ptype_contributes = */ lightcone_map_gas_only,
        /* .baseline_func = */ NULL,
        /* .units = */ UNIT_CONV_INV_AREA,
        /* .smoothing = */ map_smoothed,
        /* .compression = */ compression_write_lossless,
        /* .buffer_scale_factor = */ 3.40367719e-68, /* convert 1/Mpc^2 to
                                                  pc/cm^3 so value fits in a
                                                  float */
    },
    {
        /* NULL functions indicate end of array */
        /* .name = */ "",
        /* .update_map = */ NULL,
        /* .ptype_contributes = */ NULL,
        /* .baseline_func = */ NULL,
        /* .units = */ UNIT_CONV_NO_UNITS,
        /* .smoothing = */ map_smoothed,
        /* .compression = */ compression_write_lossless,
        /* .buffer_scale_factor = */ 1.0,
    },
};

#endif /* SWIFT_EXTRA_IO_EAGLE_H */
