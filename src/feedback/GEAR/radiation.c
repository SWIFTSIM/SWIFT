/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
 * @file src/feedback/GEAR/radiation.c
 * @brief Lifecycle of the #radiation structure for GEAR: printing,
 * initialization (angular HII pixel setup), restart dump/restore, and
 * cleanup. Per-particle and per-star physics live in radiation_gas.c and
 * radiation_pressure.c; the star-emission getters in radiation_getters.c;
 * the HDF5 table reading in radiation_table_io.c.
 */

/* Include header */
#include "radiation.h"

#include "engine.h"
#include "interpolation.h"

#include <math.h>
#include <string.h>

double radiation_lw_photon_energy_cgs = 0.;

/*! Sanity tripwires on lambda(b) - 1 (and on lambda_N(LW)). Below the floor is
    fatal: a units bug (e.g. E_lo^2 folded in twice) gives a tiny positive
    value. Above the ceiling is a warning only: it is the largest value seen on
    the PARSEC grid, rounded up, not a physical bound. */
#define RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_FLOOR 0.05
#define RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_SANITY_MAX 250.0

/**
 * @brief Set #feedback_props.band_edge_weight_pe/lw/photon_weight_lw from the
 * radiation table.
 *
 * These are the coefficients for photons redshifting through the 6/11.2/13.6
 * eV band edges. They are evaluated once over the whole IMF (at the reference
 * metallicity for 2D tables), as for #radiation_set_lw_photon_energy_cgs.
 * The defaults are kept while radiation is inactive or a denominator
 * vanishes. Call at start-up only, for the main stellar model.
 *
 * @param fb_props (output) The #feedback_props to set.
 * @param rad The main stellar model's #radiation.
 * @param sm The main #stellar_model, for its IMF mass range.
 */
void radiation_set_band_edge_coefficients(struct feedback_props *fb_props,
                                          const struct radiation *rad,
                                          const struct stellar_model *sm) {

  fb_props->band_edge_weight_pe = RADIATION_BAND_EDGE_WEIGHT_PE_DEFAULT;
  fb_props->band_edge_weight_lw = RADIATION_BAND_EDGE_WEIGHT_LW_DEFAULT;
  fb_props->band_edge_photon_weight_lw =
      RADIATION_BAND_EDGE_PHOTON_WEIGHT_LW_DEFAULT;

  if (!rad->is_active || !rad->with_ISRF) return;

  const float log_m1 = log10f(sm->imf.mass_min);
  const float log_m2 = log10f(sm->imf.mass_max);
  const float log_z =
      rad->is_2d ? radiation_get_log_metallicity(
                       RADIATION_LW_PHOTON_ENERGY_REFERENCE_METALLICITY)
                 : 0.f;

  const double l_edge_pe =
      rad->is_2d
          ? radiation_get_luminosity_edge_pe_from_integral_2d(rad, log_z,
                                                              log_m1, log_m2)
          : radiation_get_luminosity_edge_pe_from_integral(rad, log_m1, log_m2);
  const double l_pe = rad->is_2d ? radiation_get_luminosity_pe_from_integral_2d(
                                       rad, log_z, log_m1, log_m2)
                                 : radiation_get_luminosity_pe_from_integral(
                                       rad, log_m1, log_m2);
  const double l_edge_lw =
      rad->is_2d
          ? radiation_get_luminosity_edge_lw_from_integral_2d(rad, log_z,
                                                              log_m1, log_m2)
          : radiation_get_luminosity_edge_lw_from_integral(rad, log_m1, log_m2);
  const double l_lw = rad->is_2d ? radiation_get_luminosity_lw_from_integral_2d(
                                       rad, log_z, log_m1, log_m2)
                                 : radiation_get_luminosity_lw_from_integral(
                                       rad, log_m1, log_m2);
  const double mean_e_lw_cgs =
      rad->is_2d
          ? radiation_get_mean_photon_energy_lw_from_integral_2d(rad, log_z,
                                                                 log_m2)
          : radiation_get_mean_photon_energy_lw_from_integral(rad, log_m2);

  /* Denominator guard: an IMF entirely below the table's mass floor gives
   * L_PE = L_LW = 0. Exact comparison, so a legitimately zero numerator is
   * still reported below. */
  if (l_pe <= 0. || l_lw <= 0.) return;

  /* E_lo^2 is already folded into l_edge_pe/lw at read time: do not
   * reapply it. */
  const double lambda_e_pe_minus_one = l_edge_pe / l_pe;
  const double lambda_e_lw_minus_one = l_edge_lw / l_lw;

  /* Numerator guard: a nonzero L_b with a zero edge term means the spectral
   * library does not reach the band edge. The ratio is then grey, which is
   * correct arithmetic but worth reporting. */
  if (engine_rank == 0 && lambda_e_pe_minus_one <= 0.)
    message(
        "WARNING: Data/Radiation's SpectralPhotonRateAtPEEdge integrates to "
        "0 over [%.4g, %.4g] Msun at Z=%.4g while Integrated_L_PE does not: "
        "the PE band-edge weight is GREY (lambda_E(PE)=1) for this run.",
        (double)sm->imf.mass_min, (double)sm->imf.mass_max,
        (double)exp10(log_z));
  if (engine_rank == 0 && lambda_e_lw_minus_one <= 0.)
    message(
        "WARNING: Data/Radiation's SpectralPhotonRateAtLWEdge integrates to "
        "0 over [%.4g, %.4g] Msun at Z=%.4g while Integrated_L_LW does not: "
        "the LW band-edge weight is GREY (lambda_E(LW)=1) for this run.",
        (double)sm->imf.mass_min, (double)sm->imf.mass_max,
        (double)exp10(log_z));

  /* Fatal tripwire for a positive but tiny lambda(b) - 1, the signature of a
   * units bug. Every rank aborts. */
  if (lambda_e_pe_minus_one > 0. &&
      lambda_e_pe_minus_one < RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_FLOOR)
    error(
        "lambda_E(PE) - 1 = %.6g is below the sanity floor %.3g: this "
        "looks like a units bug (e.g. E_lo^2 folded in twice), not a "
        "genuine table value.",
        lambda_e_pe_minus_one, RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_FLOOR);
  if (lambda_e_lw_minus_one > 0. &&
      lambda_e_lw_minus_one < RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_FLOOR)
    error(
        "lambda_E(LW) - 1 = %.6g is below the sanity floor %.3g: this "
        "looks like a units bug (e.g. E_lo^2 folded in twice), not a "
        "genuine table value.",
        lambda_e_lw_minus_one, RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_FLOOR);

  fb_props->band_edge_weight_pe = 1. + lambda_e_pe_minus_one;
  fb_props->band_edge_weight_lw = 1. + lambda_e_lw_minus_one;

  /* Lambda_LW = (lambda_E(LW) - 1) * <E>_LW / E_lo(LW). */
  if (mean_e_lw_cgs > 0.)
    fb_props->band_edge_photon_weight_lw = lambda_e_lw_minus_one *
                                           mean_e_lw_cgs /
                                           RADIATION_LW_BAND_LOWER_EDGE_CGS;

  /* No "+1" floor for lambda_N(LW), so the tripwire applies to the value. */
  if (fb_props->band_edge_photon_weight_lw > 0. &&
      fb_props->band_edge_photon_weight_lw <
          RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_FLOOR)
    error(
        "lambda_N(LW) = %.6g is below the sanity floor %.3g: this looks "
        "like a units bug (e.g. E_lo^2 folded in twice), not a genuine "
        "table value.",
        fb_props->band_edge_photon_weight_lw,
        RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_FLOOR);

  if (engine_rank == 0) {
    message(
        "Band-edge weights (redshift transfer across the 6/11.2/13.6 eV "
        "band edges), IMF-integrated over [%.4g, %.4g] Msun at Z=%.4g: "
        "lambda_E(PE)=%.5g, lambda_E(LW)=%.5g, lambda_N(LW)=%.5g",
        (double)sm->imf.mass_min, (double)sm->imf.mass_max,
        (double)exp10(log_z), fb_props->band_edge_weight_pe,
        fb_props->band_edge_weight_lw, fb_props->band_edge_photon_weight_lw);
    if (lambda_e_pe_minus_one > RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_SANITY_MAX)
      message(
          "WARNING: lambda_E(PE) - 1 = %.5g exceeds the observed grid "
          "maximum %.1f; check the table and its unit conversion before "
          "trusting this run's PE band-edge physics.",
          lambda_e_pe_minus_one,
          RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_SANITY_MAX);
    if (lambda_e_lw_minus_one > RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_SANITY_MAX)
      message(
          "WARNING: lambda_E(LW) - 1 = %.5g exceeds the observed grid "
          "maximum %.1f; check the table and its unit conversion before "
          "trusting this run's LW band-edge physics.",
          lambda_e_lw_minus_one,
          RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_SANITY_MAX);
    if (fb_props->band_edge_photon_weight_lw >
        RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_SANITY_MAX)
      message(
          "WARNING: lambda_N(LW) = %.5g exceeds the observed grid maximum "
          "%.1f; check the table and its unit conversion before trusting "
          "this run's LW photon-number band-edge physics.",
          fb_props->band_edge_photon_weight_lw,
          RADIATION_BAND_EDGE_WEIGHT_MINUS_ONE_SANITY_MAX);
  }
}

/**
 * @brief Report the H2 photodissociation coefficient and set
 * #radiation_lw_photon_energy_cgs from the radiation table.
 *
 * The mean LW photon energy is a diagnostic only: the IMF-integrated
 * Integrated_MeanPhotonEnergyLW up to mass_max, at the reference metallicity
 * for 2D tables. Left at 0 while radiation is inactive. Call for the main
 * stellar model only.
 *
 * @param rad The main stellar model's #radiation.
 * @param sm The main #stellar_model, for its IMF mass range.
 */
void radiation_set_lw_photon_energy_cgs(const struct radiation *rad,
                                        const struct stellar_model *sm) {

  radiation_lw_photon_energy_cgs = 0.;

  if (!rad->is_active) return;

  if (engine_rank == 0)
    message(
        "H2 photodissociation sigma_H2/E_LW = %g cm^2 erg^-1, the "
        "Sternberg-anchored constant the rate uses",
        RADIATION_SIGMA_H2_OVER_E_LW_CGS);

  /* Mean over the IMF from mass_min up to mass_max. */
  const float log_m = log10f(sm->imf.mass_max);
  const double E_LW_cgs =
      rad->is_2d
          ? radiation_get_mean_photon_energy_lw_from_integral_2d(
                rad,
                radiation_get_log_metallicity(
                    RADIATION_LW_PHOTON_ENERGY_REFERENCE_METALLICITY),
                log_m)
          : radiation_get_mean_photon_energy_lw_from_integral(rad, log_m);

  if (E_LW_cgs <= 0.) return;

  radiation_lw_photon_energy_cgs = E_LW_cgs;

  if (engine_rank == 0)
    message(
        "Mean Lyman-Werner photon energy from the table = %g erg, reported "
        "only: no rate reads it",
        radiation_lw_photon_energy_cgs);
}

/**
 * @brief Print the radiation model.
 *
 * @param rad The #radiation.
 */
void radiation_print(const struct radiation *rad) {

  /* Only the master print */
  if (engine_rank != 0) {
    return;
  }

  message("Angular pixels for HII ionization = %d", rad->n_HII_pixels);
  message("Number of masses interpolated onto = %d", rad->interpolation_size);

  /* The fields below come from the table. */
  if (!rad->is_active) return;

  int n_sources = 0;
  for (int i = 0; i < RADIATION_TABLE_SOURCE_COUNT; i++) {
    if (rad->table_source[i][0] == '\0') continue;
    message("Table %s = %s", radiation_table_source_keys[i],
            rad->table_source[i]);
    n_sources++;
  }
  if (n_sources == 0)
    message("Table sources = none (no provenance attribute in the table)");

  message("Number of masses in the table = %i", rad->table_n_mass);
  message("Mass range of the table (Msun) = [%g, %g]",
          (double)rad->table_mass_min, (double)rad->table_mass_max);

  message("Metallicity dependent table? %i", rad->is_2d);

  if (rad->is_2d) {
    message("Number of metallicities in the table = %i",
            rad->table_n_metallicity);
    message("Metallicity range of the table (mass fraction) = [%g, %g]",
            (double)rad->table_metallicity_min,
            (double)rad->table_metallicity_max);
  }
}

/**
 * @brief Initialize the #radiation structure.
 *
 * @param rad The #radiation model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param us The unit system.
 * @param phys_const The physical constants.
 */
void radiation_init(struct radiation *rad, struct swift_params *params,
                    const struct stellar_model *sm,
                    const struct unit_system *us,
                    const struct phys_const *phys_const) {

  /* Before radiation_read_data(), which needs it. */
  rad->with_ISRF = (char)parser_get_opt_param_int(
      params, "GEARFeedback:with_interstellar_radiation_field", 0);

  /* Read the data */
  radiation_read_data(rad, params, sm, us, phys_const, /* restart */ 0);

  /* HEALPix split of the HII budget: nside=0 is spherical (1 pixel), else
     12*nside^2 pixels. The ceiling is the per-star dot_N_ion_pix array, sized
     by --with-number-of-hii-angular-pixels. */
  const int nside =
      parser_get_opt_param_int(params, "GEARFeedback:HII_angular_nside", 0);
  if (nside < 0) {
    error("GEARFeedback:HII_angular_nside must be >= 0; got %d.", nside);
  }
  const int n_HII_pixels_requested = (nside == 0) ? 1 : 12 * nside * nside;
  if (n_HII_pixels_requested > HII_MAX_ANGULAR_PIXELS) {
    error(
        "GEARFeedback:HII_angular_nside=%d requires %d HealPix pixels, but "
        "this build only supports up to HII_MAX_ANGULAR_PIXELS=%d. "
        "Reconfigure with "
        "--with-number-of-hii-angular-pixels=%d (or higher) and rebuild, "
        "or lower HII_angular_nside.",
        nside, n_HII_pixels_requested, HII_MAX_ANGULAR_PIXELS,
        n_HII_pixels_requested);
  }
#ifndef HAVE_CHEALPIX
  if (nside != 0) {
    error(
        "GEARFeedback:HII_angular_nside > 0 requires the HEALPix C API "
        "(chealpix). Reconfigure with --with-chealpix, or set nside=0.");
  }
#endif
  rad->n_HII_pixels = n_HII_pixels_requested;
}

/**
 * @brief Write a radiation struct to the given FILE as a stream of bytes.
 *
 * @param rad the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
void radiation_dump(const struct radiation *rad, FILE *stream,
                    const struct stellar_model *sm) {

  restart_write_blocks((void *)rad, sizeof(struct radiation), 1, stream,
                       "radiation", "radiation");
  message("Dumping GEAR radiation...");
}

/**
 * @brief Restore a radiation struct from the given FILE as a stream of bytes.
 *
 * The restored interpolation pointers are stale, so radiation_read_data()
 * rebuilds the tables from sm->yields_table. That path must still resolve and
 * the file must be unchanged since the run started.
 *
 * @param rad the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 * @param us The unit system.
 * @param phys_const The physical constants in internal units.
 * @param with_radiation Are we restoring with photoionization and/or
 * radiation pressure?
 */
void radiation_restore(struct radiation *rad, FILE *stream,
                       const struct stellar_model *sm,
                       const struct unit_system *us,
                       const struct phys_const *phys_const,
                       const char with_radiation) {

  restart_read_blocks((void *)rad, sizeof(struct radiation), 1, stream, NULL,
                      "radiation");

  if (!with_radiation) {
    /* The restored bytes hold another process's heap addresses: never
       dereference or free them. radiation_zero_pointers() also resets
       is_active. */
    radiation_zero_pointers(rad);
    return;
  }

  radiation_read_data(rad, NULL, sm, us, phys_const, /*restart=*/1);
  message("Restoring GEAR radiation struct...");
}

/**
 * @brief Free the interpolation tables.
 *
 * @param rad the #radiation.
 */
void radiation_clean(struct radiation *rad) {

  /* is_2d selects the live union member, which must be freed with the
     matching interpolate_*d_free(). */
  if (rad->is_2d) {
    interpolate_2d_free(&rad->raw.luminosities_2d);
    interpolate_2d_free(&rad->raw.dot_N_ion_2d);
    interpolate_2d_free(&rad->raw.dot_E_excess_2d);
    interpolate_2d_free(&rad->raw.teff_2d);
    interpolate_2d_free(&rad->raw.l_pe_2d);
    interpolate_2d_free(&rad->raw.l_lw_2d);
    interpolate_2d_free(&rad->raw.l_edge_pe_2d);
    interpolate_2d_free(&rad->raw.l_edge_lw_2d);
    interpolate_2d_free(&rad->raw.mean_photon_energy_lw_2d);
    interpolate_2d_free(&rad->integrated.luminosities_2d);
    interpolate_2d_free(&rad->integrated.dot_N_ion_2d);
    interpolate_2d_free(&rad->integrated.dot_E_excess_2d);
    interpolate_2d_free(&rad->integrated.l_pe_2d);
    interpolate_2d_free(&rad->integrated.l_lw_2d);
    interpolate_2d_free(&rad->integrated.l_edge_pe_2d);
    interpolate_2d_free(&rad->integrated.l_edge_lw_2d);
    interpolate_2d_free(&rad->integrated.mean_photon_energy_lw_2d);
  } else {
    interpolate_1d_free(&rad->raw.luminosities);
    interpolate_1d_free(&rad->raw.dot_N_ion);
    interpolate_1d_free(&rad->raw.dot_E_excess);
    interpolate_1d_free(&rad->raw.teff);
    interpolate_1d_free(&rad->raw.l_pe);
    interpolate_1d_free(&rad->raw.l_lw);
    interpolate_1d_free(&rad->raw.l_edge_pe);
    interpolate_1d_free(&rad->raw.l_edge_lw);
    interpolate_1d_free(&rad->raw.mean_photon_energy_lw);
    interpolate_1d_free(&rad->integrated.luminosities);
    interpolate_1d_free(&rad->integrated.dot_N_ion);
    interpolate_1d_free(&rad->integrated.dot_E_excess);
    interpolate_1d_free(&rad->integrated.l_pe);
    interpolate_1d_free(&rad->integrated.l_lw);
    interpolate_1d_free(&rad->integrated.l_edge_pe);
    interpolate_1d_free(&rad->integrated.l_edge_lw);
    interpolate_1d_free(&rad->integrated.mean_photon_energy_lw);
  }

  /* No 1D counterpart: always freed. */
  interpolate_2d_free(&rad->raw.main_sequence_lifetime_2d);
  interpolate_2d_free(&rad->raw.main_sequence_lifetime_inverse_2d);
}

/**
 * @brief Zero a #radiation struct so radiation_clean() and the printers are
 * safe. Getters must check #is_active first.
 *
 * @param rad The #radiation.
 */
void radiation_zero_pointers(struct radiation *rad) {

  rad->is_active = 0;
  rad->is_2d = 0;
  rad->interpolation_size = 0;
  rad->n_HII_pixels = 0;
  rad->age_max_myr = 0.f;
  rad->with_ISRF = 0;
  rad->has_teff = 0;
  for (int i = 0; i < RADIATION_TABLE_SOURCE_COUNT; i++)
    rad->table_source[i][0] = '\0';
  rad->table_mass_min = 0.f;
  rad->table_mass_max = 0.f;
  rad->table_n_mass = 0;
  rad->table_n_metallicity = 0;
  rad->table_metallicity_min = 0.f;
  rad->table_metallicity_max = 0.f;

  /* All-bits-zero nulls every pointer and zeroes both union members, and
     boundary_condition_error must be the first enumerator. Nothing is freed:
     after a restart the bytes are another process's heap addresses. */
  memset(&rad->raw, 0, sizeof(rad->raw));
  memset(&rad->integrated, 0, sizeof(rad->integrated));
}
