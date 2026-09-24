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

/**
 * @brief Report the H2 photodissociation coefficient and set
 * #radiation_lw_photon_energy_cgs from the radiation table.
 *
 * The rate itself reads #RADIATION_SIGMA_H2_OVER_E_LW_CGS and nothing else,
 * so the coefficient in force is announced here unconditionally: a run's own
 * log is then the record of which calibration it used.
 *
 * The population mean LW photon energy is a diagnostic beside it. A gas
 * particle's LW band sums emission from many stars and keeps no record of
 * which star contributed what, so no emitter's own mean photon energy is
 * recoverable at the consumer; one representative population value is read
 * instead, pychem's "Integrated_MeanPhotonEnergyLW" over the IMF's whole
 * mass range, at #RADIATION_LW_PHOTON_ENERGY_REFERENCE_METALLICITY for a 2D
 * table. Left at 0 while radiation is inactive.
 *
 * Call this for the main stellar model only. A run with a first-stars table
 * reads two models, and the gas-side consumer is source-anonymous.
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

  /* The population getters take a single upper mass bound and average from
     the IMF's own mass_min up to it; see
     radiation_get_mean_photon_energy_lw_from_integral(). */
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

  /* Every field below is read from the table itself, so it is meaningless
     for a run that never opened one. */
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
 */
void radiation_init(struct radiation *rad, struct swift_params *params,
                    const struct stellar_model *sm,
                    const struct unit_system *us,
                    const struct phys_const *phys_const) {

  /* Set before radiation_read_data() below, which requires the L_PE/
     L_LW datasets whenever this is on (see #radiation.with_ISRF's own
     doxygen). */
  rad->with_ISRF = (char)parser_get_opt_param_int(
      params, "GEARFeedback:with_interstellar_radiation_field", 0);

  /* Read the data */
  radiation_read_data(rad, params, sm, us, phys_const, /* restart */ 0);

  /* Angular (HEALPix) splitting of the HII ionization budget.
     - nside=0 means spherical (HEALPix disabled, today's behaviour,
     n_HII_pixels=1);
     - any nside>=1 means the standard HEALPix RING-scheme tessellation
     (n_HII_pixels=12*nside^2).
     Note that the practical ceiling is memory, not geometry: every star
     carries a fixed-size dot_N_ion_pix[HII_MAX_ANGULAR_PIXELS] array sized by
     ./configure --with-number-of-hii-angular-pixels (default 12, i.e.
     nside<=1). */
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
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
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
 * The flat restore below copies the interpolation tables' internal data
 * pointers as raw bytes, meaningless in the new process, since they held
 * the old process's heap addresses. radiation_read_data() re-derives those
 * tables from scratch instead of trying to serialize them, avoiding ever
 * leaving a dangling pointer for radiation_clean() to free(). Re-derivation
 * reads sm->yields_table again rather than recomputing from mass/Z alone, so
 * it is exact only if that path still resolves and  the file is unchanged
 * since the run started, the same uncanonicalized-path caveat already noted
 * for GEARFeedback:yields_table in general. A restart resubmitted from a
 * different working directory with a relative path can fail here.
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
    /* The raw struct bytes are always read back (radiation_dump() always
       writes them, unlike e.g. stellar_wind_dump()), but the tables, and
       #radiation.is_active, are only re-derived, and sm->yields_table only
       re-opened, when this is set; otherwise radiation_zero_pointers
       overwrites whatever stale value the raw restore above just wrote into
       is_active.
       The bytes just restored are another process's heap addresses; never
       dereference or free them. */
    radiation_zero_pointers(rad);
    return;
  }

  radiation_read_data(rad, NULL, sm, us, phys_const, /*restart=*/1);
  message("Restoring GEAR radiation struct...");
}

/**
 * @brief Clean the allocated memory.
 *
 * #raw/#integrated's luminosities/dot_N_ion/dot_E_excess/teff/l_pe/l_lw/
 * mean_photon_energy_lw
 * fields are each an anonymous union of a #interpolation_1d and a
 * #interpolation_2d variant.
 *
 * @param rad the #radiation.
 */
void radiation_clean(struct radiation *rad) {

  /* is_2d selects which one is actually live and must be freed via the matching
     interpolate_*d_free(). Freeing through the other union member's helper on
     aliased memory would be wrong. */
  if (rad->is_2d) {
    interpolate_2d_free(&rad->raw.luminosities_2d);
    interpolate_2d_free(&rad->raw.dot_N_ion_2d);
    interpolate_2d_free(&rad->raw.dot_E_excess_2d);
    interpolate_2d_free(&rad->raw.teff_2d);
    interpolate_2d_free(&rad->raw.l_pe_2d);
    interpolate_2d_free(&rad->raw.l_lw_2d);
    interpolate_2d_free(&rad->raw.mean_photon_energy_lw_2d);
    interpolate_2d_free(&rad->integrated.luminosities_2d);
    interpolate_2d_free(&rad->integrated.dot_N_ion_2d);
    interpolate_2d_free(&rad->integrated.dot_E_excess_2d);
    interpolate_2d_free(&rad->integrated.l_pe_2d);
    interpolate_2d_free(&rad->integrated.l_lw_2d);
    interpolate_2d_free(&rad->integrated.mean_photon_energy_lw_2d);
  } else {
    interpolate_1d_free(&rad->raw.luminosities);
    interpolate_1d_free(&rad->raw.dot_N_ion);
    interpolate_1d_free(&rad->raw.dot_E_excess);
    interpolate_1d_free(&rad->raw.teff);
    interpolate_1d_free(&rad->raw.l_pe);
    interpolate_1d_free(&rad->raw.l_lw);
    interpolate_1d_free(&rad->raw.mean_photon_energy_lw);
    interpolate_1d_free(&rad->integrated.luminosities);
    interpolate_1d_free(&rad->integrated.dot_N_ion);
    interpolate_1d_free(&rad->integrated.dot_E_excess);
    interpolate_1d_free(&rad->integrated.l_pe);
    interpolate_1d_free(&rad->integrated.l_lw);
    interpolate_1d_free(&rad->integrated.mean_photon_energy_lw);
  }

  /* main_sequence_lifetime_2d/main_sequence_lifetime_inverse_2d have no 1D
     counterpart, so they are always freed unconditionally. */
  interpolate_2d_free(&rad->raw.main_sequence_lifetime_2d);
  interpolate_2d_free(&rad->raw.main_sequence_lifetime_inverse_2d);
}

/**
 * @brief Zero a #radiation struct so it can be safely passed to
 * #radiation_clean, printed, or read by any getter (which must check
 * #is_active first; the zeroed pointers below are not otherwise guarded).
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

  /* All-bits-zero nulls every table pointer and zeroes every scalar in
     both members of every union, whichever one #is_2d selects, and
     boundary_condition_error is the first enumerator so the boundary
     policies land on it too. No table is freed here: the bytes may be
     another process's heap addresses on the restart path (see
     radiation_restore()), and the live-table case goes through
     radiation_clean(). */
  memset(&rad->raw, 0, sizeof(rad->raw));
  memset(&rad->integrated, 0, sizeof(rad->integrated));
}
