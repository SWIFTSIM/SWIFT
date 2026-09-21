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
  message("Interpolation table size (mass) = %d", rad->interpolation_size);
  if (rad->is_2d) {
    message("Interpolation table size (metallicity) = %d",
            rad->interpolation_size_metallicity);
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

  /* Set before radiation_read_data() below, which requires the L_FUV/
     L_LW datasets whenever this is on (see #radiation.with_ISRF's own
     doxygen). */
  rad->with_ISRF = (char)parser_get_opt_param_int(
      params, "GEARFeedback:with_photoelectric_heating", 0);

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
 * #raw/#integrated's luminosities/dot_N_ion/dot_E_excess/teff/l_pe/l_lw
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
    interpolate_2d_free(&rad->integrated.luminosities_2d);
    interpolate_2d_free(&rad->integrated.dot_N_ion_2d);
    interpolate_2d_free(&rad->integrated.dot_E_excess_2d);
    interpolate_2d_free(&rad->integrated.l_pe_2d);
    interpolate_2d_free(&rad->integrated.l_lw_2d);
  } else {
    interpolate_1d_free(&rad->raw.luminosities);
    interpolate_1d_free(&rad->raw.dot_N_ion);
    interpolate_1d_free(&rad->raw.dot_E_excess);
    interpolate_1d_free(&rad->raw.teff);
    interpolate_1d_free(&rad->raw.l_pe);
    interpolate_1d_free(&rad->raw.l_lw);
    interpolate_1d_free(&rad->integrated.luminosities);
    interpolate_1d_free(&rad->integrated.dot_N_ion);
    interpolate_1d_free(&rad->integrated.dot_E_excess);
    interpolate_1d_free(&rad->integrated.l_pe);
    interpolate_1d_free(&rad->integrated.l_lw);
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

  /* Read before being overwritten below: dispatches which union member
     (1D or 2D) gets zero-pointered, matching the INCOMING dimensionality.
     Both members' `data` pointer sits at offset 0 (verified with
     offsetof()), so the wrong dispatch would still null the right
     pointer, just leave the other member's differently-shaped tail
     fields (xmin/dx/N/... vs xmin/dx/ymin/dy/Nx/Ny/...) only partially
     cleared.

     stellar_evolution_props_init()'s `!with_radiation` branch calls this
     on a fresh, never-yet-built #radiation whose is_2d may be garbage
     (part of a non-zeroed `struct feedback_props` on swift.c's stack):
     not a memory-safety hazard even then, since either dispatch branch
     still nulls the one pointer #radiation_clean's later
     interpolate_*d_free() call can dereference, and the tail fields are
     inert scalars always overwritten by the next
     #radiation_build_tables() call before anything reads them. */
  const char was_2d = rad->is_2d;

  rad->is_active = 0;
  rad->is_2d = 0;
  rad->interpolation_size = 0;
  rad->interpolation_size_metallicity = 0;
  rad->n_HII_pixels = 0;
  rad->age_max_myr = 0.f;
  rad->ms_lifetime_inverse_log_z_min = 0.f;
  rad->ms_lifetime_inverse_log_z_step = 0.f;
  rad->ms_lifetime_inverse_n_metallicity = 0;
  rad->with_ISRF = 0;
  rad->has_teff = 0;

  if (was_2d) {
    interpolate_2d_zero_pointers(&rad->raw.luminosities_2d);
    interpolate_2d_zero_pointers(&rad->raw.dot_N_ion_2d);
    interpolate_2d_zero_pointers(&rad->raw.dot_E_excess_2d);
    interpolate_2d_zero_pointers(&rad->raw.teff_2d);
    interpolate_2d_zero_pointers(&rad->raw.l_pe_2d);
    interpolate_2d_zero_pointers(&rad->raw.l_lw_2d);
    interpolate_2d_zero_pointers(&rad->integrated.luminosities_2d);
    interpolate_2d_zero_pointers(&rad->integrated.dot_N_ion_2d);
    interpolate_2d_zero_pointers(&rad->integrated.dot_E_excess_2d);
    interpolate_2d_zero_pointers(&rad->integrated.l_pe_2d);
    interpolate_2d_zero_pointers(&rad->integrated.l_lw_2d);
  } else {
    interpolate_1d_zero_pointers(&rad->raw.luminosities);
    interpolate_1d_zero_pointers(&rad->raw.dot_N_ion);
    interpolate_1d_zero_pointers(&rad->raw.dot_E_excess);
    interpolate_1d_zero_pointers(&rad->raw.teff);
    interpolate_1d_zero_pointers(&rad->raw.l_pe);
    interpolate_1d_zero_pointers(&rad->raw.l_lw);
    interpolate_1d_zero_pointers(&rad->integrated.luminosities);
    interpolate_1d_zero_pointers(&rad->integrated.dot_N_ion);
    interpolate_1d_zero_pointers(&rad->integrated.dot_E_excess);
    interpolate_1d_zero_pointers(&rad->integrated.l_pe);
    interpolate_1d_zero_pointers(&rad->integrated.l_lw);
  }

  interpolate_2d_zero_pointers(&rad->raw.main_sequence_lifetime_2d);
  interpolate_2d_zero_pointers(&rad->raw.main_sequence_lifetime_inverse_2d);
}
