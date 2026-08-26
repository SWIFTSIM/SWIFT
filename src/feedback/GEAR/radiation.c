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

  /* Read the data */
  radiation_read_data(rad, params, sm, us, phys_const, /* restart */ 0);

  /* Angular (HEALPix) splitting of the HII ionization budget. nside=0 means
     spherical (HEALPix disabled, today's behaviour, n_HII_pixels=1); any
     nside>=1 means the standard HEALPix RING-scheme tessellation
     (n_HII_pixels=12*nside^2). RING, unlike NEST, has no power-of-2
     restriction on nside, so any positive integer is mathematically valid
     here (see /usr/include/chealpix.h). The practical ceiling is memory,
     not geometry: every star carries a fixed-size
     dot_N_ion_pix[HII_MAX_ANGULAR_PIXELS] array
     (src/feedback/GEAR_thermal/feedback_struct.h) sized by
     ./configure --with-number-of-hii-angular-pixels (default 12, i.e.
     nside<=1); a run requesting more pixels than that build was
     configured for errors clearly below rather than overflowing the
     array. */
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
 * @brief Restore a radiation struct from the given FILE as a stream of
 * bytes.
 *
 * The flat restore below copies the interpolation tables' internal data
 * pointers as raw bytes, meaningless in the new process, since they held
 * the old process's heap addresses. radiation_read_data() re-derives those
 * tables from scratch instead of trying to serialize them, avoiding ever
 * leaving a dangling pointer for radiation_clean() to free(). Re-derivation
 * reads sm->yields_table again (Data/Radiation) rather than recomputing
 * from mass/Z alone, so it is exact only if that path still resolves and
 * the file is unchanged since the run started, the same uncanonicalized-
 * path caveat already noted for GEARFeedback:yields_table in general
 * (feedback_properties.h); a restart resubmitted from a different working
 * directory with a relative path can fail here.
 *
 * @param rad the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 * @param us The unit system.
 * @param phys_const The physical constants in internal units.
 * @param with_radiation Are we restoring with photoionization and/or
 * radiation pressure? The raw struct bytes are always read back
 * (radiation_dump() always writes them, unlike e.g. stellar_wind_dump()),
 * but the tables, and #radiation.is_active, are only re-derived, and
 * sm->yields_table only re-opened, when this is set; otherwise
 * #radiation_zero_pointers overwrites whatever stale value the raw restore
 * above just wrote into is_active.
 */
void radiation_restore(struct radiation *rad, FILE *stream,
                       const struct stellar_model *sm,
                       const struct unit_system *us,
                       const struct phys_const *phys_const,
                       const char with_radiation) {

  restart_read_blocks((void *)rad, sizeof(struct radiation), 1, stream, NULL,
                      "radiation");

  if (!with_radiation) {
    /* The bytes just restored are another process's heap addresses (see the
       function's own doxygen above); never dereference or free them. */
    radiation_zero_pointers(rad);
    return;
  }

  radiation_read_data(rad, NULL, sm, us, phys_const, /*restart=*/1);
  message("Restoring GEAR radiation struct...");
}

/**
 * @brief Clean the allocated memory.
 *
 * #raw/#integrated's luminosities/dot_N_ion/dot_E_excess fields are each an
 * anonymous union of a #interpolation_1d and a #interpolation_2d variant
 * (see #radiation's own doxygen); #is_2d selects which one is actually live
 * and must be freed via the matching interpolate_*d_free(). Freeing
 * through the other union member's helper on aliased memory would be wrong
 * (though not unsafe here: both #interpolation_1d and #interpolation_2d
 * have their `data` pointer as their first member, verified with
 * `offsetof()`, so either helper happens to null/free the right pointer
 * regardless of #is_2d. Not relied upon, since the non-pointer tail
 * fields still differ, but noted as the reason this is not a
 * memory-safety hazard even if #is_2d were ever wrong here).
 * #main_sequence_lifetime_2d/#main_sequence_lifetime_inverse_2d have no 1D
 * counterpart, so they are always freed unconditionally.
 *
 * @param rad the #radiation.
 */
void radiation_clean(struct radiation *rad) {

  if (rad->is_2d) {
    interpolate_2d_free(&rad->raw.luminosities_2d);
    interpolate_2d_free(&rad->raw.dot_N_ion_2d);
    interpolate_2d_free(&rad->raw.dot_E_excess_2d);
    interpolate_2d_free(&rad->integrated.luminosities_2d);
    interpolate_2d_free(&rad->integrated.dot_N_ion_2d);
    interpolate_2d_free(&rad->integrated.dot_E_excess_2d);
  } else {
    interpolate_1d_free(&rad->raw.luminosities);
    interpolate_1d_free(&rad->raw.dot_N_ion);
    interpolate_1d_free(&rad->raw.dot_E_excess);
    interpolate_1d_free(&rad->integrated.luminosities);
    interpolate_1d_free(&rad->integrated.dot_N_ion);
    interpolate_1d_free(&rad->integrated.dot_E_excess);
  }

  interpolate_2d_free(&rad->raw.main_sequence_lifetime_2d);
  interpolate_2d_free(&rad->raw.main_sequence_lifetime_inverse_2d);
}

/**
 * @brief Zero a #radiation struct (pointers, dimensions and #is_active)
 * so a struct that was never (or is not yet) initialized can be safely
 * passed to #radiation_clean, printed, or read by any of the getters (which
 * must check #is_active first; the getters themselves do not, and would
 * dereference the zeroed pointers below). Also called internally by
 * #radiation_read_data before it (re)builds the tables; see the scalar
 * save/restore around that call for why is_2d is the only field here it can
 * leave zeroed going in.
 *
 * #raw/#integrated's union'd fields (see #radiation's own doxygen) are
 * zero-pointered via the union member matching the INCOMING #is_2d (read
 * into @c was_2d before it is overwritten below), not both. Zero-
 * pointering through the other member's helper would target the same
 * `data` pointer (offset 0 in both #interpolation_1d and #interpolation_2d,
 * verified with `offsetof()`) but leave the two structs' differently-shaped
 * tail fields (xmin/dx/N/boundary_condition vs xmin/dx/ymin/dy/Nx/Ny/
 * boundary_condition_x/_y) only partially cleared.
 *
 * One caller reads @c was_2d from memory this function cannot prove was
 * ever initialised: #stellar_evolution_props_init's `!with_radiation`
 * branch calls this directly on a fresh, never-yet-built #radiation (part
 * of a #stellar_model that is itself part of a non-`bzero`'d
 * `struct feedback_props` on the stack in swift.c): @c was_2d there may
 * be garbage, not a real prior dimensionality. This is not a
 * memory-safety hazard even so: whichever dispatch branch runs still nulls
 * the SAME `data` pointer (the offset-0 property above), which is the only
 * field #radiation_clean's later interpolate_*d_free() call can actually
 * dereference; the non-pointer tail fields left in whatever pre-call state
 * they were in are inert scalars, always overwritten by the next
 * #radiation_build_tables() call before anything ever reads them. Confirmed
 * by tracing every call site (radiation_table_io.c's radiation_read_data(),
 * this file's own radiation_restore(), and stellar_evolution.c's
 * stellar_evolution_props_init()), not assumed.
 *
 * @param rad The #radiation.
 */
void radiation_zero_pointers(struct radiation *rad) {

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

  if (was_2d) {
    interpolate_2d_zero_pointers(&rad->raw.luminosities_2d);
    interpolate_2d_zero_pointers(&rad->raw.dot_N_ion_2d);
    interpolate_2d_zero_pointers(&rad->raw.dot_E_excess_2d);
    interpolate_2d_zero_pointers(&rad->integrated.luminosities_2d);
    interpolate_2d_zero_pointers(&rad->integrated.dot_N_ion_2d);
    interpolate_2d_zero_pointers(&rad->integrated.dot_E_excess_2d);
  } else {
    interpolate_1d_zero_pointers(&rad->raw.luminosities);
    interpolate_1d_zero_pointers(&rad->raw.dot_N_ion);
    interpolate_1d_zero_pointers(&rad->raw.dot_E_excess);
    interpolate_1d_zero_pointers(&rad->integrated.luminosities);
    interpolate_1d_zero_pointers(&rad->integrated.dot_N_ion);
    interpolate_1d_zero_pointers(&rad->integrated.dot_E_excess);
  }

  interpolate_2d_zero_pointers(&rad->raw.main_sequence_lifetime_2d);
  interpolate_2d_zero_pointers(&rad->raw.main_sequence_lifetime_inverse_2d);
}
