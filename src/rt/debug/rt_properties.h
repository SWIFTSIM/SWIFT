/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_PROPERTIES_DEBUG_H
#define SWIFT_RT_PROPERTIES_DEBUG_H

/**
 * @file src/rt/debug/rt_properties.h
 * @brief Main header file for the debug radiative transfer scheme properties.
 */

/**
 * @brief Properties of the debug radiative transfer model
 */
struct rt_props {
  /* Do extended tests where we assume that all parts
   * have spart neighbours? */
  int debug_do_all_parts_have_stars_checks;

  /* Are we running with hydro or star controlled injection?
   * This is added to avoid #ifdef macros as far as possible */
  int hydro_controlled_injection;

  /* radiation emitted by stars this step. This is not really a property,
   * but a placeholder to sum up a global variable */
  int debug_radiation_emitted_this_step;

  /* total radiation emitted by stars. This is not really a property,
   * but a placeholder to sum up a global variable */
  unsigned long long debug_radiation_emitted_tot;

  /* radiation absorbed by gas this step. This is not really a property,
   * but a placeholder to sum up a global variable */
  int debug_radiation_absorbed_this_step;

  /* total radiation absorbed by gas. This is not really a property,
   * but a placeholder to sum up a global variable */
  unsigned long long debug_radiation_absorbed_tot;
};

/**
 * @brief Print the RT model.
 *
 * @param rtp The #rt_props
 */
__attribute__((always_inline)) INLINE static void rt_props_print(
    const struct rt_props* rtp) {

  /* Only the master print */
  if (engine_rank != 0) return;

  message("Radiative transfer scheme: '%s'", RT_IMPLEMENTATION);

  /* Print the RT properties */
  if (rtp->debug_do_all_parts_have_stars_checks)
    message("Doing extra checks assuming all parts have spart neighbours");
  else
    message("Skipping extra checks assuming all parts have spart neighbours");
}

/**
 * @brief Initialize the global properties of the RT scheme.
 *
 * @param rtp The #rt_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void rt_props_init(
    struct rt_props* rtp, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    struct cosmology* cosmo) {

  rtp->debug_do_all_parts_have_stars_checks =
      parser_get_opt_param_int(params, "DebugRT:all_parts_have_stars", 0);

#ifdef RT_HYDRO_CONTROLLED_INJECTION
  rtp->hydro_controlled_injection = 1;
#else
  rtp->hydro_controlled_injection = 0;
#endif

  rtp->debug_radiation_emitted_tot = 0ULL;
  rtp->debug_radiation_absorbed_tot = 0ULL;

  /* After initialisation, print params to screen */
  rt_props_print(rtp);

  /* Print a final message. */
  if (engine_rank == 0) {
    message("Radiative transfer initialized");
  }
}

/**
 * @brief Write an RT properties struct to the given FILE as a
 * stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
__attribute__((always_inline)) INLINE static void rt_struct_dump(
    const struct rt_props* props, FILE* stream) {

  restart_write_blocks((void*)props, sizeof(struct rt_props), 1, stream,
                       "RT props", "RT properties struct");
}

/**
 * @brief Restore an RT properties struct from the given FILE as
 * a stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
__attribute__((always_inline)) INLINE static void rt_struct_restore(
    struct rt_props* props, FILE* stream) {

  restart_read_blocks((void*)props, sizeof(struct rt_props), 1, stream, NULL,
                      "RT properties struct");
}

#endif /* SWIFT_RT_PROPERTIES_DEBUG_H */
