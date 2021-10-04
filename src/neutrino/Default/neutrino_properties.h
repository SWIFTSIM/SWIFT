/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Willem Elbers (willem.h.elbers@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_NEUTRINO_PROPERTIES_H
#define SWIFT_DEFAULT_NEUTRINO_PROPERTIES_H

#include "../../restart.h"

/**
 * @brief Properties of neutrinos
 */
struct neutrino_props {

  /* Whether to run with the delta-f method for neutrino weighting */
  char use_delta_f;

  /* Whether to generate random neutrino velocities in the initial conditions */
  char generate_ics;

  /* Random seed for the neutrino weighting task */
  long long neutrino_seed;

  /* Whether to use the linear respose method */
  char use_linear_response;
};

/**
 * @brief Initialise the neutrino properties from the parameter file.
 *
 * @param np The #neutrino_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
INLINE static void neutrino_props_init(struct neutrino_props *np,
                                       const struct phys_const *phys_const,
                                       const struct unit_system *us,
                                       struct swift_params *params,
                                       const struct cosmology *cosmo) {

  np->use_delta_f = parser_get_param_int(params, "Neutrino:use_delta_f");
  np->generate_ics = parser_get_param_int(params, "Neutrino:generate_ics");
  np->neutrino_seed =
      parser_get_opt_param_longlong(params, "Neutrino:neutrino_seed", 0);
  np->use_linear_response =
      parser_get_opt_param_int(params, "Neutrino:use_linear_response", 0);
}

/**
 * @brief Write a neutrino_props struct to the given FILE as a stream of
 * bytes.
 *
 * @param props the neutrino properties struct
 * @param stream the file stream
 */
INLINE static void neutrino_struct_dump(const struct neutrino_props *props,
                                        FILE *stream) {

  restart_write_blocks((void *)props, sizeof(struct neutrino_props), 1, stream,
                       "neutrino props", "Neutrino props");
}

/**
 * @brief Restore a neutrino_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the neutrino properties struct
 * @param stream the file stream
 */
INLINE static void neutrino_struct_restore(const struct neutrino_props *props,
                                           FILE *stream) {
  restart_read_blocks((void *)props, sizeof(struct neutrino_props), 1, stream,
                      NULL, "Neutrino props");
}

#endif /* SWIFT_DEFAULT_NEUTRINO_PROPERTIES_H */
