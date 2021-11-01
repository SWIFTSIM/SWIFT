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
#ifndef SWIFT_NONE_NEUTRINO_H
#define SWIFT_NONE_NEUTRINO_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "../../engine.h"
#include "neutrino_response.h"

struct neutrino_model {};

INLINE static void gather_neutrino_consts(const struct space *s,
                                          struct neutrino_model *nm) {}

INLINE static void gpart_neutrino_weight_mesh_only(
    const struct gpart *gp, const struct neutrino_model *nm, double *weight) {}

INLINE static void gpart_neutrino_mass_weight(const struct gpart *gp,
                                              const struct neutrino_model *nm,
                                              double *mass, double *weight) {
  *mass = gp->mass;
  *weight = 1.0;
}

INLINE static double neutrino_mass_factor(
    const struct cosmology *cosmo, const struct unit_system *internal_units,
    const struct phys_const *physical_constants, double volume,
    double nr_nuparts) {

  return 1.0;
}

INLINE static float relativistic_drift_factor(float v[3], float a, float c) {
  return 1.0;
}

INLINE static void neutrino_check_cosmology(
    const struct space *s, const struct cosmology *cosmo,
    const struct phys_const *physical_constants, struct swift_params *params,
    const struct neutrino_props *neutrino_props, const int rank,
    const int verbose) {

  if (cosmo->Omega_nu_0 > 0.)
    error("Not compiled with neutrinos. Configure with with-neutrinos flag.");
}

__attribute__((always_inline)) INLINE static void gravity_first_init_neutrino(
    struct gpart *gp, const struct engine *e) {}

#endif /* SWIFT_NONE_NEUTRINO_H */
