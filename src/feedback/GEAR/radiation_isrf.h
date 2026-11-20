/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_RADIATION_ISRF_GEAR_H
#define SWIFT_RADIATION_ISRF_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_isrf.h
 * @brief Receiver-side LW/PE dust extinction and per-step hooks of the packet
 * transport of the ISRF for GEAR: opacity, extinction and absorption.
 */

#include "feedback_struct.h"
#include "inline.h"

struct part;
struct xpart;
struct cosmology;
struct unit_system;
struct hydro_props;
struct engine;
struct cooling_function_data;
struct feedback_props;
struct phys_const;
struct radiation;
struct stellar_model;

void radiation_set_lw_photon_energy_cgs(const struct radiation *rad,
                                        const struct stellar_model *sm);

/**
 * @brief Set the band-edge weights of @p fb_props from the radiation table.
 *
 * @param fb_props (output) The #feedback_props to set.
 * @param rad The main stellar model's #radiation.
 * @param sm The main #stellar_model, for its IMF mass range.
 */
void radiation_set_band_edge_coefficients(struct feedback_props *fb_props,
                                          const struct radiation *rad,
                                          const struct stellar_model *sm);

struct swift_params;

#ifdef GEAR_SUBGRID_RADIATION_ISRF

void radiation_packet_build_octant_frame(const float angles_deg[3],
                                         const int seed_offset, float frame[9]);
int radiation_packet_read_test_beams(
    struct swift_params *params,
    float beams[RADIATION_PACKET_TEST_BEAMS]
               [RADIATION_PACKET_TEST_BEAM_FIELDS]);

/* The ISRF gradient pass runs in the extra hydro loop. Configure refuses a
   hydro scheme without it; this is the backstop. */
#include "part.h"
#ifndef EXTRA_HYDRO_LOOP
#error "The GEAR subgrid ISRF needs a hydro scheme with EXTRA_HYDRO_LOOP."
#endif

void radiation_first_init_part(struct part *restrict p);
void radiation_snapshot_part_propagation(struct part *p, const struct xpart *xp,
                                         const struct engine *e);
void radiation_end_gradient_propagation(struct part *p, const struct engine *e);
void radiation_end_force_propagation(struct part *p, const struct engine *e);
float radiation_isrf_part_timestep(const struct part *restrict p,
                                   const struct engine *e);

#else /* GEAR_SUBGRID_RADIATION_ISRF */

#include <float.h>

/* Without the part: no ISRF state, so every hook does nothing. */

__attribute__((always_inline)) INLINE static void
radiation_packet_build_octant_frame(const float angles_deg[3],
                                    const int seed_offset, float frame[9]) {}
__attribute__((always_inline)) INLINE static int
radiation_packet_read_test_beams(
    struct swift_params *params,
    float beams[RADIATION_PACKET_TEST_BEAMS]
               [RADIATION_PACKET_TEST_BEAM_FIELDS]) {
  return 0;
}

__attribute__((always_inline)) INLINE static void radiation_first_init_part(
    struct part *restrict p) {}
__attribute__((always_inline)) INLINE static void
radiation_snapshot_part_propagation(struct part *p, const struct xpart *xp,
                                    const struct engine *e) {}
__attribute__((always_inline)) INLINE static void
radiation_end_gradient_propagation(struct part *p, const struct engine *e) {}
__attribute__((always_inline)) INLINE static void
radiation_end_force_propagation(struct part *p, const struct engine *e) {}
__attribute__((always_inline)) INLINE static float radiation_isrf_part_timestep(
    const struct part *restrict p, const struct engine *e) {
  return FLT_MAX;
}

#endif /* GEAR_SUBGRID_RADIATION_ISRF */

#endif /* SWIFT_RADIATION_ISRF_GEAR_H */
