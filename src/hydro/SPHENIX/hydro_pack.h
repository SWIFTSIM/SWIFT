/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_SPHENIX_HYDRO_PACK_H
#define SWIFT_SPHENIX_HYDRO_PACK_H

/**
 * @file SPHENIX/hydro_pack.h
 * @brief SPHENIX specific hydro_pack operations.
 */

#include "black_holes_struct.h"
#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "csds.h"
#include "feedback_struct.h"
#include "part.h"
#include "particle_splitting_struct.h"
#include "pressure_floor_struct.h"
#include "rt_struct.h"
#include "star_formation_struct.h"
#include "timestep_limiter_struct.h"
#include "tracers_struct.h"

struct xv_message_part {
  struct part p;
} SWIFT_STRUCT_ALIGN;

__attribute__((always_inline)) INLINE static void hydro_pack_xv(
    struct xv_message_part *mp, const struct part *p) {
  memcpy(&mp->p, p, sizeof(struct part));
}

__attribute__((always_inline)) INLINE static void hydro_unpack_xv(
    const struct xv_message_part *mp, struct part *p) {
  memcpy(p, &mp->p, sizeof(struct part));
}

struct rho_message_part {

  float h;

  float rho;

  float viscosity_alpha;

  float soundspeed;

  /*! Chemistry information */
  struct chemistry_part_data chemistry_data;

  /*! Cooling information */
  struct cooling_part_data cooling_data;

  /*! Additional data used by the feedback */
  struct feedback_part_data feedback_data;

  /*! Black holes information (e.g. swallowing ID) */
  struct black_holes_part_data black_holes_data;

  /*! Additional data used by the pressure floor */
  struct pressure_floor_part_data pressure_floor_data;

  /*! Additional Radiative Transfer Data */
  struct rt_part_data rt_data;
} SWIFT_STRUCT_ALIGN;

__attribute__((always_inline)) INLINE static void hydro_pack_rho(
    struct rho_message_part *mp, const struct part *p) {
  mp->h = p->h;
  mp->rho = p->rho;
  mp->viscosity_alpha = p->viscosity.alpha;
  mp->soundspeed = p->force.soundspeed;
  memcpy(&mp->chemistry_data, &p->chemistry_data,
         sizeof(struct chemistry_part_data));
  memcpy(&mp->cooling_data, &p->cooling_data, sizeof(struct cooling_part_data));
  memcpy(&mp->feedback_data, &p->feedback_data,
         sizeof(struct feedback_part_data));
  memcpy(&mp->black_holes_data, &p->black_holes_data,
         sizeof(struct black_holes_part_data));
  memcpy(&mp->pressure_floor_data, &p->pressure_floor_data,
         sizeof(struct pressure_floor_part_data));
  memcpy(&mp->rt_data, &p->rt_data, sizeof(struct rt_part_data));
}
__attribute__((always_inline)) INLINE static void hydro_unpack_rho(
    const struct rho_message_part *mp, struct part *p) {
  p->h = mp->h;
  p->rho = mp->rho;
  p->viscosity.alpha = mp->viscosity_alpha;
  p->force.soundspeed = mp->soundspeed;
  memcpy(&p->chemistry_data, &mp->chemistry_data,
         sizeof(struct chemistry_part_data));
  memcpy(&p->cooling_data, &mp->cooling_data, sizeof(struct cooling_part_data));
  memcpy(&p->feedback_data, &mp->feedback_data,
         sizeof(struct feedback_part_data));
  memcpy(&p->black_holes_data, &mp->black_holes_data,
         sizeof(struct black_holes_part_data));
  memcpy(&p->pressure_floor_data, &mp->pressure_floor_data,
         sizeof(struct pressure_floor_part_data));
  memcpy(&p->rt_data, &mp->rt_data, sizeof(struct rt_part_data));
}

/**
 * @brief Reduced particle sent before the force loop.
 */
struct gradient_message_part {

  /* Viscosity alpha. */
  float viscosity_alpha;

  /* Diffusion alpha. */
  float diffusion_alpha;

  /* "Grad h" term */
  float force_f;

  /*! Particle pressure. */
  float pressure;

  /*! Balsara switch */
  float balsara;

  /*! Chemistry information */
  struct chemistry_part_data chemistry_data;

  /*! Cooling information */
  struct cooling_part_data cooling_data;

  /*! Additional data used by the feedback */
  struct feedback_part_data feedback_data;

  /*! Black holes information (e.g. swallowing ID) */
  struct black_holes_part_data black_holes_data;

  /*! Additional data used by the pressure floor */
  struct pressure_floor_part_data pressure_floor_data;

  /*! Additional Radiative Transfer Data */
  struct rt_part_data rt_data;

} SWIFT_STRUCT_ALIGN;

__attribute__((always_inline)) INLINE static void hydro_pack_gradient(
    struct gradient_message_part *mp, const struct part *p) {
  mp->viscosity_alpha = p->viscosity.alpha;
  mp->diffusion_alpha = p->diffusion.alpha;
  mp->force_f = p->force.f;
  mp->pressure = p->force.pressure;
  mp->balsara = p->force.balsara;
  memcpy(&mp->chemistry_data, &p->chemistry_data,
         sizeof(struct chemistry_part_data));
  memcpy(&mp->cooling_data, &p->cooling_data, sizeof(struct cooling_part_data));
  memcpy(&mp->feedback_data, &p->feedback_data,
         sizeof(struct feedback_part_data));
  memcpy(&mp->black_holes_data, &p->black_holes_data,
         sizeof(struct black_holes_part_data));
  memcpy(&mp->pressure_floor_data, &p->pressure_floor_data,
         sizeof(struct pressure_floor_part_data));
  memcpy(&mp->rt_data, &p->rt_data, sizeof(struct rt_part_data));
}

__attribute__((always_inline)) INLINE static void hydro_unpack_gradient(
    const struct gradient_message_part *mp, struct part *p) {
  p->viscosity.alpha = mp->viscosity_alpha;
  p->diffusion.alpha = mp->diffusion_alpha;
  p->force.f = mp->force_f;
  p->force.pressure = mp->pressure;
  p->force.balsara = mp->balsara;
  memcpy(&p->chemistry_data, &mp->chemistry_data,
         sizeof(struct chemistry_part_data));
  memcpy(&p->cooling_data, &mp->cooling_data, sizeof(struct cooling_part_data));
  memcpy(&p->feedback_data, &mp->feedback_data,
         sizeof(struct feedback_part_data));
  memcpy(&p->black_holes_data, &mp->black_holes_data,
         sizeof(struct black_holes_part_data));
  memcpy(&p->pressure_floor_data, &mp->pressure_floor_data,
         sizeof(struct pressure_floor_part_data));
  memcpy(&p->rt_data, &mp->rt_data, sizeof(struct rt_part_data));
}

#endif /* SWIFT_SPHENIX_HYDRO_PACK_H */
