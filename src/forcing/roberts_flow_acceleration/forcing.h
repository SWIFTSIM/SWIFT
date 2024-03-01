/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FORCING_ROBERTS_FLOW_ACCELERATION_H
#define SWIFT_FORCING_ROBERTS_FLOW_ACCELERATION_H

/* Config parameters. */
#include <config.h>

/* Standard includes. */
#include <float.h>

/* Local includes. */
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/* Type of flow */
enum flow {
  Brandenburg_flow,
  Roberts_flow_1,
  Roberts_flow_2,
  Roberts_flow_3,
  Roberts_flow_4
};

/**
 * @brief Forcing Term Properties
 */
struct forcing_terms {

  /*! Reference velocity (internal units) */
  float u0;

  /*! 'viscosity' to convert from velocity to acceleration (internal units) */
  float nu;

  /*! Velocity scaling along the z direction */
  float Vz_factor;

  /*! Kind of RobertsFlow */
  enum flow Flow_kind;

  /*! Wavenumber of the flow*/
  float kv;
};

/**
 * @brief Computes the forcing terms.
 *
 * Based on Tilgner & Brandenburg, 2008, MNRAS, 391, 1477
 *
 * @param time The current time.
 * @param terms The properties of the forcing terms.
 * @param s The #space we act on.
 * @param phys_const The physical constants in internal units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_terms_apply(
    const double time, const struct forcing_terms* terms, const struct space* s,
    const struct phys_const* phys_const, struct part* p, struct xpart* xp) {

  enum flow Flow_kind = terms->Flow_kind;
  const double L = s->dim[0];
  const float u0 = terms->u0;
  const float c_s = hydro_get_comoving_soundspeed(p);

  /* Effective viscosity from artificial viscosity, as in eq. 100 from
   * arXiv:1012.1885 */
  const float nu = terms->nu * p->viscosity.alpha * c_s * p->h /
                   (2.f * (hydro_dimension + 2.f));

  const float Vz_factor = terms->Vz_factor;
  const double k0 = (2. * M_PI / L) * terms->kv;
  const double kf = M_SQRT2 * k0;
  double v_Rob[3];
  double Psi;

  switch (Flow_kind) {

    case Brandenburg_flow:
      /* Eq. 8 of Tilgner & Brandenburg, 2008, MNRAS, 391, 1477 */
      // Psi = (u0 / k0) * cos(k0 * p->x[0]) * cos(k0 * p->x[1]);

      /* Eq. 7 of Tilgner & Brandenburg, 2008, MNRAS, 391, 1477 */
      // v_Rob[0] = u0 * cos(k0 * p->x[0]) * sin(k0 * p->x[1]);
      // v_Rob[1] = -u0 * sin(k0 * p->x[0]) * cos(k0 * p->x[1]);
      // v_Rob[2] = kf * Psi;

      // Velocity used to compare with A.B. runs (from overleaf)
      Psi = (u0 / k0) * sin(k0 * p->x[0]) * sin(k0 * p->x[1]);
      v_Rob[0] = u0 * sin(k0 * p->x[0]) * cos(k0 * p->x[1]);
      v_Rob[1] = -u0 * cos(k0 * p->x[0]) * sin(k0 * p->x[1]);
      v_Rob[2] = kf * Psi;

      break;

    case Roberts_flow_1:
      /* Eq. 5.1 of Roberts, Feb. 3, 1972, Vol. 271, No. 1216 (Feb. 3, 1972),
       * pp. 411-454.*/
      v_Rob[0] = u0 * sin(k0 * p->x[0]);
      v_Rob[1] = u0 * sin(k0 * p->x[1]);
      v_Rob[2] = u0 * (cos(k0 * p->x[0]) - cos(k0 * p->x[1]));
      break;

    case Roberts_flow_2:
      /* Eq. 6.1 of Roberts, Feb. 3, 1972, Vol. 271, No. 1216 (Feb. 3, 1972),
       * pp. 411-454.*/
      v_Rob[0] = u0 * sin(k0 * p->x[0]);
      v_Rob[1] = u0 * sin(k0 * p->x[1]);
      v_Rob[2] = u0 * (cos(k0 * p->x[0]) + cos(k0 * p->x[1]));
      break;

    case Roberts_flow_3:
      /* Eq. 6.2 of Roberts, Feb. 3, 1972, Vol. 271, No. 1216 (Feb. 3, 1972),
       * pp. 411-454.*/
      v_Rob[0] = u0 * sin(k0 * p->x[0]);
      v_Rob[1] = u0 * sin(k0 * p->x[1]);
      v_Rob[2] = u0 * 2 * cos(k0 * p->x[0]) * cos(k0 * p->x[1]);
      break;

    case Roberts_flow_4:
      /* Eq. 6.3 of Roberts, Feb. 3, 1972, Vol. 271, No. 1216 (Feb. 3, 1972),
       * pp. 411-454.*/
      v_Rob[0] = u0 * sin(k0 * p->x[0]);
      v_Rob[1] = u0 * sin(k0 * p->x[1]);
      v_Rob[2] = u0 * sin(k0 * (p->x[0] + p->x[1]));
      break;

    default:

      v_Rob[0] = 0.f;
      v_Rob[1] = 0.f;
      v_Rob[2] = 0.f;
  }

  /* Eq. 6 of Tilgner & Brandenburg, 2008, MNRAS, 391, 1477 */
  const double f[3] = {nu * kf * kf * v_Rob[0], nu * kf * kf * v_Rob[1],
                       nu * kf * kf * v_Rob[2]};

  /* Update the accelerations */
  p->a_hydro[0] += f[0];
  p->a_hydro[1] += f[1];
  p->a_hydro[2] += f[2] * Vz_factor;
}

/**
 * @brief Computes the time-step condition due to the forcing terms.
 *
 * Nothing to do here. --> Return FLT_MAX.
 *
 * @param time The current time.
 * @param terms The properties of the forcing terms.
 * @param phys_const The physical constants in internal units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static float forcing_terms_timestep(
    double time, const struct forcing_terms* terms,
    const struct phys_const* phys_const, const struct part* p,
    const struct xpart* xp) {

  return FLT_MAX;
}

/**
 * @brief Prints the properties of the forcing terms to stdout.
 *
 * @param terms The #forcing_terms properties of the run.
 */
static INLINE void forcing_terms_print(const struct forcing_terms* terms) {

  message(
      "Forcing terms is 'Roberts flow using accelerations'. U0: %.5f. nu: "
      "%.5f. Vz factor: %.5f.",
      terms->u0, terms->nu, terms->Vz_factor);
  message("Forcing 'Roberts flow' Kind: %i .", terms->Flow_kind);
}

/**
 * @brief Initialises the forcing term properties
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param s The #space object.
 * @param terms The forcing term properties to initialize
 */
static INLINE void forcing_terms_init(struct swift_params* parameter_file,
                                      const struct phys_const* phys_const,
                                      const struct unit_system* us,
                                      const struct space* s,
                                      struct forcing_terms* terms) {

  terms->u0 = parser_get_param_double(parameter_file,
                                      "RobertsFlowAccelerationForcing:u0");
  terms->nu = parser_get_param_double(parameter_file,
                                      "RobertsFlowAccelerationForcing:nu");
  terms->Vz_factor = parser_get_opt_param_float(
      parameter_file, "RobertsFlowAccelerationForcing:Vz_factor", 1.f);
  terms->Flow_kind =
      parser_get_param_int(parameter_file, "RobertsFlowForcing:Flow_kind");
  terms->kv = parser_get_param_double(parameter_file, "RobertsFlowForcing:kv");

  if (terms->Flow_kind > 4 || terms->Flow_kind < 0)
    error(
        "Error: Flow_kind variable can take integer values from [0,4] "
        "interval.");
}

#endif /* SWIFT_FORCING_ROBERTS_FLOW_ACCELERATION_H */
