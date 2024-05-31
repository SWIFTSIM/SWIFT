/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 * Copyright (c) 2024 Nikyta Shchutksyi (shchutskyi@lorentz.leidenuniv.nl)
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
#ifndef SWIFT_FORCING_ABC_FLOW_H
#define SWIFT_FORCING_ABC_FLOW_H

/* Config parameters. */
#include <config.h>

/* Standard includes. */
#include <float.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"

/**
 * @brief Forcing Term Properties
 */
struct forcing_terms {

  /*! Reference velocity (internal units) */
  float u0;

  /*! Velocity scaling along the z direction */
  float Vz_factor;

  /*! Wavenumber of the flow */
  float kv;

  /*! ABC flow coefficients */
  float A, B, C;
};

/**
 * @brief Computes the forcing terms.
 *
 * Based on David Galloway (2012) ABC flows then and now, Geophysical &
 * Astrophysical Fluid Dynamics, 106:4-5, 450-467 This version differs from the
 * paper by imposing normalized velocity
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

  const double L = s->dim[0];
  const float u0 = terms->u0;
  const float A = terms->A;
  const float B = terms->B;
  const float C = terms->C;
  const float Norm = 1 / sqrt(A * A + B * B + C * C);
  const float Vz_factor = terms->Vz_factor;
  const double k0 = (2. * M_PI / L) * terms->kv;
  double v_ABC[3];

  /* Eq. 2 of David Galloway (2012) ABC flows then and now, Geophysical &
   * Astrophysical Fluid Dynamics, 106:4-5, 450-467 */
  // Velocity normalized such that <v>rms = u0
  v_ABC[0] = u0 * Norm * (A * sin(k0 * p->x[2]) + C * cos(k0 * p->x[1]));
  v_ABC[1] = u0 * Norm * (B * sin(k0 * p->x[0]) + A * cos(k0 * p->x[2]));
  v_ABC[2] = u0 * Norm * (C * sin(k0 * p->x[1]) + B * cos(k0 * p->x[0]));

  /* Force the velocity and possibly scale the z-direction */
  xp->v_full[0] = v_ABC[0];
  xp->v_full[1] = v_ABC[1];
  xp->v_full[2] = v_ABC[2] * Vz_factor;
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

  message("Forcing terms is 'ABC flow'. U0: %.5f / Vz factor: %.5f.", terms->u0,
          terms->Vz_factor);
  message("Run using ABC parameters: A = %.5f, B = %.5f, C = %.5f", terms->A,
          terms->B, terms->C);
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

  terms->u0 = parser_get_param_double(parameter_file, "ABC_Flow_Forcing:u0");
  terms->Vz_factor = parser_get_opt_param_float(
      parameter_file, "ABC_Flow_Forcing:Vz_factor", 1.f);
  terms->kv = parser_get_param_double(parameter_file, "ABC_Flow_Forcing:kv");

  terms->A = parser_get_param_double(parameter_file, "ABC_Flow_Forcing:A");
  terms->B = parser_get_param_double(parameter_file, "ABC_Flow_Forcing:B");
  terms->C = parser_get_param_double(parameter_file, "ABC_Flow_Forcing:C");
}

#endif /* SWIFT_FORCING_ABC_FLOW_H */
