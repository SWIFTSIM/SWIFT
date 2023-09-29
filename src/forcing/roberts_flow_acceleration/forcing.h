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

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"
#include "hydro.h"

/**
 * @brief Forcing Term Properties
 */

#define forcing_propos_Vz_factor 1.f

/**
 * @brief Constants of the model
 */
struct forcing_terms {

  /*! Reference velocity (internal units) */
  float u0;

  /*! 'viscosity' to convert from velocity to acceleration */
  float nu;

  /*! Velocity scaling along the z direction */
  float Vz_factor;
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
    const struct phys_const* restrict phys_const, struct part* p,
    struct xpart* xp) {

  const float v_sig = hydro_get_signal_velocity(p);
  const double L = s->dim[0];
  const float u0 = terms->u0;
  const float nu = terms->nu*p->viscosity.alpha*v_sig*p->h; //CFL COndition?
  const float Vz_factor = terms->Vz_factor;
  const double k0 = 2. * M_PI / L;
  const double kf = M_SQRT2 * k0;
  
  /* Eq. 8 */
  const double Psi = (u0 / k0) * cos(k0 * p->x[0]) * cos(k0 * p->x[1]);

  /* Eq. 7 (with optional scaling along z) */
  const double v_Rob[3] = {u0 * cos(k0 * p->x[0]) * sin(k0 * p->x[1]),
                           -u0 * sin(k0 * p->x[0]) * cos(k0 * p->x[1]),
                           kf * Psi * Vz_factor};

  /* Eq. 6 */
  const double f[3] = {nu * kf * kf * v_Rob[0], nu * kf * kf * v_Rob[1],
                       nu * kf * kf * v_Rob[2]};

  /* Update the accelerations */
  p->a_hydro[0] += f[0];
  p->a_hydro[1] += f[1];
  p->a_hydro[2] += f[2];
  
  if(time == 0.f){
  /* Force the velocity */
  xp->v_full[0] = v_Rob[0];
  xp->v_full[1] = v_Rob[1];
  xp->v_full[2] = v_Rob[2];}
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
static INLINE void forcing_terms_print(const struct forcing_terms* terms) {

  message(
      "Forcing terms is 'Roberts flow using accelerations'. U0: %.5f. nu: "
      "%.5f. Vz factor: %.5f.",
      terms->u0, terms->nu, terms->Vz_factor);
}

/**
 * @brief Initialises the forcing term properties
 *
 * Nothing to do here.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
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
      parameter_file, "RobertsFlowAccelerationForcing:Vz_factor",
      forcing_propos_Vz_factor);
}

#endif /* SWIFT_FORCING_ROBERTS_FLOW_ACCELERATION_H */