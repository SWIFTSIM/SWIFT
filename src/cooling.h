/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Richard Bower (r.g.bower@durham.ac.uk)
 *                    Stefan Arridge  (stefan.arridge@durham.ac.uk)
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

#ifndef SWIFT_COOLING_H
#define SWIFT_COOLING_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local includes. */
#include "const.h"
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* Cooling Properties */
struct cooling_data {

#ifdef CONST_COOLING
  struct {
    double lambda;
    double min_energy;
    double cooling_tstep_mult;
  } const_cooling;
#endif
};

/* Include Cooling */
#ifdef CONST_COOLING

/**
 * @brief Computes the time-step due to cooling
 *
 * @param cooling The #cooling_data used in the run.
 * @param phys_const The physical constants in internal units.
 * @param  Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
cooling_timestep(const struct cooling_data* cooling,
		 const struct phys_const* const phys_const, 
		 const struct part* const p) {

  const float cooling_rate = get_cooling_rate( p->density, p->internal_energy, cooling );
  return  cooling->const_cooling.cooling_tstep_mult * p->internal_energy / cooling_rate;
}

/* /\** */
/*  * @brief Updates the internal energy of a particle due to cooling. */
/*  * */
/*  * @param cooling The #cooling_data used in the run. */
/*  * @param phys_const The physical constants in internal units. */
/*  * @param p Pointer to the particle data. */
/*  *\/ */
/* __attribute__((always_inline)) INLINE static float */
/* cooling_update_entropy(const struct cooling_data* cooling, */
/*                                      const struct phys_const* const phys_const, */
/* 										       const dt, */
/*                                      struct part* p) { */
/*   const float old_entropy = p->Entropy */
/*   const float cooling_rate = get_cooling_rate( p->density, p->internal_energy, cooling ); */
/*   // do other quanitities need to be updated as well?? */
/*   p->internal_energy -= dt * cooling_rate; */
/* } */
/* #endif /\* CONST_COOLING *\/ */


/* Now, some generic functions, defined in the source file */
void cooling_init(const struct swift_params* parameter_file,
                    struct UnitSystem* us,
                    struct cooling_data* cooling);

void cooling_print(const struct cooling_data* cooling);
float calculate_new_thermal_energy(float u_old, float dt, const struct cooling_data* cooling);
void update_entropy(const struct cooling_data* cooling,
		   const struct phys_const* const phys_const, struct part* p, 
		    double dt);
#endif /* SWIFT_COOLING_H */
