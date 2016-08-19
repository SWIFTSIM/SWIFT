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
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* Cooling Properties */
struct cooling_data {

#ifdef CONST_COOLING
  struct {
    float lambda;
    float min_energy;
    float cooling_tstep_mult;
  } const_cooling;
#endif

#ifdef CREASEY_COOLING
  struct {
    float lambda;
    float min_temperature;
    float hydrogen_mass_abundance;
    float mean_molecular_weight;
    float min_internal_energy;
    float min_internal_energy_cgs;
    float cooling_tstep_mult;
  } creasey_cooling;
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
__attribute__((always_inline)) INLINE static double
cooling_timestep(const struct cooling_data* cooling,
		 const struct phys_const* const phys_const,
		 const struct part* const p) {

  const double cooling_rate = cooling->const_cooling.lambda;
  const double internal_energy = hydro_get_internal_energy(p,0);// dt = 0 because using current entropy
  return  cooling->const_cooling.cooling_tstep_mult * internal_energy / cooling_rate;
}

#endif /* CONST_COOLING */


/* Now, some generic functions, defined in the source file */
void cooling_init(const struct swift_params* parameter_file,
                    struct UnitSystem* us,
		  const struct phys_const* const phys_const,
                    struct cooling_data* cooling);

void cooling_print(const struct cooling_data* cooling);
void update_entropy(const struct phys_const* const phys_const, const struct UnitSystem* us,
		    const struct cooling_data* cooling, struct part* p, float dt);
float calculate_new_thermal_energy(float u_old, float rho, float dt, 
				   const struct cooling_data* cooling,
				   const struct phys_const* const phys_const,
				   const struct UnitSystem* us);
#endif /* SWIFT_COOLING_H */
