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
#ifndef SWIFT_SN_FEEDBACK_H
#define SWIFT_SN_FEEDBACK_H
#include <float.h>
#include "equation_of_state.h"
#include "hydro.h"
// void do_supernova_feedback(const struct source_terms* source_terms, const struct phys_const* const constants, const struct part* p);


/* determine whether location is in cell */
__attribute__((always_inline)) INLINE static int is_in_cell(const double cell_min[], const double cell_width[], const double location[], const int dimen)
{
  for(int i=0; i<dimen; i++){
	 if (cell_min[i] > location[i])
		  return 0;
	 if ( (cell_min[i] + cell_width[i]) <= location[i])
		return 0;
  }
  return 1;
};

__attribute__((always_inline)) INLINE static void do_supernova_feedback(const struct sourceterms* sourceterms, struct part* p)
{
  const float density     = p->rho;
  const float u_old       = hydro_get_internal_energy(p,0);
  const float u_new       = u_old + sourceterms->supernova.energy / p->mass;
  const float entropy_old = gas_entropy_from_internal_energy(density, u_old);
  const float entropy_new = calculate_entropy(entropy_old, density, u_old, u_new);
  hydro_set_internal_energy(p, u_new);
  hydro_set_entropy(p,entropy_new);
  message(" u_old= %e u_new= %e S_old= %e S_new= %e rho= %e", u_old, u_new, entropy_old, entropy_new, density);
  message(" p->s %e", p->entropy);
};
#endif /* SWIFT_SN_FEEDBACK_H */


  


