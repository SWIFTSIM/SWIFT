/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2013 Matthieu Schaller (matthieu.schaller@durham.ac.uk),
 *                    Pedro Gonnet (pedro.gonnet@durham.ac.uk).
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


#include <stdio.h>

#include "part.h"


/**
 * @brief Looks for the particle with the given id and prints its information to the standard output.
 * 
 * @param parts The array of particles.
 * @param id The id too look for.
 * @param N The size of the array of particles.
 *
 * (Should be used for debugging only as it runs in O(N).)
 */
void printParticle ( struct part *parts , long long int id, int N ) {

    int i;

    /* Look for the particle. */
    for ( i = 0 ; i < N && parts[i].id != id; i++ );

    if ( i < N )
        printf("## Particle[%d]: id=%lld, x=[%e,%e,%e], v=[%.3e,%.3e,%.3e], a=[%.3e,%.3e,%.3e], h=%.3e, h_dt=%.3e, wcount=%.3e, m=%.3e, rho=%.3e, rho_dh=%.3e, div_v=%.3e, u=%.3e, dudt=%.3e, bals=%.3e, POrho2=%.3e, v_sig=%.3e, dt=%.3e\n",
            i,
            parts[i].id,
            parts[i].x[0], parts[i].x[1], parts[i].x[2],
            parts[i].v[0], parts[i].v[1], parts[i].v[2],
            parts[i].a[0], parts[i].a[1], parts[i].a[2],
            parts[i].h,
            parts[i].force.h_dt,
            parts[i].wcount,
            parts[i].mass,
            parts[i].rho, parts[i].rho_dh,
            parts[i].density.div_v,
            parts[i].u,
            parts[i].force.u_dt,
            parts[i].force.balsara,
            parts[i].force.POrho2,
            parts[i].force.v_sig,
            parts[i].dt
            );
    else
        printf("## Particles[???] id=%lld not found\n", id);
    
    }

