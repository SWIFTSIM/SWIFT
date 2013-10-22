/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

/* MPI headers. */
#ifdef WITH_MPI
    #include <mpi.h>
#endif

/* Local headers. */
#include "error.h"
#include "const.h"
#include "cycle.h"
#include "atomic.h"
#include "lock.h"
#include "space.h"
#include "part.h"
#include "multipole.h"
#include "cell.h"


/**
 * @brief Merge two multipoles.
 *
 * @param ma The #multipole which will contain the merged result.
 * @param mb The other #multipole.
 */
 
void multipole_merge ( struct multipole *ma , struct multipole *mb ) {

    #if multipole_order == 1
    
        /* Correct the position. */
        float mma = ma->coeffs[0], mmb = mb->coeffs[0];
        float w = 1.0f / ( mma + mmb );
        for ( int k = 0 ; k < 3 ; k++ )
            ma->x[k] = ( ma->x[k]*mma + mb->x[k]*mmb ) * w;
            
        /* Add the particle to the moments. */
        ma->coeffs[0] = mma + mmb;
    
    #else
        #error( "Multipoles of order %i not yet implemented." , multipole_order )
    #endif

    }


/**
 * @brief Add a particle to the given multipole.
 *
 * @param m The #multipole.
 * @param p The #gpart.
 */
 
void multipole_addpart ( struct multipole *m , struct gpart *p ) {
    
    #if multipole_order == 1

        /* Correct the position. */
        float mm = m->coeffs[0], mp = p->mass;
        float w = 1.0f / ( mm + mp );
        for ( int k = 0 ; k < 3 ; k++ )
            m->x[k] = ( m->x[k]*mm + p->x[k]*mp ) * w;
            
        /* Add the particle to the moments. */
        m->coeffs[0] = mm + mp;
        
    #else
        #error( "Multipoles of order %i not yet implemented." , multipole_order )
    #endif

    }


/**
 * @brief Add a group of particles to the given multipole.
 *
 * @param m The #multipole.
 * @param p The #gpart array.
 * @param N Number of parts to add.
 */
 
void multipole_addparts ( struct multipole *m , struct gpart *p , int N ) {
    
    #if multipole_order == 1
    
        /* Get the combined mass and positions. */
        double xp[3] = { 0.0 , 0.0 , 0.0 };
        float mp = 0.0f, w;
        for ( int k = 0 ; k < N ; k++ ) {
            w = p[k].mass;
            mp += w;
            xp[0] += p[k].x[0] * w;
            xp[1] += p[k].x[1] * w;
            xp[2] += p[k].x[2] * w;
            }

        /* Correct the position. */
        float mm = m->coeffs[0];
        w = 1.0f / ( mm + mp );
        for ( int k = 0 ; k < 3 ; k++ )
            m->x[k] = ( m->x[k]*mm + xp[k] ) * w;
            
        /* Add the particle to the moments. */
        m->coeffs[0] = mm + mp;
        
    #else
        #error( "Multipoles of order %i not yet implemented." , multipole_order )
    #endif

    }


/**
 * @brief Init a multipole from a set of particles.
 *
 * @param m The #multipole.
 * @param parts The #gparts.
 * @param N The number of particles.
 */
 
void multipole_init ( struct multipole *m , struct gpart *parts , int N ) {
    
    #if multipole_order == 1

        float mass = 0.0f, w;
        double x[3] = { 0.0 , 0.0 , 0.0 };
        int k;
        
        /* Collect the particle data. */
        for ( k = 0 ; k < N ; k++ ) {
            w = parts[k].mass;
            mass += w;
            x[0] += parts[k].x[0] * w;
            x[1] += parts[k].x[1] * w;
            x[2] += parts[k].x[2] * w;
            }
            
        /* Store the data on the multipole. */
        m->coeffs[0] = mass;
        m->x[0] = x[0] / mass;
        m->x[1] = x[1] / mass;
        m->x[2] = x[2] / mass;
        
    #else
        #error( "Multipoles of order %i not yet implemented." , multipole_order )
    #endif

    }


/**
 * @brief Reset the data of a #multipole.
 *
 * @param m The #multipole.
 */
 
void multipole_reset ( struct multipole *m ) {

    /* Just bzero the struct. */
    bzero( m , sizeof(struct multipole) );
    
    }
