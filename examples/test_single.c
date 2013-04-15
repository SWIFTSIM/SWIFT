/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <fenv.h>
#include <omp.h>

/* Conditional headers. */
#ifdef HAVE_LIBZ
    #include <zlib.h>
#endif

/* Local headers. */
#include "swift.h"

/* Ticks per second on this machine. */
#ifndef CPU_TPS
    #define CPU_TPS 2.67e9
#endif

/* Engine policy flags. */
#ifndef ENGINE_POLICY
    #define ENGINE_POLICY engine_policy_none
#endif

/* Error macro. */
#define error(s) { printf( "%s:%s:%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }


/**
 * @brief Main routine that loads a few particles and generates some output.
 *
 */
 
int main ( int argc , char *argv[] ) {

    int k, N = 100;
    struct part p1, p2;
    float r2, dx[3] = { 0.0f , 0.0f , 0.0f };
    
    /* Init the particles. */
    for ( k = 0 ; k < 3 ; k++ ) {
        p1.a[k] = 0.0f; p1.v[k] = 0.0f; p1.x[k] = 0.0;
        p2.a[k] = 0.0f; p2.v[k] = 0.0f; p2.x[k] = 0.0;
        }
    p1.rho = 1.0f; p1.mass = 9.7059e-4; p1.h = 0.222871287;
    p2.rho = 1.0f; p2.mass = 9.7059e-4; p2.h = 0.222871287;
    p1.force.c = 0.0f; p1.force.balsara = 0.0f;
    p2.force.c = 0.0f; p2.force.balsara = 0.0f;
    p1.u = 1.e-5 / ((const_gamma - 1.)*p1.rho);
    p2.u = 1.e-5 / ((const_gamma - 1.)*p2.rho) + 100.0f / ( 33 * p2.mass );
    p1.force.POrho2 = p1.u * ( const_gamma - 1.0f ) / p1.rho;
    p2.force.POrho2 = p2.u * ( const_gamma - 1.0f ) / p2.rho;
    
    /* Dump a header. */
    printf( "# r a_1 udt_1 a_2 udt_2\n" );
    
    /* Loop over the different radii. */
    for ( k = 1 ; k <= N ; k++ ) {
    
        /* Set the distance/radius. */
        dx[0] = -((float)k)/N * fmaxf( p1.h , p2.h );
        r2 = dx[0]*dx[0];
        
        /* Clear the particle fields. */
        p1.a[0] = 0.0f; p1.force.u_dt = 0.0f;
        p2.a[0] = 0.0f; p2.force.u_dt = 0.0f;
        
        /* Interact the particles. */
        runner_iact_force( r2 , dx , p1.h , p2.h , &p1 , &p2 );
        
        /* Output the results. */
        printf( "%.3e %.3e %.3e %.3e %.3e\n" ,
            -dx[0] , p1.a[0] , p1.force.u_dt , p2.a[0] , p2.force.u_dt );
    
        } /* loop over radii. */
    
    /* All is calm, all is bright. */
    return 0;
    
    }
