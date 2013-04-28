/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef KERNEL_H
#define KERNEL_H 


/**
 * @file kernel.h
 * @brief SPH kernel functions. Compute W(x,h) and the gradient of W(x,h).
 */


#include "vector.h"

/* Coefficients for the kernel. */ 
#define kernel_degree 3
#define kernel_ivals 2
#define kernel_gamma 2.0f
#define kernel_gamma2 4.0f
#define kernel_gamma3 8.0f
#define kernel_igamma 0.5f
#define kernel_nwneigh 4.0/3.0*M_PI*kernel_gamma3*const_eta_kernel*const_eta_kernel*const_eta_kernel
static float kernel_coeffs[ (kernel_degree + 1) * (kernel_ivals + 1) ] __attribute__ ((aligned (16))) =
    { 3.0/4.0*M_1_PI , -3.0/2.0*M_1_PI , 0.0 , M_1_PI , 
      -0.25*M_1_PI , 3.0/2.0*M_1_PI , -3.0*M_1_PI , M_2_PI , 
      0.0 , 0.0 , 0.0 , 0.0 };
#define kernel_root ( kernel_coeffs[ kernel_degree ] )
#define kernel_wroot ( 4.0/3.0*M_PI*kernel_coeffs[ kernel_degree ] )
      
      
/**
 * @brief Computes the kernel and its derivative for a given distance x. Gives a sensible answer only if x<1.
 */

__attribute__ ((always_inline)) INLINE static void kernel_deval ( float x , float *W , float *dW_dx ) {
    int ind = fmin( x , kernel_ivals );
    float *coeffs = &kernel_coeffs[ ind*(kernel_degree + 1) ];
    float w = coeffs[0]*x + coeffs[1];
    float dw_dx = coeffs[0];
    for ( int k = 2 ; k <= kernel_degree ; k++ ) {
        dw_dx = dw_dx*x + w;
        w = x*w + coeffs[k];
        }
    *W = w;
    *dW_dx = dw_dx;
    }


#ifdef VECTORIZE

/**
 * @brief Computes the kernel and its derivative for a given distance x (Vectorized version). Gives a sensible answer only if x<1.
 */

__attribute__ ((always_inline)) INLINE static void kernel_deval_vec ( vector *x , vector *w , vector *dw_dx ) {
    
    vector ind, c[kernel_degree+1];
    int j, k;
    
    /* Load x and get the interval id. */
    ind.m = vec_ftoi( vec_fmin( x->v , vec_set1( (float)kernel_ivals ) ) );
    
    /* load the coefficients. */
    for ( k = 0 ; k < VEC_SIZE ; k++ )
        for ( j = 0 ; j < kernel_degree+1 ; j++ )
            c[j].f[k] = kernel_coeffs[ ind.i[k]*(kernel_degree + 1) + j ];

    /* Init the iteration for Horner's scheme. */
    w->v = ( c[0].v * x->v ) + c[1].v;
    dw_dx->v = c[0].v;
    
    /* And we're off! */
    for ( int k = 2 ; k <= kernel_degree ; k++ ) {
        dw_dx->v = ( dw_dx->v * x->v ) + w->v;
        w->v = ( x->v * w->v ) + c[k].v;
        }
        
    }
    
#endif

/**
 * @brief Computes the kernel for a given distance x. Gives a sensible answer only if x<1.
 */

__attribute__ ((always_inline)) INLINE static void kernel_eval ( float x , float *W ) {
    int ind = fmin( x , kernel_ivals );
    float *coeffs = &kernel_coeffs[ ind*(kernel_degree + 1) ];
    float w = coeffs[0]*x + coeffs[1];
    for ( int k = 2 ; k <= kernel_degree ; k++ )
        w = x*w + coeffs[k];
    *W = w;
    }


#endif //KERNEL_H
