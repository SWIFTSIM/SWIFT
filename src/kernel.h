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
 * @brief SPH kernel functions. Compute W(x,h) and the gradient of W(x,h),
 *        as well as the blending function used for gravity.
 */

#include "vector.h"

/* Gravity kernel stuff ----------------------------------------------------------------------------------------------- */

/* The gravity kernel is defined as a degree 6 polynomial in the distance
   r. The resulting value should be post-multiplied with r^-3, resulting
   in a polynomial with terms ranging from r^-3 to r^3, which are
   sufficient to model both the direct potential as well as the splines
   near the origin. */
   
/* Coefficients for the gravity kernel. */
#define kernel_grav_degree 6
#define kernel_grav_ivals 2
#define kernel_grav_scale (2*const_iepsilon)
static float kernel_grav_coeffs[ (kernel_grav_degree+1) * (kernel_grav_ivals+1) ] =
    { 32.0f*const_iepsilon6 , -192.0f/5.0f*const_iepsilon5 , 0.0f , 32.0f/3.0f*const_iepsilon3 , 0.0f , 0.0f , 0.0f ,
      -32.0f/3.0f*const_iepsilon6 , 192.0f/5.0f*const_iepsilon5 , -48.0f*const_iepsilon4 , 64.0f/3.0f*const_iepsilon3 , 0.0f , 0.0f , -1.0f/15.0f ,
      0.0f , 0.0f , 0.0f , 0.0f , 0.0f , 0.0f , 1.0f };


/**
 * @brief Computes the gravity cubic spline for a given distance x.
 */

__attribute__ ((always_inline)) INLINE static void kernel_grav_eval ( float x , float *W ) {
    int ind = fmin( x*kernel_grav_scale , kernel_grav_ivals );
    float *coeffs = &kernel_grav_coeffs[ ind*(kernel_grav_degree + 1) ];
    float w = coeffs[0]*x + coeffs[1];
    for ( int k = 2 ; k <= kernel_grav_degree ; k++ )
        w = x*w + coeffs[k];
    *W = w;
    }


#ifdef VECTORIZE

/**
 * @brief Computes the gravity cubic spline for a given distance x (Vectorized version).
 */

__attribute__ ((always_inline)) INLINE static void kernel_grav_eval_vec ( vector *x , vector *w ) {
    
    vector ind, c[kernel_grav_degree+1];
    int j, k;
    
    /* Load x and get the interval id. */
    ind.m = vec_ftoi( vec_fmin( x->v*vec_set1( kernel_grav_scale ) , vec_set1( (float)kernel_grav_ivals ) ) );
    
    /* load the coefficients. */
    for ( k = 0 ; k < VEC_SIZE ; k++ )
        for ( j = 0 ; j < kernel_grav_degree+1 ; j++ )
            c[j].f[k] = kernel_grav_coeffs[ ind.i[k]*(kernel_grav_degree + 1) + j ];

    /* Init the iteration for Horner's scheme. */
    w->v = ( c[0].v * x->v ) + c[1].v;
    
    /* And we're off! */
    for ( int k = 2 ; k <= kernel_grav_degree ; k++ )
        w->v = ( x->v * w->v ) + c[k].v;
        
    }
    
    
#endif


/* Blending function stuff -------------------------------------------------------------------------------------------- */

/* Coefficients for the blending function. */
#define blender_degree 3
#define blender_ivals 3
#define blender_scale 4.0f
static float blender_coeffs[ (blender_degree+1) * (blender_ivals+1) ] =
    { 0.0f , 0.0f , 0.0f , 1.0f ,
      -32.0f , 24.0f , -6.0f , 1.5f , 
      -32.0f , 72.0f , -54.0f , 13.5f ,
      0.0f , 0.0f , 0.0f , 0.0f };
      
      
/**
 * @brief Computes the cubic spline blender for a given distance x.
 */

__attribute__ ((always_inline)) INLINE static void blender_eval ( float x , float *W ) {
    int ind = fmin( x*blender_scale , blender_ivals );
    float *coeffs = &blender_coeffs[ ind*(blender_degree + 1) ];
    float w = coeffs[0]*x + coeffs[1];
    for ( int k = 2 ; k <= blender_degree ; k++ )
        w = x*w + coeffs[k];
    *W = w;
    }


/**
 * @brief Computes the cubic spline blender and its derivative for a given distance x.
 */

__attribute__ ((always_inline)) INLINE static void blender_deval ( float x , float *W , float *dW_dx ) {
    int ind = fminf( x*blender_scale , blender_ivals );
    float *coeffs = &blender_coeffs[ ind*(blender_degree + 1) ];
    float w = coeffs[0]*x + coeffs[1];
    float dw_dx = coeffs[0];
    for ( int k = 2 ; k <= blender_degree ; k++ ) {
        dw_dx = dw_dx*x + w;
        w = x*w + coeffs[k];
        }
    *W = w;
    *dW_dx = dw_dx;
    }


#ifdef VECTORIZE

/**
 * @brief Computes the cubic spline blender and its derivative for a given distance x (Vectorized version). Gives a sensible answer only if x<2.
 */

__attribute__ ((always_inline)) INLINE static void blender_eval_vec ( vector *x , vector *w ) {
    
    vector ind, c[blender_degree+1];
    int j, k;
    
    /* Load x and get the interval id. */
    ind.m = vec_ftoi( vec_fmin( x->v*vec_set1( blender_scale ) , vec_set1( (float)blender_ivals ) ) );
    
    /* load the coefficients. */
    for ( k = 0 ; k < VEC_SIZE ; k++ )
        for ( j = 0 ; j < blender_degree+1 ; j++ )
            c[j].f[k] = blender_coeffs[ ind.i[k]*(blender_degree + 1) + j ];

    /* Init the iteration for Horner's scheme. */
    w->v = ( c[0].v * x->v ) + c[1].v;
    
    /* And we're off! */
    for ( int k = 2 ; k <= blender_degree ; k++ )
        w->v = ( x->v * w->v ) + c[k].v;
        
    }
    
    
/**
 * @brief Computes the cubic spline blender and its derivative for a given distance x (Vectorized version). Gives a sensible answer only if x<2.
 */

__attribute__ ((always_inline)) INLINE static void blender_deval_vec ( vector *x , vector *w , vector *dw_dx ) {
    
    vector ind, c[blender_degree+1];
    int j, k;
    
    /* Load x and get the interval id. */
    ind.m = vec_ftoi( vec_fmin( x->v*vec_set1( blender_scale ) , vec_set1( (float)blender_ivals ) ) );
    
    /* load the coefficients. */
    for ( k = 0 ; k < VEC_SIZE ; k++ )
        for ( j = 0 ; j < blender_degree+1 ; j++ )
            c[j].f[k] = blender_coeffs[ ind.i[k]*(blender_degree + 1) + j ];

    /* Init the iteration for Horner's scheme. */
    w->v = ( c[0].v * x->v ) + c[1].v;
    dw_dx->v = c[0].v;
    
    /* And we're off! */
    for ( int k = 2 ; k <= blender_degree ; k++ ) {
        dw_dx->v = ( dw_dx->v * x->v ) + w->v;
        w->v = ( x->v * w->v ) + c[k].v;
        }
        
    }
    
#endif


/* -------------------------------------------------------------------------------------------------------------------- */

#if defined(CUBIC_SPLINE_KERNEL)

/* -------------------------------------------------------------------------------------------------------------------- */

/* Coefficients for the kernel. */ 
#define kernel_name "Cubic spline"
#define kernel_degree 3
#define kernel_ivals 2
#define kernel_gamma 2.0f
#define kernel_gamma2 4.0f
#define kernel_gamma3 8.0f
#define kernel_igamma 0.5f
#define kernel_nwneigh ( 4.0/3.0*M_PI*const_eta_kernel*const_eta_kernel*const_eta_kernel*6.0858f ) 
static float kernel_coeffs[ (kernel_degree + 1) * (kernel_ivals + 1) ] __attribute__ ((aligned (16))) =
    { 3.0/4.0*M_1_PI , -3.0/2.0*M_1_PI , 0.0 , M_1_PI , 
      -0.25*M_1_PI , 3.0/2.0*M_1_PI , -3.0*M_1_PI , M_2_PI , 
      0.0 , 0.0 , 0.0 , 0.0 };
#define kernel_root ( kernel_coeffs[ kernel_degree ] )
#define kernel_wroot ( 4.0/3.0*M_PI*kernel_coeffs[ kernel_degree ] )

      
/**
 * @brief Computes the cubic spline kernel and its derivative for a given distance x. Gives a sensible answer only if x<2.
 */

__attribute__ ((always_inline)) INLINE static void kernel_deval ( float x , float *W , float *dW_dx ) {
    int ind = fminf( x , kernel_ivals );
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
 * @brief Computes the cubic spline kernel and its derivative for a given distance x (Vectorized version). Gives a sensible answer only if x<2.
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
 * @brief Computes the cubic spline kernel for a given distance x. Gives a sensible answer only if x<2.
 */

__attribute__ ((always_inline)) INLINE static void kernel_eval ( float x , float *W ) {
    int ind = fmin( x , kernel_ivals );
    float *coeffs = &kernel_coeffs[ ind*(kernel_degree + 1) ];
    float w = coeffs[0]*x + coeffs[1];
    for ( int k = 2 ; k <= kernel_degree ; k++ )
        w = x*w + coeffs[k];
    *W = w;
    }



/* -------------------------------------------------------------------------------------------------------------------- */

#elif defined(QUARTIC_SPLINE_KERNEL)

/* -------------------------------------------------------------------------------------------------------------------- */

/* Coefficients for the kernel. */ 
#define kernel_name "Quartic spline"
#define kernel_degree 4
#define kernel_ivals 3
#define kernel_gamma 2.5f
#define kernel_gamma2 6.25f
#define kernel_gamma3 15.625f
#define kernel_igamma 0.4f
#define kernel_nwneigh ( 4.0/3.0*M_PI*const_eta_kernel*const_eta_kernel*const_eta_kernel*8.2293f )
static float kernel_coeffs[ (kernel_degree + 1) * (kernel_ivals + 1) ] __attribute__ ((aligned (16))) =
  { 3.0/10.0*M_1_PI , 0.0  , -3.0/4.0*M_1_PI , 0.0 , 23.0/32.0*M_1_PI , 
    -1.0/5.0*M_1_PI , M_1_PI , -3.0/2.0*M_1_PI , 0.25*M_1_PI , 11.0/16.0*M_1_PI ,
    1.0/20.0*M_1_PI , -0.5*M_1_PI , 15.0/8.0*M_1_PI , -25.0/8.0*M_1_PI , 125.0/64.0*M_1_PI ,
    0.0 , 0.0 , 0.0 , 0.0 , 0.0 };
#define kernel_root ( kernel_coeffs[ kernel_degree ] )
#define kernel_wroot ( 4.0/3.0*M_PI*kernel_coeffs[ kernel_degree ] )
      
      
/**
 * @brief Computes the quartic spline kernel and its derivative for a given distance x. Gives a sensible answer only if x<2.5
 */

__attribute__ ((always_inline)) INLINE static void kernel_deval ( float x , float *W , float *dW_dx ) {
    int ind = fminf( x + 0.5, kernel_ivals);
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
 * @brief Computes the quartic spline kernel and its derivative for a given distance x (Vectorized version). Gives a sensible answer only if x<2.5
 */

__attribute__ ((always_inline)) INLINE static void kernel_deval_vec ( vector *x , vector *w , vector *dw_dx ) {
    
    vector ind, c[kernel_degree+1];
    int j, k;
    
    /* Load x and get the interval id. */
    ind.m = vec_ftoi( vec_fmin( x->v + 0.5f, vec_set1( (float)kernel_ivals ) ) );
    
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
 * @brief Computes the quartic spline kernel for a given distance x. Gives a sensible answer only if x<2.5
 */

__attribute__ ((always_inline)) INLINE static void kernel_eval ( float x , float *W ) {
    int ind = fmin( x + 0.5f, kernel_ivals );
    float *coeffs = &kernel_coeffs[ ind*(kernel_degree + 1) ];
    float w = coeffs[0]*x + coeffs[1];
    for ( int k = 2 ; k <= kernel_degree ; k++ )
        w = x*w + coeffs[k];
    *W = w;
    }







/* -------------------------------------------------------------------------------------------------------------------- */

#elif defined(QUINTIC_SPLINE_KERNEL)

/* -------------------------------------------------------------------------------------------------------------------- */

/* Coefficients for the kernel. */ 
#define kernel_name "Quintic spline"
#define kernel_degree 5
#define kernel_ivals 3
#define kernel_gamma 3.f
#define kernel_gamma2 9.f
#define kernel_gamma3 27.f
#define kernel_igamma 1.0f/3.0f
#define kernel_nwneigh ( 4.0/3.0*M_PI*const_eta_kernel*const_eta_kernel*const_eta_kernel*10.5868f )
static float kernel_coeffs[ (kernel_degree + 1) * (kernel_ivals + 1) ] __attribute__ ((aligned (16))) =
{ -1.0/12.0*M_1_PI  ,  1.0/4.0*M_1_PI ,  0.0            , -1.0/2.0*M_1_PI ,  0.0             , 11.0/20.0*M_1_PI,
  1.0/24.0*M_1_PI  , -3.0/8.0*M_1_PI ,  5.0/4.0*M_1_PI , -7.0/4.0*M_1_PI ,   5.0/8.0*M_1_PI , 17.0/40.0*M_1_PI ,
  -1.0/120.0*M_1_PI ,  1.0/8.0*M_1_PI , -3.0/4.0*M_1_PI ,  9.0/4.0*M_1_PI , -27.0/8.0*M_1_PI , 81.0/40.0*M_1_PI,
  0.0              , 0.0             , 0.0             , 0.0             ,   0.0            , 0.0};
#define kernel_root ( kernel_coeffs[ kernel_degree ] )
#define kernel_wroot ( 4.0/3.0*M_PI*kernel_coeffs[ kernel_degree ] )
      
      
/**
 * @brief Computes the quintic spline kernel and its derivative for a given distance x. Gives a sensible answer only if x<3.
 */

__attribute__ ((always_inline)) INLINE static void kernel_deval ( float x , float *W , float *dW_dx ) {
    int ind = fminf( x, kernel_ivals);
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
 * @brief Computes the quintic spline kernel and its derivative for a given distance x (Vectorized version). Gives a sensible answer only if x<3.
 */

__attribute__ ((always_inline)) INLINE static void kernel_deval_vec ( vector *x , vector *w , vector *dw_dx ) {
    
    vector ind, c[kernel_degree+1];
    int j, k;
    
    /* Load x and get the interval id. */
    ind.m = vec_ftoi( vec_fmin( x->v, vec_set1( (float)kernel_ivals ) ) );
    
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
 * @brief Computes the quintic spline kernel for a given distance x. Gives a sensible answer only if x<3.
 */

__attribute__ ((always_inline)) INLINE static void kernel_eval ( float x , float *W ) {
    int ind = fmin( x, kernel_ivals );
    float *coeffs = &kernel_coeffs[ ind*(kernel_degree + 1) ];
    float w = coeffs[0]*x + coeffs[1];
    for ( int k = 2 ; k <= kernel_degree ; k++ )
        w = x*w + coeffs[k];
    *W = w;
    }






/* -------------------------------------------------------------------------------------------------------------------- */

#elif defined(WENDLAND_C2_KERNEL)

/* -------------------------------------------------------------------------------------------------------------------- */

/* Coefficients for the kernel. */ 
#define kernel_name "Wendland C2"
#define kernel_degree 5
#define kernel_ivals 1
#define kernel_gamma 1.f
#define kernel_gamma2 1.f
#define kernel_gamma3 1.f
#define kernel_igamma 1.f
#define kernel_nwneigh ( 4.0/3.0*M_PI*const_eta_kernel*const_eta_kernel*const_eta_kernel*7.261825f )
static float kernel_coeffs[ (kernel_degree + 1) * (kernel_ivals + 1) ] __attribute__ ((aligned (16))) =
{  4.0f             , -15.0f           , 20.0f            , -10.0f           , 0.0f            , 1.0f,
  0.0f             , 0.0f             , 0.0f             , 0.0f            , 0.0f            , 0.0f};
#define kernel_root ( kernel_coeffs[ kernel_degree ] )
#define kernel_wroot ( 4.0/3.0*M_PI*kernel_coeffs[ kernel_degree ] )
      
      
/**
 * @brief Computes the quintic spline kernel and its derivative for a given distance x. Gives a sensible answer only if x<1.
 */

__attribute__ ((always_inline)) INLINE static void kernel_deval ( float x , float *W , float *dW_dx ) {
    int ind = fminf( x, kernel_ivals);
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
 * @brief Computes the Wendland C2 kernel and its derivative for a given distance x (Vectorized version). Gives a sensible answer only if x<1.
 */

__attribute__ ((always_inline)) INLINE static void kernel_deval_vec ( vector *x , vector *w , vector *dw_dx ) {
    
    vector ind, c[kernel_degree+1];
    int j, k;
    
    /* Load x and get the interval id. */
    ind.m = vec_ftoi( vec_fmin( x->v, vec_set1( (float)kernel_ivals ) ) );
    
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
 * @brief Computes the Wendland C2 kernel for a given distance x. Gives a sensible answer only if x<1.
 */

__attribute__ ((always_inline)) INLINE static void kernel_eval ( float x , float *W ) {
    int ind = fmin( x, kernel_ivals );
    float *coeffs = &kernel_coeffs[ ind*(kernel_degree + 1) ];
    float w = coeffs[0]*x + coeffs[1];
    for ( int k = 2 ; k <= kernel_degree ; k++ )
        w = x*w + coeffs[k];
    *W = w;
    }







/* -------------------------------------------------------------------------------------------------------------------- */

#else

/* -------------------------------------------------------------------------------------------------------------------- */

#error "A kernel function must be chosen in const.h !!"

#endif // Kernel choice

#endif //KERNEL_H
