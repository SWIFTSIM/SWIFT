/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
 

/* Load the vector stuff. */
#ifdef __SSE2__
    #define VECTORIZE
    #include <immintrin.h>
    #ifdef __AVX__
        typedef union {
            __m256 v;
            __m256i m;
            float f[8];
            int i[8];
            } vector;
        #define VEC_SIZE 8
        #define vec_load(a) _mm256_load_ps(a)
        #define vec_set1(a) _mm256_set1_ps(a)
        #define vec_sqrt(a) _mm256_sqrt_ps(a)
        #define vec_rcp(a) _mm256_rcp_ps(a)
        #define vec_rsqrt(a) _mm256_rsqrt_ps(a)
        #define vec_ftoi(a) _mm256_cvttps_epi32(a)
        #define vec_fmin(a,b) _mm256_min_ps(a,b)
    #else
        typedef union {
            __m128 v;
            __m128i m;
            float f[4];
            int i[4];
            } vector;
        #define VEC_SIZE 4
        #define vec_load(a) _mm_load_ps(a)
        #define vec_set1(a) _mm_set1_ps(a)
        #define vec_sqrt(a) _mm_sqrt_ps(a)
        #define vec_rcp(a) _mm_rcp_ps(a)
        #define vec_rsqrt(a) _mm_rsqrt_ps(a)
        #define vec_ftoi(a) _mm_cvttps_epi32(a)
        #define vec_fmin(a,b) _mm_min_ps(a,b)
    #endif
#else
    #define VEC_SIZE 4
#endif

/* Coefficients for the kernel. */ 
#define kernel_degree 3
#define kernel_ivals 2
#define kernel_gamma 0.5f
#define kernel_igamma 2.0f
#define kernel_igamma3 8.0f
#define kernel_igamma4 16.0f
static float kernel_coeffs[ (kernel_degree + 1) * (kernel_ivals + 1) ] __attribute__ ((aligned (16))) =
    { 3.0/4.0*M_1_PI , -3.0/2.0*M_1_PI , 0.0 , M_1_PI , 
      -0.25*M_1_PI , 3.0/2.0*M_1_PI , -3.0*M_1_PI , M_2_PI , 
      0.0 , 0.0 , 0.0 , 0.0 };
#define kernel_root ( 4.0/3.0*M_PI*kernel_igamma3*kernel_coeffs[ kernel_degree ] )
      
      
/**
 * @brief Helper function to evaluate the kernel at a given x.
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


__attribute__ ((always_inline)) INLINE static void kernel_eval ( float x , float *W ) {
    int ind = fmin( x , kernel_ivals );
    float *coeffs = &kernel_coeffs[ ind*(kernel_degree + 1) ];
    float w = coeffs[0]*x + coeffs[1];
    for ( int k = 2 ; k <= kernel_degree ; k++ )
        w = x*w + coeffs[k];
    *W = w;
    }


/**
 * @brief Density kernel
 */

__attribute__ ((always_inline)) INLINE static void runner_iact_density ( float r2 , float *dx , float hi , float hj , struct part *pi , struct part *pj ) {

    float r = sqrtf( r2 );
    float xi, xj;
    float h_inv, hg_inv;
    float wi, wj, wi_dx, wj_dx;
    
    if ( r2 < hi*hi && pi != NULL ) {
        
        h_inv = 1.0 / hi;
        hg_inv = kernel_igamma * h_inv;
        xi = r * hg_inv;
        kernel_deval( xi , &wi , &wi_dx );
        
        pi->rho += pj->mass * wi;
        pi->rho_dh += -pj->mass * ( 3.0*wi + xi*wi_dx );
        pi->wcount += wi * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        pi->wcount_dh -= xi * h_inv * wi_dx * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        // pi->icount += 1;
        
        }

    if ( r2 < hj*hj && pj != NULL ) {
        
        h_inv = 1.0 / hj;
        hg_inv = kernel_igamma * h_inv;
        xj = r * hg_inv;
        kernel_deval( xj , &wj , &wj_dx );
        
        pj->rho += pi->mass * wj;
        pj->rho_dh += -pi->mass * ( 3.0*wj + xj*wj_dx );
        pj->wcount += wj * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        pj->wcount_dh -= xj * h_inv * wj_dx * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        // pj->icount += 1;
        
        }
        
    #ifdef HIST
    if ( hi > hj )
        runner_hist_hit( hi / hj );
    else
        runner_hist_hit( hj / hi );
    #endif
    
    }
    

__attribute__ ((always_inline)) INLINE static void runner_iact_vec_density ( float *R2 , float *Dx , float *Hi , float *Hj , struct part **pi , struct part **pj ) {

#ifdef VECTORIZE

    vector r, xi, xj, hi, hj, hi_inv, hj_inv, hig_inv, hjg_inv, wi, wj, wi_dx, wj_dx;
    vector rhoi, rhoj, rhoi_dh, rhoj_dh, wcounti, wcountj, wcounti_dh, wcountj_dh;
    vector mi, mj, wscale;
    int k;
    
    r.v = vec_sqrt( vec_load( R2 ) );
    #if VEC_SIZE==8
        mi.v = _mm256_set_ps( pi[7]->mass , pi[6]->mass , pi[5]->mass , pi[4]->mass , pi[3]->mass , pi[2]->mass , pi[1]->mass , pi[0]->mass );
        mj.v = _mm256_set_ps( pj[7]->mass , pj[6]->mass , pj[5]->mass , pj[4]->mass , pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
    #elif VEC_SIZE==4
        mi.v = _mm_set_ps( pi[3]->mass , pi[2]->mass , pi[1]->mass , pi[0]->mass );
        mj.v = _mm_set_ps( pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
    #endif
    wscale.v = vec_set1( 4.0f * M_PI / 3.0f * kernel_igamma3 );
    
    hi.v = vec_load( Hi );
    hi_inv.v = vec_rcp( hi.v );
    hi_inv.v = hi_inv.v - hi_inv.v * ( hi_inv.v * hi.v  - vec_set1( 1.0f ) );
    hig_inv.v = hi_inv.v * vec_set1( kernel_igamma );
    xi.v = r.v * hig_inv.v;

    hj.v = vec_load( Hj );
    hj_inv.v = vec_rcp( hj.v );
    hj_inv.v = hj_inv.v - hj_inv.v * ( hj_inv.v * hj.v  - vec_set1( 1.0f ) );
    hjg_inv.v = hj_inv.v * vec_set1( kernel_igamma );
    xj.v = r.v * hjg_inv.v;
    
    kernel_deval_vec( &xi , &wi , &wi_dx );
    kernel_deval_vec( &xj , &wj , &wj_dx );
    
    rhoi.v = mj.v * wi.v;
    rhoi_dh.v = mj.v * ( vec_set1( 3.0f ) * wi.v + xi.v * wi_dx.v );
    wcounti.v = wi.v * wscale.v;
    wcounti_dh.v = xi.v * hi_inv.v * wi_dx.v * wscale.v;
        
    rhoj.v = mi.v * wj.v;
    rhoj_dh.v = mi.v * ( vec_set1( 3.0f ) * wj.v + xj.v * wj_dx.v );
    wcountj.v = wj.v * wscale.v;
    wcountj_dh.v = xj.v * hj_inv.v * wj_dx.v * wscale.v;
        
    for ( k = 0 ; k < VEC_SIZE ; k++ ) {
        pi[k]->rho += rhoi.f[k];
        pi[k]->rho_dh -= rhoi_dh.f[k];
        pi[k]->wcount += wcounti.f[k];
        pi[k]->wcount_dh -= wcounti_dh.f[k];
        pj[k]->rho += rhoj.f[k];
        pj[k]->rho_dh -= rhoj_dh.f[k];
        pj[k]->wcount += wcountj.f[k];
        pj[k]->wcount_dh -= wcountj_dh.f[k];
        }
        
#else

    for ( int k = 0 ; k < VEC_SIZE ; k++ )
        runner_iact_density( R2[k] , &Dx[3*k] , Hi[k] , Hj[k] , pi[k] , pj[k] );
        
#endif
    
    }
    


/**
 * @brief Density kernel
 */

__attribute__ ((always_inline)) INLINE static void runner_iact_nonsym_density ( float r2 , float *dx , float hi , float hj , struct part *pi , struct part *pj ) {

    float r = sqrtf( r2 );
    float xi;
    float h_inv, hg_inv;
    float wi, wi_dx;
    
    if ( r2 < hi*hi && pi != NULL ) {
        
        h_inv = 1.0 / hi;
        hg_inv = kernel_igamma * h_inv;
        xi = r * hg_inv;
        kernel_deval( xi , &wi , &wi_dx );
        
        pi->rho += pj->mass * wi;
        pi->rho_dh += -pj->mass * ( 3.0*wi + xi*wi_dx );
        pi->wcount += wi * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        pi->wcount_dh -= xi * h_inv * wi_dx * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        // pi->icount += 1;
        
        }

    }
    

__attribute__ ((always_inline)) INLINE static void runner_iact_nonsym_vec_density ( float *R2 , float *Dx , float *Hi , float *Hj , struct part **pi , struct part **pj ) {

#ifdef VECTORIZE

    vector r, xi, hi, hi_inv, hig_inv, wi, wi_dx;
    vector rhoi, rhoi_dh, wcounti, wcounti_dh;
    vector mj, wscale;
    int k;
    
    r.v = vec_sqrt( vec_load( R2 ) );
    #if VEC_SIZE==8
        mj.v = _mm256_set_ps( pj[7]->mass , pj[6]->mass , pj[5]->mass , pj[4]->mass , pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
    #elif VEC_SIZE==4
        mj.v = _mm_set_ps( pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
    #endif
    wscale.v = vec_set1( 4.0f * M_PI / 3.0f * kernel_igamma3 );
    
    hi.v = vec_load( Hi );
    hi_inv.v = vec_rcp( hi.v );
    hi_inv.v = hi_inv.v - hi_inv.v * ( hi_inv.v * hi.v  - vec_set1( 1.0f ) );
    hig_inv.v = hi_inv.v * vec_set1( kernel_igamma );
    xi.v = r.v * hig_inv.v;

    kernel_deval_vec( &xi , &wi , &wi_dx );
    
    rhoi.v = mj.v * wi.v;
    rhoi_dh.v = mj.v * ( vec_set1( 3.0f ) * wi.v + xi.v * wi_dx.v );
    wcounti.v = wi.v * wscale.v;
    wcounti_dh.v = xi.v * hi_inv.v * wi_dx.v * wscale.v;
        
    for ( k = 0 ; k < VEC_SIZE ; k++ ) {
        pi[k]->rho += rhoi.f[k];
        pi[k]->rho_dh -= rhoi_dh.f[k];
        pi[k]->wcount += wcounti.f[k];
        pi[k]->wcount_dh -= wcounti_dh.f[k];
        }
        
#else

    for ( int k = 0 ; k < VEC_SIZE ; k++ )
        runner_iact_nonsym_density( R2[k] , &Dx[3*k] , Hi[k] , Hj[k] , pi[k] , pj[k] );

#endif
        
    }
    


__attribute__ ((always_inline)) INLINE static void runner_iact_force ( float r2 , float *dx , float hi , float hj , struct part *pi , struct part *pj ) {

    float r = sqrtf( r2 ), ri = 1.0f / r;
    float xi, xj;
    float hig_inv, hig2_inv;
    float hjg_inv, hjg2_inv;
    float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr;
    float f;
    int k;
    
    /* Get the kernel for hi. */
    hig_inv = kernel_igamma / hi;
    hig2_inv = hig_inv * hig_inv;
    xi = r * hig_inv;
    kernel_deval( xi , &wi , &wi_dx );
    wi_dr = hig2_inv * hig2_inv * wi_dx;
        
    /* Get the kernel for hj. */
    hjg_inv = kernel_igamma / hj;
    hjg2_inv = hjg_inv * hjg_inv;
    xj = r * hjg_inv;
    kernel_deval( xj , &wj , &wj_dx );
    wj_dr = hjg2_inv * hjg2_inv * wj_dx;
    
    /* Get the common factor out. */
    w = ri * ( pi->POrho2 * wi_dr + pj->POrho2 * wj_dr );
        
    /* Use the force, Luke! */
    for ( k = 0 ; k < 3 ; k++ ) {
        f = dx[k] * w;
        pi->a[k] += pj->mass * f;
        pj->a[k] += -pi->mass * f;
        }
        
    /* Compute dv dot r. */
    dvdr = ( pi->v[0] - pj->v[0] ) * dx[0] + ( pi->v[1] - pj->v[1] ) * dx[1] + ( pi->v[2] - pj->v[2] ) * dx[2];
    dvdr *= ri;
        
    /* Get the time derivative for u. */
    pi->u_dt += pj->mass * dvdr * wi_dr;
    pj->u_dt += pi->mass * dvdr * wj_dr;
    
    /* Get the time derivative for h. */
    pi->h_dt += pj->mass / pj->rho * dvdr * wi_dr;
    pj->h_dt += pi->mass / pi->rho * dvdr * wj_dr;
        
    #ifdef HIST
    if ( hi > hj )
        runner_hist_hit( hi / hj );
    else
        runner_hist_hit( hj / hi );
    #endif
    
    }
    

__attribute__ ((always_inline)) INLINE static void runner_iact_vec_force ( float *R2 , float *Dx , float *Hi , float *Hj , struct part **pi , struct part **pj ) {

#ifdef VECTORIZE

    vector r, r2, ri;
    vector xi, xj;
    vector hi, hj, hi_inv, hj_inv;
    vector hig_inv, hig2_inv;
    vector hjg_inv, hjg2_inv;
    vector wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
    vector w;
    vector piPOrho2, pjPOrho2, pirho, pjrho;
    vector mi, mj;
    vector f;
    vector dx[3];
    vector vi[3], vj[3];
    vector pia[3], pja[3];
    vector piu_dt, pju_dt;
    vector pih_dt, pjh_dt;
    int j, k;

    /* Load stuff. */
    #if VEC_SIZE==8
        mi.v = _mm256_set_ps( pi[7]->mass , pi[6]->mass , pi[5]->mass , pi[4]->mass , pi[3]->mass , pi[2]->mass , pi[1]->mass , pi[0]->mass );
        mj.v = _mm256_set_ps( pj[7]->mass , pj[6]->mass , pj[5]->mass , pj[4]->mass , pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        piPOrho2.v = _mm256_set_ps( pi[7]->POrho2 , pi[6]->POrho2 , pi[5]->POrho2 , pi[4]->POrho2 , pi[3]->POrho2 , pi[2]->POrho2 , pi[1]->POrho2 , pi[0]->POrho2 );
        pjPOrho2.v = _mm256_set_ps( pj[7]->POrho2 , pj[6]->POrho2 , pj[5]->POrho2 , pj[4]->POrho2 , pj[3]->POrho2 , pj[2]->POrho2 , pj[1]->POrho2 , pj[0]->POrho2 );
        pirho.v = _mm256_set_ps( pi[7]->rho , pi[6]->rho , pi[5]->rho , pi[4]->rho , pi[3]->rho , pi[2]->rho , pi[1]->rho , pi[0]->rho );
        pjrho.v = _mm256_set_ps( pj[7]->rho , pj[6]->rho , pj[5]->rho , pj[4]->rho , pj[3]->rho , pj[2]->rho , pj[1]->rho , pj[0]->rho );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm256_set_ps( pi[7]->v[k] , pi[6]->v[k] , pi[5]->v[k] , pi[4]->v[k] , pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm256_set_ps( pj[7]->v[k] , pj[6]->v[k] , pj[5]->v[k] , pj[4]->v[k] , pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm256_set_ps( Dx[21+k] , Dx[18+k] , Dx[15+k] , Dx[12+k] , Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
    #elif VEC_SIZE==4
        mi.v = _mm_set_ps( pi[3]->mass , pi[2]->mass , pi[1]->mass , pi[0]->mass );
        mj.v = _mm_set_ps( pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        piPOrho2.v = _mm_set_ps( pi[3]->POrho2 , pi[2]->POrho2 , pi[1]->POrho2 , pi[0]->POrho2 );
        pjPOrho2.v = _mm_set_ps( pj[3]->POrho2 , pj[2]->POrho2 , pj[1]->POrho2 , pj[0]->POrho2 );
        pirho.v = _mm_set_ps( pi[3]->rho , pi[2]->rho , pi[1]->rho , pi[0]->rho );
        pjrho.v = _mm_set_ps( pj[3]->rho , pj[2]->rho , pj[1]->rho , pj[0]->rho );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm_set_ps( pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm_set_ps( pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm_set_ps( Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
    #else
        #error
    #endif
    
    /* Get the radius and inverse radius. */
    r2.v = vec_load( R2 );
    ri.v = vec_rsqrt( r2.v );
    ri.v = ri.v - vec_set1( 0.5f ) * ri.v * ( r2.v * ri.v * ri.v - vec_set1( 1.0f ) );
    r.v = r2.v * ri.v;
    
    /* Get the kernel for hi. */
    hi.v = vec_load( Hi );
    hi_inv.v = vec_rcp( hi.v );
    hi_inv.v = hi_inv.v - hi_inv.v * ( hi.v * hi_inv.v - vec_set1( 1.0f ) );
    hig_inv.v = vec_set1( kernel_igamma ) * hi_inv.v;
    hig2_inv.v = hig_inv.v * hig_inv.v;
    xi.v = r.v * hig_inv.v;
    kernel_deval_vec( &xi , &wi , &wi_dx );
    wi_dr.v = hig2_inv.v * hig2_inv.v * wi_dx.v;
        
    /* Get the kernel for hj. */
    hj.v = vec_load( Hj );
    hj_inv.v = vec_rcp( hj.v );
    hj_inv.v = hj_inv.v - hj_inv.v * ( hj.v * hj_inv.v - vec_set1( 1.0f ) );
    hjg_inv.v = vec_set1( kernel_igamma ) * hj_inv.v;
    hjg2_inv.v = hjg_inv.v * hjg_inv.v;
    xj.v = r.v * hjg_inv.v;
    kernel_deval_vec( &xj , &wj , &wj_dx );
    wj_dr.v = hjg2_inv.v * hjg2_inv.v * wj_dx.v;
        
    /* Get the common factor out. */
    w.v = ri.v * ( piPOrho2.v * wi_dr.v + pjPOrho2.v * wj_dr.v );
        
    /* Use the force, Luke! */
    for ( k = 0 ; k < 3 ; k++ ) {
        f.v = dx[k].v * w.v;
        pia[k].v = mj.v * f.v;
        pja[k].v = mi.v * f.v;
        }
        
    /* Compute dv dot r. */
    dvdr.v = ( (vi[0].v - vj[0].v) * dx[0].v ) + ( (vi[1].v - vj[1].v) * dx[1].v ) + ( (vi[2].v - vj[2].v) * dx[2].v );
    dvdr.v = dvdr.v * ri.v;
        
    /* Get the time derivative for u. */
    piu_dt.v = mj.v * dvdr.v * wi_dr.v;
    pju_dt.v = mi.v * dvdr.v * wj_dr.v;
    
    /* Get the time derivative for h. */
    pih_dt.v = mj.v / pjrho.v * dvdr.v * wi_dr.v;
    pjh_dt.v = mi.v / pirho.v * dvdr.v * wj_dr.v;
    
    /* Store the forces back on the particles. */
    for ( k = 0 ; k < VEC_SIZE ; k++ ) {
        pi[k]->u_dt += piu_dt.f[k];
        pj[k]->u_dt += pju_dt.f[k];
        pi[k]->h_dt += pih_dt.f[k];
        pj[k]->h_dt += pjh_dt.f[k];
        for ( j = 0 ; j < 3 ; j++ ) {
            pi[k]->a[j] += pia[j].f[k];
            pj[k]->a[j] -= pja[j].f[k];
            }
        }
        
#else

    for ( int k = 0 ; k < VEC_SIZE ; k++ )
        runner_iact_force( R2[k] , &Dx[3*k] , Hi[k] , Hj[k] , pi[k] , pj[k] );

#endif
        
    }
    


__attribute__ ((always_inline)) INLINE static void runner_iact_nonsym_force ( float r2 , float *dx , float hi , float hj , struct part *pi , struct part *pj ) {

    float r = sqrtf( r2 ), ri = 1.0f / r;
    float xi, xj;
    float hig_inv, hig2_inv;
    float hjg_inv, hjg2_inv;
    float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr;
    float f;
    int k;
    
    /* Get the kernel for hi. */
    hig_inv = kernel_igamma / hi;
    hig2_inv = hig_inv * hig_inv;
    xi = r * hig_inv;
    kernel_deval( xi , &wi , &wi_dx );
    wi_dr = hig2_inv * hig2_inv * wi_dx;
        
    /* Get the kernel for hj. */
    hjg_inv = kernel_igamma / hj;
    hjg2_inv = hjg_inv * hjg_inv;
    xj = r * hjg_inv;
    kernel_deval( xj , &wj , &wj_dx );
    wj_dr = hjg2_inv * hjg2_inv * wj_dx;
    
    /* Get the common factor out. */
    w = ri * ( pi->POrho2 * wi_dr + pj->POrho2 * wj_dr );
        
    /* Use the force, Luke! */
    for ( k = 0 ; k < 3 ; k++ ) {
        f = dx[k] * w;
        pi->a[k] += pj->mass * f;
        }
        
    /* Compute dv dot r. */
    dvdr = ( pi->v[0] - pj->v[0] ) * dx[0] + ( pi->v[1] - pj->v[1] ) * dx[1] + ( pi->v[2] - pj->v[2] ) * dx[2];
    dvdr *= ri;
        
    /* Get the time derivative for u. */
    pi->u_dt += pj->mass * dvdr * wi_dr;
    
    /* Get the time derivative for h. */
    pi->h_dt += pj->mass / pj->rho * dvdr * wi_dr;
        
    }
    

__attribute__ ((always_inline)) INLINE static void runner_iact_nonsym_vec_force ( float *R2 , float *Dx , float *Hi , float *Hj , struct part **pi , struct part **pj ) {

#ifdef VECTORIZE

    vector r, r2, ri;
    vector xi, xj;
    vector hi, hj, hi_inv, hj_inv;
    vector hig_inv, hig2_inv;
    vector hjg_inv, hjg2_inv;
    vector wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
    vector w;
    vector piPOrho2, pjPOrho2, pjrho;
    vector mj;
    vector f;
    vector dx[3];
    vector vi[3], vj[3];
    vector pia[3];
    vector piu_dt;
    vector pih_dt;
    int j, k;

    /* Load stuff. */
    #if VEC_SIZE==8
        mj.v = _mm256_set_ps( pj[7]->mass , pj[6]->mass , pj[5]->mass , pj[4]->mass , pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        piPOrho2.v = _mm256_set_ps( pi[7]->POrho2 , pi[6]->POrho2 , pi[5]->POrho2 , pi[4]->POrho2 , pi[3]->POrho2 , pi[2]->POrho2 , pi[1]->POrho2 , pi[0]->POrho2 );
        pjPOrho2.v = _mm256_set_ps( pj[7]->POrho2 , pj[6]->POrho2 , pj[5]->POrho2 , pj[4]->POrho2 , pj[3]->POrho2 , pj[2]->POrho2 , pj[1]->POrho2 , pj[0]->POrho2 );
        pjrho.v = _mm256_set_ps( pj[7]->rho , pj[6]->rho , pj[5]->rho , pj[4]->rho , pj[3]->rho , pj[2]->rho , pj[1]->rho , pj[0]->rho );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm256_set_ps( pi[7]->v[k] , pi[6]->v[k] , pi[5]->v[k] , pi[4]->v[k] , pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm256_set_ps( pj[7]->v[k] , pj[6]->v[k] , pj[5]->v[k] , pj[4]->v[k] , pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm256_set_ps( Dx[21+k] , Dx[18+k] , Dx[15+k] , Dx[12+k] , Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
    #elif VEC_SIZE==4
        mj.v = _mm_set_ps( pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        piPOrho2.v = _mm_set_ps( pi[3]->POrho2 , pi[2]->POrho2 , pi[1]->POrho2 , pi[0]->POrho2 );
        pjPOrho2.v = _mm_set_ps( pj[3]->POrho2 , pj[2]->POrho2 , pj[1]->POrho2 , pj[0]->POrho2 );
        pjrho.v = _mm_set_ps( pj[3]->rho , pj[2]->rho , pj[1]->rho , pj[0]->rho );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm_set_ps( pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm_set_ps( pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm_set_ps( Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
    #else
        #error
    #endif
    
    /* Get the radius and inverse radius. */
    r2.v = vec_load( R2 );
    ri.v = vec_rsqrt( r2.v );
    ri.v = ri.v - vec_set1( 0.5f ) * ri.v * ( r2.v * ri.v * ri.v - vec_set1( 1.0f ) );
    r.v = r2.v * ri.v;
    
    /* Get the kernel for hi. */
    hi.v = vec_load( Hi );
    hi_inv.v = vec_rcp( hi.v );
    hi_inv.v = hi_inv.v - hi_inv.v * ( hi.v * hi_inv.v - vec_set1( 1.0f ) );
    hig_inv.v = vec_set1( kernel_igamma ) * hi_inv.v;
    hig2_inv.v = hig_inv.v * hig_inv.v;
    xi.v = r.v * hig_inv.v;
    kernel_deval_vec( &xi , &wi , &wi_dx );
    wi_dr.v = hig2_inv.v * hig2_inv.v * wi_dx.v;
        
    /* Get the kernel for hj. */
    hj.v = vec_load( Hj );
    hj_inv.v = vec_rcp( hj.v );
    hj_inv.v = hj_inv.v - hj_inv.v * ( hj.v * hj_inv.v - vec_set1( 1.0f ) );
    hjg_inv.v = vec_set1( kernel_igamma ) * hj_inv.v;
    hjg2_inv.v = hjg_inv.v * hjg_inv.v;
    xj.v = r.v * hjg_inv.v;
    kernel_deval_vec( &xj , &wj , &wj_dx );
    wj_dr.v = hjg2_inv.v * hjg2_inv.v * wj_dx.v;
        
    /* Get the common factor out. */
    w.v = ri.v * ( piPOrho2.v * wi_dr.v + pjPOrho2.v * wj_dr.v );
        
    /* Use the force, Luke! */
    for ( k = 0 ; k < 3 ; k++ ) {
        f.v = dx[k].v * w.v;
        pia[k].v = mj.v * f.v;
        }
        
    /* Compute dv dot r. */
    dvdr.v = ( (vi[0].v - vj[0].v) * dx[0].v ) + ( (vi[1].v - vj[1].v) * dx[1].v ) + ( (vi[2].v - vj[2].v) * dx[2].v );
    dvdr.v = dvdr.v * ri.v;
        
    /* Get the time derivative for u. */
    piu_dt.v = mj.v * dvdr.v * wi_dr.v;
    
    /* Get the time derivative for h. */
    pih_dt.v = mj.v / pjrho.v * dvdr.v * wi_dr.v;
    
    /* Store the forces back on the particles. */
    for ( k = 0 ; k < VEC_SIZE ; k++ ) {
        pi[k]->u_dt += piu_dt.f[k];
        pi[k]->h_dt += pih_dt.f[k];
        for ( j = 0 ; j < 3 ; j++ )
            pi[k]->a[j] += pia[j].f[k];
        }

#else

    for ( int k = 0 ; k < VEC_SIZE ; k++ )
        runner_iact_nonsym_force( R2[k] , &Dx[3*k] , Hi[k] , Hj[k] , pi[k] , pj[k] );

#endif
        
    }
    



