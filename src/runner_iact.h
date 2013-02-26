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

#include "kernel.h"
#include "vector.h"

/**
 * @file  runner_iact.h
 * @brief SPH interaction functions following the Gadget-2 version of SPH.
 *
 * The interactions computed here are the ones presented in the Gadget-2 paper and use the same 
 * numerical coefficients as the Gadget-2 code. When used with the Spline-3 kernel, the results
 * should be equivalent to the ones obtained with Gadget-2 up to the rounding errors and interactions
 * missed by the Gadget-2 tree-code neighbours search.
 *
 * The code uses internal energy instead of entropy as a thermodynamical variable. 
 */


/**
 * @brief Density loop
 */

__attribute__ ((always_inline)) INLINE static void runner_iact_density ( float r2 , float *dx , float hi , float hj , struct part *pi , struct part *pj ) {

    float r = sqrtf( r2 ), ri = 1.0f / r;
    float xi, xj;
    float h_inv, hg_inv;
    float wi, wj, wi_dx, wj_dx;
    float dvdr;
    float curlvr[3];
    int k;
    
    /* Compute dv dot r */
    dvdr = ( pi->v[0] - pj->v[0] ) * dx[0] + ( pi->v[1] - pj->v[1] ) * dx[1] + ( pi->v[2] - pj->v[2] ) * dx[2];
    dvdr *= ri;

    /* Compute dv cross r */
    curlvr[0] = ( pi->v[1] - pj->v[1] ) * dx[2] - ( pi->v[2] - pj->v[2] ) * dx[1];
    curlvr[1] = ( pi->v[2] - pj->v[2] ) * dx[0] - ( pi->v[0] - pj->v[0] ) * dx[2];
    curlvr[2] = ( pi->v[0] - pj->v[0] ) * dx[1] - ( pi->v[1] - pj->v[1] ) * dx[0];
    for ( k = 0 ; k < 3 ; k++ )
        curlvr[k] *= ri;
            
    if ( r2 < hi*hi && pi != NULL ) {
        
        h_inv = 1.0 / hi;
        hg_inv = kernel_igamma * h_inv;
        xi = r * hg_inv;
        kernel_deval( xi , &wi , &wi_dx );
        
        pi->rho += pj->mass * wi;
        pi->rho_dh -= pj->mass * ( 3.0*wi + xi*wi_dx );
        pi->wcount += wi * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        pi->wcount_dh -= xi * h_inv * wi_dx * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        // pi->icount += 1;

	    pi->div_v += pj->mass * dvdr * wi_dx;
	    for ( k = 0 ; k < 3 ; k++ )
	        pi->curl_v[k] += pj->mass * curlvr[k] * wi_dx;
            
        }

    if ( r2 < hj*hj && pj != NULL ) {
        
        h_inv = 1.0 / hj;
        hg_inv = kernel_igamma * h_inv;
        xj = r * hg_inv;
        kernel_deval( xj , &wj , &wj_dx );
        
        pj->rho += pi->mass * wj;
        pj->rho_dh -= pi->mass * ( 3.0*wj + xj*wj_dx );
        pj->wcount += wj * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        pj->wcount_dh -= xj * h_inv * wj_dx * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        // pj->icount += 1;
        
	    pj->div_v = pi->mass * dvdr * wj_dx;
	    for ( k = 0 ; k < 3 ; k++ )
	        pj->curl_v[k] += pi->mass * curlvr[k] * wj_dx;
            
        }
        
    #ifdef HIST
    if ( hi > hj )
        runner_hist_hit( hi / hj );
    else
        runner_hist_hit( hj / hi );
    #endif
    
    }
    
/**
 * @brief Density loop (Vectorized version)
 */
__attribute__ ((always_inline)) INLINE static void runner_iact_vec_density ( float *R2 , float *Dx , float *Hi , float *Hj , struct part **pi , struct part **pj ) {

#ifdef VECTORIZE

    vector r, ri, xi, xj, hi, hj, hi_inv, hj_inv, hig_inv, hjg_inv, wi, wj, wi_dx, wj_dx;
    vector rhoi, rhoj, rhoi_dh, rhoj_dh, wcounti, wcountj, wcounti_dh, wcountj_dh;
    vector mi, mj, wscale;
    vector dx[3];
    vector vi[3], vj[3];    
    vector dvdr, div_vi, div_vj;
    vector curlvr[3], curl_vi[3], curl_vj[3];
    int k, j;
    
    r.v = vec_sqrt( vec_load( R2 ) );
    ri.v = vec_rcp( r.v );
    #if VEC_SIZE==8
        mi.v = _mm256_set_ps( pi[7]->mass , pi[6]->mass , pi[5]->mass , pi[4]->mass , pi[3]->mass , pi[2]->mass , pi[1]->mass , pi[0]->mass );
        mj.v = _mm256_set_ps( pj[7]->mass , pj[6]->mass , pj[5]->mass , pj[4]->mass , pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm256_set_ps( pi[7]->v[k] , pi[6]->v[k] , pi[5]->v[k] , pi[4]->v[k] , pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm256_set_ps( pj[7]->v[k] , pj[6]->v[k] , pj[5]->v[k] , pj[4]->v[k] , pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm256_set_ps( Dx[21+k] , Dx[18+k] , Dx[15+k] , Dx[12+k] , Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
    #elif VEC_SIZE==4
        mi.v = _mm_set_ps( pi[3]->mass , pi[2]->mass , pi[1]->mass , pi[0]->mass );
        mj.v = _mm_set_ps( pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm_set_ps( pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm_set_ps( pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm_set_ps( Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
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

    /* Compute dv dot r */
    dvdr.v = ( (vi[0].v - vj[0].v) * dx[0].v ) + ( (vi[1].v - vj[1].v) * dx[1].v ) + ( (vi[2].v - vj[2].v) * dx[2].v );
    dvdr.v = dvdr.v * ri.v;

    /* Compute dv cross r */
    curlvr[0].v = ( vi[1].v - vj[1].v ) * dx[2].v - ( vi[2].v - vj[2].v ) * dx[1].v;
    curlvr[1].v = ( vi[2].v - vj[2].v ) * dx[0].v - ( vi[0].v - vj[0].v ) * dx[2].v;
    curlvr[2].v = ( vi[0].v - vj[0].v ) * dx[1].v - ( vi[1].v - vj[1].v ) * dx[0].v;
    for ( k = 0 ; k < 3 ; k++ )
        curlvr[k].v *= ri.v;    

    rhoi.v = mj.v * wi.v;
    rhoi_dh.v = mj.v * ( vec_set1( 3.0f ) * wi.v + xi.v * wi_dx.v );
    wcounti.v = wi.v * wscale.v;
    wcounti_dh.v = xi.v * hi_inv.v * wi_dx.v * wscale.v;
    div_vi.v = mj.v * dvdr.v * wi_dx.v;
    for ( k = 0 ; k < 3 ; k++ )
        curl_vi[k].v = mj.v * curlvr[k].v * wi_dx.v;
        
    rhoj.v = mi.v * wj.v;
    rhoj_dh.v = mi.v * ( vec_set1( 3.0f ) * wj.v + xj.v * wj_dx.v );
    wcountj.v = wj.v * wscale.v;
    wcountj_dh.v = xj.v * hj_inv.v * wj_dx.v * wscale.v;
    div_vj.v = mi.v * dvdr.v * wj_dx.v;
    for ( k = 0 ; k < 3 ; k++ )
        curl_vj[k].v = mi.v * curlvr[k].v * wj_dx.v;

        
    for ( k = 0 ; k < VEC_SIZE ; k++ ) {
        pi[k]->rho += rhoi.f[k];
        pi[k]->rho_dh -= rhoi_dh.f[k];
        pi[k]->wcount += wcounti.f[k];
        pi[k]->wcount_dh -= wcounti_dh.f[k];
	    pi[k]->div_v += div_vi.f[k];
	    for( j = 0 ; j < 3 ; j++ )
   	        pi[k]->curl_v[j] += curl_vi[j].f[k];
        pj[k]->rho += rhoj.f[k];
        pj[k]->rho_dh -= rhoj_dh.f[k];
        pj[k]->wcount += wcountj.f[k];
        pj[k]->wcount_dh -= wcountj_dh.f[k];
	    pj[k]->div_v += div_vj.f[k];
	    for( j = 0 ; j < 3 ; j++ )
   	        pj[k]->curl_v[j] += curl_vj[j].f[k];
        }
        
#else

    for ( int k = 0 ; k < VEC_SIZE ; k++ )
        runner_iact_density( R2[k] , &Dx[3*k] , Hi[k] , Hj[k] , pi[k] , pj[k] );
        
#endif
    
    }
    


/**
 * @brief Density loop (non-symmetric version)
 */

__attribute__ ((always_inline)) INLINE static void runner_iact_nonsym_density ( float r2 , float *dx , float hi , float hj , struct part *pi , struct part *pj ) {

    float r = sqrtf( r2 ), ri = 1.0f / r;
    float xi;
    float h_inv, hg_inv;
    float wi, wi_dx;
    float dvdr;
    float curlvr[3];
    int k;

    /* Compute dv dot r */
    dvdr = ( pi->v[0] - pj->v[0] ) * dx[0] + ( pi->v[1] - pj->v[1] ) * dx[1] + ( pi->v[2] - pj->v[2] ) * dx[2];
    dvdr *= ri;

    /* Compute dv cross r */
    curlvr[0] = ( pi->v[1] - pj->v[1] ) * dx[2] - ( pi->v[2] - pj->v[2] ) * dx[1];
    curlvr[1] = ( pi->v[2] - pj->v[2] ) * dx[0] - ( pi->v[0] - pj->v[0] ) * dx[2];
    curlvr[2] = ( pi->v[0] - pj->v[0] ) * dx[1] - ( pi->v[1] - pj->v[1] ) * dx[0];
    for ( k = 0 ; k < 3 ; k++ )
        curlvr[k] *= ri;    

    if ( r2 < hi*hi && pi != NULL ) {
        
        h_inv = 1.0 / hi;
        hg_inv = kernel_igamma * h_inv;
        xi = r * hg_inv;
        kernel_deval( xi , &wi , &wi_dx );
        
        pi->rho += pj->mass * wi;
        pi->rho_dh -= pj->mass * ( 3.0*wi + xi*wi_dx );
        pi->wcount += wi * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        pi->wcount_dh -= xi * h_inv * wi_dx * ( 4.0f * M_PI / 3.0f * kernel_igamma3 );
        // pi->icount += 1;

	pi->div_v += pj->mass * dvdr * wi_dx;
	for ( k = 0 ; k < 3 ; k++ )
	    pi->curl_v[k] += pj->mass * curlvr[k] * wi_dx;
        }
    }
    
/**
 * @brief Density loop (non-symmetric vectorized version)
 */

__attribute__ ((always_inline)) INLINE static void runner_iact_nonsym_vec_density ( float *R2 , float *Dx , float *Hi , float *Hj , struct part **pi , struct part **pj ) {

#ifdef VECTORIZE

    vector r, ri, xi, hi, hi_inv, hig_inv, wi, wi_dx;
    vector rhoi, rhoi_dh, wcounti, wcounti_dh, div_vi;
    vector mj, wscale;
    vector dx[3];
    vector vi[3], vj[3];
    vector dvdr;
    vector curlvr[3], curl_vi[3];
    int k, j;
    
    r.v = vec_sqrt( vec_load( R2 ) );
    ri.v = vec_rcp( r.v );
    #if VEC_SIZE==8
        mj.v = _mm256_set_ps( pj[7]->mass , pj[6]->mass , pj[5]->mass , pj[4]->mass , pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm256_set_ps( pi[7]->v[k] , pi[6]->v[k] , pi[5]->v[k] , pi[4]->v[k] , pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm256_set_ps( pj[7]->v[k] , pj[6]->v[k] , pj[5]->v[k] , pj[4]->v[k] , pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm256_set_ps( Dx[21+k] , Dx[18+k] , Dx[15+k] , Dx[12+k] , Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
    #elif VEC_SIZE==4
        mj.v = _mm_set_ps( pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm_set_ps( pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm_set_ps( pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm_set_ps( Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
    #endif
    wscale.v = vec_set1( 4.0f * M_PI / 3.0f * kernel_igamma3 );
    
    hi.v = vec_load( Hi );
    hi_inv.v = vec_rcp( hi.v );
    hi_inv.v = hi_inv.v - hi_inv.v * ( hi_inv.v * hi.v  - vec_set1( 1.0f ) );
    hig_inv.v = hi_inv.v * vec_set1( kernel_igamma );
    xi.v = r.v * hig_inv.v;

    kernel_deval_vec( &xi , &wi , &wi_dx );

    /* Compute dv dot r */
    dvdr.v = ( (vi[0].v - vj[0].v) * dx[0].v ) + ( (vi[1].v - vj[1].v) * dx[1].v ) + ( (vi[2].v - vj[2].v) * dx[2].v );
    dvdr.v = dvdr.v * ri.v;

    /* Compute dv cross r */
    curlvr[0].v = ( vi[1].v - vj[1].v ) * dx[2].v - ( vi[2].v - vj[2].v ) * dx[1].v;
    curlvr[1].v = ( vi[2].v - vj[2].v ) * dx[0].v - ( vi[0].v - vj[0].v ) * dx[2].v;
    curlvr[2].v = ( vi[0].v - vj[0].v ) * dx[1].v - ( vi[1].v - vj[1].v ) * dx[0].v;
    for ( k = 0 ; k < 3 ; k++ )
        curlvr[k].v *= ri.v;    

    rhoi.v = mj.v * wi.v;
    rhoi_dh.v = mj.v * ( vec_set1( 3.0f ) * wi.v + xi.v * wi_dx.v );
    wcounti.v = wi.v * wscale.v;
    wcounti_dh.v = xi.v * hi_inv.v * wi_dx.v * wscale.v;
    div_vi.v = mj.v * dvdr.v * wi_dx.v;
    for ( k = 0 ; k < 3 ; k++ )
        curl_vi[k].v = mj.v * curlvr[k].v * wi_dx.v;
        
    for ( k = 0 ; k < VEC_SIZE ; k++ ) {
        pi[k]->rho += rhoi.f[k];
        pi[k]->rho_dh -= rhoi_dh.f[k];
        pi[k]->wcount += wcounti.f[k];
        pi[k]->wcount_dh -= wcounti_dh.f[k];
	    pi[k]->div_v += div_vi.f[k];
	    for( j = 0 ; j < 3 ; j++ )
   	        pi[k]->curl_v[j] += curl_vi[j].f[k];
        }
        
#else

    for ( int k = 0 ; k < VEC_SIZE ; k++ )
        runner_iact_nonsym_density( R2[k] , &Dx[3*k] , Hi[k] , Hj[k] , pi[k] , pj[k] );

#endif
        
    }
    

/**
 * @brief Force loop
 */

__attribute__ ((always_inline)) INLINE static void runner_iact_force ( float r2 , float *dx , float hi , float hj , struct part *pi , struct part *pj ) {

    float r = sqrtf( r2 ), ri = 1.0f / r;
    float xi, xj;
    float hi_inv, hi2_inv;
    float hj_inv, hj2_inv;
    float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr;
    float v_sig, omega_ij, Pi_ij;
    float f;
    int k;
    
    /* Get the kernel for hi. */
    hi_inv = 1.0f / hi;
    hi2_inv = hi_inv * hi_inv;
    xi = r * hi_inv * kernel_igamma;
    kernel_deval( xi , &wi , &wi_dx );
    wi_dr = hi2_inv * hi2_inv * wi_dx;
        
    /* Get the kernel for hj. */
    hj_inv = 1.0f / hj;
    hj2_inv = hj_inv * hj_inv;
    xj = r * hj_inv * kernel_igamma;
    kernel_deval( xj , &wj , &wj_dx );
    wj_dr = hj2_inv * hj2_inv * wj_dx;
                
    /* Compute dv dot r. */
    dvdr = ( pi->v[0] - pj->v[0] ) * dx[0] + ( pi->v[1] - pj->v[1] ) * dx[1] + ( pi->v[2] - pj->v[2] ) * dx[2];
    dvdr *= ri;

    /* Compute the relative velocity. (This is 0 if the particles move away from each other and negative otherwise) */
    omega_ij = fminf(dvdr, 0.f);
    
    /* Compute signal velocity */
    v_sig = pi->c + pj->c - 3.*omega_ij;

    /* Compute viscosity tensor */
    Pi_ij = -const_viscosity_alpha * v_sig * omega_ij / (pi->rho + pj->rho);

    /* Apply balsara switch */
    Pi_ij *= (pi->balsara + pj->balsara);

    /* Get the common factor out. */
    w = ri * ( ( pi->POrho2 * wi_dr + pj->POrho2 * wj_dr ) - 0.25f * Pi_ij * ( wi_dr + wj_dr ) );

    /* Use the force, Luke! */
    for ( k = 0 ; k < 3 ; k++ ) {
        f = dx[k] * w;
        pi->a[k] -= pj->mass * f;
        pj->a[k] += pi->mass * f;
        }
                
    /* Get the time derivative for u. */
    pi->u_dt += pi->POrho2 * pj->mass * dvdr * wi_dr + 0.125f * pj->mass * Pi_ij * dvdr * ( wi_dr + wj_dr );
    pj->u_dt += pj->POrho2 * pi->mass * dvdr * wj_dr + 0.125f * pi->mass * Pi_ij * dvdr * ( wi_dr + wj_dr );
    
    /* Get the time derivative for h. */
    pi->h_dt -= pj->mass / pj->rho * dvdr * wi_dr;
    pj->h_dt -= pi->mass / pi->rho * dvdr * wj_dr;
    
    /* Update the signal velocity. */
    pi->v_sig = fmaxf( pi->v_sig , v_sig );
    pj->v_sig = fmaxf( pj->v_sig , v_sig );
    
    #ifdef HIST
    if ( hi > hj )
        runner_hist_hit( hi / hj );
    else
        runner_hist_hit( hj / hi );
    #endif
    
    }
    

/**
 * @brief Force loop (Vectorized version)
 */

__attribute__ ((always_inline)) INLINE static void runner_iact_vec_force ( float *R2 , float *Dx , float *Hi , float *Hj , struct part **pi , struct part **pj ) {

#ifdef VECTORIZE

    vector r, r2, ri;
    vector xi, xj;
    vector hi, hj, hi_inv, hj_inv;
    vector hig_inv, hi2_inv;
    vector hjg_inv, hj2_inv;
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
    vector ci, cj, v_sig, vi_sig, vj_sig;
    vector omega_ij, Pi_ij, balsara;
    int j, k;

    /* Load stuff. */
    #if VEC_SIZE==8
        mi.v = _mm256_set_ps( pi[7]->mass , pi[6]->mass , pi[5]->mass , pi[4]->mass , pi[3]->mass , pi[2]->mass , pi[1]->mass , pi[0]->mass );
        mj.v = _mm256_set_ps( pj[7]->mass , pj[6]->mass , pj[5]->mass , pj[4]->mass , pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        piPOrho2.v = _mm256_set_ps( pi[7]->POrho2 , pi[6]->POrho2 , pi[5]->POrho2 , pi[4]->POrho2 , pi[3]->POrho2 , pi[2]->POrho2 , pi[1]->POrho2 , pi[0]->POrho2 );
        pjPOrho2.v = _mm256_set_ps( pj[7]->POrho2 , pj[6]->POrho2 , pj[5]->POrho2 , pj[4]->POrho2 , pj[3]->POrho2 , pj[2]->POrho2 , pj[1]->POrho2 , pj[0]->POrho2 );
        pirho.v = _mm256_set_ps( pi[7]->rho , pi[6]->rho , pi[5]->rho , pi[4]->rho , pi[3]->rho , pi[2]->rho , pi[1]->rho , pi[0]->rho );
        pjrho.v = _mm256_set_ps( pj[7]->rho , pj[6]->rho , pj[5]->rho , pj[4]->rho , pj[3]->rho , pj[2]->rho , pj[1]->rho , pj[0]->rho );
        ci.v = _mm256_set_ps( pi[7]->c , pi[6]->c , pi[5]->c , pi[4]->c , pi[3]->c , pi[2]->c , pi[1]->c , pi[0]->c );
        cj.v = _mm256_set_ps( pj[7]->c , pj[6]->c , pj[5]->c , pj[4]->c , pj[3]->c , pj[2]->c , pj[1]->c , pj[0]->c );
        vi_sig.v = _mm256_set_ps( pi[7]->v_sig , pi[6]->v_sig , pi[5]->v_sig , pi[4]->v_sig , pi[3]->v_sig , pi[2]->v_sig , pi[1]->v_sig , pi[0]->v_sig );
        vj_sig.v = _mm256_set_ps( pj[7]->v_sig , pj[6]->v_sig , pj[5]->v_sig , pj[4]->v_sig , pj[3]->v_sig , pj[2]->v_sig , pj[1]->v_sig , pj[0]->v_sig );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm256_set_ps( pi[7]->v[k] , pi[6]->v[k] , pi[5]->v[k] , pi[4]->v[k] , pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm256_set_ps( pj[7]->v[k] , pj[6]->v[k] , pj[5]->v[k] , pj[4]->v[k] , pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm256_set_ps( Dx[21+k] , Dx[18+k] , Dx[15+k] , Dx[12+k] , Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
        balsara.v = _mm256_set_ps( pi[7]->balsara , pi[6]->balsara , pi[5]->balsara , pi[4]->balsara , pi[3]->balsara , pi[2]->balsara , pi[1]->balsara , pi[0]->balsara ) +
                    _mm256_set_ps( pj[7]->balsara , pj[6]->balsara , pj[5]->balsara , pj[4]->balsara , pj[3]->balsara , pj[2]->balsara , pj[1]->balsara , pj[0]->balsara );
    #elif VEC_SIZE==4
        mi.v = _mm_set_ps( pi[3]->mass , pi[2]->mass , pi[1]->mass , pi[0]->mass );
        mj.v = _mm_set_ps( pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        piPOrho2.v = _mm_set_ps( pi[3]->POrho2 , pi[2]->POrho2 , pi[1]->POrho2 , pi[0]->POrho2 );
        pjPOrho2.v = _mm_set_ps( pj[3]->POrho2 , pj[2]->POrho2 , pj[1]->POrho2 , pj[0]->POrho2 );
        pirho.v = _mm_set_ps( pi[3]->rho , pi[2]->rho , pi[1]->rho , pi[0]->rho );
        pjrho.v = _mm_set_ps( pj[3]->rho , pj[2]->rho , pj[1]->rho , pj[0]->rho );
        ci.v = _mm_set_ps( pi[3]->c , pi[2]->c , pi[1]->c , pi[0]->c );
        cj.v = _mm_set_ps( pj[3]->c , pj[2]->c , pj[1]->c , pj[0]->c );
        vi_sig.v = _mm_set_ps( pi[3]->v_sig , pi[2]->v_sig , pi[1]->v_sig , pi[0]->v_sig );
        vj_sig.v = _mm_set_ps( pj[3]->v_sig , pj[2]->v_sig , pj[1]->v_sig , pj[0]->v_sig );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm_set_ps( pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm_set_ps( pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm_set_ps( Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
        balsara.v = _mm_set_ps( pi[3]->balsara , pi[2]->balsara , pi[1]->balsara , pi[0]->balsara ) +
                    _mm_set_ps( pj[3]->balsara , pj[2]->balsara , pj[1]->balsara , pj[0]->balsara );
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
    hi2_inv.v = hi_inv.v * hi_inv.v;
    xi.v = r.v * hig_inv.v;
    kernel_deval_vec( &xi , &wi , &wi_dx );
    wi_dr.v = hi2_inv.v * hi2_inv.v * wi_dx.v;
        
    /* Get the kernel for hj. */
    hj.v = vec_load( Hj );
    hj_inv.v = vec_rcp( hj.v );
    hj_inv.v = hj_inv.v - hj_inv.v * ( hj.v * hj_inv.v - vec_set1( 1.0f ) );
    hjg_inv.v = vec_set1( kernel_igamma ) * hj_inv.v;
    hj2_inv.v = hj_inv.v * hj_inv.v;
    xj.v = r.v * hjg_inv.v;
    kernel_deval_vec( &xj , &wj , &wj_dx );
    wj_dr.v = hj2_inv.v * hj2_inv.v * wj_dx.v;
        
    /* Compute dv dot r. */
    dvdr.v = ( (vi[0].v - vj[0].v) * dx[0].v ) + ( (vi[1].v - vj[1].v) * dx[1].v ) + ( (vi[2].v - vj[2].v) * dx[2].v );
    dvdr.v = dvdr.v * ri.v;
        
    /* Get the time derivative for h. */
    pih_dt.v = mj.v / pjrho.v * dvdr.v * wi_dr.v;
    pjh_dt.v = mi.v / pirho.v * dvdr.v * wj_dr.v;
    
    /* Compute the relative velocity. (This is 0 if the particles move away from each other and negative otherwise) */
    omega_ij.v = vec_fmin( dvdr.v , vec_set1( 0.0f ) );
    
    /* Compute signal velocity */
    v_sig.v = ci.v + cj.v - vec_set1( 3.0f )*omega_ij.v;

    /* Compute viscosity tensor */
    Pi_ij.v = -balsara.v * vec_set1( const_viscosity_alpha ) * v_sig.v * omega_ij.v / (pirho.v + pjrho.v);
    Pi_ij.v *= ( wi_dr.v + wj_dr.v );

    /* Get the common factor out. */
    w.v = ri.v * ( ( piPOrho2.v * wi_dr.v + pjPOrho2.v * wj_dr.v ) - vec_set1( 0.25f ) * Pi_ij.v );

    /* Use the force, Luke! */
    for ( k = 0 ; k < 3 ; k++ ) {
        f.v = dx[k].v * w.v;
        pia[k].v = mj.v * f.v;
        pja[k].v = mi.v * f.v;
        }
        
    /* Get the time derivative for u. */
    piu_dt.v = mj.v * dvdr.v * ( piPOrho2.v * wi_dr.v + vec_set1( 0.125f ) * Pi_ij.v );
    pju_dt.v = mi.v * dvdr.v * ( pjPOrho2.v * wj_dr.v + vec_set1( 0.125f ) * Pi_ij.v );
    
    /* compute the signal velocity (this is always symmetrical). */
    vi_sig.v = vec_fmax( vi_sig.v , v_sig.v );
    vj_sig.v = vec_fmax( vj_sig.v , v_sig.v );

    /* Store the forces back on the particles. */
    for ( k = 0 ; k < VEC_SIZE ; k++ ) {
        pi[k]->u_dt += piu_dt.f[k];
        pj[k]->u_dt += pju_dt.f[k];
        pi[k]->h_dt -= pih_dt.f[k];
        pj[k]->h_dt -= pjh_dt.f[k];
        pi[k]->v_sig = vi_sig.f[k];
        pj[k]->v_sig = vj_sig.f[k];
        for ( j = 0 ; j < 3 ; j++ ) {
            pi[k]->a[j] -= pia[j].f[k];
            pj[k]->a[j] += pja[j].f[k];
            }
        }
        
#else

    for ( int k = 0 ; k < VEC_SIZE ; k++ )
        runner_iact_force( R2[k] , &Dx[3*k] , Hi[k] , Hj[k] , pi[k] , pj[k] );

#endif
        
    }
    

/**
 * @brief Force loop (non-symmetric version)
 */

__attribute__ ((always_inline)) INLINE static void runner_iact_nonsym_force ( float r2 , float *dx , float hi , float hj , struct part *pi , struct part *pj ) {

    float r = sqrtf( r2 ), ri = 1.0f / r;
    float xi, xj;
    float hi_inv, hi2_inv;
    float hj_inv, hj2_inv;
    float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr;
    float f;
    float omega_ij, Pi_ij, v_sig;
    int k;
    
    /* Get the kernel for hi. */
    hi_inv = 1.0f / hi;
    hi2_inv = hi_inv * hi_inv;
    xi = r * hi_inv * kernel_igamma;
    kernel_deval( xi , &wi , &wi_dx );
    wi_dr = hi2_inv * hi2_inv * wi_dx;
        
    /* Get the kernel for hj. */
    hj_inv = 1.0f / hj;
    hj2_inv = hj_inv * hj_inv;
    xj = r * hj_inv * kernel_igamma;
    kernel_deval( xj , &wj , &wj_dx );
    wj_dr = hj2_inv * hj2_inv * wj_dx;

    /* Compute dv dot r. */
    dvdr = ( pi->v[0] - pj->v[0] ) * dx[0] + ( pi->v[1] - pj->v[1] ) * dx[1] + ( pi->v[2] - pj->v[2] ) * dx[2];
    dvdr *= ri;

    /* Compute the relative velocity. (This is 0 if the particles move away from each other and negative otherwise) */
    omega_ij = fminf(dvdr, 0.f);
    
    /* Compute signal velocity */
    v_sig = pi->c + pj->c - 3.*omega_ij;

    /* Compute viscosity tensor */
    Pi_ij = -const_viscosity_alpha * v_sig * omega_ij / (pi->rho + pj->rho);

    /* Apply balsara switch */
    Pi_ij *= (pi->balsara + pj->balsara);

    /* Get the common factor out. */
    w = ri * ( ( pi->POrho2 * wi_dr + pj->POrho2 * wj_dr ) - 0.25f * Pi_ij * ( wi_dr + wj_dr ) );

            
    /* Use the force, Luke! */
    for ( k = 0 ; k < 3 ; k++ ) {
        f = dx[k] * w;
        pi->a[k] -= pj->mass * f;
        }

    /* Get the time derivative for u. */
    pi->u_dt += pi->POrho2 * pj->mass * dvdr * wi_dr + 0.125f * pj->mass * Pi_ij * dvdr * ( wi_dr + wj_dr );
    
    /* Get the time derivative for h. */
    pi->h_dt -= pj->mass / pj->rho * dvdr * wi_dr;
    
    /* Update the signal velocity. */
    pi->v_sig = fmaxf( pi->v_sig , v_sig );
    pj->v_sig = fmaxf( pj->v_sig , v_sig );
    
    }
    

/**
 * @brief Force loop (Vectorized non-symmetric version)
 */

__attribute__ ((always_inline)) INLINE static void runner_iact_nonsym_vec_force ( float *R2 , float *Dx , float *Hi , float *Hj , struct part **pi , struct part **pj ) {

#ifdef VECTORIZE

    vector r, r2, ri;
    vector xi, xj;
    vector hi, hj, hi_inv, hj_inv;
    vector hig_inv, hi2_inv;
    vector hjg_inv, hj2_inv;
    vector wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
    vector w;
    vector piPOrho2, pjPOrho2, pirho, pjrho;
    vector mj;
    vector f;
    vector dx[3];
    vector vi[3], vj[3];
    vector pia[3];
    vector piu_dt;
    vector pih_dt;
    vector ci, cj, v_sig, vi_sig, vj_sig;
    vector omega_ij, Pi_ij, balsara;
    int j, k;

    /* Load stuff. */
    #if VEC_SIZE==8
        mj.v = _mm256_set_ps( pj[7]->mass , pj[6]->mass , pj[5]->mass , pj[4]->mass , pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        piPOrho2.v = _mm256_set_ps( pi[7]->POrho2 , pi[6]->POrho2 , pi[5]->POrho2 , pi[4]->POrho2 , pi[3]->POrho2 , pi[2]->POrho2 , pi[1]->POrho2 , pi[0]->POrho2 );
        pjPOrho2.v = _mm256_set_ps( pj[7]->POrho2 , pj[6]->POrho2 , pj[5]->POrho2 , pj[4]->POrho2 , pj[3]->POrho2 , pj[2]->POrho2 , pj[1]->POrho2 , pj[0]->POrho2 );
        pirho.v = _mm256_set_ps( pi[7]->rho , pi[6]->rho , pi[5]->rho , pi[4]->rho , pi[3]->rho , pi[2]->rho , pi[1]->rho , pi[0]->rho );
        pjrho.v = _mm256_set_ps( pj[7]->rho , pj[6]->rho , pj[5]->rho , pj[4]->rho , pj[3]->rho , pj[2]->rho , pj[1]->rho , pj[0]->rho );
        ci.v = _mm256_set_ps( pi[7]->c , pi[6]->c , pi[5]->c , pi[4]->c , pi[3]->c , pi[2]->c , pi[1]->c , pi[0]->c );
        cj.v = _mm256_set_ps( pj[7]->c , pj[6]->c , pj[5]->c , pj[4]->c , pj[3]->c , pj[2]->c , pj[1]->c , pj[0]->c );
        vi_sig.v = _mm256_set_ps( pi[7]->v_sig , pi[6]->v_sig , pi[5]->v_sig , pi[4]->v_sig , pi[3]->v_sig , pi[2]->v_sig , pi[1]->v_sig , pi[0]->v_sig );
        vj_sig.v = _mm256_set_ps( pj[7]->v_sig , pj[6]->v_sig , pj[5]->v_sig , pj[4]->v_sig , pj[3]->v_sig , pj[2]->v_sig , pj[1]->v_sig , pj[0]->v_sig );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm256_set_ps( pi[7]->v[k] , pi[6]->v[k] , pi[5]->v[k] , pi[4]->v[k] , pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm256_set_ps( pj[7]->v[k] , pj[6]->v[k] , pj[5]->v[k] , pj[4]->v[k] , pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm256_set_ps( Dx[21+k] , Dx[18+k] , Dx[15+k] , Dx[12+k] , Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
        balsara.v = _mm256_set_ps( pi[7]->balsara , pi[6]->balsara , pi[5]->balsara , pi[4]->balsara , pi[3]->balsara , pi[2]->balsara , pi[1]->balsara , pi[0]->balsara ) +
                    _mm256_set_ps( pj[7]->balsara , pj[6]->balsara , pj[5]->balsara , pj[4]->balsara , pj[3]->balsara , pj[2]->balsara , pj[1]->balsara , pj[0]->balsara );
    #elif VEC_SIZE==4
        mj.v = _mm_set_ps( pj[3]->mass , pj[2]->mass , pj[1]->mass , pj[0]->mass );
        piPOrho2.v = _mm_set_ps( pi[3]->POrho2 , pi[2]->POrho2 , pi[1]->POrho2 , pi[0]->POrho2 );
        pjPOrho2.v = _mm_set_ps( pj[3]->POrho2 , pj[2]->POrho2 , pj[1]->POrho2 , pj[0]->POrho2 );
        pirho.v = _mm_set_ps( pi[3]->rho , pi[2]->rho , pi[1]->rho , pi[0]->rho );
        pjrho.v = _mm_set_ps( pj[3]->rho , pj[2]->rho , pj[1]->rho , pj[0]->rho );
        ci.v = _mm_set_ps( pi[3]->c , pi[2]->c , pi[1]->c , pi[0]->c );
        cj.v = _mm_set_ps( pj[3]->c , pj[2]->c , pj[1]->c , pj[0]->c );
        vi_sig.v = _mm_set_ps( pi[3]->v_sig , pi[2]->v_sig , pi[1]->v_sig , pi[0]->v_sig );
        vj_sig.v = _mm_set_ps( pj[3]->v_sig , pj[2]->v_sig , pj[1]->v_sig , pj[0]->v_sig );
        for ( k = 0 ; k < 3 ; k++ ) {
            vi[k].v = _mm_set_ps( pi[3]->v[k] , pi[2]->v[k] , pi[1]->v[k] , pi[0]->v[k] );
            vj[k].v = _mm_set_ps( pj[3]->v[k] , pj[2]->v[k] , pj[1]->v[k] , pj[0]->v[k] );
            }
        for ( k = 0 ; k < 3 ; k++ )
            dx[k].v = _mm_set_ps( Dx[9+k] , Dx[6+k] , Dx[3+k] , Dx[0+k] );
        balsara.v = _mm_set_ps( pi[3]->balsara , pi[2]->balsara , pi[1]->balsara , pi[0]->balsara ) +
                    _mm_set_ps( pj[3]->balsara , pj[2]->balsara , pj[1]->balsara , pj[0]->balsara );
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
    hi2_inv.v = hi_inv.v * hi_inv.v;
    xi.v = r.v * hig_inv.v;
    kernel_deval_vec( &xi , &wi , &wi_dx );
    wi_dr.v = hi2_inv.v * hi2_inv.v * wi_dx.v;
        
    /* Get the kernel for hj. */
    hj.v = vec_load( Hj );
    hj_inv.v = vec_rcp( hj.v );
    hj_inv.v = hj_inv.v - hj_inv.v * ( hj.v * hj_inv.v - vec_set1( 1.0f ) );
    hjg_inv.v = vec_set1( kernel_igamma ) * hj_inv.v;
    hj2_inv.v = hj_inv.v * hj_inv.v;
    xj.v = r.v * hjg_inv.v;
    kernel_deval_vec( &xj , &wj , &wj_dx );
    wj_dr.v = hj2_inv.v * hj2_inv.v * wj_dx.v;
        
    /* Compute dv dot r. */
    dvdr.v = ( (vi[0].v - vj[0].v) * dx[0].v ) + ( (vi[1].v - vj[1].v) * dx[1].v ) + ( (vi[2].v - vj[2].v) * dx[2].v );
    dvdr.v = dvdr.v * ri.v;
        
    /* Get the time derivative for h. */
    pih_dt.v = mj.v / pjrho.v * dvdr.v * wi_dr.v;
    
    /* Compute the relative velocity. (This is 0 if the particles move away from each other and negative otherwise) */
    omega_ij.v = vec_fmin( dvdr.v , vec_set1( 0.0f ) );
    
    /* Compute signal velocity */
    v_sig.v = ci.v + cj.v - vec_set1( 3.0f )*omega_ij.v;

    /* Compute viscosity tensor */
    Pi_ij.v = -balsara.v * vec_set1( const_viscosity_alpha ) * v_sig.v * omega_ij.v / (pirho.v + pjrho.v);
    Pi_ij.v *= ( wi_dr.v + wj_dr.v );

    /* Get the common factor out. */
    w.v = ri.v * ( ( piPOrho2.v * wi_dr.v + pjPOrho2.v * wj_dr.v ) - vec_set1( 0.25f ) * Pi_ij.v );

    /* Use the force, Luke! */
    for ( k = 0 ; k < 3 ; k++ ) {
        f.v = dx[k].v * w.v;
        pia[k].v = mj.v * f.v;
        }
        
    /* Get the time derivative for u. */
    piu_dt.v = mj.v * dvdr.v * ( piPOrho2.v * wi_dr.v + vec_set1( 0.125f ) * Pi_ij.v );
    
    /* compute the signal velocity (this is always symmetrical). */
    vi_sig.v = vec_fmax( vi_sig.v , v_sig.v );
    vj_sig.v = vec_fmax( vj_sig.v , v_sig.v );

    /* Store the forces back on the particles. */
    for ( k = 0 ; k < VEC_SIZE ; k++ ) {
        pi[k]->u_dt += piu_dt.f[k];
        pi[k]->h_dt -= pih_dt.f[k];
        pi[k]->v_sig = vi_sig.f[k];
        pj[k]->v_sig = vj_sig.f[k];
        for ( j = 0 ; j < 3 ; j++ )
            pi[k]->a[j] -= pia[j].f[k];
        }

#else

    for ( int k = 0 ; k < VEC_SIZE ; k++ )
        runner_iact_nonsym_force( R2[k] , &Dx[3*k] , Hi[k] , Hj[k] , pi[k] , pj[k] );

#endif
        
    }
    



