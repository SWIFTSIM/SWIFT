/*******************************************************************************
 * This file is part of GadgetSMP.
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
 

/* Coefficients for the kernel. */
 
#define kernel_degree 3
#define kernel_ivals 2
#define kernel_gamma 0.5f
#define kernel_igamma 2.0f
#define kernel_igamma3 8.0f
#define kernel_igamma4 16.0f
static float kernel_coeffs[ (kernel_degree + 1) * (kernel_ivals + 1) ] =
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
        pi->icount += 1;
        
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
        pj->icount += 1;
        
        }
        
    #ifdef HIST
    if ( hi > hj )
        runner_hist_hit( hi / hj );
    else
        runner_hist_hit( hj / hi );
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
        pi->icount += 1;
        
        }

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
    


