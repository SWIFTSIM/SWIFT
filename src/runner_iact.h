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



__attribute__ ((always_inline)) INLINE void runner_iact_density ( float r2 , float hi , float hj , struct part *pi , struct part *pj ) {

    #define  KERNEL_COEFF_1  2.546479089470f
    #define  KERNEL_COEFF_2  15.278874536822f
    #define  KERNEL_COEFF_3  45.836623610466f
    #define  KERNEL_COEFF_4  30.557749073644f
    #define  KERNEL_COEFF_5  5.092958178941f
    #define  KERNEL_COEFF_6  (-15.278874536822f)
    #define  NORM_COEFF      4.188790204786f

    float r = sqrtf( r2 );
    float ui, uj, wi, wj;
    float ui_dh, uj_dh, wi_dh, wj_dh;
    
    if ( r2 < hi*hi && pi != NULL ) {
        
        ui = r / hi;
        ui_dh = -r / hi / hi;
        if ( ui < 0.5f ) {
            wi = KERNEL_COEFF_1 + KERNEL_COEFF_2 * (ui - 1.0f) * ui * ui;
            wi_dh = KERNEL_COEFF_2 * ui_dh * ui * ui
                  + 2 * KERNEL_COEFF_2 * (ui - 1.0f) * ui_dh * ui;
            }
        else {
            wi = KERNEL_COEFF_5 * (1.0f - ui) * (1.0f - ui) * (1.0f - ui);
            wi_dh = -3 * KERNEL_COEFF_5 * ui_dh * (1.0f - ui) * (1.0f - ui);
            }
        pi->rho += NORM_COEFF * wi;
        pi->rho_dh += NORM_COEFF * wi_dh;
        pi->icount += 1;
        
        }

    if ( r2 < hj*hj && pj != NULL ) {
        
        uj = r / hj;
        uj_dh = -r / hj / hj;
        if ( uj < 0.5f ) {
            wj = KERNEL_COEFF_1 + KERNEL_COEFF_2 * (uj - 1.0f) * uj * uj;
            wj_dh = KERNEL_COEFF_2 * uj_dh * uj * uj
                  + 2 * KERNEL_COEFF_2 * (uj - 1.0f) * uj_dh * uj;
            }
        else {
            wj = KERNEL_COEFF_5 * (1.0f - uj) * (1.0f - uj) * (1.0f - uj);
            wj_dh = -3 * KERNEL_COEFF_5 * uj_dh * (1.0f - uj) * (1.0f - uj);
            }
        pj->rho += NORM_COEFF * wj;
        pj->rho_dh += NORM_COEFF * wj_dh;
        pj->icount += 1;
            
        }
        
    #ifdef HIST
    if ( hi > hj )
        runner_hist_hit( hi / hj );
    else
        runner_hist_hit( hj / hi );
    #endif
    
    }
    


__attribute__ ((always_inline)) INLINE void runner_iact_force ( float r2 , float hi , float hj , struct part *pi , struct part *pj ) {

    #define  KERNEL_COEFF_1  2.546479089470f
    #define  KERNEL_COEFF_2  15.278874536822f
    #define  KERNEL_COEFF_3  45.836623610466f
    #define  KERNEL_COEFF_4  30.557749073644f
    #define  KERNEL_COEFF_5  5.092958178941f
    #define  KERNEL_COEFF_6  (-15.278874536822f)
    #define  NORM_COEFF      4.188790204786f
    
    }
    


