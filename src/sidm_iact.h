/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Camila Correa (camila.correa@uva.nl)
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
#ifndef SWIFT_SIDM_IACT_H
#define SWIFT_SIDM_IACT_H

/* Standard headers */
#include <float.h>

/* Local includes. */
#include "sidm.h"
#include "sidm_properties.h"

/**
 * @brief do self-interacting DM computation. Computes the probability of DM-DM interactions
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param gpi First particle.
 * @param gpj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_sidm(float r2, const float *dx, float hi, float hj,
                 struct gpart *gpi, struct gpart *gpj,
                 float a, float H, const double dt_Gyr,
                 const struct sidm_props* sidm_props) {
        
    /* Calculate probability of gparticles i & j of scattering within the next time step */
    float dv[3] = {gpi->v_full[0] - gpj->v_full[0], gpi->v_full[1] - gpj->v_full[1], gpi->v_full[2] - gpj->v_full[2]};
    const float v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
    double sigma = sidm_props->sigma_phys;
    float Pij = sqrtf(v2) * dt_Gyr * sigma;
    printf("Pij %f",Pij);
    
}


#endif /* SWIFT_SIDM_IACT_H */
