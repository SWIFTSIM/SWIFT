/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef SWIFT_PLANETARY_HYDRO_MISC_UTILS_H
#define SWIFT_PLANETARY_HYDRO_MISC_UTILS_H

/**
 * @file Planetary/hydro_misc_utils.h
 * @brief Miscellaneous hydro utilities for various scheme options.
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"


/**
 * @brief Sets stress tensor.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_set_sigma(
    struct part *restrict p, float deviatoric_stress_tensor_S[3][3], const float pressure) {

 for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        p->stress_tensor_sigma[i][j] = deviatoric_stress_tensor_S[i][j];
        if (i == j){
            // see schafer 2016 for this
            float pressure_hat = pressure;
            #ifdef MATERIAL_STRENGTH
            if (pressure < 0.f){
                pressure_hat *= (1 - p->damage_D);
            }
            #endif /* MATERIAL_STRENGTH */

            p->stress_tensor_sigma[i][j] -= pressure_hat;
            // For cylinders instead do
            //p->stress_tensor_sigma[i][j] -= 1.f*(p->rho - 1.f);
            
            // negative pressure cap (had to grep through miluphCUDA (see Schafer et al 2016) code for this. Maybe better ways)
            // This kind of thing should eventually be in EoS files
            if(pressure_hat < -p->Y_0){
              pressure_hat = -p->Y_0;   
            }
        }
    }
  }

}



/**
 * @brief Calculate eigenvalues of 3x3 matrix
 *
 */
__attribute__((always_inline)) INLINE static void compute_eigenvalues_symmetric_3x3(
    float *eigen_val1, float *eigen_val2, float *eigen_val3, float M[3][3]) {


    //Note: see calculate_all_eigenvalues in https://github.com/christophmschaefer/miluphcuda/blob/main/linalg.cu

    int i, j;
    // Current iteration towards diagonalisation
    float M_iter[3][3];
    // The largest (absolute value) off-diagonal element ...
    float max_offdiag;
    // ... and its indices
    int e, f;
    // trig functions corresponding to angle of rotation
    float sin_theta, cos_theta, tan_theta, one_over_tan_2theta;
    // rotated M_iter
    float M_iter_rotated[3][3];

    // init M_iter
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            M_iter[i][j] = M[i][j];
        }
    }

    // Note Schafer has limit of 5
    int max_iter = 5;
    float tol = 1e-10;
    // iterate until diagonalised within tolerance or until max_iter is reached
    for (int iter = 0; iter < max_iter; iter++){


        max_offdiag = 0.f;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                // init M_iter_rotated
                M_iter_rotated[i][j] = M_iter[i][j];
                // Find max off-diagonal element and its indices
                if (i != j){
                    if (fabs(M_iter[i][j]) >= max_offdiag){
                        max_offdiag = fabs(M_iter[i][j]);
                        e = i;
                        f = j;
                    }
                }
            }
        }

        // Has matrix has been diagonalised?
        if (max_offdiag < tol)
            break;

        // Get sin_theta and cos_theta

        one_over_tan_2theta = (M_iter[f][f] - M_iter[e][e])/(2 * M_iter[e][f]);
        //check this
        if (one_over_tan_2theta < 0)
            tan_theta = -1.f/(fabs(one_over_tan_2theta) + sqrtf(one_over_tan_2theta * one_over_tan_2theta + 1.f));
        else
            tan_theta = 1.f/(fabs(one_over_tan_2theta) + sqrtf(one_over_tan_2theta * one_over_tan_2theta + 1.f));

        cos_theta = 1.f/(sqrtf(tan_theta * tan_theta + 1));
        sin_theta = tan_theta * cos_theta;



        // Get M_iter_rotated by rotating M_iter
        M_iter_rotated[e][e] = cos_theta * cos_theta * M_iter[e][e] + sin_theta * sin_theta * M_iter[f][f] - 2.f * sin_theta * cos_theta * M_iter[e][f];
        M_iter_rotated[f][f] = cos_theta * cos_theta * M_iter[f][f] + sin_theta * sin_theta * M_iter[e][e] + 2.f * sin_theta * cos_theta * M_iter[e][f];
        M_iter_rotated[e][f] = (cos_theta * cos_theta - sin_theta * sin_theta) * M_iter[e][f] + sin_theta * cos_theta * (M_iter[e][e] - M_iter[f][f]);
        M_iter_rotated[f][e] = M_iter_rotated[e][f];

        /* the other element in column e and row f*/
        for (i = 0; i < 3; i++) {
            if (i != f && i != e){
                M_iter_rotated[e][i] = cos_theta * M_iter[i][e] - sin_theta * M_iter[i][f];
                M_iter_rotated[i][e] = M_iter_rotated[e][i];
                M_iter_rotated[f][i] = cos_theta * M_iter[i][f] + sin_theta * M_iter[i][e];
                M_iter_rotated[i][f] = M_iter_rotated[f][i];
            }
        }

        /* update  M_iter to M_iter_rotated for next iter or output */
        for (i = 0; i < 3; i++){
            for (j = 0; j < 3; j++){
                M_iter[i][j] = M_iter_rotated[i][j];
            }
        }
    }

    *eigen_val1 = M_iter[0][0];
    *eigen_val2 = M_iter[1][1];
    *eigen_val3 = M_iter[2][2];


}

#endif /* SWIFT_PLANETARY_HYDRO_MISC_UTILS_H */
