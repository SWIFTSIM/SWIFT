/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Zachary Rey (zachary.rey@epfl.ch)
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

/* Include Header*/
#include "stellar_wind.h"

#include <math.h>
#include <stdio.h>
// #include "unit.h"

/* convert the quantities calculated here in erg */
double convertToErg(float energyToConvert){
    return energyToConvert * 1e5 * 1e5 * 1.989e33;
}

float calculate_b_parameter(const float log_Z, const float a[]){
    return a[2] + a[1] * log_Z + a[0] * log_Z * log_Z;
}

/**
 * @brief Compute the mass-loss using quadratic relation in Z and a polynomial degree 3 relation in M_ini.
 *
 * @param log_z The log_10 of the metallicity of the star in Z_sol
 * @param log_m The log_10 of the initial mass of the star in M_sol
 * 
 * @return The mass-loss due to stellar winds in M_sol/yr
 */
float stellar_wind_get_ejected_mass(const float log_Z,const float log_m){
    float coeffs[4][3]= {
        {-1.79317697e-02, -4.21353945e-01, -2.15261087e+01},
        {0.01857143,  0.58428571, 17.99571429},
        {-0.01448827, -0.14765458, -6.89331557},
        {0.00265458, 0.00591684, 0.99311301}
    };

    float mass_ejected = 0;

    /* If the star is lower than a limit mass, calculate the function*/
    if(log_m < stellar_wind_x0){
        for (int i=0; i < 4; i++){
            mass_ejected += calculate_b_parameter(log_Z,coeffs[i]) * pow(log_m,i);
        }
    /*  Else, the function will be only a linear relation between the initial mass and the function calculated for the limit mass */
    } else{
        float A0 = 0;
        float dLogA0 = 0;

        for (int i=0; i < 4; i++){
            A0 += calculate_b_parameter(log_Z,coeffs[i]) * pow(stellar_wind_x0,i);
            // The derivative take one less step in the loop
            if (i == 0){
                continue;
            }
            dLogA0 += i * (calculate_b_parameter(log_Z,coeffs[i]) * pow(stellar_wind_x0,(i - 1)));
        }

        mass_ejected += A0 + dLogA0 * (log_m - stellar_wind_x0);
    }

    return exp10(mass_ejected);
}

/**
 * @brief Compute the terminal wind velocity using quadratic relation in Z and a polynomial degree 3 relation in M_ini.
 *
 * @param log_z The log_10 of the metallicity of the star in Z_sol
 * @param log_m The log_10 of the initial mass of the star in M_sol
 * 
 * @return The wind velocity of stellar winds in km/s
 */
float stellar_wind_get_wind_velocity(const float log_Z,const float log_m){
    float coeffs[4][3]= {
        {-0.00513859,  0.14728145,  3.24600213},
        {0.00834222, 0.0259435,  0.33506397},
        {-0.00451493, -0.01548507, -0.06350746},
        { 0.00076972,  0.00254456, -0.004458 }
    };

    float wind_velocity = 0;

    /* If the star is lower than a limit mass, calculate the function*/
    if(log_m < stellar_wind_x0){
        for (int i=0; i < 4; i++){
            message("Inside wind velocity, calculate parameter %i is = %f", i, calculate_b_parameter(log_Z,coeffs[i]));
            wind_velocity += calculate_b_parameter(log_Z,coeffs[i]) * pow(log_m,i);
        }
    /*  Else, the function will be only a linear relation between the initial mass and the function calculated for the limit mass */
    } else{
        float A0 = 0;
        float dLogA0 = 0;

        for (int i=0; i < 4; i++){
            A0 += calculate_b_parameter(log_Z,coeffs[i]) * pow(stellar_wind_x0,i);
            message("Inside wind velocity, calculate parameter %i is = %f", i, calculate_b_parameter(log_Z,coeffs[i]));
            // The derivative take one less step in the loop
            if (i == 0){
                continue;
            }
            dLogA0 += i * (calculate_b_parameter(log_Z,coeffs[i]) * pow(stellar_wind_x0,(i - 1)));
        }

        wind_velocity += A0 + dLogA0 * (log_m - stellar_wind_x0);
    }

    message("Inside wind velocity computation : V_infty = %f", wind_velocity);

    return exp10(wind_velocity);
}

/**
 * @brief Compute the wind energy  using quadratic relation in Z and a polynomial degree 3 relation in M_ini.
 *
 * @param spart The stellar particle which has the wind properties
 * @return The wind energy in [units ? /yr]
 */
double stellar_wind_get_energy_dot(const float mass_loss, const float v_infinity){
    return convertToErg(0.5 * mass_loss * v_infinity * v_infinity);
}