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

#include "interpolation.h"
#include <math.h>
#include <stdio.h>
// #include "unit.h"


/**
 * @brief Print the stellar wind model.
 *
 * @param sw The #stellarwind.
 */
// void stellar_wind_print(const struct stellar_wind *sw) {

    /* Only the master print */
    // if (engine_rank != 0) {
    //   return;
    // }
  
    // message("Mass range for SW = [%g, %g]", sw->mass_min, sw->mass_max);
    //TODO anything else ?
//}


// float stellar_wind_get_energy_discrete(const struct stellar_wind *sw, float mass, float metallicity) {
//     return interpolate_2d();//TODO


/**
 * @brief Read an array of SW from the table.
 *
 * @param sw The #stellar_wind model.
 * @param interp_raw Interpolation data to initialize (raw).
 * @param interp_int Interpolation data to initialize (integrated).
 * @param sm * The #stellar_model.
 * @param group_id The HDF5 group id where to read from.
 * @param hdf5_dataset_name The dataset name to read.
 * @param previous_count Number of element in the previous array read.
 * @param interpolation_size_m Number of element to keep in the mass interpolation data.
 * @param interpolation_size_z Number of element to keep in the metallicity interpolation data.
 */
void stellar_wind_read_yields_array(
    struct stellar_wind *sw, struct interpolation_2d *interp,
    const struct stellar_model *sm,
    hid_t group_id, const char *hdf5_dataset_name, hsize_t *previous_count,
    int interpolation_size_m, int interpolation_size_z) {

    /* Now let's get the number of elements */
    /* Open attribute */
    const hid_t h_dataset = H5Dopen(group_id, hdf5_dataset_name, H5P_DEFAULT);
    if (h_dataset < 0)
        error("Error while opening attribute '%s'", hdf5_dataset_name);

    /* Get the number of elements */
    hsize_t count = io_get_number_element_in_dataset(h_dataset);

    // message("count = %d",(int)count);

    /* Check that all the arrays have the same size */
    if (*previous_count != 0 && count != *previous_count) {
        error("The code is not able to deal with yields arrays of different size");
    }
    *previous_count = count;

    /* Read the minimal mass (in log) */
    float log_mass_min = 0;
    io_read_attribute(h_dataset, "m0", FLOAT, &log_mass_min);

    /* Read the step size in mass array (log step) */
    float step_size_m = 0;
    io_read_attribute(h_dataset, "dm", FLOAT, &step_size_m);

    int Nm = 0;
    io_read_attribute(h_dataset, "nm", INT, &Nm);
    // message("nm = %d",Nm);

    /* Read the minimal metallicity (in log) */
    float log_metallicity_min = 0;
    io_read_attribute(h_dataset, "z0", FLOAT, &log_metallicity_min);

    /* Read the step size in metallicity array (log step) */
    float step_size_z = 0;
    io_read_attribute(h_dataset, "dz", FLOAT, &step_size_z);

    int Nz = 0;
    io_read_attribute(h_dataset, "nz", INT, &Nz);
    

    /* Close the attribute */
    H5Dclose(h_dataset);
    // message("nm = %d  nz = %d",Nm,Nz);
    // message("m0 = %f , z0 = %f",log_mass_min, log_metallicity_min);
    // message("dm = %f, dz = %f", step_size_m, step_size_z);
    

    /* Allocate the memory */
    
    double *data = (double *)malloc(sizeof(double) * Nz * Nm);
    if (data == NULL)
        error("Failed to allocate the SW yields for %s.", hdf5_dataset_name);

    //   for (int i = 0; i < Nz; i++) {
    //     data[i] = (float *)malloc(Nm * sizeof(float));
    //     if (data[i] == NULL){
    //         error("Failed to allocate the SW array for %s.", hdf5_dataset_name);
    //     }
    //   }

    /* Read the dataset */
    io_read_array_dataset(group_id, hdf5_dataset_name, DOUBLE, data, count);
    //   for (int i = 0; i < (int)count; i++){
    //     message("step %d -> %g", i, data[i]);
    //   }
    

    /* Initialize the raw interpolation */
    //   interpolate_2d_init(interp_raw, log10(sw->mass_min), log10(sw->mass_max),
    //                       interpolation_size_m, log10(sw->metallicity_min), log10(sw->metallicity_max),
    //                       interpolation_size_z, log_mass_min, log_metallicity_min, step_size_m,
    //                       step_size_z, count_x, count_y, data,
    //                       boundary_condition_zero);
    float log_mass_max = log_mass_min + step_size_m * Nm;
    float log_metallicity_max = log_metallicity_min + step_size_z * Nz;
    // message("log_mass_max = %g, log_metallicity_max = %g",log_mass_max, log_metallicity_max);
    interpolate_2d_init(interp, 
                        log_metallicity_min, log_metallicity_max, Nz,
                        log_mass_min, log_mass_max, Nm, 
                        log_metallicity_min, log_mass_min, 
                        step_size_z, step_size_m, Nz, Nm, data,
                        boundary_condition_const);

    /* Cleanup the memory */
    free(data);
}

/**
 * @brief Read the SW yields from the table.
 *
 * The tables are in internal units at the end of this function.
 *
 * @param sw The #stellar_wind model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param restart Are we restarting the simulation? (Is params NULL?)
 */
void stellar_wind_read_yields(struct stellar_wind *sw,
    //struct swift_params *params,
    const struct stellar_model *sm,
    const int restart) {

    hid_t file_id, group_id;

    hsize_t previous_count = 0;

    // if (!restart) {
    // sw->interpolation_size = parser_get_opt_param_int(
    // params, "GEARSupernovaeII:interpolation_size", 200);
    // }

    /* Open IMF group */
    h5_open_group(sm->yields_table, "Data/SW/MetallicityDependent", &file_id, &group_id);
    
    /* Read the energy*/
    stellar_wind_read_yields_array(sw, &sw->raw.ejected_energy,
        sm, group_id, "Energy", &previous_count,
        sw->interpolation_size_m, sw->interpolation_size_z);

    /* Read the integrated energy*/
    stellar_wind_read_yields_array(sw, &sw->integrated.ejected_energy_per_progenitor_mass, 
        sm, group_id, "Integrated_Energy", &previous_count,
        sw->interpolation_size_m, sw->interpolation_size_z);

    /* Cleanup everything */
    h5_close_group(file_id, group_id);
};


/**
 * @brief Initialize the #supernovae_ii structure.
 *
 * @param snii The #supernovae_ii model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void stellar_wind_init(struct stellar_wind *sw, struct swift_params *params,
    const struct stellar_model *sm,
    const struct unit_system *us) {


    /* Read the supernovae yields */
    stellar_wind_read_yields(sw, /*params,*/ sm, /* restart */ 0);

}


/**
 * @brief Get the ejected energy given a discrete mass.
 *
 * @param sw The #stellar_wind model.
 * @param log_m The upper mass in log.
 * @param log_z The metallicity in log. 
 *
 * @return energy per unit time in [erg/yr].
 */
double stellar_wind_get_ejected_energy(const struct stellar_wind *sw, double log_m, float log_z) {
//   message("inside get_energy");
  return pow(10,interpolate_2d(&sw->raw.ejected_energy, log_z, log_m));
};

/**
 * @brief Get the ejected energy per progenitor mass.
 *
 * @param sw The #stellar_wind model.
 * @param log_m The upper mass in log.
 * @param log_z The metallicity in log. 
 *
 * @return energy per progenitor mass per unit time in [erg/yr].
 */
double stellar_wind_get_ejected_energy_IMF(const struct stellar_wind *sw, double log_m, float log_z) {
//   message("inside get_energy_IMF");
  return pow(10,interpolate_2d(&sw->integrated.ejected_energy_per_progenitor_mass, log_z, log_m));
};


/**
 * @brief Clean the allocated memory.
 *
 * @param sw the #stellar_wind.
 */
void stellar_wind_clean(struct stellar_wind *sw) {

  interpolate_2d_free(&sw->integrated.ejected_energy_per_progenitor_mass);
  interpolate_2d_free(&sw->raw.ejected_energy);
}









/********************************************************************************
 *  Calculus method 
 *******************************************************************************/


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
            // message("Inside wind velocity, calculate parameter %i is = %f", i, calculate_b_parameter(log_Z,coeffs[i]));
            wind_velocity += calculate_b_parameter(log_Z,coeffs[i]) * pow(log_m,i);
        }
    /*  Else, the function will be only a linear relation between the initial mass and the function calculated for the limit mass */
    } else{
        float A0 = 0;
        float dLogA0 = 0;

        for (int i=0; i < 4; i++){
            A0 += calculate_b_parameter(log_Z,coeffs[i]) * pow(stellar_wind_x0,i);
            // message("Inside wind velocity, calculate parameter %i is = %f", i, calculate_b_parameter(log_Z,coeffs[i]));
            // The derivative take one less step in the loop
            if (i == 0){
                continue;
            }
            dLogA0 += i * (calculate_b_parameter(log_Z,coeffs[i]) * pow(stellar_wind_x0,(i - 1)));
        }

        wind_velocity += A0 + dLogA0 * (log_m - stellar_wind_x0);
    }

    // message("Inside wind velocity computation : V_infty = %f", wind_velocity);

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