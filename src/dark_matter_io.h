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
#ifndef SWIFT_DARK_MATTER_IO_H
#define SWIFT_DARK_MATTER_IO_H

#include "io_properties.h"



/**
 * @brief Specifies which g-particle fields to read from a dataset
 *
 * @param gparts The g-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void darkmatter_read_as_dmparticles(struct dmpart* dmparts,
                                             struct io_props* list,
                                             int* num_fields) {
    
    /* Say how much we want to read */
    *num_fields = 6;
    
    /* List what we want to read */
    list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                  UNIT_CONV_LENGTH, dmparts, x);
    list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                  UNIT_CONV_SPEED, dmparts, v_full);
    list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                  dmparts, mass);
    list[3] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                  UNIT_CONV_NO_UNITS, dmparts, id_or_neg_offset);
    list[4] = io_make_input_field("SmoothingLength", FLOAT, 1, OPTIONAL,
                                  UNIT_CONV_LENGTH, dmparts, h);
    list[5] = io_make_input_field("Density", FLOAT, 1, OPTIONAL,
                                  UNIT_CONV_DENSITY, dmparts, rho);

}

/**
 * @brief Specifies which dm particle fields to write to a dataset
 *
 * @param gparts The dm particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static int sidm_write_dmparts(const struct dmpart* dmparts,
                                              struct io_props* list) {
    
    
    /* List what we want to write */
    list[0] = io_make_output_field("SIDM_events", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
                                   dmparts, sidm_data.number_of_sidm_events, "Total number of DM-DM collisions the particle has had");
    
    list[1] = io_make_output_field("SIDM_search_radius", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, dmparts, h,
                                   "Co-moving smoothing lengths (FWHM of the kernel) of the DM particles");
    
    /*list[2] = io_make_output_field("Densities", FLOAT, 1, UNIT_CONV_DENSITY, -3.f, dmparts, rho,
                                   "Co-moving mass densities of the particles");*/

    list[2] = io_make_output_field("DM_ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, dmparts,
                                   id_or_neg_offset, "Unique ID of the particles");
    
    list[3] = io_make_output_field("SIDM_probability", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, dmparts,
                                   sidm_probability, "Average probability of particle scattering with another (internal units)");
    
    list[4] = io_make_output_field("Number_of_neighbours", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, dmparts,
                                   num_neighbours, "Number of neighbours");
    
    list[5] = io_make_output_field("Time_step_size", FLOAT, 1, UNIT_CONV_TIME, 1.f, dmparts,
                                   time_step_size, "Average probability of particle scattering with another (internal units)");
    
    list[6] = io_make_output_field("Max_SIDM_events_per_timestep", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
                                   dmparts, sidm_data.max_sidm_events_per_timestep, "Maximum number of DM-DM collisions the particle has had in a single timestep");


    return 7;
    
}


#endif /* SWIFT_DARK_MATTER_IO_H */
