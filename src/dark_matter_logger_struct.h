/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Camila Correa (camila.correa@uva.nl)
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
 *******************************************************************************/
#ifndef SWIFT_DARK_MATTER_LOGGER_STRUCT_H
#define SWIFT_DARK_MATTER_LOGGER_STRUCT_H

/* Local includes */

/* SIDM history struct */
struct sidm_history {
    
    /*! Total kinetic energy in the simulation at the beginning and end of time-step */
    double K_before;
    double K_after;
    
    /*! Number of SIDM events */
    int num_kicks;
};

/* SIDM history struct for the engine
   Allows to integrate in time some values */
struct sidm_history_accumulator {

    /*! Total kinetic energy in the simulation at the beginning and end of time-step */
    double K_before;
    double K_after;
    
    /*! Number of SIDM events */
    int num_kicks;

    /* Number of activate particles per timestep */
    int n_parts_active;


};

#endif
