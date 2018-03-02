/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *                    Angus Lepper (angus.lepper@ed.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "velociraptor_interface.h"

/**
 * @brief Initialise VELOCIraptor with input and output file names along with cosmological info needed to run.
 *
 * @param e The #engine.
 *
 */
void init_velociraptor(struct engine *e) {
    struct cosmoinfo cosmo_info;
    struct unitinfo unit_info;
    struct siminfo sim_info;
    
    cosmo_info.atime = 1.0;
    unit_info.lengthtokpc = 1.0;
    sim_info.period = 1.0;

    InitVelociraptor("stf_input.cfg", "stf_ouput.out", cosmo_info, unit_info, sim_info);
}

/**
 * @brief Run VELOCIraptor with current particle data.
 *
 * @param e The #engine.
 *
 */
void invoke_velociraptor(struct engine *e) {

    struct gpart *gparts = e->s->gparts;
    const int nr_gparts = e->s->nr_gparts;

    InvokeVelociraptor(nr_gparts, gparts);
}
