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
    

    /* Set cosmological constants. */
    cosmo_info.atime = e->cosmology->time; /* TODO: double check this. */
    cosmo_info.littleh = e->cosmology->h;
    cosmo_info.Omega_m = e->cosmology->Omega_m;
    cosmo_info.Omega_b = e->cosmology->Omega_b;
    cosmo_info.Omega_Lambda = e->cosmology->Omega_lambda;
    cosmo_info.Omega_cdm = e->cosmology->Omega_m - e->cosmology->Omega_b; /* TODO: double check this. */
    cosmo_info.w_de = e->cosmology->w_0; /* TODO: double check this. */

    /* Set unit conversions. */
    unit_info.lengthtokpc = 3.2404407e-22 * e->internal_units->UnitLength_in_cgs; /* 1kpc <=> 3.086e21cm */
    unit_info.velocitytokms = 1e-4 * parser_get_param_double(e->parameter_file, "InternalUnitSystem:UnitVelocity_in_cgs"); /* 1km/s <=> 1e4cm/s */
    unit_info.masstosolarmass = 5.0251256e-34 * e->internal_units->UnitMass_in_cgs; /* 1M_sol <=> 1.99e33g */
    unit_info.gravity = e->physical_constants->const_newton_G; /* G = 6.67408e-8 (cgs) */
    unit_info.hubbleunit = e->cosmology->H; /* TODO: double check this. */

    /* Set simulation information. */
    sim_info.period = 1.0; /* Placeholder. */
    sim_info.zoomhigresolutionmass = 1.0; /* Placeholder. */
    sim_info.interparticlespacing = 1.0; /* Placeholder. */
    sim_info.icosmologicalsim = 1; /* Placeholder. */

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
