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
void velociraptor_init(struct engine *e) {
    struct space *s = e->s;
    struct cosmoinfo cosmo_info;
    struct unitinfo unit_info;
    struct siminfo sim_info;
    
    struct unit_system vel_us;
    struct phys_const vel_const;

    /* Initialize velociraptor unit system and constants */
    units_init(&vel_us, e->parameter_file, "VelociraptorUnitSystem");
    phys_const_init(&vel_us, e->parameter_file, &vel_const);

    /* Set cosmological constants. */
    cosmo_info.atime = e->cosmology->a;
    cosmo_info.littleh = e->cosmology->h;
    cosmo_info.Omega_m = e->cosmology->Omega_m;
    cosmo_info.Omega_b = e->cosmology->Omega_b;
    cosmo_info.Omega_Lambda = e->cosmology->Omega_lambda;
    cosmo_info.Omega_cdm = e->cosmology->Omega_m - e->cosmology->Omega_b;
    cosmo_info.w_de = e->cosmology->w_0;

    message("Scale factor: %e", cosmo_info.atime);
    message("Little h: %e", cosmo_info.littleh);
    message("Omega_m: %e", cosmo_info.Omega_m);
    message("Omega_b: %e", cosmo_info.Omega_b);
    message("Omega_Lambda: %e", cosmo_info.Omega_Lambda);
    message("Omega_cdm: %e", cosmo_info.Omega_cdm);
    message("w_de: %e", cosmo_info.w_de);

    /* Set unit conversions. */
    unit_info.lengthtokpc = units_conversion_factor(e->internal_units, &vel_us, UNIT_CONV_LENGTH); /* 1kpc <=> 3.086e21cm */
    unit_info.velocitytokms = units_conversion_factor(e->internal_units, &vel_us, UNIT_CONV_SPEED); /* 1km/s <=> 1e5cm/s */
    unit_info.masstosolarmass = units_conversion_factor(e->internal_units, &vel_us, UNIT_CONV_MASS); /* 1M_sol <=> 1.99e33g */
    unit_info.gravity = vel_const.const_newton_G; /* TODO: G = 6.67408e-8 (cgs) */
    unit_info.hubbleunit = e->cosmology->H; /* TODO: double check this. */

    message("Length conversion factor: %e", unit_info.lengthtokpc);
    message("Velocity conversion factor: %e", unit_info.velocitytokms);
    message("Mass conversion factor: %e", unit_info.masstosolarmass);
    message("G: %e", unit_info.gravity);
    message("H: %e", unit_info.hubbleunit);

    const int total_nr_gparts = e->total_nr_gparts;
    
    /* Set simulation information. */
    if(e->s->periodic) {
        sim_info.period = unit_info.lengthtokpc * s->dim[0]; /* Physical size of box in VELOCIraptor units (kpc). */
    }
    else sim_info.period = 0.0;
    sim_info.zoomhigresolutionmass = -1.0; /* Placeholder. */
    sim_info.interparticlespacing = sim_info.period / pow(total_nr_gparts, 1./3.);
    sim_info.icosmologicalsim = (e->policy & engine_policy_cosmology);
    sim_info.spacedimension[0] = unit_info.lengthtokpc * s->dim[0];
    sim_info.spacedimension[1] = unit_info.lengthtokpc * s->dim[1];
    sim_info.spacedimension[2] = unit_info.lengthtokpc * s->dim[2];
    sim_info.numcells = s->nr_cells;
 
    sim_info.cellwidth[0] = unit_info.lengthtokpc * s->cells_top[0].width[0];
    sim_info.cellwidth[1] = unit_info.lengthtokpc * s->cells_top[0].width[1];
    sim_info.cellwidth[2] = unit_info.lengthtokpc * s->cells_top[0].width[2];

    sim_info.icellwidth[0] = s->iwidth[0] / unit_info.lengthtokpc;
    sim_info.icellwidth[1] = s->iwidth[1] / unit_info.lengthtokpc;
    sim_info.icellwidth[2] = s->iwidth[2] / unit_info.lengthtokpc;
   
    /* Allocate and populate top-level cell locations. */
    /* JSW TODO: Remember to free at the end of the simulation. */
    if (posix_memalign((void **)&(sim_info.cellloc), 32,
                       s->nr_cells * sizeof(struct cell_loc)) != 0)
        error("Failed to allocate top-level cell locations for VELOCIraptor.");


    //double cell_loc_min[3] = {sim_info.spacedimension[0] - sim_info.cellwidth[0],
    //                          sim_info.spacedimension[1] - sim_info.cellwidth[1],
    //                          sim_info.spacedimension[2] - sim_info.cellwidth[2]};
    //double cell_loc_max[3] = {0.0, 0.0, 0.0};

    for(int i=0; i<s->nr_cells; i++) {
        sim_info.cellloc[i].loc[0] = unit_info.lengthtokpc * s->cells_top[i].loc[0];
        sim_info.cellloc[i].loc[1] = unit_info.lengthtokpc * s->cells_top[i].loc[1];
        sim_info.cellloc[i].loc[2] = unit_info.lengthtokpc * s->cells_top[i].loc[2];

        //for(int k=0; k<3; k++) {
        //    cell_loc_min[k] = min(cell_loc_min[k], sim_info.cellloc[i].loc[k]);
        //    cell_loc_max[k] = max(cell_loc_max[k], sim_info.cellloc[i].loc[k] + sim_info.cellwidth[k]);
        //}
    }

    //FILE *file = NULL;

    ///* Name of output file. */
    //char fname[200];
    //sprintf(fname, "cell_locs_rank_%d.dat", e->nodeID);
    //file = fopen(fname, "w");

    ///* Header. */
    //fprintf(file, "# %6s %6s %6s %6s %6s %6s %6s\n", "x", "y", "z", "xw", "yw",
    //        "zw", "rank");

    ///* Output */
    //for (int i = 0; i < s->nr_cells; i++) {
    //    struct cell *c = &s->cells_top[i];
    //    fprintf(file, "  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6d\n", sim_info.cellloc[i].loc[0] / unit_info.lengthtokpc,
    //            sim_info.cellloc[i].loc[1]/ unit_info.lengthtokpc, sim_info.cellloc[i].loc[2]/ unit_info.lengthtokpc, sim_info.cellwidth[0]/ unit_info.lengthtokpc, sim_info.cellwidth[1]/ unit_info.lengthtokpc, sim_info.cellwidth[2]/ unit_info.lengthtokpc,
    //            c->nodeID);
    //}

    //fclose(file);

    message("Period: %e", sim_info.period);
    message("Zoom high res mass: %e", sim_info.zoomhigresolutionmass);
    message("Inter-particle spacing: %e", sim_info.interparticlespacing);
    message("Cosmological: %d", sim_info.icosmologicalsim);
    message("Space dimensions: (%e,%e,%e)", sim_info.spacedimension[0], sim_info.spacedimension[1], sim_info.spacedimension[2]);
    message("No. of top-level cells: %d", sim_info.numcells);
    //message("Local top-level cell locations range: (%e,%e,%e) -> (%e,%e,%e)", cell_loc_min[0], cell_loc_min[1], cell_loc_min[2], cell_loc_max[0], cell_loc_max[1], cell_loc_max[2]);
    message("Top-level cell locations range: (%e,%e,%e) -> (%e,%e,%e)", sim_info.cellloc[0].loc[0], sim_info.cellloc[0].loc[1], sim_info.cellloc[0].loc[2], sim_info.cellloc[sim_info.numcells - 1].loc[0], sim_info.cellloc[sim_info.numcells - 1].loc[1], sim_info.cellloc[sim_info.numcells - 1].loc[2]);

    InitVelociraptor("stf_input.cfg", "stf_output.out", cosmo_info, unit_info, sim_info);

    /* Free cell locations after VELOCIraptor has copied them. */
    //free(sim_info.cellloc);
}

/**
 * @brief Run VELOCIraptor with current particle data.
 *
 * @param e The #engine.
 *
 */
void velociraptor_invoke(struct engine *e) {

    struct space *s = e->s;
    struct gpart *gparts = s->gparts;
    const int nr_gparts = s->nr_gparts;
    const int nr_cells = s->nr_cells;
    int *cell_node_ids;
    
    /* Allocate and populate array of cell node IDs. */
    /* JSW TODO: Remember to free at the end of the simulation. */
    if (posix_memalign((void **)&cell_node_ids, 32,
                       nr_cells * sizeof(int)) != 0)
        error("Failed to allocate list of cells node IDs for VELOCIraptor.");
    
    for(int i=0; i<nr_cells; i++) cell_node_ids[i] = s->cells_top[i].nodeID;    

    message("MPI rank %d sending %d gparts to VELOCIraptor.", e->nodeID, nr_gparts);
    
    //for(int i=0; i<nr_gparts; i++) message("Potential: %f", gparts[i].potential);

    InvokeVelociraptor(nr_gparts, gparts, cell_node_ids);
    
    /* Free cell node ids after VELOCIraptor has copied them. */
    //free(cell_node_ids);
}
