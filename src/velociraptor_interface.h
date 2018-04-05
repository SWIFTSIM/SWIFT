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
#ifndef SWIFT_VELOCIRAPTOR_INTERFACE_H
#define SWIFT_VELOCIRAPTOR_INTERFACE_H

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include "engine.h"

/* Structure for passing cosmological information to VELOCIraptor. */
struct cosmoinfo {
    double atime, littleh, Omega_m, Omega_b, Omega_Lambda, Omega_cdm, w_de;
};

/* Structure for passing unit information to VELOCIraptor. */
struct unitinfo {
    double lengthtokpc,velocitytokms,masstosolarmass,gravity,hubbleunit;
};

/* Structure to hold the location of a top-level cell. */
struct cell_loc {
    
    /* Coordinates x,y,z */
    double loc[3];

};

/* Structure for passing simulation information to VELOCIraptor. */
struct siminfo {
    double period, zoomhigresolutionmass, interparticlespacing, spacedimension[3];
    
    /* Number of top-cells. */
    int numcells;

    /*! Locations of top-level cells. */
    struct cell_loc *cellloc;
    
    /*! Top-level cell width. */
    double cellwidth[3];
    
    /*! Inverse of the top-level cell width. */
    double icellwidth[3];

    int icosmologicalsim;
};

/* VELOCIraptor interface. */
void InitVelociraptor(char* config_name, char* output_name, struct cosmoinfo cosmo_info, struct unitinfo unit_info, struct siminfo sim_info);
void InvokeVelociraptor(const int num_gravity_parts, struct gpart *gravity_parts, const int *cell_node_ids, char* output_name);

/* VELOCIraptor wrapper functions. */
void velociraptor_init(struct engine *e);
void velociraptor_invoke(struct engine *e);

#endif /* SWIFT_VELOCIRAPTOR_INTERFACE_H */
