/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Healpix C API header */
#ifdef HAVE_CHEALPIX
#include <chealpix.h>
#endif

#include <hdf5.h>

/* This object's header. */
#include "healpix_map.h"

/* Local headers. */
#include "engine.h"
#include "gravity_io.h"
#include "hydro_io.h"
#include "stars_io.h"
#include "black_holes_io.h"

/* Healpix map resolution */
#define NSIDE 256


void make_healpix_map(struct engine *e) {
#ifdef HAVE_CHEALPIX  

  /* Find the particle data */
  const struct space *s = e->s;
  const size_t nr_gparts = s->nr_gparts;
  const struct gpart *gparts = s->gparts;
  const struct part  *parts  = s->parts;
  const struct xpart *xparts = s->xparts;
  const struct spart *sparts = s->sparts;
  const struct bpart *bparts = s->bparts;

  /* Get the box size in each dimension */
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  /* Decide on centre and radius of spherical shell:
     For testing, we'll use radius=0.25*BoxSize[0] to 0.5*Boxsize[0]
     and centre on the centre of the box. */
  double centre[3];
  for(int i=0; i<3; i+=1)
    centre[i] = 0.5*dim[i];
  double rmin2 = pow(0.25*dim[0], 2.0);
  double rmax2 = pow(0.5*dim[0], 2.0);
  
  /* Allocate storage for the map */
  size_t npix = (size_t) nside2npix(NSIDE);
  double *map = malloc(sizeof(double)*npix);

  /* Zero out the map */
  for(size_t i=0; i<npix; i+=1)
    map[i] = 0.0;

  /* Loop over particles (should really use threadpool to parallelize)*/
  for(size_t i=0; i<nr_gparts; i+=1) {

    /* Get position of this particle and store in pos */
    double pos[3];
    switch (gparts[i].type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gparts[i].id_or_neg_offset];
      const struct xpart *xp = &xparts[-gparts[i].id_or_neg_offset];
      convert_part_pos(e, p, xp, pos);
    } break;
    case swift_type_stars: {
      const struct spart *sp = &sparts[-gparts[i].id_or_neg_offset];
      convert_spart_pos(e, sp, pos);
    } break;
    case swift_type_black_hole: {
      const struct bpart *bp = &bparts[-gparts[i].id_or_neg_offset];
      convert_bpart_pos(e, bp, pos);
    } break;
    case swift_type_dark_matter: {
      convert_gpart_pos(e, &gparts[i], pos);
    } break;
    case swift_type_dark_matter_background: {
      convert_gpart_pos(e, &gparts[i], pos);
    } break;
    default:
      error("Particle type not handled in healpix map.");
    }

    /* Get position relative to centre we chose */
    for(int j=0; j<3; j+=1)
      pos[j] -= centre[j];
    
    /* Find radius */
    double r2 = 0.0;
    for(int j=0; j<3; j+=1)
      r2 += (pos[j]-centre[j])*(pos[j]-centre[j]);

    /* Check particle is in shell */
    if(r2 > rmin2 && r2 < rmax2) {

      /* Convert position into a healpix pixel index */
      long ipring;
      vec2pix_ring((long) NSIDE, pos, &ipring);
      
      /* Accumulate mass to the map */
      map[ipring] += gparts[i].mass;
    }

    /* Next particle */
  }

  /* Whether this MPI rank should write the output file */
  int save_map = 1;
  
#ifdef WITH_MPI
  /* In MPI mode we need to sum the map over MPI ranks */
  if(e->nodeID==0) {
    MPI_Reduce(MPI_IN_PLACE, map, (int) npix, MPI_DOUBLE,
	       MPI_SUM, 0, MPI_COMM_WORLD);
    save_map = 1;
  } else {
    MPI_Reduce(map, map, (int) npix, MPI_DOUBLE,
	       MPI_SUM, 0, MPI_COMM_WORLD);
    save_map = 0;
  }
#endif
  
  /* Save the map to disk */
  if(save_map) {
    
    /* Determine name of the output file */
    const int snapnum = e->snapshot_output_count;
    char fname[100];
    sprintf(fname, "mass_map_%04d.fits", snapnum);

    /* Create the file */
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(file_id < 0)error("Unable to create HDF5 file for healpix map");

    /* Create the dataspace */
    hsize_t dims[1] = { (hsize_t) npix };
    hid_t space_id = H5Screate_simple(1, dims, NULL);

    /* Create the dataset */
    hid_t dset_id = H5Dcreate(file_id, "mass_map", H5T_NATIVE_DOUBLE, space_id,
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(dset_id<0)error("Unable to create dataset for healpix map");

    /* Write the data */
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, space_id, space_id, H5P_DEFAULT, map); 
    
    /* Tidy up */
    H5Dclose(dset_id);
    H5Sclose(space_id);
    H5Fclose(file_id);
  }
  
  /* Deallocate the map */
  free(map);

#else
  error("Can't make healpix map because healpix library was not found!");
#endif
}
