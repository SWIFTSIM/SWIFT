/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 James Willis (james.s.willis@durham.ac.uk)
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

/* Forward declaration */
struct engine;

/* Structure for passing cosmological information to VELOCIraptor. */
struct cosmoinfo {

  /*! Current expansion factor of the Universe. (cosmology.a) */
  double atime;

  /*! Reduced Hubble constant (H0 / (100km/s/Mpc) (cosmology.h) */
  double littleh;

  /*! Matter density parameter (cosmology.Omega_m) */
  double Omega_m;

  /*! Baryon density parameter (cosmology.Omega_b) */
  double Omega_b;

  /*! Radiation constant density parameter (cosmology.Omega_lambda) */
  double Omega_Lambda;

  /*! Dark matter density parameter (cosmology.Omega_m - cosmology.Omega_b) */
  double Omega_cdm;

  /*! Dark-energy equation of state at the current time (cosmology.w)*/
  double w_de;
};

/* Structure for passing unit information to VELOCIraptor. */
struct unitinfo {

  /* Length conversion factor to kpc. */
  double lengthtokpc;

  /* Velocity conversion factor to km/s. */
  double velocitytokms;

  /* Mass conversion factor to solar masses. */
  double masstosolarmass;

  /* Potential conversion factor. */
  double energyperunitmass;

  /*! Newton's gravitationl constant (phys_const.const_newton_G)*/
  double gravity;

  /*! Hubble constant at the current redshift (cosmology.H) */
  double hubbleunit;
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
  struct cell_loc *cell_loc;

  /*! Top-level cell width. */
  double cellwidth[3];

  /*! Inverse of the top-level cell width. */
  double icellwidth[3];

  int icosmologicalsim;
};

/* VELOCIraptor wrapper functions. */
void velociraptor_init(struct engine *e);
void velociraptor_invoke(struct engine *e);

#endif /* SWIFT_VELOCIRAPTOR_INTERFACE_H */
