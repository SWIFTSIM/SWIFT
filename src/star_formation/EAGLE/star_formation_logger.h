/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_SCHAYE_STARFORMATION_LOGGER_H
#define SWIFT_SCHAYE_STARFORMATION_LOGGER_H

/* Some standard headers */
#include <stdlib.h>

/* Local includes */
#include "adiabatic_index.h"
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "equation_of_state.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "stars.h"
#include "units.h"

/* Starformation history struct */
struct star_formation_history {
  /*! Numb of stars */
  unsigned long int N_stars;

  /*! Total new stellar mass */
  float new_stellar_mass;

  /*! Time union */
  union {
    /*! Time */
    float time;

    /*! Scale factor */
    float scale_factor;
  };
};

INLINE static void starformation_update_SFH(struct spart* sp, struct star_formation_history* sf, const struct cosmology* cosmo, 
    const int with_cosmology){ 
  /* Add mass of created sparticle to the total stellar mass in this cell*/
  sf->new_stellar_mass = sf->new_stellar_mass + sp->mass;

  /* Increase the counter */
  sf->N_stars++;

}

INLINE static void starformation_init_SFH(struct star_formation_history* sf, const struct cosmology* cosmo, 
    const int with_cosmology){ 
  /* Initialize the stellar mass to zero*/
  sf->new_stellar_mass = 0.f;

  /* Initialize the counter at zero */
  sf->N_stars=0;

}

INLINE static void starformation_add_progeny_SFH(struct star_formation_history* sf, 
    const struct star_formation_history* sfprogeny, const struct cosmology* cosmo, 
    const int with_cosmology){
  /* Add the new stellar mass from the progeny */
  sf->new_stellar_mass = sf->new_stellar_mass + sfprogeny->new_stellar_mass;

  /* Increase amount of new stars formed */
  sf->N_stars = sf->N_stars + sfprogeny->N_stars;
}

INLINE static void star_formation_get_total_cell(const struct cell *c, struct star_formation_history *sf, 
    const struct cosmology* cosmo, const int with_cosmology){
  /* Get the star formation history from the cell */
  struct star_formation_history *sfcell = c->stars.sfh;
  sf->new_stellar_mass += sfcell->new_stellar_mass;
  
  sf->N_stars += sfcell->new_stellar_mass;
}

#endif /* SWIFT_SCHAYE_STARFORMATION_LOGGER_H */
