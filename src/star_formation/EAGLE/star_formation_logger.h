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
#include "cosmology.h"
#include "cell.h"
#include "hydro.h"
#include "part.h"

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

INLINE static void star_formation_add_progeny_SFH(struct star_formation_history* sf, 
    const struct star_formation_history* sfprogeny){
  /* Add the new stellar mass from the progeny */
  sf->new_stellar_mass = sf->new_stellar_mass + sfprogeny->new_stellar_mass;

  /* Increase amount of new stars formed */
  sf->N_stars = sf->N_stars + sfprogeny->N_stars;
}

/**
 * @brief Get the total star formation in this cell and add it to the star 
 * formation history struct
 *
 * @param c the cell of which we want to know the star formation
 * @param sf the star formation structure to which we want to add the star 
 * formation
 * @param cosmo the cosmology struct
 * @param with_cosmology if we run with cosmology
 */
INLINE static void star_formation_get_total_cell(const struct cell *c, struct star_formation_history *sf){
  /* Get the star formation history from the cell */
  struct star_formation_history *sfcell = c->stars.sfh;
  sf->new_stellar_mass += sfcell->new_stellar_mass;
  
  sf->N_stars += sfcell->new_stellar_mass;
}

/**
 * @brief Clear the total star formation in this cell 
 * 
 * @param c the cell of which we want to know the star formation
 */
INLINE static void star_formation_clear_total_cell(const struct cell *c){
  /* Get the star formation history from the cell */
  struct star_formation_history *sfcell = c->stars.sfh;
  sfcell->new_stellar_mass = 0.f;
  
  sfcell->N_stars = 0;
}

/**
 * @brief add the star formation to the parent cell 
 * 
 * @param c the cell for which we want to add the star formation
 * @param sf the combined star formation history of the progeny
 */
INLINE static void star_formation_add_to_parent_cell(const struct cell *c, struct star_formation_history *sf){
  /* Get the star formation history from the cell */
  struct star_formation_history *sfcell = c->stars.sfh;
  sfcell->new_stellar_mass = sf->new_stellar_mass;
  
  sfcell->N_stars = sf->N_stars;
}

/** 
 * @brief Initialize the star formation history structure
 *
 * @param The pointer to the star formation history structure
 * */
INLINE static void star_formation_init_SFH_engine(struct star_formation_history *sfh){
  sfh->new_stellar_mass = 0.f;

  sfh->N_stars = 0;
}

#endif /* SWIFT_SCHAYE_STARFORMATION_LOGGER_H */
