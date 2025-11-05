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
#include <config.h>

/* Local includes */
#include "black_holes.h"
#include "engine.h"
#include "gravity.h"
#include "hydro.h"
#include "lightcone/lightcone_map.h"
#include "neutrino.h"
#include "part.h"
#include "star_formation.h"
#include "stars.h"

/* This object's header */
#include "lightcone/lightcone_map_types.h"

/* Required for the xrays */
#include "extra_io.h"
#include "io_properties.h"

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_gas_only(int ptype) {

  switch (ptype) {
    case swift_type_gas:
      return 1;
    default:
      return 0;
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_total_mass_type_contributes(int ptype) {

  switch (ptype) {
    case swift_type_gas:
    case swift_type_stars:
    case swift_type_black_hole:
    case swift_type_dark_matter:
    case swift_type_dark_matter_background:
    case swift_type_neutrino:
      return 1;
    default:
      return 0;
  }
}

/**
 * @brief Make a healpix map of projected mass in each pixel
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_total_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  /* const struct xpart *xparts = s->xparts; */ /* Currently not used */
  const struct spart *sparts = s->sparts;
  const struct bpart *bparts = s->bparts;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      return hydro_get_mass(p);
    } break;
    case swift_type_stars: {
      const struct spart *sp = &sparts[-gp->id_or_neg_offset];
      return sp->mass;
    } break;
    case swift_type_black_hole: {
      const struct bpart *bp = &bparts[-gp->id_or_neg_offset];
      return bp->mass;
    } break;
    case swift_type_dark_matter:
    case swift_type_dark_matter_background: {
      return gp->mass;
    } break;
    case swift_type_neutrino: {
      struct neutrino_model nu_model;
      bzero(&nu_model, sizeof(struct neutrino_model));
      if (e->neutrino_properties->use_delta_f_mesh_only)
        gather_neutrino_consts(e->s, &nu_model);
      double weight = 1.0;
      gpart_neutrino_weight_mesh_only(gp, &nu_model, &weight);
      return gp->mass * weight;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/**
 * Return background neutrino density to add to the total mass map
 *
 * The neutrino particles trace density perturbations relative
 * to a constant background so we need to add in the background
 * to get the total mass.
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param map The lightcone map
 *
 */
double lightcone_map_total_mass_baseline_value(
    const struct cosmology *c, const struct lightcone_props *lightcone_props,
    const struct lightcone_map *map) {
  return lightcone_map_neutrino_baseline_value(c, lightcone_props, map);
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_gas_mass_type_contributes(int ptype) {

  switch (ptype) {
    case swift_type_gas:
      return 1;
    default:
      return 0;
  }
}

/**
 * @brief Make a healpix map of projected gas mass in each pixel
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_gas_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      return hydro_get_mass(p);
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_dark_matter_mass_type_contributes(int ptype) {

  switch (ptype) {
    case swift_type_dark_matter:
    case swift_type_dark_matter_background:
      return 1;
    default:
      return 0;
  }
}

/**
 * @brief Make a healpix map of projected dark matter mass in each pixel
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_dark_matter_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {
  switch (gp->type) {
    case swift_type_dark_matter:
    case swift_type_dark_matter_background: {
      return gp->mass;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_stellar_mass_type_contributes(int ptype) {

  switch (ptype) {
    case swift_type_stars:
      return 1;
    default:
      return 0;
  }
}

/**
 * @brief Make a healpix map of stellar mass in each pixel
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_stellar_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct spart *sparts = s->sparts;

  switch (gp->type) {
    case swift_type_stars: {
      const struct spart *sp = &sparts[-gp->id_or_neg_offset];
      return sp->mass;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_black_hole_mass_type_contributes(int ptype) {

  switch (ptype) {
    case swift_type_black_hole:
      return 1;
    default:
      return 0;
  }
}

/**
 * @brief Make a healpix map of black hole mass in each pixel
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_black_hole_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct bpart *bparts = s->bparts;

  switch (gp->type) {
    case swift_type_black_hole: {
      const struct bpart *bp = &bparts[-gp->id_or_neg_offset];
      return bp->mass;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_sfr_type_contributes(int ptype) {

  switch (ptype) {
    case swift_type_gas:
      return 1;
    default:
      return 0;
  }
}

/**
 * @brief Make a healpix map of star formation rate in each pixel
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_sfr_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart *xp = &xparts[-gp->id_or_neg_offset];
      return star_formation_get_SFR(p, xp);
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}
