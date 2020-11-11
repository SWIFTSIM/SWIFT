/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 John Helly (j.c.helly@durham.ac.uk)
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

/* Some standard headers. */
#include <stdlib.h>
#include <stdio.h>

/* This object's header. */
#include "hbt_interface.h"

/* Local includes. */
#include "engine.h"

#ifdef HAVE_HBT
#include "libHBT.h"


/* Data needed by the hbt_get_particle function */
struct hbt_particle_data
{
  size_t nr_gparts;
  struct gpart  *gparts;
  struct part   *parts;
  struct xpart  *xparts;
  struct spart  *sparts;
  struct bpart  *bparts;
};


/**
 * @brief Callback function to pass particles to HBT
 *
 * @param data Pointer to data needed by the function
 * @param index Index of the particle to look up
 * @param fofid Returns FoF group ID of the particle
 * @param type Returns particle type (HBT's enum ParticleType_t)
 * @param pos Returns coordinates of the particle
 * @param vel Returns velocity of the particle
 * @param id Returns ID of the particle
 * @param mass Returns mass of the particle
 * @param u Returns internal energy of the particle
 *
 */
void hbt_get_particle(void *data, size_t index, HBTInt *fofid,
                      HBTInt *type, HBTReal *pos, HBTReal *vel,
                      HBTInt *id, HBTReal *mass, HBTReal *u) {
  
  /* Unpack input data */
  const struct hbt_particle_data *pdata = (struct hbt_particle_data *) data;
  const size_t nr_gparts = pdata->nr_gparts;
  const struct *gparts = pdata->gparts;
  const struct *parts  = pdata->parts;
  const struct *xparts = pdata->xparts;
  const struct *sparts = pdata->sparts;
  const struct *bparts = pdata->bparts;

  /* Look up particle properties */
  const size_t i = index;
  switch (gparts[i].type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gparts[i].id_or_neg_offset];
    const struct xpart *xp = &xparts[-gparts[i].id_or_neg_offset];
    convert_part_pos(e, p, xp, pos);
    convert_part_vel(e, p, xp, vel);
    *id   = p->id;
    *mass = gravity_get_mass(&gparts[i]);
    *type = 0;
    *fofid = gparts[i]->fof_data.group_id;
    *u = hydro_get_drifted_physical_internal_energy(p, cosmo);
  } break;
  case swift_type_stars: {
    const struct spart *sp = &sparts[-gparts[i].id_or_neg_offset];
    convert_spart_pos(e, sp, pos);
    convert_spart_vel(e, sp, vel);
    *id   = sp->id;
    *mass = gravity_get_mass(&gparts[i]);
    *type = 4;
    *fofid = gparts[i]->fof_data.group_id;
    *u = 0.0;
  } break;
  case swift_type_black_hole: {
    const struct bpart *bp = &bparts[-gparts[i].id_or_neg_offset];
    convert_bpart_pos(e, bp, pos);
    convert_bpart_vel(e, bp, vel);
    *id   = bp->id;
    *mass = gravity_get_mass(&gparts[i]);
    *type = 5;
    *fofid = gparts[i]->fof_data.group_id;
    *u = 0.0;
  } break;
  case swift_type_dark_matter: {
    convert_gpart_pos(e, &gparts[i], pos);
    convert_gpart_vel(e, &gparts[i], vel);
    *id   = gparts[i].id_or_neg_offset;
    *mass = gravity_get_mass(&gparts[i]);
    *type = 1;
    *fofid = gparts[i]->fof_data.group_id;
    *u = 0.0;
  } break;
  case swift_type_dark_matter_background: {
    convert_gpart_pos(e, &gparts[i], pos);
    convert_gpart_vel(e, &gparts[i], vel);
    *id   = gparts[i].id_or_neg_offset;
    *mass = gravity_get_mass(&gparts[i]);
    *type = 2;
    *fofid = gparts[i]->fof_data.group_id;
    *u = 0.0;
  } break;
  default:
    error("Particle type not handled by HBT.");
  }
}

/**
 * @brief Initialise HBT library
 */
void hbt_init(void) {

  /* Internal SWIFT units */
  const struct unit_system *swift_us = e->internal_units;

  /* Retrieve information needed by HBT */
  const double h              = e->cosmology->h;
  const char * config_file    = e->hbt_config_file_name; // TODO: add this to Swift!
  const int    num_threads    = e->nr_threads;
  const double omega_m0       = e->cosmology->Omega_m;
  const double omega_lambda0  = e->cosmology->Omega_lambda;
  const double MassInMsunh    = units_cgs_conversion_factor(swift_us, UNIT_CONV_MASS) /
    cgs_pc.const_solar_mass * h;
  const double LengthInMpch   = units_cgs_conversion_factor(swift_us, UNIT_CONV_LENGTH) /
    (1000000. * cgs_pc.const_parsec) * h;
  const double VelInKmS       = units_cgs_conversion_factor(swift_us, UNIT_CONV_VELOCITY) / 1.0e5;
  const long long NullGroupId = e->fof_properties->group_id_default;  
  
  /* Initialise HBT */
  libhbt_init(config_file, num_threads, omega_m0, omega_lambda0,
              MassInMsunh, LengthInMpch, VelInKmS, NullGroupId);
}


/**
 * @brief Invoke the HBT halo finder
 *
 * @param output_nr Index of the current output
 *
 */
void hbt_invoke(int output_nr) {

  /* Find the particle data */
  const struct space *s = e->s;
  const size_t nr_gparts = s->nr_gparts;
  struct hbt_particle_data hbt_data;
  hbt_data.gparts = e->s->gparts;
  hbt_data.parts  = e->s->parts;
  hbt_data.xparts = e->s->xparts;
  hbt_data.sparts = e->s->sparts;
  hbt_data.bparts = e->s->bparts;

  /* Find current scale factor */
  const double scalefactor = e->cosmology->a;

  /* Call HBT */
  libhbt_invoke_hbt(output_nr, scalefactor,
                    &hbt_data, nr_gparts, &hbt_get_particle);
}


/**
 * @brief Free resources used by HBT library
 */
void hbt_free(void) {
  libhbt_free();
}
