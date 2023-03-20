/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2022 Will Roper (w.roper@sussex.ac.uk)
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
#ifndef SWIFT_HALO_FINDER_HALO_H
#define SWIFT_HALO_FINDER_HALO_H

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "cosmology.h"
#include "part_type.h"
#include "physical_constants.h"
#include "space.h"

/* Avoid cyclic inclusions */
struct cell;
struct gpart;
struct fof_props;

/**
 * @brief What kind of halo is this?
 *
 * 0 = Not a halo.
 * 1 = A 3D FOF group.
 * 1 = A 6D Host halo.
 * 2 = A 6D Subhalo
 */
enum halo_types { no_halo, fof_group, host_halo, sub_halo };

/**
 * @brief Properties of a halo.
 */
struct halo_props {

  /*! Is this halo real? */
  int is_real;

  /*! The total number of particles in the halo. */
  size_t npart_tot;

  /*! Number of particles in halo per species. */
  size_t npart[swift_type_count];

  /*! Mass. */
  double mass_tot;

  /*! Mass per species. */
  double mass[swift_type_count];

  /*! The groups mass weighted bulk velocity vector. */
  double velocity[3];

  /*! The groups mass weighted bulk velocity vector. */
  double velocity_mag;

  /*! The mean physical acceleration vector of the halo. */
  double a_phys[3];

  /*! Start pointer into halo particle property arrays. */
  size_t part_start_index;

  /*! Unisolated gravitaitonal binding energy. */
  double grav_pot;

  /*! Boosted frame potential. (https://arxiv.org/pdf/2107.13008.pdf) */
  double grav_boost;

  /*! Kinetic energy. */
  double kinetic_nrg;

  /*! First particles position. */
  double first_position[3];

  /*! Centre of potential. */
  double centre_of_potential[3];

  /*! Centre of mass. */
  double centre_of_mass[3];

  /*! Most bound particle. */
  struct gpart *most_bound_gpart;

  /*! The extent of the halo in each dimension (used when defining tasks). */
  double extent[6];

  /*! The width of the halo in each dimension (used when defining tasks). */
  double width[3];
  
};

  
/**
 * @brief Halo substructure properties.
 */
struct halo_substructures {

  /*! Number of substructures. */
  int nsubs;

  /*! Pointers to substructures. */
  struct halo *substructures;
  
};

/**
 * @brief Temporal linking halo properties.
 */
struct halo_temporal_links {
  
  /*! Pointers to progenitor halos. */
  struct halo *progs;
  
  /*! Pointers to descendant halos. */
  struct halo *descs;
  
};

/**
 * @brief A halo used to store all information for a single halo.
 */
struct halo {

  /*! Halo ID. */
  size_t halo_id;

  /*! The total number of particles in the halo. */
  size_t size;

  /*! Parent halo pointer (used to indicate the halo this halo was extracted
   *  from). */
  struct halo *parent;

  /*! Parent halo pointer (used to indicate the parent in the
   *  FOF->host->subhalo heirarchy). */
  struct halo *host;

  /*! Halo type flag. */
  enum halo_types type;

  /*! The current velocity space linking length coefficient. */
  double alpha_vel;

  /*! Halo properties (only initialised if halo is found to be real). */
  struct halo_props *props;

  /*! Halo substructures. */
  struct halo_substructures substructure;

  /*! Halo progenitor and descendant data. */
  struct halo_temporal_links temporal_links;
  
};

/* Prototypes */
void fof_to_halo_finder_mapper(void *map_data, int num_elements,
                               void *extra_data);
void halo_search_tree(struct fof_props *props,
                      const struct phys_const *constants,
                      const struct cosmology *cosmo, struct space *s,
                      const int dump_results, const int dump_debug_results);

#endif /* SWIFT_HALO_FINDER_HALO_H */
