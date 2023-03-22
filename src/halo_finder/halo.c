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

/* Config parameters. */
#include <config.h>

/* Includes. */
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

/* Local includes. */
#include "cell.h"
#include "engine.h"
#include "fof.h"
#include "halo_finder/halo.h"
#include "hashmap.h"

#define FOF_COMPRESS_PATHS_MIN_LENGTH (2)

#ifdef WITH_MPI

/* MPI types used for communications */
MPI_Datatype fof_mpi_type;
MPI_Datatype group_length_mpi_type;
MPI_Datatype fof_final_index_type;
MPI_Datatype fof_final_mass_type;

/*! Offset between the first particle on this MPI rank and the first particle in
 * the global order */
size_t node_offset;
#endif

/**
 * @brief Finds the local root ID of the group a particle exists in.
 *
 * We follow the group_index array until reaching the root of the group.
 *
 * Also performs path compression if the path is long.
 *
 * @param i The index of the particle.
 * @param group_index Array of group root indices.
 */
__attribute__((always_inline)) INLINE static size_t fof_find(
    const size_t i, size_t *group_index) {

  size_t root = i;
  int tree_depth = 0;

  while (root != group_index[root]) {
#ifdef PATH_HALVING
    atomic_cas(&group_index[root], group_index[root],
               group_index[group_index[root]]);
#endif
    root = group_index[root];
    tree_depth++;
  }

  /* Only perform path compression on trees with a depth of
   * FOF_COMPRESS_PATHS_MIN_LENGTH or higher. */
  if (tree_depth >= FOF_COMPRESS_PATHS_MIN_LENGTH)
    atomic_cas(&group_index[i], group_index[i], root);

  return root;
}

/**
 * @brief   Finds the local root ID of the group a particle exists in
 * when group_index contains globally unique identifiers -
 * i.e. we stop *before* we advance to a foreign root.
 *
 * Here we assume that the input i is a local index and we
 * return the local index of the root.
 *
 * @param i Index of the particle.
 * @param nr_gparts The number of g-particles on this node.
 * @param group_index Array of group root indices.
 */
__attribute__((always_inline)) INLINE static size_t fof_find_local(
    const size_t i, const size_t nr_gparts, const size_t *group_index) {
#ifdef WITH_MPI
  size_t root = node_offset + i;

  while ((group_index[root - node_offset] != root) &&
         (group_index[root - node_offset] >= node_offset) &&
         (group_index[root - node_offset] < node_offset + nr_gparts)) {
    root = group_index[root - node_offset];
  }

  return root - node_offset;
#else
  size_t root = i;

  while ((group_index[root] != root) && (group_index[root] < nr_gparts)) {
    root = group_index[root];
  }

  return root;
#endif
}

/**
 * @brief Mapper function to calculate the group sizes.
 *
 * @param map_data An array of #gpart%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void fof_to_halo_finder_mapper(void *map_data, int num_elements,
                                void *extra_data) {

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  struct gpart *gparts = (struct gpart *)map_data;

  /* Get the indexes and default id. */
  size_t *group_index = s->e->fof_properties->group_index;
  size_t group_id_default = s->e->fof_properties->group_id_default;

  /* Offset into gparts array. */
  ptrdiff_t gparts_offset = (ptrdiff_t)(gparts - s->gparts);
  size_t *const group_index_offset = group_index + gparts_offset;

  /* Loop over particles and find which cells are in range of each other to
   * perform the FOF search. */
  for (int ind = 0; ind < num_elements; ind++) {

    hashmap_key_t root =
        (hashmap_key_t)fof_find(group_index_offset[ind], group_index);
    const size_t gpart_index = gparts_offset + ind;

    /* Skip the root, it's already attatched. */
    if (root == gpart_index || root == group_id_default) continue;

    /* Attatch the root group to this particle. */
    gparts[ind].fof_data.group = gparts[root - gparts_offset].fof_data.group;
  }
}

/**
 * @brief Mapper function to set the initial group IDs.
 *
 * @param map_data The array of #gpart%s.
 * @param num_elements Chunk size.
 * @param extra_data FOF properties.
 */
void set_initial_halo_id_mapper(void *map_data, int num_elements,
                                void *extra_data) {

  /* Unpack the information */
  struct gpart *gparts = (struct gpart *)map_data;
  struct space *s = (struct space *)extra_data;
  struct fof_props *props = s->e->fof_properties;

  /* Extract FOF properties */
  const size_t group_id_default = props->group_id_default;
  enum halo_types current_level = props->current_level;

  /* Local copy of the halo arrays. */
  struct halo *groups;
  if (current_level == fof_group) {
    groups = props->groups;
  } else if (current_level == host_halo) {
    groups = props->hosts;
  } else {
    groups = props->subhalos;
  }

  /* Offset into gparts array. */
  ptrdiff_t gparts_offset = (ptrdiff_t)(gparts - s->gparts);

  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the index of this particle. */
    const size_t gpart_index = gparts_offset + ind;

    /* Get this particle's halo. */
    struct halo *halo = &groups[gpart_index];

    /* Set halo ID */
    halo->halo_id = group_id_default;

    /* And set the halo type. */
    halo->type = no_halo;
    
  }
}

/**
 * @brief Mapper function to calculate the numer of particles (size) in halos.
 *
 * @param map_data An array of #gpart%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void halo_size_mapper(void *map_data, int num_elements,
                      void *extra_data) {

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  struct gpart *gparts = (struct gpart *)map_data;
  enum halo_types current_level = s->e->fof_properties->current_level;
  
  /* Local copy of the halo arrays. */
  struct halo *groups;
  if (current_level == fof_group) {
    groups = s->e->fof_properties->groups;
  } else if (current_level == host_halo) {
    groups = s->e->fof_properties->hosts;
  } else {
    groups = s->e->fof_properties->subhalos;
  }

  /* Offset into gparts array. */
  ptrdiff_t gparts_offset = (ptrdiff_t)(gparts - s->gparts);

  /* Loop over particles and find which cells are in range of each other to
   * perform the FOF search. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the index of this particle. */
    const size_t gpart_ind = gparts_offset + ind;

    /* Get this particle's halo. */
    struct halo *halo = &groups[gpart_ind];

    /* Skip particles not in a halo. */
    if (halo->type == no_halo) continue;

    /* Count this particle. */
    atomic_inc(&halo->size);
    atomic_inc(&halo->props->npart_tot);
    atomic_inc(&halo->props->npart[gparts[gpart_ind].type]);
  }
}

/**
 * @brief Mapper function to calculate the group masses.
 *
 * @param map_data An array of #gpart%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void calc_halo_props_mapper(void *map_data, int num_elements,
                           void *extra_data) {

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  struct engine *e = s->e;
  const struct cosmology* cosmo = e->cosmology;
  struct gpart *gparts = (struct gpart *)map_data;
  enum halo_types current_level = s->e->fof_properties->current_level;

  /* Local copy of the halo arrays. */
  struct halo *groups;
  if (current_level == fof_group) {
    groups = s->e->fof_properties->groups;
  } else if (current_level == host_halo) {
    groups = s->e->fof_properties->hosts;
  } else {
    groups = s->e->fof_properties->subhalos;
  }

  /* Offset into gparts array. */
  ptrdiff_t gparts_offset = (ptrdiff_t)(gparts - s->gparts);

  /* Loop over particles and increment the group mass for groups above
   * min_group_size. */
  for (int ind = 0; ind < num_elements; ind++) {
    
    /* Get the index of this particle. */
    const size_t gpart_ind = gparts_offset + ind;

    /* Get this particle's halo. */
    struct halo *halo = &groups[gpart_ind];
    
    /* Only check actual halos above the minimum size. */
    if (halo->type == no_halo) continue;

    /* Get the particle. */
    struct gpart *gp = &gparts[gpart_ind];

    /* Get the particle properties. */
    const double mass = gp->mass;
    double x[3] = {gp->x[0],  gp->x[1], gp->x[2]};
    double vel[3] = {gp->v_full[0],  gp->v_full[1], gp->v_full[2]};
    double a_phys[3] = {gp->a_grav[0] * cosmo->a_factor_grav_accel,
                      gp->a_grav[1] * cosmo->a_factor_grav_accel,
                      gp->a_grav[2] * cosmo->a_factor_grav_accel};
    a_phys[0] += gp->a_grav_mesh[0] * cosmo->a_factor_grav_accel;
    a_phys[1] += gp->a_grav_mesh[1] * cosmo->a_factor_grav_accel;
    a_phys[2] += gp->a_grav_mesh[2] * cosmo->a_factor_grav_accel;
    double grav_pot = gp->potential;

    /* Add this particle's mass. */
    atomic_add_d(&halo->props->mass_tot, mass);
    atomic_add_d(&halo->props->mass[gp->type], mass);

    /* Set the reference frame for CoM calcualtion. */
    if (halo->props->first_position[0] == (double)(-FLT_MAX)) {
      halo->props->first_position[0] = x[0];
      halo->props->first_position[1] = x[1];
      halo->props->first_position[2] = x[2];

      /* Set this particle as the initial most bound while we are at it. */
      halo->props->most_bound_gpart = gp;
    }

    /* Wrap coordinates. */
    if (s->periodic) {
      x[0] = nearest(x[0] - halo->props->first_position[0], s->dim[0]);
      x[1] = nearest(x[1] - halo->props->first_position[1], s->dim[1]);
      x[2] = nearest(x[2] - halo->props->first_position[2], s->dim[2]);
    }

    /* Include particle in centre of mass. */
    atomic_add_d(&halo->props->centre_of_mass[0], mass * x[0]);
    atomic_add_d(&halo->props->centre_of_mass[1], mass * x[1]);
    atomic_add_d(&halo->props->centre_of_mass[2], mass * x[2]);

    /* Include this particle's contribution to the halo's velocity. */
    atomic_add_d(&halo->props->velocity[0], mass * vel[0]);
    atomic_add_d(&halo->props->velocity[1], mass * vel[1]);
    atomic_add_d(&halo->props->velocity[2], mass * vel[2]);

    /* Include this particles potential. */
    atomic_add_d(&halo->props->grav_pot, grav_pot);

    /* Include this particles kintetic energy. */
    double v2 = 0.0f;
    for (int k = 0; k < 3; k++)
      v2 = vel[k] * vel[k];
    atomic_add_d(&halo->props->kinetic_nrg,
                 0.5 * cosmo->a2_inv * mass * v2);

    /* Is this particle the most bound? */
    if (halo->props->most_bound_gpart->potential < grav_pot) {

      /* We have a new most bound particle. */
      halo->props->most_bound_gpart = gp;

      /* Set the centre of potential. */
      halo->props->centre_of_potential[0] = gp->x[0];
      halo->props->centre_of_potential[1] = gp->x[1];
      halo->props->centre_of_potential[2] = gp->x[2];
    }

    /* Add this particles contribution to the boosted potential. */
    atomic_add_d(&halo->props->a_phys[0], a_phys[0] * mass);
    atomic_add_d(&halo->props->a_phys[1], a_phys[1] * mass);
    atomic_add_d(&halo->props->a_phys[2], a_phys[2] * mass);

    /* Update the extent of the halo. */
    if (gp->x[0] < halo->props->extent[0])
      halo->props->extent[0] = gp->x[0];
    if (gp->x[0] > halo->props->extent[1])
      halo->props->extent[0] = gp->x[0];
    if (gp->x[1] < halo->props->extent[2])
      halo->props->extent[0] = gp->x[1];
    if (gp->x[1] > halo->props->extent[3])
      halo->props->extent[0] = gp->x[1];
    if (gp->x[2] < halo->props->extent[4])
      halo->props->extent[0] = gp->x[2];
    if (gp->x[2] > halo->props->extent[5])
      halo->props->extent[0] = gp->x[2];

    /* Store the host halo of this halo. */
    if (halo->type == host_halo) {
      halo->host = gp->fof_data.group;
    } else if (halo->type == sub_halo) {
      halo->host = gp->fof_data.host;
    } else {
      halo->host = NULL;
    }
  }
}

/**
 * @brief Function to clean up halo properties after the mapper over particles.
 *
 * @param s The #space.
 * @param props The properties of the FOF.
 * @param halo_props An array of pointers to halo property structs.
 * @param num_halos The number of halos to loop over.
 */
void finalise_halo_data(const struct space *s, struct fof_props *props,
                        struct halo_props *halo_props,
                        const int num_halos) {

  /* Define some useful variables. */
  const int periodic = s->periodic;
  const double *dim = s->dim;
  int nr_real_halos = 0;

  /* Loop over halos. */
  for (int i = 0; i < num_halos; i++) {

    /* Get this halo. */
    struct halo_props *halo = &halo_props[i];

    /* Calculate the width of the halo in each dimension (simply the distance
     * between minimum and maximum coordinate along each axis). */
    halo->width[0] = halo->extent[1] - halo->extent[0];
    halo->width[1] = halo->extent[3] - halo->extent[2];
    halo->width[2] = halo->extent[5] - halo->extent[4];

    /* Finish the Centre of mass, including possible box wrapping */
    double CoM[3] = {halo->centre_of_mass[0] / halo->mass_tot,
                     halo->centre_of_mass[1] / halo->mass_tot,
                     halo->centre_of_mass[2] / halo->mass_tot};
    if (periodic) {
      CoM[0] =
          box_wrap(CoM[0] + halo->first_position[i * 3 + 0], 0., dim[0]);
      CoM[1] =
          box_wrap(CoM[1] + halo->first_position[i * 3 + 1], 0., dim[1]);
      CoM[2] =
          box_wrap(CoM[2] + halo->first_position[i * 3 + 2], 0., dim[2]);
    }
    halo->centre_of_mass[i * 3 + 0] = CoM[0];
    halo->centre_of_mass[i * 3 + 1] = CoM[1];
    halo->centre_of_mass[i * 3 + 2] = CoM[2];

    /* Finish the velocity calculation. */
    halo->velocity[0] /= halo->mass_tot;
    halo->velocity[1] /= halo->mass_tot;
    halo->velocity[2] /= halo->mass_tot;
    halo->velocity_mag = sqrt(halo->velocity[0] * halo->velocity[0] +
                              halo->velocity[1] * halo->velocity[1] +
                              halo->velocity[2] * halo->velocity[2]);

    /* Finish the acceleration calculation. */
    halo->a_phys[0] /= halo->mass_tot;
    halo->a_phys[1] /= halo->mass_tot;
    halo->a_phys[2] /= halo->mass_tot;

    /* Compute the boosted potential. */
    halo->grav_boost = (halo->a_phys[0] * halo->centre_of_potential[0] +
                        halo->a_phys[1] * halo->centre_of_potential[1] +
                        halo->a_phys[2] * halo->centre_of_potential[2]);

    /* Is this halo real? */
    double grav_nrg_boost = halo->grav_pot + halo->grav_boost;
    double kin_nrg_boost = halo->kinetic_nrg + halo->grav_boost;
    if (grav_nrg_boost / kin_nrg_boost >= 1) {
      halo->is_real = 1;
      nr_real_halos++;
    }

    message("Found %d real halos (of %d)", nr_real_halos, num_halos);
    
  }
}


/**
 * @brief Perform a FOF search on gravity particles using the cells and applying
 * the Union-Find algorithm.
 *
 * @param props The properties of the FOF scheme.
 * @param constants The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param s The #space containing the particles.
 * @param dump_debug_results Are we writing txt-file debug catalogues including
 * BH-seeding info?
 * @param dump_results Do we want to write the group catalogue to a hdf5 file?
 */
void halo_search_tree(struct fof_props *props,
                      const struct phys_const *constants,
                      const struct cosmology *cosmo, struct space *s,
                      const int dump_results, const int dump_debug_results) {

  const size_t nr_gparts = s->nr_gparts;
  const size_t min_group_size = props->min_group_size;
  const size_t group_id_offset = props->group_id_offset;
  struct gpart *gparts = s->gparts;
  size_t *group_index;
  struct halo *groups;
  long long num_groups = 0, num_parts_in_groups = 0, max_group_size = 0;
  const int verbose = s->e->verbose;
  const ticks tic_total = getticks();

  /* What level of the search are we doing? */
  enum halo_types current_level= props->current_level;

  char output_file_name[PARSER_MAX_LINE_SIZE];
  snprintf(output_file_name, PARSER_MAX_LINE_SIZE, "%s", props->base_name);

  if (verbose)
    message("Searching %zu gravity particles for links with l_x: %lf",
            nr_gparts, sqrt(props->l_x2));

  if (engine_rank == 0 && verbose)
    message("Size of hash table element: %ld", sizeof(hashmap_element_t));

  /* Local copy of the arrays. */
  if (current_level == fof_group) {
    group_index = props->group_index;
    groups = props->groups;
  } else if (current_level == host_halo) {
    group_index = props->host_index;
    groups = props->hosts;
  } else {
    group_index = props->subhalo_index;
    groups = props->subhalos;
  }
  
  const ticks tic_calc_group_size = getticks();

  threadpool_map(&s->e->threadpool, halo_size_mapper, gparts,
                 nr_gparts, sizeof(struct gpart), threadpool_auto_chunk_size,
                 s);
  if (verbose)
    message("FOF calc group size took (FOF SCALING): %.3f %s.",
            clocks_from_ticks(getticks() - tic_calc_group_size),
            clocks_getunit());

  size_t num_groups_local = 0;
  size_t num_parts_in_groups_local = 0;
  size_t max_group_size_local = 0;

  const ticks tic_num_groups_calc = getticks();

  for (size_t i = 0; i < nr_gparts; i++) {

    if (group_index[i] == i)
        message("Halo %ld has %ld particles", group_index[i], groups[i].size);

    /* Find the total number of groups. */
    if (group_index[i] == i && groups[i].size >= min_group_size)
      num_groups_local++;

    /* Find the total number of particles in groups. */
    if (groups[i].size >= min_group_size)
      num_parts_in_groups_local += groups[i].size;

    /* Find the largest group. */
    if (groups[i].size > max_group_size_local)
      max_group_size_local = groups[i].size;
  }

  if (verbose)
    message(
        "Calculating the total no. of local groups took: (FOF SCALING): %.3f "
        "%s.",
        clocks_from_ticks(getticks() - tic_num_groups_calc), clocks_getunit());

  /* Sort the groups in descending order based upon size and re-label their
   * IDs 0-num_groups. */
  struct group_length *high_group_sizes = NULL;
  int group_count = 0;

  if (swift_memalign("fof_high_group_sizes", (void **)&high_group_sizes, 32,
                     num_groups_local * sizeof(struct group_length)) != 0)
    error("Failed to allocate list of large groups.");

  /* Store the group_sizes and their offset. */
  for (size_t i = 0; i < nr_gparts; i++) {

    if (group_index[i] == i && groups[i].size >= min_group_size) {
      high_group_sizes[group_count].index = i;
      high_group_sizes[group_count++].size = groups[i].size;
    }
  }

  ticks tic = getticks();

  /* Find global properties. */
  num_groups = num_groups_local;
  num_parts_in_groups = num_parts_in_groups_local;
  max_group_size = max_group_size_local;

  /* Store what we've found... */
  if (current_level == fof_group) {
    props->num_groups = num_groups;
    props->num_groups_rank = num_groups_local;
  } else if (current_level == host_halo) {
    props->num_hosts = num_groups;
    props->num_hosts_rank = num_groups_local;
  } else {
    props->num_subhalos = num_groups;
    props->num_subhalos_rank = num_groups_local;
  }

  if (verbose)
    message("Finding the total no. of groups took: (FOF SCALING): %.3f %s.",
            clocks_from_ticks(getticks() - tic_num_groups_calc),
            clocks_getunit());

  /* Sort local groups into descending order of size */
  qsort(high_group_sizes, num_groups_local, sizeof(struct group_length),
        cmp_func_group_size);

  tic = getticks();

  /* Set default group ID for all particles */
  threadpool_map(&s->e->threadpool, set_initial_halo_id_mapper, s->gparts,
                 s->nr_gparts, sizeof(struct gpart), threadpool_auto_chunk_size,
                 s);

  if (verbose)
    message("Setting default group ID took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Assign the group and final group IDs to local root particles where the
   * global root is on this node and the group is large enough. Within a node
   * IDs are assigned in descending order of particle number. */
  /* Attatch the halo properties to the halos. */
  for (size_t i = 0; i < num_groups_local; i++) {
    size_t part_ind = high_group_sizes[i].index;
    if (current_level == fof_group) {
      gparts[part_ind].fof_data.group = &groups[part_ind];
    } else if (current_level == host_halo) {
      gparts[part_ind].fof_data.host = &groups[part_ind];
    } else {
      gparts[part_ind].fof_data.subhalo = &groups[part_ind];
    }
    groups[part_ind].halo_id = group_id_offset + i;
  }

  /* Assign every particle the pointer to its local root. */
  for (size_t i = 0; i < nr_gparts; i++) {
    const size_t root = fof_find_local(i, nr_gparts, group_index);
    if (current_level == fof_group) {
      gparts[i].fof_data.group = &groups[root];
    } else if (current_level == host_halo) {
      gparts[i].fof_data.host = &groups[root];
    } else {
      gparts[i].fof_data.subhalo = &groups[root];
    }
  }

  if (verbose)
    message("Group sorting took: %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* Local copy of the arrays. */
  struct halo_props *halo_props;
  if (current_level == fof_group) {
    halo_props = props->group_props;
  } else if (current_level == host_halo) {
    halo_props = props->host_props;
  } else {
    halo_props = props->subhalo_props;
  }

  /* Allocate and initialise property arrays. */
  if (swift_memalign("fof_group_props", (void **)&halo_props,
                     halo_align,
                     num_groups_local * sizeof(struct halo_props)) != 0)
    error("Failed to allocate list of group properties.");
    bzero(halo_props, num_groups_local * sizeof(struct halo_props));

  /* Attatch the halo properties to the halos. */
  for (size_t i = 0; i < num_groups_local; i++) {
    size_t part_ind = high_group_sizes[i].index;
    if (current_level == fof_group) {
      groups[part_ind].props = &props->group_props[i];
      props->group_props[i].npart_tot = groups[part_ind].size;
    } else if (current_level == host_halo) {
      groups[part_ind].props = &props->host_props[i];
      props->host_props[i].npart_tot = groups[part_ind].size;
    } else {
      groups[part_ind].props = &props->subhalo_props[i];
      props->subhalo_props[i].npart_tot = groups[part_ind].size;
    }
  }

  /* Set up halo first positions (used as a reference frame for centre of
   * mass calculation) and extents. */
  for (size_t i = 0; i < num_groups_local; i++) {
    for (int k = 0; k < 3; k++) {
    halo_props[i].first_position[k] = -FLT_MAX;
    }
    for (int k = 0; k < 6; k++) {
      if (k % 2 == 0) {
        halo_props[i].extent[k] = -FLT_MAX;
      } else {
        halo_props[i].extent[k] = FLT_MAX;
      }
    }
  }

  ticks tic_seeding = getticks();

  /* Calculate halo properties. */
  threadpool_map(&s->e->threadpool, calc_halo_props_mapper, gparts,
                 nr_gparts, sizeof(struct gpart), threadpool_auto_chunk_size,
                 (struct space *)s);

  /* Finalise the group data before dump */
  finalise_halo_data(s, props, halo_props, num_groups_local);
  
  if (verbose)
    message("Computing group properties took: %.3f %s.",
            clocks_from_ticks(getticks() - tic_seeding), clocks_getunit());
    
  tic_seeding = getticks();

  /* Get the right particle arrays. */
  size_t *group_particle_inds;
  double *group_particle_pos;
  if (current_level == fof_group) {
    group_particle_inds = props->group_particle_inds;
    group_particle_pos = props->group_particle_pos;
  } else if (current_level == host_halo) {
    group_particle_inds = props->host_particle_inds;
    group_particle_pos = props->host_particle_pos;
  } else {
    group_particle_inds = props->subhalo_particle_inds;
    group_particle_pos = props->subhalo_particle_pos;
  }
  
  /* Allocate arrays to hold particle indices and positions. */
  group_particle_inds =
    (size_t *)swift_malloc("fof_group_particle_indices",
                           num_parts_in_groups_local * sizeof(size_t));
  bzero(group_particle_inds, num_parts_in_groups_local * sizeof(size_t));
  group_particle_pos =
    (double *)swift_malloc("fof_group_particle_postions",
                           num_parts_in_groups_local * 3 * sizeof(double));
  bzero(group_particle_pos,
        num_parts_in_groups_local * 3 * sizeof(double));

  /* Allocate and initialise temporary counters for assigning particles. */
  int *part_counters  =
    (int *)swift_malloc("fof_particle_counters",
                        num_groups_local * sizeof(int));
  bzero(part_counters, num_groups_local * sizeof(int));
  
  /* Populate pointers. */
  halo_props[0].part_start_index = 0;
  for (size_t i = 1; i < num_groups_local; i++) {
    halo_props[i].part_start_index =
      halo_props[i - 1].part_start_index + halo_props[i - 1].npart_tot;
  }
  
  /* Populate particle arrays. */
  /* TODO: threadpool this. */
  for (size_t i = 0; i < nr_gparts; i++) {
    
    /* Get this particle's halo. */
    struct halo *halo = &groups[i];
    
    /* Skip particles not in a group. */
    if (halo->type == no_halo) continue;
    
    /* Get the index for this group. */
    size_t halo_ind = halo->halo_id - group_id_offset;
    
       /* Get the start pointer for this group. */
    size_t start = halo->props->part_start_index;
    
    /* Assign this particle's index to the corresponding positon. */
    group_particle_inds[start + part_counters[halo_ind]] = i;
    
    /* Assign this particles position. */
    for (int k = 0; k < 3; k++)
      group_particle_pos[(start + part_counters[halo_ind]) * 3 + k] =
           gparts[i].x[k];
    
    /* Increment halo particle counter. */
    part_counters[halo_ind]++;
    
  }
     
  swift_free("fof_particle_counters", part_counters); 

  if (verbose)
    message("Sorting group particles took: %.3f %s.",
            clocks_from_ticks(getticks() - tic_seeding), clocks_getunit());
  
  /* Assign the number of particles in groups */
  props->num_parts_in_groups = num_parts_in_groups_local;
  if (current_level == fof_group) {
    props->num_parts_in_groups = num_parts_in_groups_local;
  } else if (current_level == host_halo) {
    props->num_parts_in_hosts = num_parts_in_groups_local;
  } else {
    props->num_parts_in_subhalos = num_parts_in_groups_local;
  }

  /* Free the left-overs */
  swift_free("fof_high_group_sizes", high_group_sizes);

  if (engine_rank == 0) {
    message(
        "No. of groups: %lld. No. of particles in groups: %lld. No. of "
        "particles not in groups: %lld.",
        num_groups, num_parts_in_groups,
        s->e->total_nr_gparts - num_parts_in_groups);

    message("Largest group by size: %lld", max_group_size);
  }
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic_total),
            clocks_getunit());
}
