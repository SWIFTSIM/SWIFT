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

/* Config parameters. */
#include <config.h>

#ifdef WITH_FOF

/* Some standard headers. */
#include <errno.h>
#include <libgen.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "fof.h"

/* Local headers. */
#include "black_holes.h"
#include "common_io.h"
#include "engine.h"
#include "fof_catalogue_io.h"
#include "hashmap.h"
#include "memuse.h"
#include "proxy.h"
#include "threadpool.h"
#include "tools.h"
#include "tracers.h"

#define fof_props_default_group_id 2147483647
#define fof_props_default_group_id_offset 1
#define fof_props_default_group_link_size 20000

/* Constants. */
#define UNION_BY_SIZE_OVER_MPI (1)
#define FOF_COMPRESS_PATHS_MIN_LENGTH (2)

/* The FoF policy we are running */
int current_fof_linking_type;

/* The FoF policy for particles attached to the main type */
int current_fof_attach_type;

/* The FoF policy for particles ignored altogether */
int current_fof_ignore_type;

/* Are we timing calculating group properties in the FOF? */
// #define WITHOUT_GROUP_PROPS

/**
 * @brief Properties of a group used for black hole seeding
 */
enum fof_halo_seeding_props {
  fof_halo_has_no_gas = -1LL,
  fof_halo_has_black_hole = -2LL,
  fof_halo_has_too_low_mass = -3LL
};

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

#ifdef SWIFT_DEBUG_CHECKS
static integertime_t ti_current;
#endif

void fof_set_current_types(const struct fof_props *props) {

  /* Initialize the FoF linking mode */
  current_fof_linking_type = 0;
  for (int i = 0; i < swift_type_count; ++i)
    if (props->fof_linking_types[i]) {
      current_fof_linking_type |= (1 << (i + 1));
    }

  /* Initialize the FoF attaching mode */
  current_fof_attach_type = 0;
  for (int i = 0; i < swift_type_count; ++i)
    if (props->fof_attach_types[i]) {
      current_fof_attach_type |= (1 << (i + 1));
    }

  /* Construct the combined mask of ignored particles */
  current_fof_ignore_type =
      ~(current_fof_linking_type | current_fof_attach_type);
}

/**
 * @brief Initialise the properties of the FOF code.
 *
 * @param props the #fof_props structure to fill.
 * @param params the parameter file parser.
 * @param phys_const The physical constants in internal units.
 * @param us The internal unit system.
 * @param stand_alone_fof Are we initialising for stand-alone? (1) or
 * on-the-fly? (0)
 */
void fof_init(struct fof_props *props, struct swift_params *params,
              const struct phys_const *phys_const, const struct unit_system *us,
              const int stand_alone_fof) {

  /* Base name for the FOF output file */
  parser_get_param_string(params, "FOF:basename", props->base_name);

  /* Check that we can write outputs by testing if the output
   * directory exists and is searchable and writable. */
  char directory[PARSER_MAX_LINE_SIZE] = {0};
  sprintf(directory, "%s", props->base_name);
  const char *dirp = dirname(directory);
  if (engine_rank == 0) safe_checkdir(dirp, /*create=*/1);

  /* Read the minimum group size. */
  props->min_group_size = parser_get_param_int(params, "FOF:min_group_size");

  /* Read whether we're doing FoF calls to seed black holes. */
  props->seed_black_holes_enabled =
      parser_get_param_int(params, "FOF:seed_black_holes_enabled");

  /* Read the default group ID of particles in groups below the minimum group
   * size. */
  props->group_id_default = parser_get_opt_param_int(
      params, "FOF:group_id_default", fof_props_default_group_id);

  /* Read the starting group ID. */
  props->group_id_offset = parser_get_opt_param_int(
      params, "FOF:group_id_offset", fof_props_default_group_id_offset);

  /* Read the linking length ratio to the mean inter-particle separation. */
  props->l_x_ratio =
      parser_get_opt_param_double(params, "FOF:linking_length_ratio", -1.);

  /* Read value of absolute linking length aksed by the user */
  props->l_x_absolute =
      parser_get_opt_param_double(params, "FOF:absolute_linking_length", -1.);

  if (props->l_x_ratio == -1. && props->l_x_absolute <= 0.)
    error("The FOF linking length can't be negative!");

  if (props->l_x_ratio <= 0. && props->l_x_absolute == -1.)
    error("The FOF linking length ratio can't be negative!");

  if (!stand_alone_fof && props->seed_black_holes_enabled) {

    /* Read the minimal halo mass for black hole seeding */
    props->seed_halo_mass =
        parser_get_param_double(params, "FOF:black_hole_seed_halo_mass_Msun");

    /* Convert to internal units */
    props->seed_halo_mass *= phys_const->const_solar_mass;
  }

  /* Read what particle types we want to run FOF on */
  parser_get_param_int_array(params, "FOF:linking_types", swift_type_count,
                             props->fof_linking_types);

  /* Read what particle types we want to attach to FOF groups */
  parser_get_param_int_array(params, "FOF:attaching_types", swift_type_count,
                             props->fof_attach_types);

  /* Check that there is something to do */
  int sum = 0;
  for (int i = 0; i < swift_type_count; ++i)
    if (props->fof_linking_types[i]) sum++;
  if (sum == 0) error("FOF must run on at least one type of particles!");

  for (int i = 0; i < swift_type_count; ++i)
    if (props->fof_linking_types[i] && props->fof_attach_types[i])
      error("FOF can't use a type (%s) as both linking and attaching type!",
            part_type_names[i]);

  /* Set the current FOF types */
  fof_set_current_types(props);

  /* Report what we do */
  if (engine_rank == 0) {
    printf("FOF using the following types for linking:");
    for (int i = 0; i < swift_type_count; ++i)
      if (props->fof_linking_types[i]) printf("'%s' ", part_type_names[i]);
    printf("\n");

    printf("FOF using the following types for attaching:");
    for (int i = 0; i < swift_type_count; ++i)
      if (props->fof_attach_types[i]) printf("'%s' ", part_type_names[i]);
    printf("\n");

    printf("FOF ignoring the following types:");
    for (int i = 0; i < swift_type_count; ++i)
      if (current_fof_ignore_type & (1 << (i + 1)))
        printf("'%s' ", part_type_names[i]);
    printf("\n");
  }

#if defined(WITH_MPI) && defined(UNION_BY_SIZE_OVER_MPI)
  if (engine_rank == 0)
    message(
        "Performing FOF over MPI using union by size and union by rank "
        "locally.");
#else
  message("Performing FOF using union by rank.");
#endif
}

/**
 * @brief Registers MPI types used by FOF.
 */
void fof_create_mpi_types(void) {

#ifdef WITH_MPI
  if (MPI_Type_contiguous(sizeof(struct fof_mpi) / sizeof(unsigned char),
                          MPI_BYTE, &fof_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&fof_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for fof.");
  }
  if (MPI_Type_contiguous(sizeof(struct group_length) / sizeof(unsigned char),
                          MPI_BYTE, &group_length_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&group_length_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for group_length.");
  }
  /* Define type for sending fof_final_index struct */
  if (MPI_Type_contiguous(sizeof(struct fof_final_index), MPI_BYTE,
                          &fof_final_index_type) != MPI_SUCCESS ||
      MPI_Type_commit(&fof_final_index_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for fof_final_index.");
  }
  /* Define type for sending fof_final_mass struct */
  if (MPI_Type_contiguous(sizeof(struct fof_final_mass), MPI_BYTE,
                          &fof_final_mass_type) != MPI_SUCCESS ||
      MPI_Type_commit(&fof_final_mass_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for fof_final_mass.");
  }
#else
  error("Calling an MPI function in non-MPI code.");
#endif
}

/**
 * @brief Mapper function to set the initial group indices.
 *
 * @param map_data The array of group indices.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to first group index.
 */
void fof_set_initial_group_index_mapper(void *map_data, int num_elements,
                                        void *extra_data) {
  size_t *group_index = (size_t *)map_data;
  size_t *group_index_start = (size_t *)extra_data;

  const ptrdiff_t offset = group_index - group_index_start;

  for (int i = 0; i < num_elements; ++i) {
    group_index[i] = i + offset;
  }
}

/**
 * @brief Mapper function to set the initial attach group indices.
 *
 * @param map_data The array of attach group indices.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to first group index.
 */
void fof_set_initial_attach_index_mapper(void *map_data, int num_elements,
                                         void *extra_data) {
  size_t *attach_index = (size_t *)map_data;
  size_t *attach_index_start = (size_t *)extra_data;

  const ptrdiff_t offset = attach_index - attach_index_start;

  for (int i = 0; i < num_elements; ++i) {
    attach_index[i] = i + offset;
  }
}

/**
 * @brief Mapper function to set the initial group sizes.
 *
 * @param map_data The array of group sizes.
 * @param num_elements Chunk size.
 * @param extra_data N/A.
 */
void fof_set_initial_group_size_mapper(void *map_data, int num_elements,
                                       void *extra_data) {

  size_t *group_size = (size_t *)map_data;
  for (int i = 0; i < num_elements; ++i) {
    group_size[i] = 1;
  }
}

/**
 * @brief Mapper function to set the initial distances.
 *
 * @param map_data The array of distance.
 * @param num_elements Chunk size.
 * @param extra_data N/A.
 */
void fof_set_initial_part_distances_mapper(void *map_data, int num_elements,
                                           void *extra_data) {

  float *distance = (float *)map_data;
  for (int i = 0; i < num_elements; ++i) {
    distance[i] = FLT_MAX;
  }
}

/**
 * @brief Mapper function to set the initial group IDs.
 *
 * @param map_data The array of #gpart%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to the default group ID.
 */
void fof_set_initial_group_id_mapper(void *map_data, int num_elements,
                                     void *extra_data) {

  /* Unpack the information */
  struct gpart *gparts = (struct gpart *)map_data;
  const size_t group_id_default = *((size_t *)extra_data);

  for (int i = 0; i < num_elements; ++i) {
    gparts[i].fof_data.group_id = group_id_default;
  }
}

/**
 * @brief Allocate the memory and initialise the arrays for a FOF calculation.
 *
 * @param s The #space to act on.
 * @param props The properties of the FOF structure.
 */
void fof_allocate(const struct space *s, struct fof_props *props) {

  const int verbose = s->e->verbose;
  const ticks total_tic = getticks();

  /* Start by computing the mean inter DM particle separation */

  /* Collect the mean mass of the non-background gpart */
  double high_res_DM_mass = 0.;
  long long num_high_res_DM = 0;
  for (size_t i = 0; i < s->nr_gparts; ++i) {
    const struct gpart *gp = &s->gparts[i];
    if (gp->type == swift_type_dark_matter &&
        gp->time_bin != time_bin_inhibited &&
        gp->time_bin != time_bin_not_created) {
      high_res_DM_mass += gp->mass;
      num_high_res_DM++;
    }
  }

#ifdef WITH_MPI
  /* Gather the information from all ranks */
  MPI_Allreduce(MPI_IN_PLACE, &high_res_DM_mass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &num_high_res_DM, 1, MPI_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  high_res_DM_mass /= (double)num_high_res_DM;

  /* Are we using the aboslute value or the one derived from the mean
     inter-particle sepration? */
  if (props->l_x_absolute != -1.) {
    props->l_x2 = props->l_x_absolute * props->l_x_absolute;

    if (s->e->nodeID == 0) {
      message("Linking length is set to %e [internal units].",
              props->l_x_absolute);
    }
  } else {

    /* Safety check */
    if (!(s->e->policy & engine_policy_cosmology))
      error(
          "Attempting to run FoF on a simulation using cosmological "
          "information but cosmology was not initialised");

    /* Calculate the mean inter-particle separation as if we were in
       a scenario where the entire box was filled with high-resolution
         particles */
    const double Omega_cdm = s->e->cosmology->Omega_cdm;
    const double Omega_b = s->e->cosmology->Omega_b;
    const double Omega_m = Omega_cdm + Omega_b;
    const double critical_density_0 = s->e->cosmology->critical_density_0;
    double mean_matter_density;
    if (s->with_hydro)
      mean_matter_density = Omega_cdm * critical_density_0;
    else
      mean_matter_density = Omega_m * critical_density_0;

    /* Mean inter-particle separation of the DM particles */
    const double mean_inter_particle_sep =
        cbrt(high_res_DM_mass / mean_matter_density);

    /* Calculate the particle linking length based upon the mean inter-particle
     * spacing of the DM particles. */
    const double l_x = props->l_x_ratio * mean_inter_particle_sep;

    props->l_x2 = l_x * l_x;

    if (s->e->nodeID == 0) {
      message(
          "Linking length is set to %e [internal units] (%f of mean "
          "inter-DM-particle separation).",
          l_x, props->l_x_ratio);
    }
  }

#ifdef WITH_MPI
  /* Check size of linking length against the top-level cell dimensions. */
  if (props->l_x2 > s->width[0] * s->width[0])
    error(
        "Linking length greater than the width of a top-level cell. Need to "
        "check more than one layer of top-level cells for links.");
#endif

  /* Allocate and initialise a group index array. */
  if (swift_memalign("fof_group_index", (void **)&props->group_index, 64,
                     s->nr_gparts * sizeof(size_t)) != 0)
    error("Failed to allocate list of particle group indices for FOF search.");

  /* Allocate and initialise a group index array for attachables. */
  if (swift_memalign("fof_attach_index", (void **)&props->attach_index, 64,
                     s->nr_gparts * sizeof(size_t)) != 0)
    error(
        "Failed to allocate list of particle distances array for FOF search.");

  /* Allocate and initialise a group index array for attachables. */
  if (swift_memalign("fof_found_attach", (void **)&props->found_attachable_link,
                     64, s->nr_gparts * sizeof(char)) != 0)
    error(
        "Failed to allocate list of particle distances array for FOF search.");

  /* Allocate and initialise the closest distance array. */
  if (swift_memalign("fof_distance", (void **)&props->distance_to_link, 64,
                     s->nr_gparts * sizeof(float)) != 0)
    error(
        "Failed to allocate list of particle distances array for FOF search.");

  /* Allocate and initialise a group size array. */
  if (swift_memalign("fof_group_size", (void **)&props->group_size, 64,
                     s->nr_gparts * sizeof(size_t)) != 0)
    error("Failed to allocate list of group size for FOF search.");

  ticks tic = getticks();

  /* Set initial group index */
  threadpool_map(&s->e->threadpool, fof_set_initial_group_index_mapper,
                 props->group_index, s->nr_gparts, sizeof(size_t),
                 threadpool_auto_chunk_size, props->group_index);

  if (verbose)
    message("Setting initial group index took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Set initial attach index */
  threadpool_map(&s->e->threadpool, fof_set_initial_attach_index_mapper,
                 props->attach_index, s->nr_gparts, sizeof(size_t),
                 threadpool_auto_chunk_size, props->attach_index);

  bzero(props->found_attachable_link, s->nr_gparts * sizeof(char));

  if (verbose)
    message("Setting initial attach index took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Set initial distances */
  threadpool_map(&s->e->threadpool, fof_set_initial_part_distances_mapper,
                 props->distance_to_link, s->nr_gparts, sizeof(float),
                 threadpool_auto_chunk_size, NULL);

  if (verbose)
    message("Setting initial distances took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Set initial group sizes */
  threadpool_map(&s->e->threadpool, fof_set_initial_group_size_mapper,
                 props->group_size, s->nr_gparts, sizeof(size_t),
                 threadpool_auto_chunk_size, NULL);

  if (verbose)
    message("Setting initial group sizes took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#ifdef SWIFT_DEBUG_CHECKS
  ti_current = s->e->ti_current;
#endif

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - total_tic),
            clocks_getunit());
}

/**
 * @brief Comparison function for qsort call comparing group sizes.
 *
 * @param a The first #group_length object.
 * @param b The second #group_length object.
 * @return 1 if the size of the group b is larger than the size of group a, -1
 * if a is the largest and 0 if they are equal.
 */
int cmp_func_group_size(const void *a, const void *b) {
  struct group_length *a_group_size = (struct group_length *)a;
  struct group_length *b_group_size = (struct group_length *)b;
  if (b_group_size->size > a_group_size->size)
    return 1;
  else if (b_group_size->size < a_group_size->size)
    return -1;
  else
    return 0;
}

#ifdef WITH_MPI

/**
 * @brief Comparison function for qsort call comparing group global roots.
 *
 * @param a The first #fof_final_index object.
 * @param b The second #fof_final_index object.
 * @return 1 if the global of the group b is *smaller* than the global group of
 * group a, -1 if a is the smaller one and 0 if they are equal.
 */
int compare_fof_final_index_global_root(const void *a, const void *b) {
  struct fof_final_index *fof_final_index_a = (struct fof_final_index *)a;
  struct fof_final_index *fof_final_index_b = (struct fof_final_index *)b;
  if (fof_final_index_b->global_root < fof_final_index_a->global_root)
    return 1;
  else if (fof_final_index_b->global_root > fof_final_index_a->global_root)
    return -1;
  else
    return 0;
}

/**
 * @brief Comparison function for qsort call comparing group global roots
 *
 * @param a The first #fof_final_mass object.
 * @param b The second #fof_final_mass object.
 * @return 1 if the global of the group b is *smaller* than the global group of
 * group a, -1 if a is the smaller one and 0 if they are equal.
 */
int compare_fof_final_mass_global_root(const void *a, const void *b) {
  struct fof_final_mass *fof_final_mass_a = (struct fof_final_mass *)a;
  struct fof_final_mass *fof_final_mass_b = (struct fof_final_mass *)b;
  if (fof_final_mass_b->global_root < fof_final_mass_a->global_root)
    return 1;
  else if (fof_final_mass_b->global_root > fof_final_mass_a->global_root)
    return -1;
  else
    return 0;
}

/**
 * @brief Check whether a given group ID is on the local node.
 *
 * This function only makes sense in MPI mode.
 *
 * @param group_id The ID to check.
 * @param nr_gparts The number of gparts on this node.
 */
__attribute__((always_inline)) INLINE static int is_local(
    const size_t group_id, const size_t nr_gparts) {
#ifdef WITH_MPI
  return (group_id >= node_offset && group_id < node_offset + nr_gparts);
#else
  error("Calling MPI function in non-MPI mode");
  return 1;
#endif
}

/**
 * @brief Find the global root ID of a given particle
 *
 * This function only makes sense in MPI mode.
 *
 * @param i Index of the particle.
 * @param group_index Array of group root indices.
 * @param nr_gparts The number of g-particles on this node.
 */
__attribute__((always_inline)) INLINE static size_t fof_find_global(
    const size_t i, const size_t *group_index, const size_t nr_gparts) {

#ifdef WITH_MPI
  size_t root = node_offset + i;
  if (!is_local(root, nr_gparts)) {

    /* Non local --> This is the root */
    return root;
  } else {

    /* Local --> Follow the links until we find the root */
    while (root != group_index[root - node_offset]) {
      root = group_index[root - node_offset];
      if (!is_local(root, nr_gparts)) break;
    }
  }

  /* Perform path compression. */
  // int index = i;
  // while(index != root) {
  //  int next = group_index[index];
  //  group_index[index] = root;
  //  index = next;
  //}

  return root;
#else
  error("Calling MPI function in non-MPI mode");
  return -1;
#endif
}

#endif /* WITH_MPI */

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
 * @brief Returns whether a #gpart is of the 'attachable' kind.
 */
__attribute__((always_inline)) INLINE static int gpart_is_attachable(
    const struct gpart *gp) {

  return current_fof_attach_type & (1 << (gp->type + 1));
}

/**
 * @brief Returns whether a #gpart is of the 'linkable' kind.
 */
__attribute__((always_inline)) INLINE static int gpart_is_linkable(
    const struct gpart *gp) {

  return current_fof_linking_type & (1 << (gp->type + 1));
}

/**
 * @brief Returns whether a #gpart is to be ignored by FOF.
 */
__attribute__((always_inline)) INLINE static int gpart_is_ignorable(
    const struct gpart *gp) {

  return current_fof_ignore_type & (1 << (gp->type + 1));
}

/**
 * @brief Returns whether a foreign #gpart is of the 'attachable' kind.
 */
__attribute__((always_inline)) INLINE static int gpart_foreign_is_attachable(
    const struct gpart_fof_foreign *gp) {

  return current_fof_attach_type & (1 << (gp->type + 1));
}

/**
 * @brief Returns whether a foreign #gpart is of the 'linkable' kind.
 */
__attribute__((always_inline)) INLINE static int gpart_foreign_is_linkable(
    const struct gpart_fof_foreign *gp) {

  return current_fof_linking_type & (1 << (gp->type + 1));
}

/**
 * @brief Returns whether a foreign #gpart is to be ignored by FOF.
 */
__attribute__((always_inline)) INLINE static int gpart_foreign_is_ignorable(
    const struct gpart_fof_foreign *gp) {

  return current_fof_ignore_type & (1 << (gp->type + 1));
}

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
 * @brief Atomically update the root of a group
 *
 * @param address The address of the value to update.
 * @param y The new value to write.
 *
 * @return 1 If successful, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int atomic_update_root(
    volatile size_t *address, const size_t y) {

  size_t *size_t_ptr = (size_t *)address;

  size_t old_val = *address;
  size_t test_val = old_val;
  size_t new_val = y;

  /* atomic_cas returns old_val if *size_t_ptr has not changed since being
   * read.*/
  old_val = atomic_cas(size_t_ptr, test_val, new_val);

  if (test_val == old_val)
    return 1;
  else
    return 0;
}

/**
 * @brief Unifies two groups by setting them to the same root.
 *
 * @param root_i The root of the first group. Will be updated.
 * @param root_j The root of the second group.
 * @param group_index The list of group roots.
 */
__attribute__((always_inline)) INLINE static void fof_union(
    size_t *restrict root_i, const size_t root_j,
    size_t *restrict group_index) {

  int result = 0;

  /* Loop until the root can be set to a new value. */
  do {
    size_t root_i_new = fof_find(*root_i, group_index);
    const size_t root_j_new = fof_find(root_j, group_index);

    /* Skip particles in the same group. */
    if (root_i_new == root_j_new) return;

    /* If the root ID of pj is lower than pi's root ID set pi's root to point to
     * pj's. Otherwise set pj's root to point to pi's.*/
    if (root_j_new < root_i_new) {

      /* Updates the root and checks that its value has not been changed since
       * being read. */
      result = atomic_update_root(&group_index[root_i_new], root_j_new);

      /* Update root_i on the fly. */
      *root_i = root_j_new;
    } else {

      /* Updates the root and checks that its value has not been changed since
       * being read. */
      result = atomic_update_root(&group_index[root_j_new], root_i_new);

      /* Update root_i on the fly. */
      *root_i = root_i_new;
    }
  } while (result != 1);
}

/**
 * @brief Compute th minimal distance between any two points in two cells.
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param dim The size of the simulation domain.
 */
__attribute__((always_inline)) INLINE static double cell_min_dist(
    const struct cell *restrict ci, const struct cell *restrict cj,
    const double dim[3]) {

  /* Get cell locations. */
  const double cix_min = ci->loc[0];
  const double ciy_min = ci->loc[1];
  const double ciz_min = ci->loc[2];
  const double cjx_min = cj->loc[0];
  const double cjy_min = cj->loc[1];
  const double cjz_min = cj->loc[2];

  const double cix_max = ci->loc[0] + ci->width[0];
  const double ciy_max = ci->loc[1] + ci->width[1];
  const double ciz_max = ci->loc[2] + ci->width[2];
  const double cjx_max = cj->loc[0] + cj->width[0];
  const double cjy_max = cj->loc[1] + cj->width[1];
  const double cjz_max = cj->loc[2] + cj->width[2];

  double not_same_range[3];

  /* If two cells are in the same range of coordinates along
     any of the 3 axis, the distance along this axis is 0 */
  if (ci->width[0] > cj->width[0]) {
    if ((cix_min <= cjx_min) && (cjx_max <= cix_max))
      not_same_range[0] = 0.;
    else
      not_same_range[0] = 1.;
  } else {
    if ((cjx_min <= cix_min) && (cix_max <= cjx_max))
      not_same_range[0] = 0.;
    else
      not_same_range[0] = 1.;
  }
  if (ci->width[1] > cj->width[1]) {
    if ((ciy_min <= cjy_min) && (cjy_max <= ciy_max))
      not_same_range[1] = 0.;
    else
      not_same_range[1] = 1.;
  } else {
    if ((cjy_min <= ciy_min) && (ciy_max <= cjy_max))
      not_same_range[1] = 0.;
    else
      not_same_range[1] = 1.;
  }
  if (ci->width[2] > cj->width[2]) {
    if ((ciz_min <= cjz_min) && (cjz_max <= ciz_max))
      not_same_range[2] = 0.;
    else
      not_same_range[2] = 1.;
  } else {
    if ((cjz_min <= ciz_min) && (ciz_max <= cjz_max))
      not_same_range[2] = 0.;
    else
      not_same_range[2] = 1.;
  }

  /* Find the shortest distance between cells, remembering to account for
   * periodic boundary conditions. */
  double dx[3];
  dx[0] = min4(fabs(nearest(cix_min - cjx_min, dim[0])),
               fabs(nearest(cix_min - cjx_max, dim[0])),
               fabs(nearest(cix_max - cjx_min, dim[0])),
               fabs(nearest(cix_max - cjx_max, dim[0])));

  dx[1] = min4(fabs(nearest(ciy_min - cjy_min, dim[1])),
               fabs(nearest(ciy_min - cjy_max, dim[1])),
               fabs(nearest(ciy_max - cjy_min, dim[1])),
               fabs(nearest(ciy_max - cjy_max, dim[1])));

  dx[2] = min4(fabs(nearest(ciz_min - cjz_min, dim[2])),
               fabs(nearest(ciz_min - cjz_max, dim[2])),
               fabs(nearest(ciz_max - cjz_min, dim[2])),
               fabs(nearest(ciz_max - cjz_max, dim[2])));

  double r2 = 0.;
  for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k] * not_same_range[k];

  return r2;
}

#ifdef WITH_MPI

/* Add a group to the hash table. */
__attribute__((always_inline)) INLINE static void hashmap_add_group(
    const size_t group_id, const size_t group_offset, hashmap_t *map) {

  int created_new_element = 0;
  hashmap_value_t *offset =
      hashmap_get_new(map, group_id, &created_new_element);

  if (offset != NULL) {

    /* If the element is a new entry set its value. */
    if (created_new_element) {
      (*offset).value_st = group_offset;
    }
  } else
    error("Couldn't find key (%zu) or create new one.", group_id);
}

/* Find a group in the hash table. */
__attribute__((always_inline)) INLINE static size_t hashmap_find_group_offset(
    const size_t group_id, hashmap_t *map) {

  hashmap_value_t *group_offset = hashmap_get(map, group_id);

  if (group_offset == NULL)
    error("Couldn't find key (%zu) or create new one.", group_id);

  return (size_t)(*group_offset).value_st;
}

/* Compute send/recv offsets for MPI communication. */
__attribute__((always_inline)) INLINE static void fof_compute_send_recv_offsets(
    const int nr_nodes, int *sendcount, int **recvcount, int **sendoffset,
    int **recvoffset, size_t *nrecv) {

  /* Determine number of entries to receive */
  *recvcount = (int *)malloc(nr_nodes * sizeof(int));
  MPI_Alltoall(sendcount, 1, MPI_INT, *recvcount, 1, MPI_INT, MPI_COMM_WORLD);

  /* Compute send/recv offsets */
  *sendoffset = (int *)malloc(nr_nodes * sizeof(int));

  (*sendoffset)[0] = 0;
  for (int i = 1; i < nr_nodes; i++)
    (*sendoffset)[i] = (*sendoffset)[i - 1] + sendcount[i - 1];

  *recvoffset = (int *)malloc(nr_nodes * sizeof(int));

  (*recvoffset)[0] = 0;
  for (int i = 1; i < nr_nodes; i++)
    (*recvoffset)[i] = (*recvoffset)[i - 1] + (*recvcount)[i - 1];

  /* Allocate receive buffer */
  *nrecv = 0;
  for (int i = 0; i < nr_nodes; i++) (*nrecv) += (*recvcount)[i];
}

#endif /* WITH_MPI */

/**
 * @brief Perform a FOF search using union-find on a given leaf-cell
 *
 * @param props The properties fof the FOF scheme.
 * @param l_x2 The square of the FOF linking length.
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param c The #cell in which to perform FOF.
 */
void fof_search_self_cell(const struct fof_props *props, const double l_x2,
                          const struct gpart *const space_gparts,
                          const struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->split) error("Performing the FOF search at a non-leaf level!");
#endif

  const size_t count = c->grav.count;
  const struct gpart *gparts = c->grav.parts;

  /* Index of particles in the global group list */
  size_t *const group_index = props->group_index;

  /* Make a list of particle offsets into the global gparts array. */
  size_t *const offset = group_index + (ptrdiff_t)(gparts - space_gparts);

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank)
    error("Performing self FOF search on foreign cell.");
#endif

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count; i++) {

    const struct gpart *pi = &gparts[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

    /* Check whether we ignore this particle type altogether */
    if (gpart_is_ignorable(pi)) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    /* Find the root of pi. */
    size_t root_i = fof_find(offset[i], group_index);

    /* Get the nature of the linking */
    const int is_link_i = gpart_is_linkable(pi);

    for (size_t j = i + 1; j < count; j++) {

      const struct gpart *pj = &gparts[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Check whether we ignore this particle type altogether */
      if (gpart_is_ignorable(pj)) continue;

      /* Get the nature of the linking */
      const int is_link_j = gpart_is_linkable(pj);

      /* Both particles must be of the linking kind */
      if (!(is_link_i && is_link_j)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      /* Find the root of pj. */
      const size_t root_j = fof_find(offset[j], group_index);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute the pairwise distance */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) {

        /* Merge the groups` */
        fof_union(&root_i, root_j, group_index);
      }
    }
  }
}

/**
 * @brief Perform a FOF search using union-find between two cells
 *
 * @param props The properties fof the FOF scheme.
 * @param dim The dimension of the simulation volume.
 * @param l_x2 The square of the FOF linking length.
 * @param periodic Are we using periodic BCs?
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param ci The first #cell in which to perform FOF.
 * @param cj The second #cell in which to perform FOF.
 */
void fof_search_pair_cells(const struct fof_props *props, const double dim[3],
                           const double l_x2, const int periodic,
                           const struct gpart *const space_gparts,
                           const struct cell *restrict ci,
                           const struct cell *restrict cj) {

  const size_t count_i = ci->grav.count;
  const size_t count_j = cj->grav.count;
  const struct gpart *gparts_i = ci->grav.parts;
  const struct gpart *gparts_j = cj->grav.parts;

  /* Index of particles in the global group list */
  size_t *const group_index = props->group_index;

  /* Make a list of particle offsets into the global gparts array. */
  size_t *const offset_i = group_index + (ptrdiff_t)(gparts_i - space_gparts);
  size_t *const offset_j = group_index + (ptrdiff_t)(gparts_j - space_gparts);

#ifdef SWIFT_DEBUG_CHECKS
  if (offset_j > offset_i && (offset_j < offset_i + count_i))
    error("Overlapping cells");
  if (offset_i > offset_j && (offset_i < offset_j + count_j))
    error("Overlapping cells");
  if (ci->nodeID != cj->nodeID) error("Searching foreign cells!");
#endif

  /* Account for boundary conditions.*/
  double shift[3] = {0.0, 0.0, 0.0};

  /* Get the relative distance between the pairs, wrapping. */
  double diff[3];
  for (int k = 0; k < 3; k++) {
    diff[k] = cj->loc[k] - ci->loc[k];
    if (periodic && diff[k] < -dim[k] * 0.5)
      shift[k] = dim[k];
    else if (periodic && diff[k] > dim[k] * 0.5)
      shift[k] = -dim[k];
    else
      shift[k] = 0.0;
    diff[k] += shift[k];
  }

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count_i; i++) {

    const struct gpart *restrict pi = &gparts_i[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

    /* Check whether we ignore this particle type altogether */
    if (gpart_is_ignorable(pi)) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    const double pix = pi->x[0] - shift[0];
    const double piy = pi->x[1] - shift[1];
    const double piz = pi->x[2] - shift[2];

    /* Find the root of pi. */
    size_t root_i = fof_find(offset_i[i], group_index);

    /* Get the nature of the linking */
    const int is_link_i = gpart_is_linkable(pi);

    for (size_t j = 0; j < count_j; j++) {

      const struct gpart *restrict pj = &gparts_j[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Check whether we ignore this particle type altogether */
      if (gpart_is_ignorable(pj)) continue;

      /* Get the nature of the linking */
      const int is_link_j = gpart_is_linkable(pj);

      /* At least one of the particles has to be of linking type */
      if (!(is_link_i && is_link_j)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      /* Find the root of pj. */
      const size_t root_j = fof_find(offset_j[j], group_index);

      /* Skip particles in the same group. */
      if (root_i == root_j) continue;

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute pairwise distance (periodic BCs were accounted
       for by the shift vector) */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) {

        /* Merge the groups */
        fof_union(&root_i, root_j, group_index);
      }
    }
  }
}

#ifdef WITH_MPI

/**
 * @brief Add a local<->foreign pair in range to the list of links
 *
 * Possibly reallocates the local_group_links if we run out of space.
 */
static INLINE void add_foreign_link_to_list(
    int *local_link_count, int *group_links_size, struct fof_mpi **group_links,
    struct fof_mpi **local_group_links, const size_t root_i,
    const size_t root_j, const size_t size_i, const size_t size_j) {

  /* If the group_links array is not big enough re-allocate it. */
  if (*local_link_count + 1 > *group_links_size) {

    const int new_size = 2 * (*group_links_size);

    *group_links_size = new_size;

    (*group_links) = (struct fof_mpi *)realloc(
        *group_links, new_size * sizeof(struct fof_mpi));

    /* Reset the local pointer */
    (*local_group_links) = *group_links;

    message("Re-allocating local group links from %d to %d elements.",
            *local_link_count, new_size);

    if (new_size < 0) error("Overflow in size of list of foreign links");
  }

  /* Store the particle group properties for communication. */

  (*local_group_links)[*local_link_count].group_i = root_i;
  (*local_group_links)[*local_link_count].group_i_size = size_i;

  (*local_group_links)[*local_link_count].group_j = root_j;
  (*local_group_links)[*local_link_count].group_j_size = size_j;

  (*local_link_count)++;
}

#endif

/* Perform a FOF search between a local and foreign cell using the Union-Find
 * algorithm. Store any links found between particles.*/
void fof_search_pair_cells_foreign(
    const struct fof_props *props, const double dim[3], const double l_x2,
    const int periodic, const struct gpart *const space_gparts,
    const size_t nr_gparts, const struct cell *restrict ci,
    const struct cell *restrict cj, int *restrict link_count,
    struct fof_mpi **group_links, int *restrict group_links_size) {

#ifdef WITH_MPI
  const size_t count_i = ci->grav.count;
  const size_t count_j = cj->grav.count;
  const struct gpart *gparts_i = ci->grav.parts;
  const struct gpart_fof_foreign *gparts_j = cj->grav.parts_fof_foreign;

  /* Get local pointers */
  const size_t *restrict group_index = props->group_index;
  const size_t *restrict group_size = props->group_size;

  /* Values local to this function to avoid dereferencing */
  struct fof_mpi *local_group_links = *group_links;
  int local_link_count = *link_count;

  /* Make a list of particle offsets into the global gparts array. */
  const size_t *const offset_i =
      group_index + (ptrdiff_t)(gparts_i - space_gparts);

#ifdef SWIFT_DEBUG_CHECKS

  /* Check whether cells are local to the node. */
  const int ci_local = (ci->nodeID == engine_rank);
  const int cj_local = (cj->nodeID == engine_rank);

  if ((ci_local && cj_local) || (!ci_local && !cj_local))
    error(
        "FOF search of foreign cells called on two local cells or two foreign "
        "cells.");

  if (!ci_local) {
    error("Cell ci is not local!");
  }
  if (cj_local) {
    error("Cell cj is local!");
  }
#endif

  double shift[3] = {0.0, 0.0, 0.0};

  /* Get the relative distance between the pairs, wrapping. */
  for (int k = 0; k < 3; k++) {
    const double diff = cj->loc[k] - ci->loc[k];
    if (periodic && diff < -dim[k] / 2)
      shift[k] = dim[k];
    else if (periodic && diff > dim[k] / 2)
      shift[k] = -dim[k];
    else
      shift[k] = 0.0;
  }

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count_i; i++) {

    const struct gpart *pi = &gparts_i[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

    /* Check whether we ignore this particle type altogether */
    if (gpart_is_ignorable(pi)) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    const double pix = pi->x[0] - shift[0];
    const double piy = pi->x[1] - shift[1];
    const double piz = pi->x[2] - shift[2];

    /* Find the root of pi. */
    const size_t root_i =
        fof_find_global(offset_i[i] - node_offset, group_index, nr_gparts);

    /* Get the nature of the linking */
    const int is_link_i = gpart_is_linkable(pi);

    for (size_t j = 0; j < count_j; j++) {

      const struct gpart_fof_foreign *pj = &gparts_j[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Check whether we ignore this particle type altogether */
      if (gpart_foreign_is_ignorable(pj)) continue;

      /* Get the nature of the linking */
      const int is_link_j = gpart_foreign_is_linkable(pj);

      /* Only consider linkable<->linkable pairs */
      if (!(is_link_i && is_link_j)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute pairwise distance (periodic BCs were accounted
       for by the shift vector) */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) {

        /* Check that the links have not already been added to the list. */
        for (int l = 0; l < local_link_count; l++) {
          if (local_group_links[l].group_i == root_i &&
              local_group_links[l].group_j == pj->fof_data.group_id) {
            continue;
          }
        }

        /* Add a possible link to the list */
        add_foreign_link_to_list(
            &local_link_count, group_links_size, group_links,
            &local_group_links, root_i, pj->fof_data.group_id,
            group_size[root_i - node_offset], pj->fof_data.group_size);
      }
    }
  }

  /* Update the returned values */
  *link_count = local_link_count;

#else
  error("Calling MPI function in non-MPI mode.");
#endif
}

/**
 * @brief Recursively perform a union-find FOF between two cells.
 *
 * If cells are more distant than the linking length, we abort early.
 *
 * @param props The properties fof the FOF scheme.
 * @param dim The dimension of the space.
 * @param search_r2 the square of the FOF linking length.
 * @param periodic Are we using periodic BCs?
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param ci The first #cell in which to perform FOF.
 * @param cj The second #cell in which to perform FOF.
 */
void rec_fof_search_pair(const struct fof_props *props, const double dim[3],
                         const double search_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         struct cell *restrict ci, struct cell *restrict cj) {

  /* Find the shortest distance between cells, remembering to account for
   * boundary conditions. */
  const double r2 = cell_min_dist(ci, cj, dim);

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == cj) error("Pair FOF called on same cell!!!");
#endif

  /* Return if cells are out of range of each other. */
  if (r2 > search_r2) return;

  /* Recurse on both cells if they are both split. */
  if (ci->split && cj->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {

        for (int l = 0; l < 8; l++)
          if (cj->progeny[l] != NULL)
            rec_fof_search_pair(props, dim, search_r2, periodic, space_gparts,
                                ci->progeny[k], cj->progeny[l]);
      }
    }
  } else if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL)
        rec_fof_search_pair(props, dim, search_r2, periodic, space_gparts,
                            ci->progeny[k], cj);
    }
  } else if (cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL)
        rec_fof_search_pair(props, dim, search_r2, periodic, space_gparts, ci,
                            cj->progeny[k]);
    }
  } else {
    /* Perform FOF search between pairs of cells that are within the linking
     * length and not the same cell. */
    fof_search_pair_cells(props, dim, search_r2, periodic, space_gparts, ci,
                          cj);
  }
}
#ifdef WITH_MPI

/* Recurse on a pair of cells (one local, one foreign) and perform a FOF search
 * between cells that are within range. */
void rec_fof_search_pair_foreign(
    const struct fof_props *props, const double dim[3], const double search_r2,
    const int periodic, const struct gpart *const space_gparts,
    const size_t nr_gparts, const struct cell *ci, const struct cell *cj,
    int *restrict link_count, struct fof_mpi **group_links,
    int *restrict group_links_size) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == cj) error("Pair FOF called on same cell!!!");
  if (ci->nodeID == cj->nodeID) error("Fully local pair!");
#endif

  /* Find the shortest distance between cells, remembering to account for
   * boundary conditions. */
  const double r2 = cell_min_dist(ci, cj, dim);

  /* Return if cells are out of range of each other. */
  if (r2 > search_r2) return;

  /* Recurse on both cells if they are both split. */
  if (ci->split && cj->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {

        for (int l = 0; l < 8; l++)
          if (cj->progeny[l] != NULL)
            rec_fof_search_pair_foreign(props, dim, search_r2, periodic,
                                        space_gparts, nr_gparts, ci->progeny[k],
                                        cj->progeny[l], link_count, group_links,
                                        group_links_size);
      }
    }
  } else if (ci->split) {

    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL)
        rec_fof_search_pair_foreign(props, dim, search_r2, periodic,
                                    space_gparts, nr_gparts, ci->progeny[k], cj,
                                    link_count, group_links, group_links_size);
    }
  } else if (cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL)
        rec_fof_search_pair_foreign(props, dim, search_r2, periodic,
                                    space_gparts, nr_gparts, ci, cj->progeny[k],
                                    link_count, group_links, group_links_size);
    }
  } else {
    /* Perform FOF search between pairs of cells that are within the linking
     * length and not the same cell. */
    fof_search_pair_cells_foreign(props, dim, search_r2, periodic, space_gparts,
                                  nr_gparts, ci, cj, link_count, group_links,
                                  group_links_size);
  }
}

#endif /* WITH_MPI */

/**
 * @brief Recursively perform a union-find FOF on a cell.
 *
 * @param props The properties fof the FOF scheme.
 * @param dim The dimension of the space.
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param search_r2 the square of the FOF linking length.
 * @param periodic Are we using periodic BCs?
 * @param c The #cell in which to perform FOF.
 */
void rec_fof_search_self(const struct fof_props *props, const double dim[3],
                         const double search_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         struct cell *c) {

  /* Recurse? */
  if (c->split) {

    /* Loop over all progeny. Perform pair and self recursion on progenies.*/
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {

        rec_fof_search_self(props, dim, search_r2, periodic, space_gparts,
                            c->progeny[k]);

        for (int l = k + 1; l < 8; l++)
          if (c->progeny[l] != NULL)
            rec_fof_search_pair(props, dim, search_r2, periodic, space_gparts,
                                c->progeny[k], c->progeny[l]);
      }
    }
  }
  /* Otherwise, compute self-interaction. */
  else
    fof_search_self_cell(props, search_r2, space_gparts, c);
}

/**
 * @brief Perform the attaching operation using union-find on a given leaf-cell
 *
 * @param props The properties fof the FOF scheme.
 * @param l_x2 The square of the FOF linking length.
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param nr_gparts The number of #gpart in the local #space structure.
 * @param c The #cell in which to perform FOF.
 */
void fof_attach_self_cell(const struct fof_props *props, const double l_x2,
                          const struct gpart *const space_gparts,
                          const size_t nr_gparts, const struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->split) error("Performing the FOF search at a non-leaf level!");
#endif

  const size_t count = c->grav.count;
  struct gpart *gparts = (struct gpart *)c->grav.parts;

  /* Distances of particles in the global list */
  float *const offset_dist =
      props->distance_to_link + (ptrdiff_t)(gparts - space_gparts);

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank)
    error("Performing self FOF search on foreign cell.");
#endif

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count; i++) {

    struct gpart *pi = &gparts[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

    /* Check whether we ignore this particle type altogether */
    if (gpart_is_ignorable(pi)) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    /* Get the nature of the linking */
    const int is_link_i = gpart_is_linkable(pi);
    const int is_attach_i = gpart_is_attachable(pi);

#ifdef SWIFT_DEBUG_CHECKS
    if (is_link_i && is_attach_i)
      error("Particle cannot be both linkable and attachable!");
#endif

    for (size_t j = i + 1; j < count; j++) {

      struct gpart *pj = &gparts[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Check whether we ignore this particle type altogether */
      if (gpart_is_ignorable(pj)) continue;

      /* Get the nature of the linking */
      const int is_link_j = gpart_is_linkable(pj);
      const int is_attach_j = gpart_is_attachable(pj);

#ifdef SWIFT_DEBUG_CHECKS
      if (is_link_j && is_attach_j)
        error("Particle cannot be both linkable and attachable!");
#endif

      /* We only want link<->attach pairs */
      if (is_attach_i && is_attach_j) continue;
      if (is_link_i && is_link_j) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute the pairwise distance */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) {

        /* Now that we are within the linking length,
         * decide what to do based on linking types */

        if (is_link_i && is_link_j) {
#ifdef SWIFT_DEBUG_CHECKS
          error("Fundamental logic error!");
#endif
        } else if (is_link_i && is_attach_j) {

          /* We got a linkable and an attachable.
           * See whether it is closer and if so re-link.
           * This is safe to do as the attachables are never roots and
           * nothing is attached to them */
          const float dist = sqrtf(r2);
          if (dist < offset_dist[j]) {

            /* Store the new min dist */
            offset_dist[j] = dist;

            /* Store the current best root */
            pj->fof_data.group_id = pi->fof_data.group_id;
          }

        } else if (is_link_j && is_attach_i) {

          /* We got a linkable and an attachable.
           * See whether it is closer and if so re-link.
           * This is safe to do as the attachables are never roots and
           * nothing is attached to them */
          const float dist = sqrtf(r2);
          if (dist < offset_dist[i]) {

            /* Store the new min dist */
            offset_dist[i] = dist;

            /* Store the current best root */
            pi->fof_data.group_id = pj->fof_data.group_id;
          }

        } else {
#ifdef SWIFT_DEBUG_CHECKS
          error("Fundamental logic error!");
#endif
        }
      }
    }
  }
}

/**
 * @brief Perform the attaching operation using union-find between two cells
 *
 * @param props The properties fof the FOF scheme.
 * @param dim The dimension of the simulation volume.
 * @param l_x2 The square of the FOF linking length.
 * @param periodic Are we using periodic BCs?
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param nr_gparts The number of #gpart in the local #space structure.
 * @param ci The first #cell in which to perform FOF.
 * @param cj The second #cell in which to perform FOF.
 */
void fof_attach_pair_cells_both_local(const struct fof_props *props,
                                      const double dim[3], const double l_x2,
                                      const int periodic,
                                      const struct gpart *const space_gparts,
                                      const size_t nr_gparts,
                                      const struct cell *restrict ci,
                                      const struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("ci not local!");
  if (cj->nodeID != engine_rank) error("cj not local!");
#endif

  const size_t count_i = ci->grav.count;
  const size_t count_j = cj->grav.count;
  struct gpart *gparts_i = (struct gpart *)ci->grav.parts;
  struct gpart *gparts_j = (struct gpart *)cj->grav.parts;

  /* Distances of particles in the global list */
  float *const offset_dist_i =
      props->distance_to_link + (ptrdiff_t)(gparts_i - space_gparts);
  float *const offset_dist_j =
      props->distance_to_link + (ptrdiff_t)(gparts_j - space_gparts);

  /* Account for boundary conditions.*/
  double shift[3] = {0.0, 0.0, 0.0};

  /* Get the relative distance between the pairs, wrapping. */
  double diff[3];
  for (int k = 0; k < 3; k++) {
    diff[k] = cj->loc[k] - ci->loc[k];
    if (periodic && diff[k] < -dim[k] * 0.5)
      shift[k] = dim[k];
    else if (periodic && diff[k] > dim[k] * 0.5)
      shift[k] = -dim[k];
    else
      shift[k] = 0.0;
    diff[k] += shift[k];
  }

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count_i; i++) {

    struct gpart *restrict pi = &gparts_i[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

    /* Check whether we ignore this particle type altogether */
    if (gpart_is_ignorable(pi)) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    const double pix = pi->x[0] - shift[0];
    const double piy = pi->x[1] - shift[1];
    const double piz = pi->x[2] - shift[2];

    /* Get the nature of the linking */
    const int is_link_i = gpart_is_linkable(pi);
    const int is_attach_i = gpart_is_attachable(pi);

#ifdef SWIFT_DEBUG_CHECKS
    if (is_link_i && is_attach_i)
      error("Particle cannot be both linkable and attachable!");
#endif

    for (size_t j = 0; j < count_j; j++) {

      struct gpart *restrict pj = &gparts_j[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Check whether we ignore this particle type altogether */
      if (gpart_is_ignorable(pj)) continue;

      /* Get the nature of the linking */
      const int is_link_j = gpart_is_linkable(pj);
      const int is_attach_j = gpart_is_attachable(pj);

#ifdef SWIFT_DEBUG_CHECKS
      if (is_link_j && is_attach_j)
        error("Particle cannot be both linkable and attachable!");
#endif

      /* We only want link<->attach pairs */
      if (is_attach_i && is_attach_j) continue;
      if (is_link_i && is_link_j) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute pairwise distance (periodic BCs were accounted
       for by the shift vector) */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) {

        /* Now that we are within the linking length,
         * decide what to do based on linking types */

        if (is_link_i && is_link_j) {

#ifdef SWIFT_DEBUG_CHECKS
          error("Fundamental logic error!");
#endif

        } else if (is_link_i && is_attach_j) {

          /* We got a linkable and an attachable.
           * See whether it is closer and if so re-link.
           * This is safe to do as the attachables are never roots and
           * nothing is attached to them */
          const float dist = sqrtf(r2);

          if (dist < offset_dist_j[j]) {

            /* Store the new min dist */
            offset_dist_j[j] = dist;

            /* Store the current best root */
            pj->fof_data.group_id = pi->fof_data.group_id;
          }

        } else if (is_link_j && is_attach_i) {

          /* We got a linkable and an attachable.
           * See whether it is closer and if so re-link.
           * This is safe to do as the attachables are never roots and
           * nothing is attached to them */
          const float dist = sqrtf(r2);

          if (dist < offset_dist_i[i]) {

            /* Store the new min dist */
            offset_dist_i[i] = dist;

            /* Store the current best root */
            pi->fof_data.group_id = pj->fof_data.group_id;
          }

        } else {
#ifdef SWIFT_DEBUG_CHECKS
          error("Fundamental logic error!");
#endif
        }
      }
    }
  }
}

void fof_attach_pair_cells_ci_local(const struct fof_props *props,
                                    const double dim[3], const double l_x2,
                                    const int periodic,
                                    const struct gpart *const space_gparts,
                                    const size_t nr_gparts,
                                    const struct cell *restrict ci,
                                    const struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("ci not local!");
  if (cj->nodeID == engine_rank) error("cj local!");
#endif

  const size_t count_i = ci->grav.count;
  const size_t count_j = cj->grav.count;
  struct gpart *gparts_i = ci->grav.parts;
  struct gpart_fof_foreign *gparts_j = cj->grav.parts_fof_foreign;

  /* Distances of particles in the global list */
  float *const offset_dist_i =
      props->distance_to_link + (ptrdiff_t)(gparts_i - space_gparts);

  /* Account for boundary conditions.*/
  double shift[3] = {0.0, 0.0, 0.0};

  /* Get the relative distance between the pairs, wrapping. */
  double diff[3];
  for (int k = 0; k < 3; k++) {
    diff[k] = cj->loc[k] - ci->loc[k];
    if (periodic && diff[k] < -dim[k] * 0.5)
      shift[k] = dim[k];
    else if (periodic && diff[k] > dim[k] * 0.5)
      shift[k] = -dim[k];
    else
      shift[k] = 0.0;
    diff[k] += shift[k];
  }

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count_i; i++) {

    struct gpart *restrict pi = &gparts_i[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

    /* Check whether we ignore this particle type altogether */
    if (gpart_is_ignorable(pi)) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    const double pix = pi->x[0] - shift[0];
    const double piy = pi->x[1] - shift[1];
    const double piz = pi->x[2] - shift[2];

    /* Get the nature of the linking */
    const int is_link_i = gpart_is_linkable(pi);
    const int is_attach_i = gpart_is_attachable(pi);

#ifdef SWIFT_DEBUG_CHECKS
    if (is_link_i && is_attach_i)
      error("Particle cannot be both linkable and attachable!");
#endif

    for (size_t j = 0; j < count_j; j++) {

      struct gpart_fof_foreign *restrict pj = &gparts_j[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Check whether we ignore this particle type altogether */
      if (gpart_foreign_is_ignorable(pj)) continue;

      /* Get the nature of the linking */
      const int is_link_j = gpart_foreign_is_linkable(pj);
      const int is_attach_j = gpart_foreign_is_attachable(pj);

#ifdef SWIFT_DEBUG_CHECKS
      if (is_link_j && is_attach_j)
        error("Particle cannot be both linkable and attachable!");
#endif

      /* We only want link<->attach pairs */
      if (is_attach_i && is_attach_j) continue;
      if (is_link_i && is_link_j) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute pairwise distance (periodic BCs were accounted
       for by the shift vector) */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) {

        /* Now that we are within the linking length,
         * decide what to do based on linking types */

        if (is_link_i && is_link_j) {

#ifdef SWIFT_DEBUG_CHECKS
          error("Fundamental logic error!");
#endif

        } else if (is_link_i && is_attach_j) {

          /* Nothing to do here. The reverse action will be done
             in the converse call on the other node. */

        } else if (is_link_j && is_attach_i) {

          /* We got a linkable and an attachable.
           * See whether it is closer and if so re-link.
           * This is safe to do as the attachables are never roots and
           * nothing is attached to them */
          const float dist = sqrtf(r2);

          if (dist < offset_dist_i[i]) {

            /* Store the new min dist */
            offset_dist_i[i] = dist;

            /* Store the current best root */
            pi->fof_data.group_id = pj->fof_data.group_id;
          }
        } else {
#ifdef SWIFT_DEBUG_CHECKS
          error("Fundamental logic error!");
#endif
        }
      }
    }
  }
}

void fof_attach_pair_cells_cj_local(const struct fof_props *props,
                                    const double dim[3], const double l_x2,
                                    const int periodic,
                                    const struct gpart *const space_gparts,
                                    const size_t nr_gparts,
                                    const struct cell *restrict ci,
                                    const struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID == engine_rank) error("ci local!");
  if (cj->nodeID != engine_rank) error("cj not local!");
#endif

  const size_t count_i = ci->grav.count;
  const size_t count_j = cj->grav.count;
  struct gpart_fof_foreign *gparts_i = ci->grav.parts_fof_foreign;
  struct gpart *gparts_j = cj->grav.parts;

  /* Distances of particles in the global list */
  float *const offset_dist_j =
      props->distance_to_link + (ptrdiff_t)(gparts_j - space_gparts);

  /* Account for boundary conditions.*/
  double shift[3] = {0.0, 0.0, 0.0};

  /* Get the relative distance between the pairs, wrapping. */
  double diff[3];
  for (int k = 0; k < 3; k++) {
    diff[k] = cj->loc[k] - ci->loc[k];
    if (periodic && diff[k] < -dim[k] * 0.5)
      shift[k] = dim[k];
    else if (periodic && diff[k] > dim[k] * 0.5)
      shift[k] = -dim[k];
    else
      shift[k] = 0.0;
    diff[k] += shift[k];
  }

  /* Loop over particles and find which particles belong in the same group. */
  for (size_t i = 0; i < count_i; i++) {

    struct gpart_fof_foreign *restrict pi = &gparts_i[i];

    /* Ignore inhibited particles */
    if (pi->time_bin >= time_bin_inhibited) continue;

    /* Check whether we ignore this particle type altogether */
    if (gpart_foreign_is_ignorable(pi)) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (pi->ti_drift != ti_current)
      error("Running FOF on an un-drifted particle!");
#endif

    const double pix = pi->x[0] - shift[0];
    const double piy = pi->x[1] - shift[1];
    const double piz = pi->x[2] - shift[2];

    /* Get the nature of the linking */
    const int is_link_i = gpart_foreign_is_linkable(pi);
    const int is_attach_i = gpart_foreign_is_attachable(pi);

#ifdef SWIFT_DEBUG_CHECKS
    if (is_link_i && is_attach_i)
      error("Particle cannot be both linkable and attachable!");
#endif

    for (size_t j = 0; j < count_j; j++) {

      struct gpart *restrict pj = &gparts_j[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Check whether we ignore this particle type altogether */
      if (gpart_is_ignorable(pj)) continue;

      /* Get the nature of the linking */
      const int is_link_j = gpart_is_linkable(pj);
      const int is_attach_j = gpart_is_attachable(pj);

#ifdef SWIFT_DEBUG_CHECKS
      if (is_link_j && is_attach_j)
        error("Particle cannot be both linkable and attachable!");
#endif

      /* We only want link<->attach pairs */
      if (is_attach_i && is_attach_j) continue;
      if (is_link_i && is_link_j) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != ti_current)
        error("Running FOF on an un-drifted particle!");
#endif

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute pairwise distance (periodic BCs were accounted
       for by the shift vector) */
      float dx[3], r2 = 0.0f;
      dx[0] = pix - pjx;
      dx[1] = piy - pjy;
      dx[2] = piz - pjz;

      for (int k = 0; k < 3; k++) r2 += dx[k] * dx[k];

      /* Hit or miss? */
      if (r2 < l_x2) {

        /* Now that we are within the linking length,
         * decide what to do based on linking types */

        if (is_link_i && is_link_j) {

#ifdef SWIFT_DEBUG_CHECKS
          error("Fundamental logic error!");
#endif

        } else if (is_link_i && is_attach_j) {

          /* We got a linkable and an attachable.
           * See whether it is closer and if so re-link.
           * This is safe to do as the attachables are never roots and
           * nothing is attached to them */

          const float dist = sqrtf(r2);

          if (dist < offset_dist_j[j]) {

            /* Store the new min dist */
            offset_dist_j[j] = dist;

            /* Store the current best root */
            pj->fof_data.group_id = pi->fof_data.group_id;
          }

        } else if (is_link_j && is_attach_i) {

          /* Nothing to do here. The reverse action will be done
             in the converse call on the other node. */

        } else {
#ifdef SWIFT_DEBUG_CHECKS
          error("Fundamental logic error!");
#endif
        }
      }
    }
  }
}

/**
 * @brief Recursively perform a union-find attaching between two cells.
 *
 * If cells are more distant than the linking length, we abort early.
 *
 * @param props The properties fof the FOF scheme.
 * @param dim The dimension of the space.
 * @param attach_r2 the square of the FOF linking length.
 * @param periodic Are we using periodic BCs?
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param nr_gparts The number of #gpart in the local #space structure.
 * @param ci The first #cell in which to perform FOF.
 * @param cj The second #cell in which to perform FOF.
 * @param ci_local Is the #cell ci on the local MPI rank?
 * @param cj_local Is the #cell cj on the local MPI rank?
 */
void rec_fof_attach_pair(const struct fof_props *props, const double dim[3],
                         const double attach_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         const size_t nr_gparts, struct cell *restrict ci,
                         struct cell *restrict cj, const int ci_local,
                         const int cj_local) {

  /* Find the shortest distance between cells, remembering to account for
   * boundary conditions. */
  const double r2 = cell_min_dist(ci, cj, dim);

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == cj) error("Pair FOF called on same cell!!!");
#endif

  /* Return if cells are out of range of each other. */
  if (r2 > attach_r2) return;

  /* Recurse on both cells if they are both split. */
  if (ci->split && cj->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {

        for (int l = 0; l < 8; l++)
          if (cj->progeny[l] != NULL)
            rec_fof_attach_pair(props, dim, attach_r2, periodic, space_gparts,
                                nr_gparts, ci->progeny[k], cj->progeny[l],
                                ci_local, cj_local);
      }
    }
  } else if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL)
        rec_fof_attach_pair(props, dim, attach_r2, periodic, space_gparts,
                            nr_gparts, ci->progeny[k], cj, ci_local, cj_local);
    }
  } else if (cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL)
        rec_fof_attach_pair(props, dim, attach_r2, periodic, space_gparts,
                            nr_gparts, ci, cj->progeny[k], ci_local, cj_local);
    }
  } else {
    /* Perform FOF attach between pairs of cells that are within the linking
     * length and not the same cell. */
    if (ci_local && cj_local) {

      fof_attach_pair_cells_both_local(props, dim, attach_r2, periodic,
                                       space_gparts, nr_gparts, ci, cj);
    } else if (ci_local && !cj_local) {
      fof_attach_pair_cells_ci_local(props, dim, attach_r2, periodic,
                                     space_gparts, nr_gparts, ci, cj);

    } else if (!ci_local && cj_local) {
      fof_attach_pair_cells_cj_local(props, dim, attach_r2, periodic,
                                     space_gparts, nr_gparts, ci, cj);
    } else {
      error("Error in the recursion logic");
    }
  }
}

/**
 * @brief Recursively perform a the attaching operation on a cell.
 *
 * @param props The properties fof the FOF scheme.
 * @param dim The dimension of the space.
 * @param attach_r2 the square of the FOF linking length.
 * @param periodic Are we using periodic BCs?
 * @param space_gparts The start of the #gpart array in the #space structure.
 * @param nr_gparts The number of #gpart in the local #space structure.
 * @param c The #cell in which to perform FOF.
 */
void rec_fof_attach_self(const struct fof_props *props, const double dim[3],
                         const double attach_r2, const int periodic,
                         const struct gpart *const space_gparts,
                         const size_t nr_gparts, struct cell *c) {

  /* Recurse? */
  if (c->split) {

    /* Loop over all progeny. Perform pair and self recursion on progenies.*/
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {

        rec_fof_attach_self(props, dim, attach_r2, periodic, space_gparts,
                            nr_gparts, c->progeny[k]);

        for (int l = k + 1; l < 8; l++)
          if (c->progeny[l] != NULL)
            rec_fof_attach_pair(props, dim, attach_r2, periodic, space_gparts,
                                nr_gparts, c->progeny[k], c->progeny[l],
                                /*ci_local=*/1,
                                /*cj_local=*/1);
      }
    }
  } else {

    /* Otherwise, compute self-interaction. */
    fof_attach_self_cell(props, attach_r2, space_gparts, nr_gparts, c);
  }
}

/* Mapper function to atomically update the group size array. */
void fof_update_group_size_mapper(hashmap_key_t key, hashmap_value_t *value,
                                  void *data) {

  size_t *group_size = (size_t *)data;

  /* Use key to index into group size array. */
  atomic_add(&group_size[key], value->value_st);
}

/**
 * @brief Mapper function to calculate the group sizes.
 *
 * @param map_data An array of #gpart%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void fof_calc_group_size_mapper(void *map_data, int num_elements,
                                void *extra_data) {

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  struct gpart *gparts = (struct gpart *)map_data;
  size_t *restrict group_index = s->e->fof_properties->group_index;
  size_t *restrict group_size = s->e->fof_properties->group_size;

  /* Offset into gparts array. */
  const ptrdiff_t gparts_offset = (ptrdiff_t)(gparts - s->gparts);
  size_t *const group_index_offset = group_index + gparts_offset;

  /* Create hash table. */
  hashmap_t map;
  hashmap_init(&map);

  for (int ind = 0; ind < num_elements; ind++) {

    const hashmap_key_t root =
        (hashmap_key_t)fof_find(group_index_offset[ind], group_index);
    const size_t gpart_index = gparts_offset + ind;

    /* Only add particles which aren't the root of a group. Stops groups of size
     * 1 being added to the hash table. */
    if (root != gpart_index) {
      hashmap_value_t *size = hashmap_get(&map, root);

      if (size != NULL)
        (*size).value_st++;
      else
        error("Couldn't find key (%zu) or create new one.", root);
    }
  }

  /* Update the group size array. */
  if (map.size > 0)
    hashmap_iterate(&map, fof_update_group_size_mapper, group_size);

  hashmap_free(&map);
}

/**
 * @brief Mapper function to perform FOF search.
 *
 * @param map_data An array of #cell pair indices.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to a #space.
 */
void fof_find_foreign_links_mapper(void *map_data, int num_elements,
                                   void *extra_data) {

#ifdef WITH_MPI

  /* Retrieve mapped data. */
  struct space *s = (struct space *)extra_data;
  const int periodic = s->periodic;
  const size_t nr_gparts = s->nr_gparts;
  const struct gpart *const gparts = s->gparts;
  const struct engine *e = s->e;
  struct fof_props *props = e->fof_properties;
  struct cell_pair_indices *cell_pairs = (struct cell_pair_indices *)map_data;

  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = props->l_x2;

  /* Store links in an array local to this thread. */
  int local_link_count = 0;
  int local_group_links_size = props->group_links_size / e->nr_threads;

  /* Init the local group links buffer. */
  struct fof_mpi *local_group_links = (struct fof_mpi *)swift_calloc(
      "fof_local_group_links", sizeof(struct fof_mpi), local_group_links_size);
  if (local_group_links == NULL)
    error("Failed to allocate temporary group links buffer.");

  /* Loop over all pairs of local and foreign cells, recurse then perform a
   * FOF search. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the local and foreign cells to recurse on. */
    const struct cell *restrict local_cell = cell_pairs[ind].local;
    const struct cell *restrict foreign_cell = cell_pairs[ind].foreign;

    rec_fof_search_pair_foreign(props, dim, search_r2, periodic, gparts,
                                nr_gparts, local_cell, foreign_cell,
                                &local_link_count, &local_group_links,
                                &local_group_links_size);
  }

  /* Add links found by this thread to the global link list. */
  /* Lock to prevent race conditions while adding to the global link list.*/
  if (lock_lock(&s->lock) == 0) {

    /* Get pointers to global arrays. */
    int *restrict group_links_size = &props->group_links_size;
    int *restrict group_link_count = &props->group_link_count;
    struct fof_mpi **group_links = &props->group_links;

    /* If the global group_links array is not big enough re-allocate it. */
    if (*group_link_count + local_link_count > *group_links_size) {

      const int old_size = *group_links_size;
      const int new_size =
          max(*group_link_count + local_link_count, 2 * old_size);

      (*group_links) = (struct fof_mpi *)realloc(
          *group_links, new_size * sizeof(struct fof_mpi));

      *group_links_size = new_size;

      message("Re-allocating global group links from %d to %d elements.",
              old_size, new_size);
    }

    /* Copy the local links to the global list. */
    for (int i = 0; i < local_link_count; i++) {

      int found = 0;

      /* Check that the links have not already been added to the list by another
       * thread. */
      for (int l = 0; l < *group_link_count; l++) {
        if ((*group_links)[l].group_i == local_group_links[i].group_i &&
            (*group_links)[l].group_j == local_group_links[i].group_j) {
          found = 1;
          break;
        }
      }

      if (!found) {

        (*group_links)[*group_link_count].group_i =
            local_group_links[i].group_i;
        (*group_links)[*group_link_count].group_i_size =
            local_group_links[i].group_i_size;

        (*group_links)[*group_link_count].group_j =
            local_group_links[i].group_j;
        (*group_links)[*group_link_count].group_j_size =
            local_group_links[i].group_j_size;

        (*group_link_count) = (*group_link_count) + 1;
      }
    }
  }

  /* Release lock. */
  if (lock_unlock(&s->lock) != 0) error("Failed to unlock the space");

  swift_free("fof_local_group_links", local_group_links);
#endif
}

/**
 * @brief Compute the group properties for all groups.
 *
 * At the end of the process, all MPI ranks have all the information
 * about every group. We compute:
 * - Group size,
 * - Group total mass,
 * - Group centre of mass,
 * - Maximal gas particle density,
 * - Whether a group has a BH particle or not.
 *
 * TODO: Possible improvement is to delay the allocation of the CoM array
 * until after the min/max density arrays have played their role.
 * TODO: Can the loop over particles be threadpoolized?
 */
void fof_calc_group_mass(struct fof_props *props, const struct space *s,
                         int *restrict number_of_local_seeds,
                         int *restrict number_of_global_seeds) {

  const size_t nr_gparts = s->nr_gparts;
  const struct gpart *gparts = s->gparts;
  const struct part *parts = s->parts;
  const size_t group_id_default = props->group_id_default;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double seed_halo_mass = props->seed_halo_mass;

  /* Direct pointers to the arrays */
  long long *final_group_size = props->final_group_size;
  double *group_mass = props->group_mass;
  double *centre_of_mass = props->group_centre_of_mass;
  char *has_black_hole = props->has_black_hole;
  float *max_part_density = props->max_part_density;

  /* Temporary arrays to help with the CoMs */
  float *max_positions, *min_positions;
  if (swift_memalign("fof_group_max_position", (void **)&max_positions, 32,
                     props->num_groups * 3 * sizeof(double)) != 0)
    error("Unable to allocate memory for the max positions");
  if (swift_memalign("fof_group_min_position", (void **)&min_positions, 32,
                     props->num_groups * 3 * sizeof(double)) != 0)
    error("Unable to allocate memory for the min positions");

  /* Initialise the min/max arrays to the limits */
  for (size_t i = 0; i < 3 * (size_t)props->num_groups; ++i) {
    min_positions[i] = DBL_MAX;
    max_positions[i] = -DBL_MAX;
  }

  /* Collect information about the local particles and update the array of
   * properties. Recall the array is as big as all the haloes accross
   * all domains */
  for (size_t i = 0; i < nr_gparts; i++) {

    /* Ignore inhibited particles */
    if (gparts[i].time_bin >= time_bin_inhibited) continue;

    /* Check whether we ignore this particle type altogether */
    if (gpart_is_ignorable(&gparts[i])) continue;

    /* Ignore particles not in groups */
    if (gparts[i].fof_data.group_id == group_id_default) continue;

    if (gparts[i].fof_data.group_id > (size_t)props->num_groups)
      error("Found an invalid group ID!");

    /* Entry into the global list of group properties */
    const size_t index = gparts[i].fof_data.group_id - 1;

    /********************
     * We know in which group this particle is: compute props
     ********************/

    /* Count the number of particles */
    final_group_size[index]++;

    /* Add to the total mass */
    group_mass[index] += gparts[i].mass;

    /* Get the min/max position along each axis */
    min_positions[index * 3 + 0] =
        fmin(min_positions[index * 3 + 0], gparts[i].x[0]);
    min_positions[index * 3 + 1] =
        fmin(min_positions[index * 3 + 1], gparts[i].x[1]);
    min_positions[index * 3 + 2] =
        fmin(min_positions[index * 3 + 2], gparts[i].x[2]);
    max_positions[index * 3 + 0] =
        fmax(max_positions[index * 3 + 0], gparts[i].x[0]);
    max_positions[index * 3 + 1] =
        fmax(max_positions[index * 3 + 1], gparts[i].x[1]);
    max_positions[index * 3 + 2] =
        fmax(max_positions[index * 3 + 2], gparts[i].x[2]);

    /* Check whether there is a black hole */
    if (gparts[i].type == swift_type_black_hole) {
      has_black_hole[index] = 1;
    }

    /* Idntify the densest gas particle in the group */
    if (gparts[i].type == swift_type_gas) {
      const size_t gas_index = -gparts[i].id_or_neg_offset;
      const float rho_com = hydro_get_comoving_density(&parts[gas_index]);
      max_part_density[index] = fmaxf(rho_com, max_part_density[index]);
    }
  }

#ifdef WITH_MPI
  /* Now all-reduce the local fragments so that everyone has a full catalog */
  MPI_Allreduce(MPI_IN_PLACE, final_group_size, props->num_groups,
                MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, group_mass, props->num_groups, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, min_positions, 3 * props->num_groups, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, max_positions, 3 * props->num_groups, MPI_DOUBLE,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, has_black_hole, props->num_groups, MPI_CHAR,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, max_part_density, props->num_groups, MPI_FLOAT,
                MPI_MAX, MPI_COMM_WORLD);
#endif

  *number_of_local_seeds = 0;
  *number_of_global_seeds = 0;

  /* We can now do a second pass to compute the centre of mass
   * Because of periodic BCs, we need to shift the positions in cases where the
   * (max-min) is larger than 1/2 of the box size */
  for (size_t i = 0; i < nr_gparts; i++) {

    /* Ignore inhibited particles */
    if (gparts[i].time_bin >= time_bin_inhibited) continue;

    /* Check whether we ignore this particle type altogether */
    if (gpart_is_ignorable(&gparts[i])) continue;

    /* Ignore particles not in groups */
    if (gparts[i].fof_data.group_id == group_id_default) continue;

    /* Entry into the global list of group properties */
    const size_t index = gparts[i].fof_data.group_id - 1;

    const int halo_is_on_edge[3] = {
        max_positions[index * 3 + 0] - min_positions[index * 3 + 0] >
            0.5 * dim[0],
        max_positions[index * 3 + 1] - min_positions[index * 3 + 1] >
            0.5 * dim[1],
        max_positions[index * 3 + 2] - min_positions[index * 3 + 2] >
            0.5 * dim[2]};

    /* Get particle position, including necessary wrapping */
    double x[3] = {gparts[i].x[0], gparts[i].x[1], gparts[i].x[2]};
    if (periodic) {
      for (int k = 0; k < 3; k++) {
        if (halo_is_on_edge[k]) {
          x[k] = box_wrap(x[k] + 0.5 * dim[k], 0., dim[k]);
        }
      }
    }

    /* Centre of mass */
    centre_of_mass[index * 3 + 0] += gparts[i].mass * x[0];
    centre_of_mass[index * 3 + 1] += gparts[i].mass * x[1];
    centre_of_mass[index * 3 + 2] += gparts[i].mass * x[2];

    /* Should we seed a BH in this group? */
    if (!has_black_hole[index] && group_mass[index] > seed_halo_mass) {

      /* Is this a gas particle? */
      if (gparts[i].type == swift_type_gas) {
        const size_t gas_index = -gparts[i].id_or_neg_offset;
        const float rho_com = hydro_get_comoving_density(&parts[gas_index]);

        /* Is this the gas paricle which is the densest? */
        if (rho_com == max_part_density[index]) {
          (*number_of_local_seeds)++;
        }
      }
    }
  }

#ifdef WITH_MPI
  /* Now all-reduce the CoMs so that everyone has a full catalog */
  MPI_Allreduce(MPI_IN_PLACE, centre_of_mass, 3 * props->num_groups, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(number_of_local_seeds, number_of_global_seeds, 1, MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);
#else
  *number_of_global_seeds = *number_of_local_seeds;
#endif

  /* Finalise the operations on the group catalog */
  for (size_t i = 0; i < (size_t)props->num_groups; ++i) {

    centre_of_mass[i * 3 + 0] /= group_mass[i];
    centre_of_mass[i * 3 + 1] /= group_mass[i];
    centre_of_mass[i * 3 + 2] /= group_mass[i];

    const int halo_is_on_edge[3] = {
        max_positions[i * 3 + 0] - min_positions[i * 3 + 0] > 0.5 * dim[0],
        max_positions[i * 3 + 1] - min_positions[i * 3 + 1] > 0.5 * dim[1],
        max_positions[i * 3 + 2] - min_positions[i * 3 + 2] > 0.5 * dim[2]};

    /* Undo the half-box periodic shift */
    if (periodic) {
      for (int k = 0; k < 3; k++) {
        if (halo_is_on_edge[k]) {
          centre_of_mass[i * 3 + k] -= 0.5 * dim[k];
        }
        centre_of_mass[i * 3 + k] =
            box_wrap(centre_of_mass[i * 3 + k], 0., dim[k]);
      }
    }

    /* Relabel the group ids */
    props->final_group_index[i] = i + 1;
  }

  /* Free temporary arrays */
  swift_free("max_positions", max_positions);
  swift_free("min_positions", min_positions);
}

/**
 * @brief Seed black holes from gas particles in the haloes on the local MPI
 * rank that passed the criteria.
 *
 * @param props The properties of the FOF scheme.
 * @param bh_props The properties of the black hole scheme.
 * @param constants The physical constants.
 * @param cosmo The cosmological model.
 * @param s The @space we act on.
 * @param number_of_local_seeds Number of BHs to create on this rank.
 * @param number_of_global_seeds Number of BHs to create in total.
 */
void fof_seed_black_holes(const struct fof_props *props,
                          const struct black_holes_props *bh_props,
                          const struct phys_const *constants,
                          const struct cosmology *cosmo, struct space *s,
                          const int number_of_local_seeds,
                          const int number_of_global_seeds) {

  if (engine_rank == 0)
    message("Seeding %d black hole(s)", number_of_global_seeds);

  /* Anything to do this time on this rank? */
  if (number_of_local_seeds == 0) return;

  /* Do we need to reallocate the black hole array for the new particles? */
  if (s->nr_bparts + number_of_local_seeds > s->size_bparts) {
    const size_t nr_bparts_new = s->nr_bparts + number_of_local_seeds;

    s->size_bparts = engine_parts_size_grow * nr_bparts_new;

    struct bpart *bparts_new = NULL;
    if (swift_memalign("bparts", (void **)&bparts_new, bpart_align,
                       sizeof(struct bpart) * s->size_bparts) != 0)
      error("Failed to allocate new bpart data.");
    memcpy(bparts_new, s->bparts, sizeof(struct bpart) * s->nr_bparts);
    swift_free("bparts", s->bparts);

    s->bparts = bparts_new;
  }

  const size_t nr_gparts = s->nr_gparts;
  struct gpart *gparts = s->gparts;
  struct part *parts = s->parts;
  struct xpart *xparts = s->xparts;
  struct bpart *bparts = s->bparts;
  const size_t group_id_default = props->group_id_default;
  const double seed_halo_mass = props->seed_halo_mass;

  /* Direct pointers to the arrays */
  double *group_mass = props->group_mass;
  char *has_black_hole = props->has_black_hole;
  float *max_part_density = props->max_part_density;

  size_t k = s->nr_bparts;

  for (size_t i = 0; i < nr_gparts; i++) {

    /* Check whether this is a gas particle */
    if (gparts[i].type != swift_type_gas) continue;

    /* Ignore inhibited particles */
    if (gparts[i].time_bin >= time_bin_inhibited) continue;

    /* Ignore particles not in groups */
    if (gparts[i].fof_data.group_id == group_id_default) continue;

    if (gparts[i].fof_data.group_id > (size_t)props->num_groups)
      error("Found an invalid group ID!");

    /* Get the density of this particle */
    const size_t gas_index = -gparts[i].id_or_neg_offset;
    const float rho_com = hydro_get_comoving_density(&parts[gas_index]);

    /* Entry into the global list of group properties of this particle */
    const size_t index = gparts[i].fof_data.group_id - 1;

    /* Should we seed a BH in this group? */
    if (!has_black_hole[index] && group_mass[index] > seed_halo_mass) {

      /* Does it match the max density for this group?
       * (i.e. is it the particle we identified as the one to convert?) */
      if (rho_com == max_part_density[index]) {

        /* Handle on the particle to convert */
        struct part *p = &parts[gas_index];
        struct xpart *xp = &xparts[gas_index];
        struct gpart *gp = p->gpart;

#ifdef SWIFT_DEBUG_CHECKS
        if (gp != &gparts[i]) error("Weird gas<->gpart linking error!");
#endif

        /* Let's destroy the gas particle */
        p->time_bin = time_bin_inhibited;
        p->gpart = NULL;

        /* Mark the gpart as black hole */
        gp->type = swift_type_black_hole;

        /* Basic properties of the black hole */
        struct bpart *bp = &bparts[k];
        bzero(bp, sizeof(struct bpart));
        bp->time_bin = gp->time_bin;

        /* Re-link things */
        bp->gpart = gp;
        gp->id_or_neg_offset = -(bp - s->bparts);

        /* Synchronize masses, positions and velocities */
        bp->mass = gp->mass;
        bp->x[0] = gp->x[0];
        bp->x[1] = gp->x[1];
        bp->x[2] = gp->x[2];
        bp->v[0] = gp->v_full[0];
        bp->v[1] = gp->v_full[1];
        bp->v[2] = gp->v_full[2];

        /* Set a smoothing length */
        bp->h = p->h;

        /* Save the ID */
        bp->id = p->id;

        /* Save the tree depth */
        bp->depth_h = p->depth_h;

#ifdef SWIFT_DEBUG_CHECKS
        bp->ti_kick = p->ti_kick;
        bp->ti_drift = p->ti_drift;
#endif

        /* Copy over all the gas properties that we want */
        black_holes_create_from_gas(bp, bh_props, constants, cosmo, p, xp,
                                    s->e->ti_current);
        tracers_first_init_bpart(bp, s->e->internal_units,
                                 s->e->physical_constants, cosmo);

        /* Move to the next BH slot */
        k++;
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (s->nr_bparts + number_of_local_seeds != k) {
    error("Seeded the wrong number of black holes!");
  }
#endif

  /* Update the count of black holes. */
  s->nr_bparts = k;
}

/* Dump FOF group data. */
void fof_dump_group_data(const struct fof_props *props, const int my_rank,
                         const int nr_nodes, const char *out_file_name,
                         struct space *s, const int num_groups) {

  FILE *file = NULL;

  long long *final_group_size = props->final_group_size;
  long long *group_index = props->final_group_index;
  double *group_mass = props->group_mass;
  double *group_centre_of_mass = props->group_centre_of_mass;

  for (int rank = 0; rank < nr_nodes; ++rank) {

    if (rank == 0) {

      const char *mode;
      if (my_rank == 0)
        mode = "w";
      else
        mode = "a";

      file = fopen(out_file_name, mode);

      if (file == NULL)
        error("Could not open the file '%s' with mode '%s'.", out_file_name,
              mode);

      if (my_rank == 0) {
        fprintf(file, "# %8s %12s %12s %12s %12s %12s %12s %24s %24s \n",
                "Group ID", "Group Size", "Group Mass", "CoM_x", "CoM_y",
                "CoM_z", "Max Density", "Max Density Local Index",
                "Particle ID");
        fprintf(file,
                "#-------------------------------------------------------------"
                "------"
                "------------------------------\n");
      }

      for (int i = 0; i < num_groups; i++) {

        fprintf(file, "  %8lld %12lld %12e %12e %12e %12e %12e %24lld %24lld\n",
                group_index[i], final_group_size[i], group_mass[i],
                group_centre_of_mass[i * 3 + 0],
                group_centre_of_mass[i * 3 + 1],
                group_centre_of_mass[i * 3 + 2], 0., -1ll, -1ll);
      }
      fclose(file);
    }
  }
}

struct mapper_data {
  size_t *group_index;
  size_t *group_size;
  float *distance_to_link;
  size_t nr_gparts;
  struct gpart *space_gparts;
};

/**
 * @brief Mapper function to set the roots of #gpart%s going to other MPI ranks.
 *
 * @param map_data The list of outgoing local cells.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to mapper data.
 */
void fof_set_outgoing_root_mapper(void *map_data, int num_elements,
                                  void *extra_data) {
#ifdef WITH_MPI

  /* Unpack the data */
  struct cell **local_cells = (struct cell **)map_data;
  const struct mapper_data *data = (struct mapper_data *)extra_data;
  const size_t *const group_index = data->group_index;
  const size_t *const group_size = data->group_size;
  const size_t nr_gparts = data->nr_gparts;
  const struct gpart *const space_gparts = data->space_gparts;

  /* Loop over the out-going local cells */
  for (int i = 0; i < num_elements; ++i) {

    /* Get the cell and its gparts */
    struct cell *local_cell = local_cells[i];
    struct gpart *gparts = local_cell->grav.parts;

    /* Make a list of particle offsets into the global gparts array. */
    const size_t *const offset =
        group_index + (ptrdiff_t)(gparts - space_gparts);

    /* Set each particle's root and group properties found in the local FOF.*/
    for (int k = 0; k < local_cell->grav.count; k++) {

      /* TODO: Can we skip ignorable particles here?
       * Likely makes no difference */

      /* Recall we did alter the group_index with a global_offset.
       * We need to remove that here as we want the *local* root */
      const size_t root =
          fof_find_global(offset[k] - node_offset, group_index, nr_gparts);

      /* TODO: Could we call fof_find() here instead?
       * Likely yes but we  don't want path compression at this stage.
       * So, probably not */

      gparts[k].fof_data.group_id = root;
      gparts[k].fof_data.group_size = group_size[root - node_offset];
    }
  }

#else
  error("Calling MPI function in non-MPI mode");
#endif
}

/**
 * @brief Search foreign cells for links and communicate any found to the
 * appropriate node.
 *
 * @param props the properties of the FOF scheme.
 * @param s Pointer to a #space.
 */
void fof_search_foreign_cells(struct fof_props *props, const struct space *s) {

#ifdef WITH_MPI
  struct engine *e = s->e;
  const int verbose = e->verbose;

  /* Abort if only one node */
  if (e->nr_nodes == 1) return;

  size_t *restrict group_index = props->group_index;
  size_t *restrict group_size = props->group_size;
  const size_t nr_gparts = s->nr_gparts;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = props->l_x2;

  const ticks tic_total = getticks();
  ticks tic = getticks();

  /* Make group IDs globally unique. */
  for (size_t i = 0; i < nr_gparts; i++) group_index[i] += node_offset;

  struct cell_pair_indices *cell_pairs = NULL;
  int cell_pair_count = 0;

  props->group_links_size = fof_props_default_group_link_size;

  int num_cells_out = 0;
  int num_cells_in = 0;

  /* Find the maximum no. of cell pairs that can communicate. */
  for (int i = 0; i < e->nr_proxies; i++) {

    for (int j = 0; j < e->proxies[i].nr_cells_out; j++) {

      /* Only include gravity cells. */
      if (e->proxies[i].cells_out_type[j] & proxy_cell_type_gravity)
        num_cells_out++;
    }

    for (int j = 0; j < e->proxies[i].nr_cells_in; j++) {

      /* Only include gravity cells. */
      if (e->proxies[i].cells_in_type[j] & proxy_cell_type_gravity)
        num_cells_in++;
    }
  }

  if (verbose)
    message(
        "Finding max no. of cells + offset IDs"
        "took: %.3f %s.",
        clocks_from_ticks(getticks() - tic), clocks_getunit());

  const int cell_pair_size = num_cells_in * num_cells_out;

  /* Allocate memory for all the possible cell links */
  if (swift_memalign("fof_group_links", (void **)&props->group_links,
                     SWIFT_STRUCT_ALIGNMENT,
                     props->group_links_size * sizeof(struct fof_mpi)) != 0)
    error("Error while allocating memory for FOF links over an MPI domain");

  if (swift_memalign("fof_cell_pairs", (void **)&cell_pairs,
                     SWIFT_STRUCT_ALIGNMENT,
                     cell_pair_size * sizeof(struct cell_pair_indices)) != 0)
    error("Error while allocating memory for FOF cell pair indices");

  ticks tic_pairs = getticks();

  /* Loop over cells_in and cells_out for each proxy and find which cells are in
   * range of each other to perform the FOF search. Store local cells that are
   * touching foreign cells in a list. */
  for (int i = 0; i < e->nr_proxies; i++) {

    /* Only find links across an MPI domain on one rank. */
    if (engine_rank == min(engine_rank, e->proxies[i].nodeID)) {

      for (int j = 0; j < e->proxies[i].nr_cells_out; j++) {

        /* Skip non-gravity cells. */
        if (!(e->proxies[i].cells_out_type[j] & proxy_cell_type_gravity))
          continue;

        struct cell *restrict local_cell = e->proxies[i].cells_out[j];

        /* Skip empty cells. */
        if (local_cell->grav.count == 0) continue;

        for (int k = 0; k < e->proxies[i].nr_cells_in; k++) {

          /* Skip non-gravity cells. */
          if (!(e->proxies[i].cells_in_type[k] & proxy_cell_type_gravity))
            continue;

          struct cell *restrict foreign_cell = e->proxies[i].cells_in[k];

          /* Skip empty cells. */
          if (foreign_cell->grav.count == 0) continue;

          /* Add candidates in range to the list of pairs of cells to treat. */
          const double r2 = cell_min_dist(local_cell, foreign_cell, dim);
          if (r2 < search_r2) {
            cell_pairs[cell_pair_count].local = local_cell;
            cell_pairs[cell_pair_count].foreign = foreign_cell;

            cell_pair_count++;
          }
        }
      }
    }
  }

  if (verbose)
    message("Finding local/foreign cell pairs took: %.3f %s.",
            clocks_from_ticks(getticks() - tic_pairs), clocks_getunit());

  const ticks tic_set_roots = getticks();

  /* Set the root of outgoing particles. */

  /* Allocate array of outgoing cells and populate it */
  struct cell **local_cells =
      (struct cell **)malloc(num_cells_out * sizeof(struct cell *));
  int count = 0;
  for (int i = 0; i < e->nr_proxies; i++) {
    for (int j = 0; j < e->proxies[i].nr_cells_out; j++) {

      /* Only include gravity cells. */
      if (e->proxies[i].cells_out_type[j] & proxy_cell_type_gravity) {

        local_cells[count] = e->proxies[i].cells_out[j];
        ++count;
      }
    }
  }

  /* Now set the *local* roots of all the gparts we are sending */
  struct mapper_data data;
  data.group_index = group_index;
  data.group_size = group_size;
  data.nr_gparts = nr_gparts;
  data.space_gparts = s->gparts;
  threadpool_map(&e->threadpool, fof_set_outgoing_root_mapper, local_cells,
                 num_cells_out, sizeof(struct cell **),
                 threadpool_auto_chunk_size, &data);

  if (verbose)
    message("Initialising particle roots took: %.3f %s.",
            clocks_from_ticks(getticks() - tic_set_roots), clocks_getunit());

  free(local_cells);

  if (verbose)
    message(
        "Finding local/foreign cell pairs and initialising particle roots "
        "took: %.3f %s.",
        clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Activate the tasks exchanging all the required gparts */
  engine_activate_gpart_comms(e);

  ticks local_fof_tic = getticks();

  /* Wait for all the communication tasks to be ready */
  MPI_Barrier(MPI_COMM_WORLD);

  if (verbose)
    message("Local FOF imbalance: %.3f %s.",
            clocks_from_ticks(getticks() - local_fof_tic), clocks_getunit());

  tic = getticks();

  /* Perform send and receive tasks. */
  engine_launch(e, "fof comms");

  if (verbose)
    message("MPI send/recv comms took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* We have now recevied the foreign particles. Each particle received
   * carries information about its own *foreign* (to us) root and the
   * size of the group fragment it belongs too its original foreign rank. */

  tic = getticks();

  props->group_link_count = 0;

  /* Perform search of group links between local and foreign cells with the
   * threadpool. */
  threadpool_map(&s->e->threadpool, fof_find_foreign_links_mapper, cell_pairs,
                 cell_pair_count, sizeof(struct cell_pair_indices), 1,
                 (struct space *)s);

  /* Clean up memory used by foreign particles. */
  swift_free("fof_cell_pairs", cell_pairs);

  if (verbose)
    message("Searching for foreign links took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  const ticks comms_tic = getticks();

  MPI_Barrier(MPI_COMM_WORLD);

  if (verbose)
    message("Imbalance took: %.3f %s.",
            clocks_from_ticks(getticks() - comms_tic), clocks_getunit());

  if (verbose)
    message("fof_search_foreign_cells() took (FOF SCALING): %.3f %s.",
            clocks_from_ticks(getticks() - tic_total), clocks_getunit());

#endif /* WITH_MPI */
}

/**
 * @brief Run all the tasks attaching the attachables to their
 * nearest linkable particle.
 *
 * @param props The properties fof the FOF scheme.
 * @param s The #space we work with.
 */
void fof_link_attachable_particles(struct fof_props *props,
                                   const struct space *s) {

  /* Is there anything to attach? */
  if (!current_fof_attach_type) return;

  const ticks tic_total = getticks();

  /* Activate the tasks attaching attachable particles to the linkable ones */
  engine_activate_fof_attach_tasks(s->e);

  /* Perform FOF tasks for attachable particles. */
  engine_launch(s->e, "fof");

  if (s->e->verbose)
    message("fof_link_attachable_particles() took (FOF SCALING): %.3f %s.",
            clocks_from_ticks(getticks() - tic_total), clocks_getunit());
}

/**
 * @brief Process all the group fragments spanning more than
 * one rank to link them.
 *
 * This is the final global union-find pass which concludes
 * the MPI-FOF-algorithm.
 *
 * @param props The properties fof the FOF scheme.
 * @param s The #space we work with.
 */
void fof_link_foreign_fragments(struct fof_props *props,
                                const struct space *s) {

#ifdef WITH_MPI

  struct engine *e = s->e;
  const int verbose = e->verbose;

  /* Abort if only one node */
  if (e->nr_nodes == 1) return;

  const size_t nr_gparts = s->nr_gparts;
  size_t *restrict group_index = props->group_index;
  size_t *restrict group_size = props->group_size;

  const ticks tic_total = getticks();
  ticks tic = getticks();
  const ticks comms_tic = getticks();

  if (verbose)
    message(
        "Searching %zu gravity particles for cross-node links with l_x: %lf",
        nr_gparts, sqrt(props->l_x2));

  /* Local copy of the variable set in the mapper */
  const int group_link_count = props->group_link_count;

  /* Sum the total number of links across MPI domains over each MPI rank. */
  int global_group_link_count = 0;
  MPI_Allreduce(&group_link_count, &global_group_link_count, 1, MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);

  if (global_group_link_count < 0)
    error("Overflow of the size of the global list of foreign links");

  struct fof_mpi *global_group_links = NULL;
  int *displ = NULL, *group_link_counts = NULL;

  if (swift_memalign("fof_global_group_links", (void **)&global_group_links,
                     SWIFT_STRUCT_ALIGNMENT,
                     global_group_link_count * sizeof(struct fof_mpi)) != 0)
    error("Error while allocating memory for the global list of group links");

  if (posix_memalign((void **)&group_link_counts, SWIFT_STRUCT_ALIGNMENT,
                     e->nr_nodes * sizeof(int)) != 0)
    error(
        "Error while allocating memory for the number of group links on each "
        "MPI rank");

  if (posix_memalign((void **)&displ, SWIFT_STRUCT_ALIGNMENT,
                     e->nr_nodes * sizeof(int)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  /* Gather the total number of links on each rank. */
  MPI_Allgather(&group_link_count, 1, MPI_INT, group_link_counts, 1, MPI_INT,
                MPI_COMM_WORLD);

  /* Set the displacements into the global link list using the link counts from
   * each rank */
  displ[0] = 0;
  for (int i = 1; i < e->nr_nodes; i++) {
    displ[i] = displ[i - 1] + group_link_counts[i - 1];
    if (displ[i] < 0) error("Number of group links overflowing!");
  }

  /* Gather the global link list on all ranks. */
  MPI_Allgatherv(props->group_links, group_link_count, fof_mpi_type,
                 global_group_links, group_link_counts, displ, fof_mpi_type,
                 MPI_COMM_WORLD);

  /* Clean up memory. */
  free(group_link_counts);
  free(displ);
  swift_free("fof_group_links", props->group_links);
  props->group_links = NULL;

  if (verbose) {
    message("Communication took: %.3f %s.",
            clocks_from_ticks(getticks() - comms_tic), clocks_getunit());

    message("Global comms took: %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
  }

  tic = getticks();

  /* Transform the group IDs to a local list going from 0-group_count so a
   * union-find can be performed.
   * Each member of a link is stored separately --> Need 2x as many entries */
  size_t *global_group_index = NULL, *global_group_id = NULL,
         *global_group_size = NULL;
  const int global_group_list_size = 2 * global_group_link_count;

  if (swift_memalign("fof_global_group_index", (void **)&global_group_index,
                     SWIFT_STRUCT_ALIGNMENT,
                     global_group_list_size * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  if (swift_memalign("fof_global_group_id", (void **)&global_group_id,
                     SWIFT_STRUCT_ALIGNMENT,
                     global_group_list_size * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  if (swift_memalign("fof_global_group_size", (void **)&global_group_size,
                     SWIFT_STRUCT_ALIGNMENT,
                     global_group_list_size * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  bzero(global_group_size, global_group_list_size * sizeof(size_t));

  /* Create hash table. */
  hashmap_t map;
  hashmap_init(&map);

  /* Store each group ID and its properties. */
  int group_count = 0;
  for (int k = 0; k < global_group_link_count; k++) {

    const size_t group_i = global_group_links[k].group_i;
    const size_t group_j = global_group_links[k].group_j;

    global_group_size[group_count] += global_group_links[k].group_i_size;
    global_group_id[group_count] = group_i;
    hashmap_add_group(group_i, group_count, &map);
    group_count++;

    global_group_size[group_count] += global_group_links[k].group_j_size;
    global_group_id[group_count] = group_j;
    hashmap_add_group(group_j, group_count, &map);
    group_count++;
  }

  if (verbose)
    message("Global list compression took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Create a global_group_index list of groups across MPI domains so that you
   * can perform a union-find locally on each node.
   * The value of which is an offset into global_group_id, which is the actual
   * root. */
  for (int i = 0; i < group_count; i++) global_group_index[i] = i;

  /* Store the original group size before incrementing in the Union-Find. */
  size_t *orig_global_group_size = NULL;

  if (swift_memalign("fof_orig_global_group_size",
                     (void **)&orig_global_group_size, SWIFT_STRUCT_ALIGNMENT,
                     group_count * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  memcpy(orig_global_group_size, global_group_size,
         group_count * sizeof(size_t));

  /* Perform a union-find on the group links. */
  for (int k = 0; k < global_group_link_count; k++) {

    /* Use the hash table to find the group offsets in the index array. */
    const size_t find_i =
        hashmap_find_group_offset(global_group_links[k].group_i, &map);
    const size_t find_j =
        hashmap_find_group_offset(global_group_links[k].group_j, &map);

    /* Use the offset to find the group's root. */
    const size_t root_i = fof_find(find_i, global_group_index);
    const size_t root_j = fof_find(find_j, global_group_index);

    const size_t group_i = global_group_id[root_i];
    const size_t group_j = global_group_id[root_j];

    if (group_i == group_j) continue;

    /* Update roots accordingly. */
    const size_t size_i = global_group_size[root_i];
    const size_t size_j = global_group_size[root_j];
#ifdef UNION_BY_SIZE_OVER_MPI
    if (size_i < size_j) {
      global_group_index[root_i] = root_j;
      global_group_size[root_j] += size_i;
    } else {
      global_group_index[root_j] = root_i;
      global_group_size[root_i] += size_j;
    }
#else
    if (group_j < group_i) {
      global_group_index[root_i] = root_j;
      global_group_size[root_j] += size_i;
    } else {
      global_group_index[root_j] = root_i;
      global_group_size[root_i] += size_j;
    }
#endif
  }

  hashmap_free(&map);

  if (verbose)
    message("global_group_index construction took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Update each group locally with new root information. */
  for (int i = 0; i < group_count; i++) {

    const size_t group_id = global_group_id[i];
    const size_t offset = fof_find(global_group_index[i], global_group_index);
    const size_t new_root = global_group_id[offset];

    /* If the group is local update its root and size. */
    if (is_local(group_id, nr_gparts) && new_root != group_id) {
      group_index[group_id - node_offset] = new_root;

      group_size[group_id - node_offset] -= orig_global_group_size[i];
    }

    /* If the group linked to a local root update its size. */
    if (is_local(new_root, nr_gparts) && new_root != group_id) {

      /* Use group sizes before Union-Find */
      group_size[new_root - node_offset] += orig_global_group_size[i];
    }
  }

  if (verbose)
    message("Updating groups locally took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Clean up memory. */
  swift_free("fof_global_group_links", global_group_links);
  swift_free("fof_global_group_index", global_group_index);
  swift_free("fof_global_group_size", global_group_size);
  swift_free("fof_global_group_id", global_group_id);
  swift_free("fof_orig_global_group_size", orig_global_group_size);

  if (verbose) {
    message("link_foreign_fragmens() took (FOF SCALING): %.3f %s.",
            clocks_from_ticks(getticks() - tic_total), clocks_getunit());
  }

#endif /* WITH_MPI */
}

/**
 * @brief Compute the local size of each FOF group fragment.
 *
 * @param props The properties of the FOF scheme.
 * @param s The #space containing the particles.
 */
void fof_compute_local_sizes(struct fof_props *props, struct space *s) {

  const int verbose = s->e->verbose;

  struct gpart *gparts = s->gparts;
  const size_t nr_gparts = s->nr_gparts;

  const ticks tic_total = getticks();

  if (engine_rank == 0 && verbose)
    message("Size of hash table element: %ld", sizeof(hashmap_element_t));

#ifdef WITH_MPI

  const ticks comms_tic = getticks();

  /* Determine number of gparts on lower numbered MPI ranks */
  const long long nr_gparts_local = s->nr_gparts;
  long long nr_gparts_cumulative;
  MPI_Scan(&nr_gparts_local, &nr_gparts_cumulative, 1, MPI_LONG_LONG, MPI_SUM,
           MPI_COMM_WORLD);

  if (verbose)
    message("MPI_Scan Imbalance took: %.3f %s.",
            clocks_from_ticks(getticks() - comms_tic), clocks_getunit());

  /* Reset global variable containing the rank particle count offset */
  node_offset = nr_gparts_cumulative - nr_gparts_local;
#endif

  /* Compute the group sizes of the local fragments
   * (in non-MPI land that is the final group size of the haloes) */
  const ticks tic_calc_group_size = getticks();

  threadpool_map(&s->e->threadpool, fof_calc_group_size_mapper, gparts,
                 nr_gparts, sizeof(struct gpart), threadpool_auto_chunk_size,
                 s);
  if (verbose)
    message("FOF calc group size took (FOF SCALING): %.3f %s.",
            clocks_from_ticks(getticks() - tic_calc_group_size),
            clocks_getunit());

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic_total),
            clocks_getunit());
}

/**
 * @brief Compute all the group properties
 *
 * @param props The properties of the FOF scheme.
 * @param bh_props The properties of the black hole scheme.
 * @param constants The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param s The #space containing the particles.
 * @param dump_debug_results Are we writing txt-file debug catalogues including
 * BH-seeding info?
 * @param dump_results Do we want to write the group catalogue to a hdf5 file?
 * @param seed_black_holes Do we want to seed black holes in haloes?
 */
void fof_assign_group_ids(struct fof_props *props, struct space *s) {

  const int verbose = s->e->verbose;
#ifdef WITH_MPI
  const int nr_nodes = s->e->nr_nodes;
#endif
  const ticks tic_total = getticks();

  struct gpart *gparts = s->gparts;
  const size_t nr_gparts = s->nr_gparts;

  const size_t min_group_size = props->min_group_size;
  const size_t group_id_offset = props->group_id_offset;
  const size_t group_id_default = props->group_id_default;

  size_t num_groups_local = 0;
  size_t num_parts_in_groups_local = 0;
  size_t max_group_size_local = 0;

  /* Local copy of the arrays */
  size_t *restrict group_index = props->group_index;
  size_t *restrict group_size = props->group_size;

  const ticks tic_num_groups_calc = getticks();

  for (size_t i = 0; i < nr_gparts; i++) {

#ifdef WITH_MPI
    /* Find the total number of groups. */
    if (group_index[i] == i + node_offset && group_size[i] >= min_group_size)
      num_groups_local++;
#else
    /* Find the total number of groups. */
    if (group_index[i] == i && group_size[i] >= min_group_size)
      num_groups_local++;
#endif

    /* Find the total number of particles in groups. */
    if (group_size[i] >= min_group_size)
      num_parts_in_groups_local += group_size[i];

    /* Find the largest group. */
    if (group_size[i] > max_group_size_local)
      max_group_size_local = group_size[i];
  }

  if (verbose)
    message(
        "Calculating the total no. of local groups took: (FOF SCALING): %.3f "
        "%s.",
        clocks_from_ticks(getticks() - tic_num_groups_calc), clocks_getunit());

  /* Sort the groups in descending order based upon size and re-label their
   * IDs 0-num_groups. */
  struct group_length *high_group_sizes = NULL;
  size_t group_count = 0;

  if (swift_memalign("fof_high_group_sizes", (void **)&high_group_sizes, 32,
                     num_groups_local * sizeof(struct group_length)) != 0)
    error("Failed to allocate list of large groups.");

  /* Store the group_sizes and their offset. */
  for (size_t i = 0; i < nr_gparts; i++) {

#ifdef WITH_MPI
    if (group_index[i] == i + node_offset && group_size[i] >= min_group_size) {
      high_group_sizes[group_count].index = node_offset + i;
      high_group_sizes[group_count++].size = group_size[i];
    }
#else
    if (group_index[i] == i && group_size[i] >= min_group_size) {
      high_group_sizes[group_count].index = i;
      high_group_sizes[group_count++].size = group_size[i];
    }
#endif
  }
  props->high_group_sizes = high_group_sizes;

  ticks tic = getticks();

  /* Find global properties. */
  long long num_groups = 0, num_parts_in_groups = 0, max_group_size = 0;
#ifdef WITH_MPI
  MPI_Allreduce(&num_groups_local, &num_groups, 1, MPI_LONG_LONG_INT, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Reduce(&num_parts_in_groups_local, &num_parts_in_groups, 1,
             MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max_group_size_local, &max_group_size, 1, MPI_LONG_LONG_INT,
             MPI_MAX, 0, MPI_COMM_WORLD);
#else
  num_groups = num_groups_local;
  num_parts_in_groups = num_parts_in_groups_local;
  max_group_size = max_group_size_local;
#endif /* WITH_MPI */

  props->num_groups = num_groups;
  // message("num_groups_local=%zd", num_groups_local);

  /* Find number of groups on lower numbered MPI ranks */
#ifdef WITH_MPI
  long long nglocal = num_groups_local;
  long long ngsum;
  MPI_Scan(&nglocal, &ngsum, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  props->num_groups_prev = (size_t)(ngsum - nglocal);
#endif /* WITH_MPI */

  if (verbose)
    message("Finding the total no. of groups took: (FOF SCALING): %.3f %s.",
            clocks_from_ticks(getticks() - tic_num_groups_calc),
            clocks_getunit());

  /* Sort local groups into descending order of size */
  qsort(high_group_sizes, num_groups_local, sizeof(struct group_length),
        cmp_func_group_size);

  tic = getticks();

  /* Set default group ID for all particles */
  threadpool_map(&s->e->threadpool, fof_set_initial_group_id_mapper, s->gparts,
                 s->nr_gparts, sizeof(struct gpart), threadpool_auto_chunk_size,
                 (void *)&group_id_default);

  if (verbose)
    message("Setting default group ID took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* Assign final group IDs to local root particles where the global root is
   * on this node and the group is large enough. Within a node IDs are
   * assigned in descending order of particle number. */
  for (size_t i = 0; i < num_groups_local; i++) {
#ifdef WITH_MPI
    gparts[high_group_sizes[i].index - node_offset].fof_data.group_id =
        group_id_offset + i + props->num_groups_prev;
#else
    gparts[high_group_sizes[i].index].fof_data.group_id = group_id_offset + i;
#endif
  }

#ifdef WITH_MPI

  /* Now, for each local root where the global root is on some other node
   * AND the total size of the group is >= min_group_size we need to
   * retrieve the gparts.group_id we just assigned to the global root.
   *
   * Will do that by sending the group_index of these local roots to the
   * node where their global root is stored and receiving back the new
   * group_id associated with that particle.
   *
   * Identify local roots with global root on another node and large enough
   * group_size. Store index of the local and global roots in these cases.
   *
   * NOTE: if group_size only contains the total FoF mass for global roots,
   * then we have to communicate ALL fragments where the global root is not
   * on this node. Hence the commented out extra conditions below.*/
  size_t nsend = 0;
  for (size_t i = 0; i < nr_gparts; i += 1) {
    if ((!is_local(group_index[i],
                   nr_gparts))) { /* && (group_size[i] >= min_group_size)) { */
      nsend += 1;
    }
  }
  struct fof_final_index *fof_index_send =
      (struct fof_final_index *)swift_malloc(
          "fof_index_send", sizeof(struct fof_final_index) * nsend);
  nsend = 0;
  for (size_t i = 0; i < nr_gparts; i += 1) {
    if ((!is_local(group_index[i],
                   nr_gparts))) { /* && (group_size[i] >= min_group_size)) { */
      fof_index_send[nsend].local_root = node_offset + i;
      fof_index_send[nsend].global_root = group_index[i];
      nsend += 1;
    }
  }

  /* Sort by global root - this puts the groups in order of which node they're
   * stored on */
  qsort(fof_index_send, nsend, sizeof(struct fof_final_index),
        compare_fof_final_index_global_root);

  /* Determine range of global indexes (i.e. particles) on each node */
  props->num_on_node =
      (size_t *)swift_malloc("fof_num_on_node", nr_nodes * sizeof(size_t));
  MPI_Allgather(&nr_gparts, sizeof(size_t), MPI_BYTE, props->num_on_node,
                sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);
  props->first_on_node =
      (size_t *)swift_malloc("fof_first_on_node", nr_nodes * sizeof(size_t));
  props->first_on_node[0] = 0;
  for (int i = 1; i < nr_nodes; i++)
    props->first_on_node[i] =
        props->first_on_node[i - 1] + props->num_on_node[i - 1];

  /* Determine how many entries go to each node */
  int *sendcount = (int *)calloc(nr_nodes, sizeof(int));
  int dest = 0;
  for (size_t i = 0; i < nsend; i += 1) {
    while ((fof_index_send[i].global_root >=
            props->first_on_node[dest] + props->num_on_node[dest]) ||
           (props->num_on_node[dest] == 0))
      dest += 1;
    if (dest >= nr_nodes) error("Node index out of range!");
    sendcount[dest] += 1;
  }

  int *recvcount = NULL, *sendoffset = NULL, *recvoffset = NULL;
  size_t nrecv = 0;

  fof_compute_send_recv_offsets(nr_nodes, sendcount, &recvcount, &sendoffset,
                                &recvoffset, &nrecv);

  struct fof_final_index *fof_index_recv =
      (struct fof_final_index *)swift_malloc(
          "fof_index_recv", nrecv * sizeof(struct fof_final_index));

  /* Exchange group indexes */
  MPI_Alltoallv(fof_index_send, sendcount, sendoffset, fof_final_index_type,
                fof_index_recv, recvcount, recvoffset, fof_final_index_type,
                MPI_COMM_WORLD);

  /* For each received global root, look up the group ID we assigned and store
   * it in the struct */
  for (size_t i = 0; i < nrecv; i += 1) {
    if ((fof_index_recv[i].global_root < node_offset) ||
        (fof_index_recv[i].global_root >= node_offset + nr_gparts)) {
      error("Received global root index out of range!");
    }
    fof_index_recv[i].global_root =
        gparts[fof_index_recv[i].global_root - node_offset].fof_data.group_id;
  }

  /* Send the result back */
  MPI_Alltoallv(fof_index_recv, recvcount, recvoffset, fof_final_index_type,
                fof_index_send, sendcount, sendoffset, fof_final_index_type,
                MPI_COMM_WORLD);

  /* Update local gparts.group_id */
  for (size_t i = 0; i < nsend; i += 1) {
    if ((fof_index_send[i].local_root < node_offset) ||
        (fof_index_send[i].local_root >= node_offset + nr_gparts)) {
      error("Sent local root index out of range!");
    }
    gparts[fof_index_send[i].local_root - node_offset].fof_data.group_id =
        fof_index_send[i].global_root;
  }

  free(sendcount);
  free(recvcount);
  free(sendoffset);
  free(recvoffset);
  swift_free("fof_index_send", fof_index_send);
  swift_free("fof_index_recv", fof_index_recv);

#endif /* WITH_MPI */

  size_t max_id = 0;

  /* Assign every particle the group_id of its local root. */
  for (size_t i = 0; i < nr_gparts; i++) {

    if (gpart_is_ignorable(&gparts[i])) continue;
    if (gpart_is_attachable(&gparts[i])) continue;

    const size_t root = fof_find_local(i, nr_gparts, group_index);
    gparts[i].fof_data.group_id = gparts[root].fof_data.group_id;

    if (gparts[i].fof_data.group_id != fof_props_default_group_id)
      max_id = max(max_id, gparts[i].fof_data.group_id);
  }

  /* Give some info */
  if (engine_rank == 0) {
    message(
        "No. of groups: %lld. No. of particles in groups: %lld. No. of "
        "particles not in groups: %lld.",
        num_groups, num_parts_in_groups,
        s->e->total_nr_gparts - num_parts_in_groups);
    message("Largest group (linkables only) by size: %lld", max_group_size);
  }

  /* Free data we are done with */
  swift_free("fof_group_index", props->group_index);
  props->group_index = NULL;

  if (verbose)
    message("took: %.3f %s.", clocks_from_ticks(getticks() - tic_total),
            clocks_getunit());
}

/**
 * @brief Compute all the group properties
 *
 * @param props The properties of the FOF scheme.
 * @param bh_props The properties of the black hole scheme.
 * @param constants The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param s The #space containing the particles.
 * @param dump_debug_results Are we writing txt-file debug catalogues including
 * BH-seeding info?
 * @param dump_results Do we want to write the group catalogue to a hdf5 file?
 * @param seed_black_holes Do we want to seed black holes in haloes?
 */
void fof_compute_group_props(struct fof_props *props,
                             const struct black_holes_props *bh_props,
                             const struct phys_const *constants,
                             const struct cosmology *cosmo, struct space *s,
                             const int dump_results,
                             const int dump_debug_results,
                             const int seed_black_holes) {

  const int verbose = s->e->verbose;
  const ticks tic_total = getticks();

  const size_t num_groups = props->num_groups;

  /* Allocate and initialise a group mass and centre of mass array
     for *all* groups. */
  if (swift_memalign("fof_group_mass", (void **)&props->group_mass, 32,
                     num_groups * sizeof(double)) != 0)
    error("Failed to allocate list of group masses for FOF search.");
  if (swift_memalign("fof_group_size", (void **)&props->final_group_size, 32,
                     num_groups * sizeof(long long)) != 0)
    error("Failed to allocate list of group masses for FOF search.");
  if (swift_memalign("fof_group_index", (void **)&props->final_group_index, 32,
                     num_groups * sizeof(long long)) != 0)
    error("Failed to allocate list of group masses for FOF search.");
  if (swift_memalign("fof_has_black_hole", (void **)&props->has_black_hole, 32,
                     num_groups * sizeof(char)) != 0)
    error("Failed to allocate list of black holes for FOF search.");
  if (swift_memalign("fof_group_centre_of_mass",
                     (void **)&props->group_centre_of_mass, 32,
                     num_groups * 3 * sizeof(double)) != 0)
    error("Failed to allocate list of group CoM for FOF search.");
  if (swift_memalign("fof_max_part_density", (void **)&props->max_part_density,
                     32, num_groups * sizeof(float)) != 0)
    error("Failed to allocate list of max group densities for FOF search.");

  bzero(props->group_mass, num_groups * sizeof(double));
  bzero(props->final_group_size, num_groups * sizeof(long long));
  bzero(props->final_group_index, num_groups * sizeof(long long));
  bzero(props->has_black_hole, num_groups * sizeof(char));
  bzero(props->group_centre_of_mass, num_groups * 3 * sizeof(double));
  bzero(props->max_part_density, num_groups * sizeof(float));

  const ticks tic_props = getticks();

  /* Now, compute all things */
  int number_of_local_seeds = 0, number_of_global_seeds = 0;
  fof_calc_group_mass(props, s, &number_of_local_seeds,
                      &number_of_global_seeds);

  if (verbose)
    message("Computing group properties took: %.3f %s.",
            clocks_from_ticks(getticks() - tic_props), clocks_getunit());

  /* All MPI ranks now have information about all haloes, even
   * the ones where there are no local particles in. */

  /* Dump group data (only rank 0 since everyone has everything anyway). */
  if (dump_results && engine_rank == 0) {
#ifdef HAVE_HDF5
    write_fof_hdf5_catalogue(props, s->e);
#else
    error("Can't dump hdf5 catalogues with hdf5 switched off!");
#endif
  }

  if (dump_debug_results && engine_rank == 0) {

    char output_file_name[PARSER_MAX_LINE_SIZE];
    snprintf(output_file_name, PARSER_MAX_LINE_SIZE, "%s", props->base_name);

#ifdef WITH_MPI
    snprintf(output_file_name + strlen(output_file_name), FILENAME_BUFFER_SIZE,
             "_mpi.dat");
#else
    snprintf(output_file_name + strlen(output_file_name), FILENAME_BUFFER_SIZE,
             ".dat");
#endif
    fof_dump_group_data(props, s->e->nodeID, s->e->nr_nodes, output_file_name,
                        s, num_groups);
  }

  /* Seed black holes */
  if (seed_black_holes) {
    fof_seed_black_holes(props, bh_props, constants, cosmo, s,
                         number_of_local_seeds, number_of_global_seeds);
  }

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic_total),
            clocks_getunit());
}

/**
 * @brief Free all the arrays we allocated on the way
 */
void fof_free_arrays(struct fof_props *props) {

  swift_free("fof_high_group_sizes", props->high_group_sizes);
  swift_free("fof_group_mass", props->group_mass);
  swift_free("fof_group_size", props->final_group_size);
  swift_free("fof_group_index", props->final_group_index);
  swift_free("fof_group_centre_of_mass", props->group_centre_of_mass);
  swift_free("fof_max_part_density", props->max_part_density);
  swift_free("fof_has_black_hole", props->has_black_hole);
  swift_free("fof_distance", props->distance_to_link);
  swift_free("fof_attach_index", props->attach_index);
  swift_free("fof_found_attach", props->found_attachable_link);
  swift_free("fof_group_size", props->group_size);
  props->group_mass = NULL;
  props->final_group_size = NULL;
  props->final_group_index = NULL;
  props->group_centre_of_mass = NULL;
  props->max_part_density = NULL;
  props->has_black_hole = NULL;
  props->group_size = NULL;

#ifdef WITH_MPI
  swift_free("fof_num_on_node", props->num_on_node);
  swift_free("fof_first_on_node", props->first_on_node);
  props->num_on_node = NULL;
  props->first_on_node = NULL;
#endif
}

void fof_struct_dump(const struct fof_props *props, FILE *stream) {

  struct fof_props temp = *props;
  temp.num_groups = 0;
  temp.group_link_count = 0;
  temp.group_links_size = 0;
  temp.group_index = NULL;
  temp.group_size = NULL;
  temp.group_mass = NULL;
  temp.final_group_size = NULL;
  temp.group_centre_of_mass = NULL;
  temp.max_part_density = NULL;
  temp.group_links = NULL;

  restart_write_blocks((void *)&temp, sizeof(struct fof_props), 1, stream,
                       "fof_props", "fof_props");
}

void fof_struct_restore(struct fof_props *props, FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct fof_props), 1, stream, NULL,
                      "fof_props");

  fof_set_current_types(props);
}

#endif /* WITH_FOF */
