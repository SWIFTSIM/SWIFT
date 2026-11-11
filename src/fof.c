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
#include <float.h>
#include <libgen.h>
#include <limits.h>
#include <stdint.h>
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
#include "custom_mpi.h"
#include "engine.h"
#include "fof_catalogue_io.h"
#include "hashmap.h"
#include "memuse.h"
#include "proxy.h"
#include "star_formation.h"
#include "threadpool.h"
#include "tools.h"
#include "tracers.h"

#define fof_props_default_group_id 2147483647
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
  const int min_group_size = parser_get_param_int(params, "FOF:min_group_size");
  if (min_group_size < 2)
    error("FOF:min_group_size must be at least 2 (got %d).", min_group_size);
  props->min_group_size = min_group_size;

  /* Read whether we're doing FoF calls to seed black holes. */
#ifdef BLACK_HOLES_NONE
  props->seed_black_holes_enabled = 0;
#else
  props->seed_black_holes_enabled =
      parser_get_param_int(params, "FOF:seed_black_holes_enabled");
#endif

  /* Read the default group ID of particles in groups below the minimum group
   * size. */
  const int group_id_default = parser_get_opt_param_int(
      params, "FOF:group_id_default", fof_props_default_group_id);
  if (group_id_default < 0)
    error("FOF:group_id_default must not be negative (got %d).",
          group_id_default);
  props->group_id_default = group_id_default;

  /* The group IDs always start at 1. The option to change this was removed but
   * we still accept parameter files that set it to that value. */
  if (parser_get_opt_param_int(params, "FOF:group_id_offset", 1) != 1)
    error(
        "FOF:group_id_offset has been removed. Group IDs always start at 1. "
        "Please remove the parameter from the parameter file.");

  /* Read the linking length ratio to the mean inter-particle separation. */
  props->l_x_ratio =
      parser_get_opt_param_double(params, "FOF:linking_length_ratio", -1.);

  /* Read value of absolute linking length aksed by the user */
  props->l_x_absolute =
      parser_get_opt_param_double(params, "FOF:absolute_linking_length", -1.);

  /* -1 means the parameter was not given. Any value given must be positive. */
  if (props->l_x_absolute != -1. && props->l_x_absolute <= 0.)
    error("FOF:absolute_linking_length must be positive (got %e).",
          props->l_x_absolute);

  if (props->l_x_ratio != -1. && props->l_x_ratio <= 0.)
    error("FOF:linking_length_ratio must be positive (got %e).",
          props->l_x_ratio);

  if (props->l_x_ratio == -1. && props->l_x_absolute == -1.)
    error(
        "One of FOF:linking_length_ratio or FOF:absolute_linking_length must "
        "be specified.");

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

  /* The seeding only checks for existing black holes among the members of a
   * group, so the black holes must be part of the groups. */
  if (!stand_alone_fof && props->seed_black_holes_enabled &&
      !props->fof_linking_types[swift_type_black_hole] &&
      !props->fof_attach_types[swift_type_black_hole])
    error(
        "Seeding black holes with FOF requires the black holes to be a "
        "linking or attaching type. Otherwise, black holes would be seeded "
        "again in groups that already contain one.");

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
        "Performing FOF over MPI using union by size and union by index "
        "locally.");
#else
  if (engine_rank == 0) message("Performing FOF using union by index.");
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

  /* The linking length can only be derived from the dark matter. */
  if (props->l_x_absolute == -1. && num_high_res_DM == 0)
    error(
        "No dark matter particles to compute the FOF linking length from the "
        "mean inter-particle separation. Please specify "
        "FOF:absolute_linking_length in the parameter file.");

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

  /* Check size of linking length against the top-level cell dimensions. The
   * FOF tasks are only constructed between neighbouring top-level cells (and
   * the split tasks rely on the same assumption). */
  const double min_width = min3(s->width[0], s->width[1], s->width[2]);
  if (props->l_x2 > min_width * min_width)
    error(
        "Linking length greater than the width of a top-level cell. Need to "
        "check more than one layer of top-level cells for links.");

  /* Are there any particles to attach? The attaching tasks and the
   * associated communications are only needed if there are. For the types
   * that only exist as gparts we have no separate count, so assume they are
   * present. Note that the counts include the extra (not yet created)
   * particles, which is conservative. */
  props->with_attachables = 0;
  for (int i = 0; i < swift_type_count; ++i) {
    if (!props->fof_attach_types[i]) continue;
    switch (i) {
      case swift_type_gas:
        if (s->nr_parts > 0) props->with_attachables = 1;
        break;
      case swift_type_sink:
        if (s->nr_sinks > 0) props->with_attachables = 1;
        break;
      case swift_type_stars:
        if (s->nr_sparts > 0) props->with_attachables = 1;
        break;
      case swift_type_black_hole:
        if (s->nr_bparts > 0) props->with_attachables = 1;
        break;
      default:
        props->with_attachables = 1;
        break;
    }
  }

#ifdef WITH_MPI
  /* All ranks must agree as the attaching requires communications */
  MPI_Allreduce(MPI_IN_PLACE, &props->with_attachables, 1, MPI_INT, MPI_LOR,
                MPI_COMM_WORLD);
#endif

  if (s->e->nodeID == 0 && !props->with_attachables && verbose)
    message("No particles to attach. Skipping the attaching phase.");

  /* Allocate and initialise a group index array. */
  if (swift_memalign("fof_group_index", (void **)&props->group_index, 64,
                     s->nr_gparts * sizeof(size_t)) != 0)
    error("Failed to allocate list of particle group indices for FOF search.");

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

/* ------------------------------------------------------------------------ */
/* TEMPORARY: instrumentation for the FOF scaling tests. Not to be merged.  */
/* ------------------------------------------------------------------------ */

ticks fof_timers[fof_timer_count];
struct fof_counter fof_counters[fof_counter_count];

static const char *fof_timer_names[fof_timer_count] = {"free_start",
                                                       "allocate",
                                                       "make_tasks",
                                                       "activate_tasks",
                                                       "local_link",
                                                       "local_sizes",
                                                       "alloc_foreign",
                                                       "foreign_setup",
                                                       "foreign_roots",
                                                       "foreign_wait",
                                                       "foreign_comms",
                                                       "foreign_search",
                                                       "foreign_merge_thr",
                                                       "foreign_wait2",
                                                       "frag_gather",
                                                       "frag_hash",
                                                       "frag_union",
                                                       "frag_update",
                                                       "ids_count",
                                                       "ids_collectives",
                                                       "ids_sort_default",
                                                       "ids_exchange",
                                                       "ids_final",
                                                       "attach_comms",
                                                       "attach",
                                                       "free_end",
                                                       "props_alloc",
                                                       "props_calc",
                                                       "props_allreduce",
                                                       "props_dump",
                                                       "props_seed",
                                                       "free_arrays",
                                                       "end_wait",
                                                       "realloc_foreign",
                                                       "total",
                                                       "task_link_self",
                                                       "task_link_pair",
                                                       "task_attach_self",
                                                       "task_attach_pair"};

static const char *fof_counter_names[fof_counter_count] = {
    "nr_gparts",
    "nr_threads",
    "tasks_link_self",
    "tasks_link_pair",
    "tasks_attach_self",
    "tasks_attach_pair",
    "leaves_link_self",
    "leaves_link_pair",
    "leaves_foreign",
    "leaves_attach_pair",
    "cand_link_self",
    "cand_link_pair",
    "cand_link_pair_bbox",
    "cand_foreign",
    "cand_foreign_bbox",
    "cand_attach_self",
    "cand_attach_pair",
    "cand_attach_pair_bbox",
    "foreign_cell_pairs",
    "foreign_cell_pairs_max",
    "gparts_foreign_recv",
    "foreign_links_found",
    "foreign_links_kept",
    "global_links",
    "num_groups",
    "ids_nsend",
    "ids_nrecv"};

/* Number of FOF calls so far (for the CSV output) */
static int fof_timers_call_count = 0;

/**
 * @brief Reset all the FOF timers and counters.
 */
void fof_timers_reset(void) {
  bzero(fof_timers, sizeof(fof_timers));
  bzero(fof_counters, sizeof(fof_counters));
}

/**
 * @brief Gather the FOF timers and counters of all ranks on rank 0, print a
 * min/mean/max summary and append the per-rank values to fof_timers.csv.
 *
 * Must be called by all the ranks.
 *
 * @param e The #engine.
 */
void fof_timers_report(const struct engine *e) {

  /* Derived quantities: the fraction of the thread time spent in the FOF
   * tasks during the two "fof" launches. */
  enum { fof_derived_eff_link, fof_derived_eff_attach, fof_derived_count };
  static const char *fof_derived_names[fof_derived_count] = {"eff_link",
                                                             "eff_attach"};

  const int nr_vals = fof_timer_count + fof_counter_count + fof_derived_count;
  double *vals = (double *)malloc(nr_vals * sizeof(double));
  if (vals == NULL) error("Failed to allocate FOF timer values");

  for (int i = 0; i < fof_timer_count; i++)
    vals[i] = clocks_from_ticks(fof_timers[i]);
  for (int i = 0; i < fof_counter_count; i++)
    vals[fof_timer_count + i] = (double)fof_counters[i].value;

  const double nr_threads = e->nr_threads;
  const double t_link = vals[fof_timer_local_link];
  const double t_attach = vals[fof_timer_attach];
  vals[fof_timer_count + fof_counter_count + fof_derived_eff_link] =
      t_link > 0.
          ? (vals[fof_timer_task_link_self] + vals[fof_timer_task_link_pair]) /
                (nr_threads * t_link)
          : 0.;
  vals[fof_timer_count + fof_counter_count + fof_derived_eff_attach] =
      t_attach > 0. ? (vals[fof_timer_task_attach_self] +
                       vals[fof_timer_task_attach_pair]) /
                          (nr_threads * t_attach)
                    : 0.;

  const int nr_nodes = e->nr_nodes;
  double *all = NULL;
  if (e->nodeID == 0) {
    all = (double *)malloc((size_t)nr_nodes * nr_vals * sizeof(double));
    if (all == NULL) error("Failed to allocate FOF timer values");
  }
#ifdef WITH_MPI
  MPI_Gather(vals, nr_vals, MPI_DOUBLE, all, nr_vals, MPI_DOUBLE, 0,
             MPI_COMM_WORLD);
#else
  memcpy(all, vals, nr_vals * sizeof(double));
#endif

  if (e->nodeID == 0) {

    message("FOF timers for call %d (step %d, %d ranks x %d threads), in %s:",
            fof_timers_call_count, e->step, nr_nodes, e->nr_threads,
            clocks_getunit());
    message("  %-22s %12s %12s %12s %6s %6s", "name", "min", "mean", "max",
            "rk_max", "max/mn");

    for (int k = 0; k < nr_vals; k++) {
      double v_min = all[k], v_max = all[k], v_sum = 0.;
      int r_max = 0;
      for (int r = 0; r < nr_nodes; r++) {
        const double v = all[(size_t)r * nr_vals + k];
        v_sum += v;
        if (v < v_min) v_min = v;
        if (v > v_max) {
          v_max = v;
          r_max = r;
        }
      }
      const double v_mean = v_sum / nr_nodes;
      const char *name =
          k < fof_timer_count ? fof_timer_names[k]
          : k < fof_timer_count + fof_counter_count
              ? fof_counter_names[k - fof_timer_count]
              : fof_derived_names[k - fof_timer_count - fof_counter_count];

      if (k == fof_timer_count)
        message("  FOF counters (per rank), sum over ranks in last column:");
      if (k == fof_timer_count + fof_counter_count)
        message("  FOF derived (per rank):");

      if (k >= fof_timer_count && k < fof_timer_count + fof_counter_count)
        message("  %-22s %12.4g %12.4g %12.4g %6d %6.2f %12.4g", name, v_min,
                v_mean, v_max, r_max, v_mean > 0. ? v_max / v_mean : 0., v_sum);
      else
        message("  %-22s %12.3f %12.3f %12.3f %6d %6.2f", name, v_min, v_mean,
                v_max, r_max, v_mean > 0. ? v_max / v_mean : 0.);
    }

    /* Append the per-rank values to the CSV file */
    const int new_file = (access("fof_timers.csv", F_OK) != 0);
    FILE *file = fopen("fof_timers.csv", "a");
    if (file == NULL) error("Could not open fof_timers.csv");
    if (new_file) {
      fprintf(file, "call,step,rank");
      for (int k = 0; k < fof_timer_count; k++)
        fprintf(file, ",%s", fof_timer_names[k]);
      for (int k = 0; k < fof_counter_count; k++)
        fprintf(file, ",%s", fof_counter_names[k]);
      for (int k = 0; k < fof_derived_count; k++)
        fprintf(file, ",%s", fof_derived_names[k]);
      fprintf(file, "\n");
    }
    for (int r = 0; r < nr_nodes; r++) {
      fprintf(file, "%d,%d,%d", fof_timers_call_count, e->step, r);
      for (int k = 0; k < nr_vals; k++)
        fprintf(file, ",%.6g", all[(size_t)r * nr_vals + k]);
      fprintf(file, "\n");
    }
    fclose(file);
    free(all);
  }

  free(vals);
  fof_timers_call_count++;
}

#if FOF_COUNT_CANDIDATES

/*! Kinds of particles for the candidate counting */
enum fof_count_kind { fof_count_linkable, fof_count_attachable };

/*! Bounding box of the particles of one kind in a leaf cell */
struct fof_count_box {
  long long n;
  double min[3], max[3];
};

/**
 * @brief Is a local #gpart of the requested kind?
 */
INLINE static int fof_count_is_kind_local(const struct gpart *gp,
                                          const enum fof_count_kind kind) {
  if (gp->time_bin >= time_bin_inhibited) return 0;
  if (gpart_is_ignorable(gp)) return 0;
  return kind == fof_count_linkable ? gpart_is_linkable(gp) != 0
                                    : gpart_is_attachable(gp) != 0;
}

/**
 * @brief Is a foreign #gpart of the requested kind?
 */
INLINE static int fof_count_is_kind_foreign(const struct gpart_fof_foreign *gp,
                                            const enum fof_count_kind kind) {
  if (gp->time_bin >= time_bin_inhibited) return 0;
  if (gpart_foreign_is_ignorable(gp)) return 0;
  return kind == fof_count_linkable ? gpart_foreign_is_linkable(gp) != 0
                                    : gpart_foreign_is_attachable(gp) != 0;
}

/**
 * @brief Square of the distance between a point and a box (0 if inside).
 */
INLINE static double fof_count_box_dist2(const struct fof_count_box *b,
                                         const double x[3]) {
  double r2 = 0.;
  for (int k = 0; k < 3; k++) {
    const double d = max3(b->min[k] - x[k], 0., x[k] - b->max[k]);
    r2 += d * d;
  }
  return r2;
}

/**
 * @brief Number and bounding box of the particles of a kind in a leaf.
 *
 * @param c The leaf #cell.
 * @param foreign Read the foreign (#gpart_fof_foreign) particles?
 * @param kind The kind of particles to consider.
 * @param shift Shift to apply to the positions.
 * @param b (return) The box.
 */
static void fof_count_make_box(const struct cell *c, const int foreign,
                               const enum fof_count_kind kind,
                               const double shift[3], struct fof_count_box *b) {
  b->n = 0;
  for (int k = 0; k < 3; k++) {
    b->min[k] = DBL_MAX;
    b->max[k] = -DBL_MAX;
  }
  for (int i = 0; i < c->grav.count; i++) {
    double x[3];
    if (foreign) {
      const struct gpart_fof_foreign *gp = &c->grav.parts_fof_foreign[i];
      if (!fof_count_is_kind_foreign(gp, kind)) continue;
      for (int k = 0; k < 3; k++) x[k] = gp->x[k] + shift[k];
    } else {
      const struct gpart *gp = &c->grav.parts[i];
      if (!fof_count_is_kind_local(gp, kind)) continue;
      for (int k = 0; k < 3; k++) x[k] = gp->x[k] + shift[k];
    }
    b->n++;
    for (int k = 0; k < 3; k++) {
      b->min[k] = min(b->min[k], x[k]);
      b->max[k] = max(b->max[k], x[k]);
    }
  }
}

/**
 * @brief Number of particles of a kind in a leaf within sqrt(r2) of a box.
 */
static long long fof_count_near_box(const struct cell *c, const int foreign,
                                    const enum fof_count_kind kind,
                                    const double shift[3],
                                    const struct fof_count_box *b,
                                    const double r2) {
  if (b->n == 0) return 0;
  long long n = 0;
  for (int i = 0; i < c->grav.count; i++) {
    double x[3];
    if (foreign) {
      const struct gpart_fof_foreign *gp = &c->grav.parts_fof_foreign[i];
      if (!fof_count_is_kind_foreign(gp, kind)) continue;
      for (int k = 0; k < 3; k++) x[k] = gp->x[k] + shift[k];
    } else {
      const struct gpart *gp = &c->grav.parts[i];
      if (!fof_count_is_kind_local(gp, kind)) continue;
      for (int k = 0; k < 3; k++) x[k] = gp->x[k] + shift[k];
    }
    if (fof_count_box_dist2(b, x) < r2) n++;
  }
  return n;
}

/**
 * @brief Count the candidate pairs of one direction (kind_i in ci against
 * kind_j in cj) of a leaf pair, and the pairs left after a bounding-box
 * pre-filter (keep only the particles within the linking length of the
 * other side's box).
 */
static void fof_count_pair_direction(const struct cell *ci, const int foreign_i,
                                     const enum fof_count_kind kind_i,
                                     const struct cell *cj, const int foreign_j,
                                     const enum fof_count_kind kind_j,
                                     const double shift_j[3], const double r2,
                                     long long *cand, long long *cand_bbox) {

  const double no_shift[3] = {0., 0., 0.};
  struct fof_count_box box_i, box_j;
  fof_count_make_box(ci, foreign_i, kind_i, no_shift, &box_i);
  fof_count_make_box(cj, foreign_j, kind_j, shift_j, &box_j);

  *cand += box_i.n * box_j.n;
  if (box_i.n == 0 || box_j.n == 0) return;

  const long long n_i =
      fof_count_near_box(ci, foreign_i, kind_i, no_shift, &box_j, r2);
  const long long n_j =
      fof_count_near_box(cj, foreign_j, kind_j, shift_j, &box_i, r2);
  *cand_bbox += n_i * n_j;
}

/**
 * @brief Periodic shift to apply to the particles of cj to bring them next to
 * ci (same convention as the pair kernels).
 */
static void fof_count_shift(const struct cell *ci, const struct cell *cj,
                            const double dim[3], const int periodic,
                            double shift_j[3]) {
  for (int k = 0; k < 3; k++) {
    const double diff = cj->loc[k] - ci->loc[k];
    if (periodic && diff < -dim[k] * 0.5)
      shift_j[k] = dim[k];
    else if (periodic && diff > dim[k] * 0.5)
      shift_j[k] = -dim[k];
    else
      shift_j[k] = 0.;
  }
}

/**
 * @brief Count the candidates of a linking leaf pair (local or foreign cj).
 */
static void fof_count_link_pair(const struct cell *ci, const struct cell *cj,
                                const int foreign_j, const double dim[3],
                                const int periodic, const double r2) {
  double shift_j[3];
  fof_count_shift(ci, cj, dim, periodic, shift_j);
  long long cand = 0, cand_bbox = 0;
  fof_count_pair_direction(ci, 0, fof_count_linkable, cj, foreign_j,
                           fof_count_linkable, shift_j, r2, &cand, &cand_bbox);
  if (foreign_j) {
    atomic_inc(&fof_counters[fof_counter_leaves_foreign].value);
    atomic_add(&fof_counters[fof_counter_cand_foreign].value, cand);
    atomic_add(&fof_counters[fof_counter_cand_foreign_bbox].value, cand_bbox);
  } else {
    atomic_inc(&fof_counters[fof_counter_leaves_link_pair].value);
    atomic_add(&fof_counters[fof_counter_cand_link_pair].value, cand);
    atomic_add(&fof_counters[fof_counter_cand_link_pair_bbox].value, cand_bbox);
  }
}

/**
 * @brief Count the candidates of an attaching leaf pair (attachables against
 * linkables, only in the directions the kernels can update).
 */
static void fof_count_attach_pair(const struct cell *ci, const struct cell *cj,
                                  const int ci_local, const int cj_local,
                                  const double dim[3], const int periodic,
                                  const double r2) {
  double shift_j[3], shift_i[3];
  fof_count_shift(ci, cj, dim, periodic, shift_j);
  fof_count_shift(cj, ci, dim, periodic, shift_i);
  long long cand = 0, cand_bbox = 0;
  if (ci_local)
    fof_count_pair_direction(ci, 0, fof_count_attachable, cj, !cj_local,
                             fof_count_linkable, shift_j, r2, &cand,
                             &cand_bbox);
  if (cj_local)
    fof_count_pair_direction(cj, 0, fof_count_attachable, ci, !ci_local,
                             fof_count_linkable, shift_i, r2, &cand,
                             &cand_bbox);
  atomic_inc(&fof_counters[fof_counter_leaves_attach_pair].value);
  atomic_add(&fof_counters[fof_counter_cand_attach_pair].value, cand);
  atomic_add(&fof_counters[fof_counter_cand_attach_pair_bbox].value, cand_bbox);
}

/**
 * @brief Count the candidates of a self leaf (linking or attaching).
 */
static void fof_count_self(const struct cell *c, const int attach) {
  long long n_link = 0, n_attach = 0;
  for (int i = 0; i < c->grav.count; i++) {
    const struct gpart *gp = &c->grav.parts[i];
    if (fof_count_is_kind_local(gp, fof_count_linkable)) n_link++;
    if (attach && fof_count_is_kind_local(gp, fof_count_attachable)) n_attach++;
  }
  if (attach) {
    atomic_add(&fof_counters[fof_counter_cand_attach_self].value,
               n_attach * n_link);
  } else {
    atomic_inc(&fof_counters[fof_counter_leaves_link_self].value);
    atomic_add(&fof_counters[fof_counter_cand_link_self].value,
               n_link * (n_link - 1) / 2);
  }
}

#endif /* FOF_COUNT_CANDIDATES */

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
 * @brief Atomically attach a root to a new root.
 *
 * The update only happens if the old root is still a root, i.e. if no other
 * thread has attached it to another group since we found it.
 *
 * @param group_index Array of group root indices.
 * @param root The current root to attach.
 * @param new_root The root to attach it to.
 *
 * @return 1 If successful, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int atomic_update_root(
    size_t *group_index, const size_t root, const size_t new_root) {

  /* atomic_cas returns root if group_index[root] was still root, i.e. if root
   * was still the root of its group. */
  return atomic_cas(&group_index[root], root, new_root) == root;
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
    const size_t root_i_new = fof_find(*root_i, group_index);
    const size_t root_j_new = fof_find(root_j, group_index);

    /* Skip particles in the same group. */
    if (root_i_new == root_j_new) return;

    /* If the root ID of pj is lower than pi's root ID set pi's root to point to
     * pj's. Otherwise set pj's root to point to pi's.*/
    if (root_j_new < root_i_new) {

      /* Updates the root if it is still a root. Otherwise another thread
       * linked it in the meantime and we try again from the new roots. */
      result = atomic_update_root(group_index, root_i_new, root_j_new);

      /* Update root_i on the fly. */
      if (result) *root_i = root_j_new;
    } else {

      /* Updates the root if it is still a root. Otherwise another thread
       * linked it in the meantime and we try again from the new roots. */
      result = atomic_update_root(group_index, root_j_new, root_i_new);

      /* Update root_i on the fly. */
      if (result) *root_i = root_i_new;
    }
  } while (result != 1);
}

/**
 * @brief Verify that the particles of a foreign cell were received before
 * reading them in a FOF operation.
 *
 * The FOF data is only exchanged for the cells (and their progeny) that carry
 * a gravity recv task. If FOF ever needs the particles of another foreign cell
 * we would read garbage, so stop here instead.
 *
 * @param c The foreign #cell.
 */
__attribute__((always_inline)) INLINE static void fof_check_foreign_cell(
    const struct cell *c) {

  if (c->grav.count > 0 && c->grav.parts_fof_foreign == NULL)
    error(
        "FOF needs the particles of a foreign cell (rank %d, depth %d, %d "
        "gparts, loc=[%e %e %e], width=[%e %e %e]) that were not received.",
        c->nodeID, c->depth, c->grav.count, c->loc[0], c->loc[1], c->loc[2],
        c->width[0], c->width[1], c->width[2]);
}

#ifdef WITH_MPI

/* Add a group to the hash table. */
__attribute__((always_inline)) INLINE static void hashmap_add_group(
    const size_t group_id, const size_t group_offset, hashmap_t *map) {

  int created_new_element = 0;
  hashmap_value_t *offset =
      hashmap_get_new(map, group_id, &created_new_element);

  if (offset != NULL) {

    /* If the element is a new entry set its value.
     * The hashmap stores a long long: fine for anything below 2^63 groups. */
    if (created_new_element) {
#ifdef SWIFT_DEBUG_CHECKS
      if (group_offset > (size_t)LLONG_MAX)
        error("Group offset %zu does not fit in the hashmap", group_offset);
#endif
      (*offset).value_st = (long long)group_offset;
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

    /* Get the nature of the linking */
    const int is_link_i = gpart_is_linkable(pi);

    /* Only linkable particles can link */
    if (!is_link_i) continue;

    /* Find the root of pi. */
    size_t root_i = fof_find(offset[i], group_index);

    for (size_t j = i + 1; j < count; j++) {

      const struct gpart *pj = &gparts[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Check whether we ignore this particle type altogether */
      if (gpart_is_ignorable(pj)) continue;

      /* Get the nature of the linking */
      const int is_link_j = gpart_is_linkable(pj);

      /* Only consider linkable<->linkable pairs */
      if (!is_link_j) continue;

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

        /* Find the root of pj (only for pairs in range). */
        const size_t root_j = fof_find(offset[j], group_index);

        /* Skip particles in the same group. */
        if (root_i == root_j) continue;

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

    /* Get the nature of the linking */
    const int is_link_i = gpart_is_linkable(pi);

    /* Only linkable particles can link */
    if (!is_link_i) continue;

    /* Find the root of pi. */
    size_t root_i = fof_find(offset_i[i], group_index);

    for (size_t j = 0; j < count_j; j++) {

      const struct gpart *restrict pj = &gparts_j[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Check whether we ignore this particle type altogether */
      if (gpart_is_ignorable(pj)) continue;

      /* Get the nature of the linking */
      const int is_link_j = gpart_is_linkable(pj);

      /* Only consider linkable<->linkable pairs */
      if (!is_link_j) continue;

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

        /* Find the root of pj (only for pairs in range). */
        const size_t root_j = fof_find(offset_j[j], group_index);

        /* Skip particles in the same group. */
        if (root_i == root_j) continue;

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
    size_t *local_link_count, size_t *group_links_size,
    struct fof_mpi **group_links, struct fof_mpi **local_group_links,
    const size_t root_i, const size_t root_j, const size_t size_i,
    const size_t size_j, const int verbose) {

  /* If the group_links array is not big enough re-allocate it. */
  if (*local_link_count + 1 > *group_links_size) {

    const size_t old_size = *group_links_size;

    /* Grow geometrically but make sure we always make progress. */
    const size_t new_size = (old_size > 0) ? 2 * old_size : 1;

    /* The doubling must not wrap around and the buffer must remain
     * addressable. */
    if (new_size <= old_size || new_size > SIZE_MAX / sizeof(struct fof_mpi))
      error("Overflow in the size of the list of foreign links (%zu elements)",
            old_size);

    struct fof_mpi *temp =
        (struct fof_mpi *)swift_realloc("fof_local_group_links", *group_links,
                                        new_size * sizeof(struct fof_mpi));
    if (temp == NULL)
      error("Failed to re-allocate the list of foreign links to %zu elements",
            new_size);

    *group_links = temp;
    *group_links_size = new_size;

    /* Reset the local pointer */
    (*local_group_links) = *group_links;

    if (verbose)
      message("Re-allocating local group links from %zu to %zu elements.",
              old_size, new_size);
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
    const struct cell *restrict cj, size_t *restrict link_count,
    struct fof_mpi **group_links, size_t *restrict group_links_size,
    const int verbose) {

#ifdef WITH_MPI
  const size_t count_i = ci->grav.count;
  const size_t count_j = cj->grav.count;
  fof_check_foreign_cell(cj);
  const struct gpart *gparts_i = ci->grav.parts;
  const struct gpart_fof_foreign *gparts_j = cj->grav.parts_fof_foreign;

  /* Get local pointers */
  const size_t *restrict group_index = props->group_index;
  const size_t *restrict group_size = props->group_size;

  /* Values local to this function to avoid dereferencing */
  struct fof_mpi *local_group_links = *group_links;
  size_t local_link_count = *link_count;

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

    /* Get the nature of the linking */
    const int is_link_i = gpart_is_linkable(pi);

    /* Only linkable particles can link */
    if (!is_link_i) continue;

    /* Find the root of pi. */
    const size_t root_i =
        fof_find_global(offset_i[i] - node_offset, group_index, nr_gparts);

    for (size_t j = 0; j < count_j; j++) {

      const struct gpart_fof_foreign *pj = &gparts_j[j];

      /* Ignore inhibited particles */
      if (pj->time_bin >= time_bin_inhibited) continue;

      /* Check whether we ignore this particle type altogether */
      if (gpart_foreign_is_ignorable(pj)) continue;

      /* Get the nature of the linking */
      const int is_link_j = gpart_foreign_is_linkable(pj);

      /* Only consider linkable<->linkable pairs */
      if (!is_link_j) continue;

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

        /* Check that the links have not already been added to the list.
         * Search backwards as the most recent links are the most likely
         * duplicates. */
        int found = 0;
        for (size_t l = local_link_count; l > 0; l--) {
          if (local_group_links[l - 1].group_i == root_i &&
              local_group_links[l - 1].group_j == pj->fof_data.group_id) {
            found = 1;
            break;
          }
        }

        /* Add a possible link to the list */
        if (!found)
          add_foreign_link_to_list(&local_link_count, group_links_size,
                                   group_links, &local_group_links, root_i,
                                   pj->fof_data.group_id,
                                   group_size[root_i - node_offset],
                                   pj->fof_data.group_size, verbose);
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
#if FOF_COUNT_CANDIDATES
    fof_count_link_pair(ci, cj, /*foreign_j=*/0, dim, periodic, search_r2);
#endif
  }
}
#ifdef WITH_MPI

/* Recurse on a pair of cells (one local, one foreign) and perform a FOF search
 * between cells that are within range. */
void rec_fof_search_pair_foreign(
    const struct fof_props *props, const double dim[3], const double search_r2,
    const int periodic, const struct gpart *const space_gparts,
    const size_t nr_gparts, const struct cell *ci, const struct cell *cj,
    size_t *restrict link_count, struct fof_mpi **group_links,
    size_t *restrict group_links_size, const int verbose) {

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
                                        group_links_size, verbose);
      }
    }
  } else if (ci->split) {

    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL)
        rec_fof_search_pair_foreign(props, dim, search_r2, periodic,
                                    space_gparts, nr_gparts, ci->progeny[k], cj,
                                    link_count, group_links, group_links_size,
                                    verbose);
    }
  } else if (cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL)
        rec_fof_search_pair_foreign(
            props, dim, search_r2, periodic, space_gparts, nr_gparts, ci,
            cj->progeny[k], link_count, group_links, group_links_size, verbose);
    }
  } else {
    /* Perform FOF search between pairs of cells that are within the linking
     * length and not the same cell. */
    fof_search_pair_cells_foreign(props, dim, search_r2, periodic, space_gparts,
                                  nr_gparts, ci, cj, link_count, group_links,
                                  group_links_size, verbose);
#if FOF_COUNT_CANDIDATES
    fof_count_link_pair(ci, cj, /*foreign_j=*/1, dim, periodic, search_r2);
#endif
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
  else {
    fof_search_self_cell(props, search_r2, space_gparts, c);
#if FOF_COUNT_CANDIDATES
    fof_count_self(c, /*attach=*/0);
#endif
  }
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
  fof_check_foreign_cell(cj);
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
  fof_check_foreign_cell(ci);
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
#if FOF_COUNT_CANDIDATES
    fof_count_attach_pair(ci, cj, ci_local, cj_local, dim, periodic, attach_r2);
#endif
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
#if FOF_COUNT_CANDIDATES
    fof_count_self(c, /*attach=*/1);
#endif
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
  size_t local_link_count = 0;
  size_t local_group_links_size = props->group_links_size / e->nr_threads;

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
                                &local_group_links_size, e->verbose);
  }

  /* TEMPORARY: FOF scaling timers. */
  const ticks tic_merge = getticks();
  atomic_add(&fof_counters[fof_counter_foreign_links_found].value,
             (long long)local_link_count);

  /* Add links found by this thread to the global link list. */
  /* Lock to prevent race conditions while adding to the global link list.*/
  if (lock_lock(&s->lock) == 0) {

    /* Get pointers to global arrays. */
    size_t *restrict group_links_size = &props->group_links_size;
    size_t *restrict group_link_count = &props->group_link_count;
    struct fof_mpi **group_links = &props->group_links;

    /* Total number of links we need to be able to store. */
    const size_t needed = *group_link_count + local_link_count;
    if (needed < *group_link_count)
      error("Overflow in the number of foreign links (%zu + %zu)",
            *group_link_count, local_link_count);

    /* If the global group_links array is not big enough re-allocate it. */
    if (needed > *group_links_size) {

      const size_t old_size = *group_links_size;
      const size_t new_size = max(needed, 2 * old_size);

      /* The growth must not wrap around and the buffer must remain
       * addressable. */
      if (new_size <= old_size || new_size > SIZE_MAX / sizeof(struct fof_mpi))
        error(
            "Overflow in the size of the list of foreign links (%zu elements)",
            old_size);

      struct fof_mpi *temp = (struct fof_mpi *)swift_realloc(
          "fof_group_links", *group_links, new_size * sizeof(struct fof_mpi));
      if (temp == NULL)
        error("Failed to re-allocate the list of foreign links to %zu elements",
              new_size);

      (*group_links) = temp;
      *group_links_size = new_size;

      if (e->verbose)
        message("Re-allocating global group links from %zu to %zu elements.",
                old_size, new_size);
    }

    /* Copy the local links to the global list. */
    for (size_t i = 0; i < local_link_count; i++) {

      int found = 0;

      /* Check that the links have not already been added to the list by another
       * thread. */
      for (size_t l = 0; l < *group_link_count; l++) {
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

  /* TEMPORARY: FOF scaling timers. */
  atomic_add(&fof_timers[fof_timer_foreign_merge_threads],
             getticks() - tic_merge);

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
 * - Group radii (maximum distance from COM to a particle),
 * - Group stellar mass,
 * - Group gas mass,
 * - Group SFR,
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
  const struct xpart *xparts = s->xparts;
  const size_t group_id_default = props->group_id_default;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double seed_halo_mass = props->seed_halo_mass;

  /* Direct pointers to the arrays */
  long long *final_group_size = props->final_group_size;
  double *group_mass = props->group_mass;
  double *centre_of_mass = props->group_centre_of_mass;
  float *radii = props->group_radii;
  char *has_black_hole = props->has_black_hole;
  float *gas_mass = props->group_gas_mass;
  float *stellar_mass = props->group_stellar_mass;
  float *star_formation_rate = props->group_star_formation_rate;
  long long *id_gas_particle_to_convert = props->id_gas_particle_to_convert;

  /* Temporary arrays to help with the CoMs */
  double *max_positions = NULL, *min_positions = NULL;
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

  float *max_part_density = NULL;
  if (swift_memalign("fof_group_max_density", (void **)&max_part_density, 32,
                     props->num_groups * sizeof(float)) != 0)
    error("Unable to allocate memory for the max densities");

  /* Initialise the max density */
  for (size_t i = 0; i < (size_t)props->num_groups; ++i) {
    max_part_density[i] = -FLT_MAX;
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

    /* Identify the densest gas particle in the group */
    if (gparts[i].type == swift_type_gas) {
      const size_t gas_index = -gparts[i].id_or_neg_offset;
      const float rho_com = hydro_get_comoving_density(&parts[gas_index]);
#ifdef SWIFT_DEBUG_CHECKS
      if (rho_com == 0.f) {
        error("Found a particle with 0-density! id=%lld", parts[gas_index].id);
      }
#endif
      max_part_density[index] = fmaxf(rho_com, max_part_density[index]);
      star_formation_rate[index] +=
          star_formation_get_SFR(&parts[gas_index], &xparts[gas_index]);
      gas_mass[index] += gparts[i].mass;
    }

    /* Add to the stellar mass */
    if (gparts[i].type == swift_type_stars) {
      stellar_mass[index] += gparts[i].mass;
    }
  }

#ifdef WITH_MPI
  /* The MPI calls below take an int number of elements (up to 3 per group) */
  if (3LL * props->num_groups > (long long)INT_MAX)
    error(
        "Too many FOF groups (%lld) for the MPI reductions of the group "
        "properties.",
        props->num_groups);

  /* TEMPORARY: FOF scaling timers. */
  const ticks tic_allreduce_1 = getticks();

  /* Now all-reduce the local fragments so that everyone has a full catalog */
  MPI_Allreduce(MPI_IN_PLACE, final_group_size, props->num_groups,
                MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, group_mass, props->num_groups, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, min_positions, 3 * props->num_groups, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, max_positions, 3 * props->num_groups, MPI_DOUBLE,
                MPI_MAX, MPI_COMM_WORLD);
  /* Note: 0/1 flags, so a byte-wise OR (MPI_CHAR is not a valid type for
   * MPI_MAX). */
  MPI_Allreduce(MPI_IN_PLACE, has_black_hole, props->num_groups, MPI_BYTE,
                MPI_BOR, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, max_part_density, props->num_groups, MPI_FLOAT,
                MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, gas_mass, props->num_groups, MPI_FLOAT, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, stellar_mass, props->num_groups, MPI_FLOAT,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, star_formation_rate, props->num_groups, MPI_FLOAT,
                MPI_SUM, MPI_COMM_WORLD);

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_props_allreduce] += getticks() - tic_allreduce_1;
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

          /* If so, is this the minimum ID we have seen with this density? */
          if (parts[gas_index].id < id_gas_particle_to_convert[index])
            id_gas_particle_to_convert[index] = parts[gas_index].id;
        }
      }
    }
  }

#ifdef WITH_MPI
  /* TEMPORARY: FOF scaling timers. */
  const ticks tic_allreduce_2 = getticks();

  /* Now all-reduce the CoMs so that everyone has a full catalog */
  MPI_Allreduce(MPI_IN_PLACE, centre_of_mass, 3 * props->num_groups, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, id_gas_particle_to_convert, props->num_groups,
                MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_props_allreduce] += getticks() - tic_allreduce_2;
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

  /* Calculate the maximum radius for each FOF & number of BHs to seed */
  for (size_t i = 0; i < nr_gparts; i++) {

    /* Ignore particles not in groups */
    if (gparts[i].fof_data.group_id == group_id_default) continue;

    /* Entry into the global list of group properties */
    const size_t index = gparts[i].fof_data.group_id - 1;

    /* Set CoM as the origin*/
    double x[3] = {gparts[i].x[0], gparts[i].x[1], gparts[i].x[2]};
    for (int k = 0; k < 3; k++) {
      if (periodic) {
        x[k] = box_wrap(x[k] + (dim[k] / 2.) - centre_of_mass[index * 3 + k],
                        0., dim[k]);
        x[k] -= dim[k] / 2.;
      } else {
        x[k] -= centre_of_mass[index * 3 + k];
      }
    }

    /* Calculate the radius*/
    const float r = sqrtf((x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]));
    radii[index] = fmaxf(radii[index], r);

    /* Should we seed a BH in this group? */
    if (!has_black_hole[index] && group_mass[index] > seed_halo_mass) {

      /* Is this a gas particle? */
      if (gparts[i].type == swift_type_gas) {
        const size_t gas_index = -gparts[i].id_or_neg_offset;

        /* Increment the local BH counter */
        if (parts[gas_index].id == id_gas_particle_to_convert[index])
          (*number_of_local_seeds)++;
      }
    }
  }

#ifdef WITH_MPI
  /* TEMPORARY: FOF scaling timers. */
  const ticks tic_allreduce_3 = getticks();

  MPI_Allreduce(MPI_IN_PLACE, radii, props->num_groups, MPI_FLOAT, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(number_of_local_seeds, number_of_global_seeds, 1, MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_props_allreduce] += getticks() - tic_allreduce_3;
#else
  *number_of_global_seeds = *number_of_local_seeds;
#endif

  /* Free temporary arrays */
  swift_free("fof_group_max_positions", max_positions);
  swift_free("fof_group_min_positions", min_positions);
  swift_free("fof_group_max_density", max_part_density);
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
  long long *id_gas_particle_to_convert = props->id_gas_particle_to_convert;

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
    const long long gas_part_id = parts[gas_index].id;

    /* Entry into the global list of group properties of this particle */
    const size_t index = gparts[i].fof_data.group_id - 1;

    /* Should we seed a BH in this group? */
    if (!has_black_hole[index] && group_mass[index] > seed_halo_mass) {

      /* Does it match the max density for this group?
       * (i.e. is it the particle we identified as the one to convert?) */
      if (gas_part_id == id_gas_particle_to_convert[index]) {

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
  float *group_radii = props->group_radii;
  float *group_gas_mass = props->group_gas_mass;
  float *group_stellar_mass = props->group_stellar_mass;
  float *group_star_formation_rate = props->group_star_formation_rate;

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
        fprintf(file,
                "# %8s %12s %12s %12s %24s %24s %12s %12s %12s %12s %12s %24s "
                "%24s \n",
                "Group ID", "Group Size", "Group Mass", "Group Radii",
                "Group Gas Mass", "Group Stellar Mass", "Group SFR", "CoM_x",
                "CoM_y", "CoM_z", "Max Density", "Max Density Local Index",
                "Particle ID");
        fprintf(file,
                "#-------------------------------------------------------------"
                "------"
                "------------------------------\n");
      }

      for (int i = 0; i < num_groups; i++) {

        fprintf(file,
                "  %8lld %12lld %12e %12e %12e %12e %12e %12e %12e %12e %12e "
                "%24lld %24lld\n",
                group_index[i], final_group_size[i], group_mass[i],
                group_radii[i], group_gas_mass[i], group_stellar_mass[i],
                group_star_formation_rate[i], group_centre_of_mass[i * 3 + 0],
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
  size_t cell_pair_size = 0;

  /* Find the maximum no. of cell pairs that can communicate. */
  for (int i = 0; i < e->nr_proxies; i++) {

    int proxy_cells_out = 0;
    for (int j = 0; j < e->proxies[i].nr_cells_out; j++) {

      /* Only include gravity cells. */
      if (e->proxies[i].cells_out_type[j] & proxy_cell_type_gravity)
        proxy_cells_out++;
    }

    int proxy_cells_in = 0;
    for (int j = 0; j < e->proxies[i].nr_cells_in; j++) {

      /* Only include gravity cells. */
      if (e->proxies[i].cells_in_type[j] & proxy_cell_type_gravity)
        proxy_cells_in++;
    }

    num_cells_out += proxy_cells_out;

    /* Pairs are only made within a proxy and only on the lower rank */
    if (engine_rank == min(engine_rank, e->proxies[i].nodeID))
      cell_pair_size += (size_t)proxy_cells_out * (size_t)proxy_cells_in;
  }

  if (cell_pair_size > INT_MAX)
    error("Too many FOF cell pairs (%zu) for an int counter.", cell_pair_size);

  if (verbose)
    message(
        "Finding max no. of cells + offset IDs"
        "took: %.3f %s.",
        clocks_from_ticks(getticks() - tic), clocks_getunit());

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_foreign_setup] = getticks() - tic_total;
  fof_counters[fof_counter_foreign_cell_pairs].value = cell_pair_count;
  fof_counters[fof_counter_foreign_cell_pairs_max].value = cell_pair_size;

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_foreign_roots] = getticks() - tic_set_roots;

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_foreign_wait] = getticks() - local_fof_tic;

  if (verbose)
    message("Local FOF imbalance: %.3f %s.",
            clocks_from_ticks(getticks() - local_fof_tic), clocks_getunit());

  tic = getticks();

  /* Perform send and receive tasks. */
  engine_launch(e, "fof comms");

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_foreign_comms] = getticks() - tic;

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_foreign_search] = getticks() - tic;
  fof_counters[fof_counter_foreign_links_kept].value = props->group_link_count;

  /* Clean up memory used by foreign particles. */
  swift_free("fof_cell_pairs", cell_pairs);

  if (verbose)
    message("Searching for foreign links took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  const ticks comms_tic = getticks();

  MPI_Barrier(MPI_COMM_WORLD);

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_foreign_wait2] = getticks() - comms_tic;

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
  if (!current_fof_attach_type || !props->with_attachables) return;

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
  const size_t group_link_count = props->group_link_count;

  /* Sum the total number of links across MPI domains over each MPI rank. */
  size_t global_group_link_count = 0;
  MPI_Allreduce(&group_link_count, &global_group_link_count, 1,
                MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  /* The global list is allocated as a single block on every rank, as are the
   * three 2x-as-long size_t arrays derived from it further down. The largest
   * of those is the list itself, so checking it covers them all. */
  if (global_group_link_count > SIZE_MAX / sizeof(struct fof_mpi))
    error(
        "Overflow of the size of the global list of foreign links (%zu "
        "elements)",
        global_group_link_count);

  struct fof_mpi *global_group_links = NULL;
  size_t *displ = NULL, *group_link_counts = NULL;

  if (swift_memalign("fof_global_group_links", (void **)&global_group_links,
                     SWIFT_STRUCT_ALIGNMENT,
                     global_group_link_count * sizeof(struct fof_mpi)) != 0)
    error("Error while allocating memory for the global list of group links");

  if (posix_memalign((void **)&group_link_counts, SWIFT_STRUCT_ALIGNMENT,
                     e->nr_nodes * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the number of group links on each "
        "MPI rank");

  if (posix_memalign((void **)&displ, SWIFT_STRUCT_ALIGNMENT,
                     e->nr_nodes * sizeof(size_t)) != 0)
    error(
        "Error while allocating memory for the displacement in memory for the "
        "global group link list");

  /* Gather the total number of links on each rank. */
  MPI_Allgather(&group_link_count, 1, MPI_UNSIGNED_LONG, group_link_counts, 1,
                MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  /* Set the displacements into the global link list using the link counts from
   * each rank */
  displ[0] = 0;
  for (int i = 1; i < e->nr_nodes; i++) {
    displ[i] = displ[i - 1] + group_link_counts[i - 1];

    /* The running sum is monotonic unless it has wrapped around. */
    if (displ[i] < displ[i - 1])
      error(
          "Overflow of the displacements into the global list of foreign "
          "links at rank %d",
          i);
  }

  /* Gather the global link list on all ranks. */
  swift_mpi_allgatherv_sizet(props->group_links, group_link_count, fof_mpi_type,
                             global_group_links, group_link_counts, displ,
                             fof_mpi_type, MPI_COMM_WORLD);

  /* Clean up memory. */
  free(group_link_counts);
  free(displ);
  swift_free("fof_group_links", props->group_links);
  props->group_links = NULL;

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_fragments_gather] = getticks() - comms_tic;
  fof_counters[fof_counter_global_links].value = global_group_link_count;

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
  const size_t global_group_list_size = 2 * global_group_link_count;

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
  size_t group_count = 0;
  for (size_t k = 0; k < global_group_link_count; k++) {

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_fragments_hash] = getticks() - tic;

  if (verbose)
    message("Global list compression took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Create a global_group_index list of groups across MPI domains so that you
   * can perform a union-find locally on each node.
   * The value of which is an offset into global_group_id, which is the actual
   * root. */
  for (size_t i = 0; i < group_count; i++) global_group_index[i] = i;

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
  for (size_t k = 0; k < global_group_link_count; k++) {

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_fragments_union] = getticks() - tic;

  if (verbose)
    message("global_group_index construction took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Update each group locally with new root information. */
  for (size_t i = 0; i < group_count; i++) {

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_fragments_update] = getticks() - tic;

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_ids_count] = getticks() - tic_num_groups_calc;

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

  /* Group IDs run from 1 to num_groups. The default ID must not be one of
   * them or the particles of that group would be treated as not in a group. */
  if (group_id_default >= 1 && group_id_default <= (size_t)num_groups)
    error(
        "FOF:group_id_default (%zu) collides with the IDs of the %lld groups "
        "found. Use 0 or a value larger than the number of groups.",
        group_id_default, num_groups);

  /* Find number of groups on lower numbered MPI ranks */
#ifdef WITH_MPI
  long long nglocal = num_groups_local;
  long long ngsum;
  MPI_Scan(&nglocal, &ngsum, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  props->num_groups_prev = (size_t)(ngsum - nglocal);
#endif /* WITH_MPI */

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_ids_collectives] = getticks() - tic;
  fof_counters[fof_counter_num_groups].value = num_groups;

  if (verbose)
    message("Finding the total no. of groups took: (FOF SCALING): %.3f %s.",
            clocks_from_ticks(getticks() - tic_num_groups_calc),
            clocks_getunit());

  /* TEMPORARY: FOF scaling timers. */
  const ticks tic_sort_default = getticks();

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_ids_sort_default] = getticks() - tic_sort_default;

  /* Assign final group IDs to local root particles where the global root is
   * on this node and the group is large enough. Within a node IDs are
   * assigned in descending order of particle number, starting from 1. */
  for (size_t i = 0; i < num_groups_local; i++) {
#ifdef WITH_MPI
    gparts[high_group_sizes[i].index - node_offset].fof_data.group_id =
        1 + i + props->num_groups_prev;
#else
    gparts[high_group_sizes[i].index].fof_data.group_id = 1 + i;
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
  /* TEMPORARY: FOF scaling timers. */
  const ticks tic_exchange = getticks();

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_ids_exchange] = getticks() - tic_exchange;
  fof_counters[fof_counter_ids_nsend].value = nsend;
  fof_counters[fof_counter_ids_nrecv].value = nrecv;

#endif /* WITH_MPI */

  /* TEMPORARY: FOF scaling timers. */
  const ticks tic_final = getticks();

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_ids_final] = getticks() - tic_final;

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
  if (swift_memalign("fof_group_radii", (void **)&props->group_radii, 32,
                     num_groups * sizeof(float)) != 0)
    error("Failed to allocate list of group radii for FOF search.");
  if (swift_memalign("fof_group_gas_mass", (void **)&props->group_gas_mass, 32,
                     num_groups * sizeof(float)) != 0)
    error("Failed to allocate list of group gas mass for FOF search.");
  if (swift_memalign("fof_group_stellar_mass",
                     (void **)&props->group_stellar_mass, 32,
                     num_groups * sizeof(float)) != 0)
    error("Failed to allocate list of group stellar mass for FOF search.");
  if (swift_memalign("fof_group_star_formation_rate",
                     (void **)&props->group_star_formation_rate, 32,
                     num_groups * sizeof(float)) != 0)
    error(
        "Failed to allocate list of group star formation rate for FOF search.");
  if (swift_memalign("fof_id_gas_particle_to_convert",
                     (void **)&props->id_gas_particle_to_convert, 32,
                     num_groups * sizeof(long long)) != 0)
    error("Failed to allocate list of particles ID to convert for FOF search.");

  bzero(props->group_mass, num_groups * sizeof(double));
  bzero(props->final_group_size, num_groups * sizeof(long long));
  bzero(props->final_group_index, num_groups * sizeof(long long));
  bzero(props->has_black_hole, num_groups * sizeof(char));
  bzero(props->group_centre_of_mass, num_groups * 3 * sizeof(double));
  bzero(props->group_radii, num_groups * sizeof(float));
  bzero(props->group_gas_mass, num_groups * sizeof(float));
  bzero(props->group_stellar_mass, num_groups * sizeof(float));
  bzero(props->group_star_formation_rate, num_groups * sizeof(float));

  /* Initialise the IDs to convert to max possible */
  for (size_t i = 0; i < (size_t)num_groups; ++i) {
    props->id_gas_particle_to_convert[i] = LLONG_MAX;
  }

  const ticks tic_props = getticks();

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_props_alloc] = tic_props - tic_total;

  /* Now, compute all things */
  int number_of_local_seeds = 0, number_of_global_seeds = 0;
  fof_calc_group_mass(props, s, &number_of_local_seeds,
                      &number_of_global_seeds);

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_props_calc] = getticks() - tic_props;

  if (verbose)
    message("Computing group properties took: %.3f %s.",
            clocks_from_ticks(getticks() - tic_props), clocks_getunit());

  /* All MPI ranks now have information about all haloes, even
   * the ones where there are no local particles in. */

  /* TEMPORARY: FOF scaling timers. */
  const ticks tic_dump = getticks();

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

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_props_dump] = getticks() - tic_dump;
  const ticks tic_seed = getticks();

  /* Seed black holes */
  if (seed_black_holes) {
    fof_seed_black_holes(props, bh_props, constants, cosmo, s,
                         number_of_local_seeds, number_of_global_seeds);
  }

  /* TEMPORARY: FOF scaling timers. */
  fof_timers[fof_timer_props_seed] = getticks() - tic_seed;

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
  swift_free("fof_group_radii", props->group_radii);
  swift_free("fof_group_gas_mass", props->group_gas_mass);
  swift_free("fof_group_stellar_mass", props->group_stellar_mass);
  swift_free("fof_group_star_formation_rate", props->group_star_formation_rate);
  swift_free("fof_id_gas_particle_to_convert",
             props->id_gas_particle_to_convert);
  swift_free("fof_has_black_hole", props->has_black_hole);
  swift_free("fof_distance", props->distance_to_link);
  swift_free("fof_group_size", props->group_size);
  props->group_mass = NULL;
  props->final_group_size = NULL;
  props->final_group_index = NULL;
  props->group_centre_of_mass = NULL;
  props->group_radii = NULL;
  props->group_gas_mass = NULL;
  props->group_stellar_mass = NULL;
  props->group_star_formation_rate = NULL;
  props->id_gas_particle_to_convert = NULL;
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
  temp.group_radii = NULL;
  temp.id_gas_particle_to_convert = NULL;
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
