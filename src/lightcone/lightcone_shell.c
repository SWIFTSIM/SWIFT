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

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>

/* HEALPix C API */
#ifdef HAVE_CHEALPIX
#include <chealpix.h>
#endif

/* Local headers */
#include "cosmology.h"
#include "engine.h"
#include "exchange_structs.h"
#include "hydro.h"
#include "lightcone/healpix_util.h"

/* This object's header. */
#include "lightcone/lightcone_shell.h"

/**
 * @brief Read in shell radii for lightcone healpix maps
 *
 * Allocates the output array, shell_out.
 *
 * @param cosmo the #cosmology structure
 * @param radius_file name of the file with shell radii
 * @param nr_shells returns the number of shells
 * @param shell_out returns the array of shells
 */
static void read_shell_radii(const struct cosmology *cosmo,
                             const char *radius_file, int *nr_shells,
                             struct lightcone_shell **shell_out) {

  /* Allow shell radii to be specified in several different units */
  enum shell_units {
    not_known = 0,
    comoving_distance = 1,
    redshift = 2,
    expansion_factor = 3
  };

  FILE *fd = fopen(radius_file, "r");
  if (!fd) error("Failed to open lightcone radius file %s", radius_file);

  /* Count number of non-zero length lines */
  size_t len = 0;
  char *line = NULL;
  int nr_lines = 0;
  while (getline(&line, &len, fd) != -1 && strlen(line) > 0) nr_lines += 1;
  rewind(fd);

  /* Allocate output array */
  struct lightcone_shell *shell = (struct lightcone_shell *)malloc(
      sizeof(struct lightcone_shell) * (nr_lines - 1));

  /* Check header */
  enum shell_units units = not_known;
  if (getline(&line, &len, fd) != -1) {
    if (strcmp(line,
               "# Minimum comoving distance, Maximum comoving distance\n") ==
        0) {
      units = comoving_distance;
    } else if (strcmp(line, "# Minimum redshift, Maximum redshift\n") == 0) {
      units = redshift;
    } else if (strcmp(
                   line,
                   "# Maximum expansion factor, Minimum expansion factor\n") ==
               0) {
      units = expansion_factor;
    } else {
      error("Unrecognized header in radius file");
    }
  } else {
    error("Unable to read header in radius file");
  }

  /* Read lines */
  for (int i = 0; i < nr_lines - 1; i += 1) {
    if (fscanf(fd, "%le, %le\n", &shell[i].rmin, &shell[i].rmax) != 2)
      error("Failed to read line from radius file");
  }
  fclose(fd);
  *nr_shells = nr_lines - 1;
  const int nr = *nr_shells;
  free(line);

  /* Convert units */
  switch (units) {
    case comoving_distance:
      /* Input is already comoving distance */
      break;
    case redshift:
      /* Convert redshift to comoving distance */
      for (int i = 0; i < nr; i += 1) {
        const double a_at_rmin = 1.0 / (1.0 + shell[i].rmin);
        shell[i].rmin = cosmology_get_comoving_distance(cosmo, a_at_rmin);
        const double a_at_rmax = 1.0 / (1.0 + shell[i].rmax);
        shell[i].rmax = cosmology_get_comoving_distance(cosmo, a_at_rmax);
      }
      break;
    case expansion_factor:
      /* Convert expansion factor to comoving distance */
      for (int i = 0; i < nr; i += 1) {
        shell[i].rmin = cosmology_get_comoving_distance(cosmo, shell[i].rmin);
        shell[i].rmax = cosmology_get_comoving_distance(cosmo, shell[i].rmax);
      }
      break;
    default:
      error("unknown unit type");
  }

  /* Do some sanity checks on the radii */
  /* All values should be monotonically increasing */
  for (int i = 1; i < nr; i += 1) {
    if (shell[i].rmin <= shell[i - 1].rmin)
      error("Minimum radii should be monotonically increasing");
    if (shell[i].rmax <= shell[i - 1].rmax)
      error("Maximum radii should be monotonically increasing");
  }

  /* Maximum radius should be greater than minimum */
  for (int i = 0; i < nr; i += 1)
    if (shell[i].rmin >= shell[i].rmax)
      error("Maximum radius should be greater than minimum");

  /* Shells should not overlap */
  for (int i = 1; i < nr; i += 1)
    if (shell[i].rmin < shell[i - 1].rmax) error("Shells should not overlap");

  /* Return pointer to array */
  *shell_out = shell;
}

/**
 * @brief Creates an array of struct lightcone_shell
 *
 * Returns a pointer to the newly allocated array. Each shell
 * contains one #lightcone_map for each healpix map to be produced
 * by this lightcone.
 *
 * @param cosmo the #cosmology structure
 * @param radius_file file with the shell radii
 * @param nr_maps number of lightcone_maps per shell
 * @param map_type specifies the types of healpix maps to make
 * @param nside healpix resolution parameter
 * @param total_nr_pix number of pixels in each map
 * @param part_type specifies which particle types update which maps
 * @param elements_per_block size of blocks used in the update buffers
 * @param nr_shells_out returns the number of lightcone shells in the array
 *
 */
struct lightcone_shell *lightcone_shell_array_init(
    const struct cosmology *cosmo, const char *radius_file, int nr_maps,
    struct lightcone_map_type *map_type, int nside, pixel_index_t total_nr_pix,
    struct lightcone_particle_type *part_type, size_t elements_per_block,
    int *nr_shells_out) {

  /* Read in the shell radii */
  int nr_shells = 0;
  struct lightcone_shell *shell = NULL;
  if (engine_rank == 0)
    read_shell_radii(cosmo, radius_file, &nr_shells, &shell);
#ifdef WITH_MPI
  MPI_Bcast(&nr_shells, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (engine_rank != 0)
    shell = (struct lightcone_shell *)malloc(sizeof(struct lightcone_shell) *
                                             nr_shells);
  MPI_Bcast(shell, sizeof(struct lightcone_shell) * nr_shells, MPI_BYTE, 0,
            MPI_COMM_WORLD);
#endif

  /* Compute expansion factor at shell edges */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    /* Inner edge of the shell */
    shell[shell_nr].amax = cosmology_scale_factor_at_comoving_distance(
        cosmo, shell[shell_nr].rmin);
    /* Outer edge of the shell */
    shell[shell_nr].amin = cosmology_scale_factor_at_comoving_distance(
        cosmo, shell[shell_nr].rmax);
  }

  /* Set initial state of the lightcone shells */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1)
    shell[shell_nr].state = shell_uninitialized;

  /* Allocate lightcone_map structs for each shell */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    shell[shell_nr].nr_maps = nr_maps;
    shell[shell_nr].map =
        (struct lightcone_map *)malloc(nr_maps * sizeof(struct lightcone_map));
  }

  int comm_rank = 0, comm_size = 1;
#ifdef WITH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif

  /* Determine how healpix maps will be distributed between MPI ranks */
  const pixel_index_t pix_per_rank = total_nr_pix / comm_size;
  if (pix_per_rank == 0) error("Must have healpix npix > number of MPI ranks!");
  const pixel_index_t local_pix_offset = comm_rank * pix_per_rank;
  pixel_index_t local_nr_pix;
  if (comm_rank < comm_size - 1)
    local_nr_pix = pix_per_rank;
  else
    local_nr_pix = total_nr_pix - (comm_size - 1) * pix_per_rank;

  /* Store this information in the shells */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    shell[shell_nr].nside = nside;
    shell[shell_nr].total_nr_pix = total_nr_pix;
    shell[shell_nr].pix_per_rank = total_nr_pix / comm_size;
    shell[shell_nr].local_nr_pix = local_nr_pix;
    shell[shell_nr].local_pix_offset = local_pix_offset;
  }

  /* Initialize lightcone_maps for each shell */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    for (int map_nr = 0; map_nr < nr_maps; map_nr += 1) {
      lightcone_map_init(&shell[shell_nr].map[map_nr], nside, total_nr_pix,
                         pix_per_rank, local_nr_pix, local_pix_offset,
                         shell[shell_nr].rmin, shell[shell_nr].rmax,
                         map_type[map_nr]);
    }
  }

  /* Initialize data buffers for map updates - one per particle type per shell
   */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
      particle_buffer_init(&shell[shell_nr].buffer[ptype],
                           part_type[ptype].buffer_element_size,
                           elements_per_block, "lightcone_map_updates");
    }
  }

  /* Return the array of shells */
  *nr_shells_out = nr_shells;
  return shell;
}

/**
 * @brief Free an array of struct lightcone_shell
 *
 * This also cleans up the lightcone_maps in the shell and the
 * update buffers.
 *
 * @param shell pointer to the array of lightcone_shells
 * @param nr_shells number of shells in the array
 */
void lightcone_shell_array_free(struct lightcone_shell *shell, int nr_shells) {

  /* Free the lightcone healpix maps for each shell */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    const int nr_maps = shell[shell_nr].nr_maps;
    for (int map_nr = 0; map_nr < nr_maps; map_nr += 1) {
      lightcone_map_clean(&shell[shell_nr].map[map_nr]);
    }
  }

  /* Free the arrays of lightcone_map structs for each shell */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    free(shell[shell_nr].map);
  }

  /* Free the buffers associated with each shell */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
      particle_buffer_free(&shell[shell_nr].buffer[ptype]);
    }
  }

  /* Free the array of shells */
  free(shell);
}

/**
 * @brief Dump the shell array to a restart file
 *
 * @param shell pointer to the array of lightcone_shells
 * @param nr_shells number of shells in the array
 * @param stream the output stream to write to
 */
void lightcone_shell_array_dump(const struct lightcone_shell *shell,
                                int nr_shells, FILE *stream) {

  /* Dump the array of shell structs  */
  restart_write_blocks((void *)shell, sizeof(struct lightcone_shell), nr_shells,
                       stream, "lightcone_shells", "lightcone_shells");

  /* Dump the lightcone maps associated with each shell */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    const int nr_maps = shell[shell_nr].nr_maps;
    for (int map_nr = 0; map_nr < nr_maps; map_nr += 1) {
      lightcone_map_struct_dump(&shell[shell_nr].map[map_nr], stream);
    }
  }
}

/**
 * @brief Restore the shell array from a restart file
 *
 * @param stream the output stream to write to
 * @param nr_shells number of shells in the array
 * @param part_type specifies which particle types update which maps
 * @param elements_per_block size of blocks used in the update buffers
 *
 */
struct lightcone_shell *lightcone_shell_array_restore(
    FILE *stream, int nr_shells, struct lightcone_particle_type *part_type,
    size_t elements_per_block) {

  /* Restore the array of lightcone_shell structs */
  struct lightcone_shell *shell = (struct lightcone_shell *)malloc(
      sizeof(struct lightcone_shell) * nr_shells);
  restart_read_blocks((void *)shell, sizeof(struct lightcone_shell), nr_shells,
                      stream, NULL, "lightcone_shells");

  /* Restore the lightcone maps associated with each shell */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    const int nr_maps = shell[shell_nr].nr_maps;
    shell[shell_nr].map =
        (struct lightcone_map *)malloc(sizeof(struct lightcone_map) * nr_maps);
    for (int map_nr = 0; map_nr < nr_maps; map_nr += 1) {
      lightcone_map_struct_restore(&shell[shell_nr].map[map_nr], stream);
    }
  }

  /* Initialise the map update buffers */
  for (int shell_nr = 0; shell_nr < nr_shells; shell_nr += 1) {
    for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
      particle_buffer_init(&shell[shell_nr].buffer[ptype],
                           part_type[ptype].buffer_element_size,
                           elements_per_block, "lightcone_map_updates");
    }
  }

  return shell;
}

struct healpix_smoothing_mapper_data {

  /*! MPI rank */
  int comm_rank, comm_size;

  /*! Pointer to the lightcone shell we're updating */
  struct lightcone_shell *shell;

  /*! Information about the particle type we're updating */
  struct lightcone_particle_type *part_type;

  /*! Pointer to the send buffer for communication */
  union lightcone_map_buffer_entry *sendbuf;

  /*! Pointer to the projected kernel table */
  struct projected_kernel_table *kernel_table;
};

#ifdef HAVE_CHEALPIX
static pixel_index_t angle_to_pixel(int nside, double theta, double phi) {
  int64_t ipring;
  ang2pix_ring64(nside, theta, phi, &ipring);
  return ipring;
}
#endif

#ifdef WITH_MPI

struct buffer_block_info {

  /*! Pointer to the buffer block */
  struct particle_buffer_block *block;

  /*! Number of elements from this block to go to each MPI rank */
  size_t *count;

  /*! Offsets at which to write elements in the send buffer */
  size_t *offset;

  /*! First destination rank each entry is to be sent to */
  int *first_dest;

  /*! Last destination rank each entry is to be sent to */
  int *last_dest;
};

#ifdef HAVE_CHEALPIX
static int pixel_to_rank(int comm_size, pixel_index_t pix_per_rank,
                         pixel_index_t pixel) {
  int rank = pixel / pix_per_rank;
  if (rank >= comm_size) rank = comm_size - 1;
  return rank;
}
#endif

/**
 * @brief Count elements to send to each rank from each buffer block
 *
 * For each buffer_block_info in the input array, this counts how
 * many lightcone map updates are to be sent to each MPI rank.
 * It also determines the range of MPI ranks which each update
 * needs to be sent to. Updates must be copied to several ranks
 * if we're smoothing the maps and the smoothing kernel overlaps parts
 * of the healpix map which are stored on different ranks.
 *
 * @param map_data Pointer to an array of buffer_block_info
 * @param num_elements Number of elements buffer_block_info array
 * @param extra_data Pointer to healpix_smoothing_mapper_data struct
 *
 */
static void count_elements_to_send_mapper(void *map_data, int num_elements,
                                          void *extra_data) {
#ifdef HAVE_CHEALPIX

  /* Unpack information about the array of blocks to process */
  struct buffer_block_info *block_info = (struct buffer_block_info *)map_data;

  /* Unpack extra input parameters we need */
  struct healpix_smoothing_mapper_data *mapper_data =
      (struct healpix_smoothing_mapper_data *)extra_data;
  struct lightcone_particle_type *part_type = mapper_data->part_type;
  struct lightcone_shell *shell = mapper_data->shell;

  /* Number of healpix maps we're updating */
  const int nr_maps = part_type->nr_maps;

  /* Number of MPI ranks we have */
  const int comm_size = mapper_data->comm_size;

  /* Maximum radius of a HEALPix pixel */
  const double max_pixrad = healpix_max_pixrad(shell->nside);

  /* Loop over buffer blocks to process */
  for (int block_nr = 0; block_nr < num_elements; block_nr += 1) {

    /* Find the count and offset for this block */
    size_t *count = block_info[block_nr].count;
    size_t *offset = block_info[block_nr].offset;
    int *first_dest = block_info[block_nr].first_dest;
    int *last_dest = block_info[block_nr].last_dest;

    /* Get a pointer to the block itself */
    struct particle_buffer_block *block = block_info[block_nr].block;

    /* Initialise count and offset into the send buffer for this block */
    for (int i = 0; i < comm_size; i += 1) {
      count[i] = 0;
      offset[i] = 0;
    }

    /* Loop over lightcone map contributions in this block */
    union lightcone_map_buffer_entry *update_data =
        (union lightcone_map_buffer_entry *)block->data;
    for (size_t i = 0; i < block->num_elements; i += 1) {

      /* Find the particle angular coordinates and size for this update */
      size_t index = i * (3 + nr_maps);
      const double theta = int_to_angle(update_data[index + 0].i);
      const double phi = int_to_angle(update_data[index + 1].i);
      /* Retrieve angular smoothing length for this particle */
      const double smoothing_radius = update_data[index + 2].f;
      /* Compute angular radius at which the projected kernel reaches zero */
      const double search_radius = smoothing_radius * kernel_gamma;

      /* Determine which MPI ranks this contribution needs to go to */
      pixel_index_t first_pixel, last_pixel;

      /* Check whether this particle updates multiple pixels */
      if (search_radius < max_pixrad) {

        /* If the radius is small, we'll just assign the contribution to one
         * pixel */
        first_pixel = last_pixel = angle_to_pixel(shell->nside, theta, phi);

      } else {

        /* If the radius is large we will update a range of pixels */
        double vec[3];
        ang2vec(theta, phi, vec);
        pixel_index_t pix_min, pix_max;
        healpix_query_disc_range(shell->nside, vec, search_radius, &pix_min,
                                 &pix_max, NULL, NULL);
        first_pixel = pix_min;
        last_pixel = pix_max;
      }

      first_dest[i] =
          pixel_to_rank(comm_size, shell->pix_per_rank, first_pixel);
      last_dest[i] = pixel_to_rank(comm_size, shell->pix_per_rank, last_pixel);

      /* Update the counts for this block */
      for (int dest = first_dest[i]; dest <= last_dest[i]; dest += 1)
        count[dest] += 1;
    }

    /* Next block */
  }
#else
  error("Need HEALPix C API for lightcone maps");
#endif
}

/**
 * @brief Store elements to send to each MPI rank from each buffer block
 *
 * This stores the updates to be sent to MPI ranks in order of which
 * rank they need to be sent to. It also duplicates updates which need
 * to go to multiple ranks.
 *
 * @param map_data Pointer to an array of buffer_block_info
 * @param num_elements Number of elements buffer_block_info array
 * @param extra_data Pointer to healpix_smoothing_mapper_data struct
 *
 */
static void store_elements_to_send_mapper(void *map_data, int num_elements,
                                          void *extra_data) {

  /* Unpack input data */
  struct buffer_block_info *block_info = (struct buffer_block_info *)map_data;
  struct healpix_smoothing_mapper_data *mapper_data =
      (struct healpix_smoothing_mapper_data *)extra_data;
  struct lightcone_particle_type *part_type = mapper_data->part_type;

  /* Find the send buffer where we will place the updates from this block */
  union lightcone_map_buffer_entry *sendbuf = mapper_data->sendbuf;

  /* Find how many elements we have per update */
  const int nr_elements_per_update = 3 + part_type->nr_maps;

  /* Loop over blocks to process on this call */
  for (int block_nr = 0; block_nr < num_elements; block_nr += 1) {

    /* Find the offset into the send buffer where we will place the
       the first element from this block to go to each MPI rank.
       Offset is in units of number of updates. */
    size_t *offset = block_info[block_nr].offset;

    /* Find range of MPI ranks to send each element in this block to */
    int *first_dest_rank = block_info[block_nr].first_dest;
    int *last_dest_rank = block_info[block_nr].last_dest;

    /* Get a pointer to the block itself */
    struct particle_buffer_block *block = block_info[block_nr].block;

    /* Loop over lightcone map updates in this block */
    union lightcone_map_buffer_entry *update_data =
        (union lightcone_map_buffer_entry *)block->data;
    for (size_t i = 0; i < block->num_elements; i += 1) {

      /* Find the data to send for this update */
      union lightcone_map_buffer_entry *block_data =
          &update_data[i * nr_elements_per_update];

      /* Store this contribution to the send buffer (possibly multiple times) */
      for (int rank = first_dest_rank[i]; rank <= last_dest_rank[i];
           rank += 1) {

        /* Find where in the send buffer to write the update */
        union lightcone_map_buffer_entry *dest =
            sendbuf + (offset[rank] * nr_elements_per_update);

        /* Copy the update to the send buffer */
        memcpy(
            dest, block_data,
            sizeof(union lightcone_map_buffer_entry) * nr_elements_per_update);
        offset[rank] += 1;
      }

      /* Next element in this block */
    }
    /* Next block */
  }
}
#endif

/**
 * @brief Mapper function for updating the healpix map
 *
 * map_data is a pointer to an array of doubles. If there are
 * N lightcone maps to update and M updates to apply then the array
 * contains (3+N)*M doubles. Each group of 3+N doubles consists of
 * (theta, phi, radius, value1, value2, ...) where theta and phi
 * are angular coordinates of the particle, radius is the angular
 * smoothing length and the values are the quantities to add to the
 * healpix maps.
 *
 * @param map_data Pointer to an array of doubles
 * @param num_elements Number of elements in map_data
 * @param extra_data Pointer to healpix_smoothing_mapper_data struct
 *
 */
void healpix_smoothing_mapper(void *map_data, int num_elements,
                              void *extra_data) {

#ifdef HAVE_CHEALPIX

  /* Unpack pointers to the lightcone shell and particle_type structs */
  struct healpix_smoothing_mapper_data *mapper_data =
      (struct healpix_smoothing_mapper_data *)extra_data;
  struct lightcone_shell *shell = mapper_data->shell;
  struct lightcone_particle_type *part_type = mapper_data->part_type;
  struct projected_kernel_table *kernel_table = mapper_data->kernel_table;

  /* Get maximum radius of any pixel in the map */
  const double max_pixrad = healpix_max_pixrad(shell->nside);

  /* Find the array of updates to apply to the healpix maps */
  union lightcone_map_buffer_entry *update_data =
      (union lightcone_map_buffer_entry *)map_data;

  /* Find range of pixel indexes stored locally. Here we assume all maps
     have the same number of pixels and distribution between MPI ranks */
  if (shell->nr_maps < 1)
    error("called on lightcone_shell which contributes to no maps");
  pixel_index_t local_pix_offset = shell->map[0].local_pix_offset;
  pixel_index_t local_nr_pix = shell->map[0].local_nr_pix;

  /* Loop over updates to apply */
  for (int i = 0; i < num_elements; i += 1) {

    /* Find the data for this update */
    size_t index = i * (3 + part_type->nr_maps);
    const double theta = int_to_angle(update_data[index + 0].i);
    const double phi = int_to_angle(update_data[index + 1].i);
    /* Retrieve angular smoothing length for this particle */
    const double smoothing_radius = update_data[index + 2].f;
    /* Compute angular radius at which the projected kernel reaches zero */
    const double search_radius = smoothing_radius * kernel_gamma;
    const union lightcone_map_buffer_entry *value = &update_data[index + 3];

    if (search_radius < max_pixrad) {

      /*
        Small particles are added to the maps directly regardless of
        whether the map is smoothed. Find the pixel index.
      */
      pixel_index_t global_pix = angle_to_pixel(shell->nside, theta, phi);

      /* Check the pixel is stored on this MPI rank */
      if ((global_pix >= local_pix_offset) &&
          (global_pix < local_pix_offset + local_nr_pix)) {

        /* Find local index of the pixel to update */
        const pixel_index_t local_pix = global_pix - local_pix_offset;

        /* Add this particle to all healpix maps */
        for (int j = 0; j < part_type->nr_maps; j += 1) {
          const int map_index = part_type->map_index[j];
          const double buffered_value = value[j].f;
          const double fac_inv = shell->map[map_index].buffer_scale_factor_inv;
          const double value_to_add = buffered_value * fac_inv;
          atomic_add_d(&shell->map[map_index].data[local_pix], value_to_add);
        }
      }

    } else {

      /*
         Large particles are SPH smoothed onto smoothed maps and just added
         to the appropriate pixel in un-smoothed maps.

         First do the smoothed maps
      */
      if (part_type->nr_smoothed_maps > 0) {

        /* Get array of ranges of pixels to update */
        double part_vec[3];
        ang2vec(theta, phi, part_vec);
        pixel_index_t pix_min, pix_max;
        int nr_ranges;
        struct pixel_range *range;
        healpix_query_disc_range(shell->nside, part_vec, search_radius,
                                 &pix_min, &pix_max, &nr_ranges, &range);

        /* Compute total weight of pixels to update */
        double total_weight = 0;
        for (int range_nr = 0; range_nr < nr_ranges; range_nr += 1) {
          for (pixel_index_t pix = range[range_nr].first;
               pix <= range[range_nr].last; pix += 1) {

            /* Get vector at the centre of this pixel */
            double pixel_vec[3];
            pix2vec_ring64(shell->nside, pix, pixel_vec);

            /* Find angle between this pixel centre and the particle.
               Dot product may be a tiny bit greater than one due to rounding
               error */
            const double dp =
                (pixel_vec[0] * part_vec[0] + pixel_vec[1] * part_vec[1] +
                 pixel_vec[2] * part_vec[2]);
            const double angle = dp < 1.0 ? acos(dp) : 0.0;

            /* Evaluate the kernel at this radius */
            total_weight +=
                projected_kernel_eval(kernel_table, angle / smoothing_radius);
          }
        }

        /* Update the pixels */
        for (int range_nr = 0; range_nr < nr_ranges; range_nr += 1) {
          for (pixel_index_t pix = range[range_nr].first;
               pix <= range[range_nr].last; pix += 1) {

            /* Check if this pixel is stored locally */
            pixel_index_t global_pix = pix;
            if ((global_pix >= local_pix_offset) &&
                (global_pix < local_pix_offset + local_nr_pix)) {

              /* Get vector at the centre of this pixel */
              double pixel_vec[3];
              pix2vec_ring64(shell->nside, pix, pixel_vec);

              /* Find angle between this pixel centre and the particle.
                 Dot product may be a tiny bit greater than one due to rounding
                 error */
              const double dp =
                  (pixel_vec[0] * part_vec[0] + pixel_vec[1] * part_vec[1] +
                   pixel_vec[2] * part_vec[2]);
              const double angle = dp < 1.0 ? acos(dp) : 0.0;

              /* Evaluate the kernel at this radius */
              const double weight =
                  projected_kernel_eval(kernel_table,
                                        angle / smoothing_radius) /
                  total_weight;

              /* Find local index of the pixel to update */
              const pixel_index_t local_pix = global_pix - local_pix_offset;

              /* Update the smoothed healpix maps */
              for (int j = 0; j < part_type->nr_smoothed_maps; j += 1) {
                const int map_index = part_type->map_index[j];
                const double buffered_value = value[j].f;
                const double fac_inv =
                    shell->map[map_index].buffer_scale_factor_inv;
                const double value_to_add = buffered_value * fac_inv;
                atomic_add_d(&shell->map[map_index].data[local_pix],
                             value_to_add * weight);
              } /* Next smoothed map */
            }
          } /* Next pixel in this range */
        }   /* Next range of pixels */

        /* Free array of pixel ranges */
        free(range);

      } /* if nr_smoothed_maps > 0*/

      /* Then do any un-smoothed maps */
      if (part_type->nr_unsmoothed_maps > 0) {

        /* Find the index of the pixel containing the particle */
        pixel_index_t global_pix = angle_to_pixel(shell->nside, theta, phi);

        /* Check the pixel is stored on this MPI rank */
        if ((global_pix >= local_pix_offset) &&
            (global_pix < local_pix_offset + local_nr_pix)) {

          /* Find local index of the pixel to update */
          const pixel_index_t local_pix = global_pix - local_pix_offset;

          /* Update the un-smoothed healpix maps */
          for (int j = part_type->nr_smoothed_maps; j < part_type->nr_maps;
               j += 1) {
            const int map_index = part_type->map_index[j];
            const double buffered_value = value[j].f;
            const double fac_inv =
                shell->map[map_index].buffer_scale_factor_inv;
            const double value_to_add = buffered_value * fac_inv;
            atomic_add_d(&shell->map[map_index].data[local_pix], value_to_add);
          }
        }
      } /* if part_type->nr_unsmoothed_maps > 0 */
    }
  } /* End loop over updates to apply */
#else
  error("Need HEALPix C API for lightcone maps");
#endif
}

/**
 * @brief Apply updates for one particle type to all lightcone maps in a shell
 *
 * When a particle of type ptype crosses the lightcone it generates an entry
 * in shell->buffer[ptype] which contains the angular position and size of
 * the particle and the values it contributes to the lightcone_maps in the
 * shell. This function applies these buffered updates to the lightcone
 * map pixel data.
 *
 * We carry out all the updates for one particle type at the same time so that
 * we avoid repeating the healpix neighbour search for every healpix map.
 *
 * Applying the updates involves copying them to a send buffer then a receive
 * buffer, so if there are a lot we process them in chunks of up to
 * max_map_update_send_size_mb megabytes to save memory.
 *
 * @param shell the #lightcone_shell to update
 * @param tp the #threadpool used to execute the updates
 * @param part_type contains information about each particle type to be updated
 * @param smoothing_info contains parameters relating to smoothing onto the
 * sphere
 * @param ptype index of the particle type to update
 * @param max_map_update_send_size_mb maximum amount of data each ranks sends
 *
 */
void lightcone_shell_flush_map_updates_for_type(
    struct lightcone_shell *shell, struct threadpool *tp,
    struct lightcone_particle_type *part_type, int ptype,
    const double max_map_update_send_size_mb,
    struct projected_kernel_table *kernel_table, int verbose) {

  int comm_rank = 0, comm_size = 1;
#ifdef WITH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif

  /* Information needed by mapper functions */
  struct healpix_smoothing_mapper_data mapper_data;
  mapper_data.shell = shell;
  mapper_data.part_type = &part_type[ptype];
  mapper_data.comm_rank = comm_rank;
  mapper_data.comm_size = comm_size;
  mapper_data.sendbuf = NULL;
  mapper_data.kernel_table = kernel_table;

#ifdef WITH_MPI

  /* Count data blocks and ensure number of elements is in range */
  int nr_blocks = 0;
  struct particle_buffer *buffer = &shell->buffer[ptype];
  struct particle_buffer_block *block = buffer->first_block;
  while (block) {
    if (block->num_elements > buffer->elements_per_block)
      block->num_elements = buffer->elements_per_block;
    nr_blocks += 1;
    block = block->next;
  }

  /* Allocate array with counts and offsets for each block */
  struct buffer_block_info *block_info = (struct buffer_block_info *)malloc(
      sizeof(struct buffer_block_info) * nr_blocks);

  /* Initialize array of blocks */
  nr_blocks = 0;
  block = buffer->first_block;
  while (block) {
    block_info[nr_blocks].block = block;
    block_info[nr_blocks].count = (size_t *)malloc(sizeof(size_t) * comm_size);
    block_info[nr_blocks].offset = (size_t *)malloc(sizeof(size_t) * comm_size);
    block_info[nr_blocks].first_dest =
        (int *)malloc(sizeof(int) * block->num_elements);
    block_info[nr_blocks].last_dest =
        (int *)malloc(sizeof(int) * block->num_elements);
    nr_blocks += 1;
    block = block->next;
  }

  /* To minimize memory usage we don't process all of the blocks at once.
     Determine the maximum number of blocks on any rank. */
  int max_nr_blocks;
  MPI_Allreduce(&nr_blocks, &max_nr_blocks, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);

  /* Determine the maximum number of blocks to process per iteration */
  size_t max_bytes = max_map_update_send_size_mb * 1024.0 * 1024.0;
  int max_blocks_per_iteration =
      max_bytes / (buffer->element_size * buffer->elements_per_block);
  if (max_blocks_per_iteration < 1)
    error("max_map_update_send_size_mb is too small to process even one block");

  /* Determine how many iterations we need */
  int nr_iterations = max_nr_blocks / max_blocks_per_iteration;
  if (max_nr_blocks % max_blocks_per_iteration != 0) nr_iterations += 1;
  if (engine_rank == 0 && nr_iterations > 0 && verbose)
    message("will require %d iterations with %d blocks per iteration",
            nr_iterations, max_blocks_per_iteration);

  /* Loop over iterations */
  int nr_blocks_done = 0;
  for (int iter = 0; iter < nr_iterations; iter += 1) {

    /* Find number of blocks to do on this iteration (may be zero) */
    int nr_blocks_iter = nr_blocks - nr_blocks_done;
    if (nr_blocks_iter > max_blocks_per_iteration)
      nr_blocks_iter = max_blocks_per_iteration;

    /* Get a pointer to the blocks to do on this iteration */
    struct buffer_block_info *block_info_iter = block_info + nr_blocks_done;

    /* For each block, count how many elements are to be sent to each MPI rank
     */
    threadpool_map(tp, count_elements_to_send_mapper, block_info_iter,
                   nr_blocks_iter, sizeof(struct buffer_block_info), 1,
                   &mapper_data);

    /* Find total number of elements to go to each rank */
    size_t *send_count = (size_t *)malloc(sizeof(size_t) * comm_size);
    for (int i = 0; i < comm_size; i += 1) send_count[i] = 0;
    for (int block_nr = 0; block_nr < nr_blocks_iter; block_nr += 1) {
      for (int i = 0; i < comm_size; i += 1)
        send_count[i] += block_info_iter[block_nr].count[i];
    }

    /* Find offset to the first element to go to each rank if we sort them by
     * destination */
    size_t *send_offset = (size_t *)malloc(sizeof(size_t) * comm_size);
    send_offset[0] = 0;
    for (int i = 1; i < comm_size; i += 1) {
      send_offset[i] = send_offset[i - 1] + send_count[i - 1];
    }

    /* For each block, find the location in the send buffer where we need to
       place the first element to go to each MPI rank */
    for (int block_nr = 0; block_nr < nr_blocks_iter; block_nr += 1) {
      for (int i = 0; i < comm_size; i += 1) {
        if (block_nr == 0) {
          /* This is the first block */
          block_info_iter[block_nr].offset[i] = send_offset[i];
        } else {
          /* Not first, so elements are written after those of the previous
           * block */
          block_info_iter[block_nr].offset[i] =
              block_info_iter[block_nr - 1].offset[i] +
              block_info_iter[block_nr - 1].count[i];
        }
      }
    }

    /* Find the total number of elements to be sent */
    size_t total_nr_send = 0;
    for (int i = 0; i < comm_size; i += 1) total_nr_send += send_count[i];

    /* Allocate the send buffer */
    union lightcone_map_buffer_entry *sendbuf =
        (union lightcone_map_buffer_entry *)malloc(
            part_type[ptype].buffer_element_size * total_nr_send);
    mapper_data.sendbuf = sendbuf;

    /* Populate the send buffer */
    threadpool_map(tp, store_elements_to_send_mapper, block_info_iter,
                   nr_blocks_iter, sizeof(struct buffer_block_info), 1,
                   &mapper_data);

    /* Determine number of elements to receive */
    size_t *recv_count = (size_t *)malloc(comm_size * sizeof(size_t));
    MPI_Alltoall(send_count, sizeof(size_t), MPI_BYTE, recv_count,
                 sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);
    size_t total_nr_recv = 0;
    for (int i = 0; i < comm_size; i += 1) total_nr_recv += recv_count[i];

    /* Allocate receive buffer */
    union lightcone_map_buffer_entry *recvbuf =
        (union lightcone_map_buffer_entry *)malloc(
            part_type[ptype].buffer_element_size * total_nr_recv);

    /* Exchange data */
    exchange_structs(send_count, sendbuf, recv_count, recvbuf,
                     part_type[ptype].buffer_element_size);

    /* Apply received updates to the healpix map */
    threadpool_map(tp, healpix_smoothing_mapper, recvbuf, total_nr_recv,
                   part_type[ptype].buffer_element_size,
                   threadpool_auto_chunk_size, &mapper_data);

    /* Tidy up */
    free(send_count);
    free(send_offset);
    free(sendbuf);
    free(recv_count);
    free(recvbuf);

    /* Advance to next set of blocks */
    nr_blocks_done += nr_blocks_iter;
  }
  if (nr_blocks_done != nr_blocks)
    error("not all map update blocks were processed");

  /* We no longer need the array of blocks */
  for (int block_nr = 0; block_nr < nr_blocks; block_nr += 1) {
    free(block_info[block_nr].count);
    free(block_info[block_nr].offset);
    free(block_info[block_nr].first_dest);
    free(block_info[block_nr].last_dest);
  }
  free(block_info);

  /* Empty the particle buffer now that we copied the data from it */
  particle_buffer_empty(buffer);

#else

  /* If not using MPI, we can update the healpix maps directly from the buffer
   */
  struct particle_buffer_block *block = NULL;
  size_t num_elements;
  double *update_data;
  do {
    particle_buffer_iterate(&shell->buffer[ptype], &block, &num_elements,
                            (void **)&update_data);
    threadpool_map(tp, healpix_smoothing_mapper, update_data, num_elements,
                   part_type[ptype].buffer_element_size,
                   threadpool_auto_chunk_size, &mapper_data);
  } while (block);
  particle_buffer_empty(&shell->buffer[ptype]);

#endif
}

/**
 * @brief Apply buffered updates to all lightcone maps in a shell
 *
 * @param shell the #lightcone_shell to update
 * @param tp the #threadpool used to execute the updates
 * @param part_type contains information about each particle type to be updated
 * sphere
 *
 */
void lightcone_shell_flush_map_updates(
    struct lightcone_shell *shell, struct threadpool *tp,
    struct lightcone_particle_type *part_type,
    const double max_map_update_send_size_mb,
    struct projected_kernel_table *kernel_table, int verbose) {

  if (shell->state != shell_current)
    error("Attempt to flush updates for non-current shell!");

  for (int ptype = 0; ptype < swift_type_count; ptype += 1) {
    if ((shell->nr_maps > 0) && (part_type[ptype].nr_maps > 0)) {
      lightcone_shell_flush_map_updates_for_type(shell, tp, part_type, ptype,
                                                 max_map_update_send_size_mb,
                                                 kernel_table, verbose);
    }
  }
}
