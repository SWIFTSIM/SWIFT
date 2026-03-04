/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Peter W. Draper (p.w.draper@durham.ac.uk)
 *                    Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/**
 *  @file partition.c
 *  @brief file of various techniques for partitioning and repartitioning
 *  a grid of cells into geometrically connected regions and distributing
 *  these around a number of MPI nodes.
 *
 *  Currently supported partitioning types: grid, vectorise and METIS/ParMETIS.
 */

/* Config parameters. */
#include <config.h>

/* Standard headers. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

/* Include int min and max values. Define these limits in C++ as well. */
#define __STDC_LIMIT_MACROS
#include <stdint.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
/* METIS/ParMETIS headers only used when MPI is also available. */
#ifdef HAVE_PARMETIS
#include <parmetis.h>
#endif
#ifdef HAVE_METIS
#include <metis.h>
#endif
#endif

/* Local headers. */
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "partition.h"
#include "restart.h"
#include "space.h"
#include "threadpool.h"
#include "tools.h"

/* Simple descriptions of initial partition types for reports. */
const char *initial_partition_name[] = {
    "axis aligned grids of cells", "vectorized point associated cells",
    "memory balanced, using particle weighted cells",
    "similar sized regions, using unweighted cells",
    "memory and edge balanced cells using particle weights"};

/* Simple descriptions of repartition types for reports. */
const char *repartition_name[] = {
    "none", "edge and vertex task cost weights", "task cost edge weights",
    "memory balanced, using particle vertex weights",
    "vertex task costs and edge delta timebin weights"};

/* Local functions, if needed. */
/* TODO: We don't actually need this? */
static int check_complete(struct space *s, int verbose, int nregions);

/*
 * Repartition fixed costs per type/subtype. These are determined from the
 * statistics output produced when running with task debugging enabled.
 */
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
double repartition_costs[task_type_count][task_subtype_count];
#endif
#if defined(WITH_MPI)
static int repart_init_fixed_costs(void);
#endif

/*  Vectorisation support */
/*  ===================== */

#if defined(WITH_MPI)
/**
 *  @brief Pick a number of cell positions from a vectorised list.
 *
 *  Vectorise the cell space and pick positions in it for the number of
 *  expected regions using a single step. Vectorisation is guaranteed
 *  to work, providing there are more cells than regions.
 *
 *  @param s the space.
 *  @param nregions the number of regions
 *  @param samplecells the list of sample cell positions, size of 3*nregions
 */
static void pick_vector(struct space *s, int nregions, int *samplecells) {

  /* Get length of space and divide up. */
  int length = s->cdim[0] * s->cdim[1] * s->cdim[2];
  if (nregions > length) {
    error("Too few cells (%d) for this number of regions (%d)", length,
          nregions);
  }

  int step = length / nregions;
  int n = 0;
  int m = 0;
  int l = 0;

  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        if (n == 0 && l < nregions) {
          samplecells[m++] = i;
          samplecells[m++] = j;
          samplecells[m++] = k;
          l++;
        }
        n++;
        if (n == step) n = 0;
      }
    }
  }
}
#endif

#if defined(WITH_MPI)
/**
 * @brief Partition the space.
 *
 * Using the sample positions as seeds pick cells that are geometrically
 * closest and apply the partition to the space.
 */
static void split_vector(struct space *s, int nregions, int *samplecells) {
  int n = 0;
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        int select = -1;
        float rsqmax = FLT_MAX;
        int m = 0;
        for (int l = 0; l < nregions; l++) {
          float dx = samplecells[m++] - i;
          float dy = samplecells[m++] - j;
          float dz = samplecells[m++] - k;
          float rsq = (dx * dx + dy * dy + dz * dz);
          if (rsq < rsqmax) {
            rsqmax = rsq;
            select = l;
          }
        }
        s->cells_top[n++].nodeID = select;
      }
    }
  }
}
#endif

/**
 * @brief Initial partition of space cells.
 *
 * Cells are assigned to a node on the basis of various schemes, all of which
 * should attempt to distribute them in geometrically close regions to
 * minimise the movement of particles.
 *
 * Note that the partition type is a suggestion and will be ignored if that
 * scheme fails. In that case we fallback to a vectorised scheme, that is
 * guaranteed to work provided we have more cells than nodes.
 *
 * @param initial_partition the type of partitioning to try.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells.
 */
void partition_initial_partition(struct partition *initial_partition,
                                 int nodeID, int nr_nodes, struct space *s) {
  ticks tic = getticks();

  /* Geometric grid partitioning. */
  if (initial_partition->type == INITPART_GRID) {
    int j, k;
    int ind[3];
    struct cell *c;

    /* If we've got the wrong number of nodes, fail. */
    if (nr_nodes != initial_partition->grid[0] * initial_partition->grid[1] *
                        initial_partition->grid[2])
      error("Grid size does not match number of nodes.");

    /* Run through the cells and set their nodeID. */
    for (k = 0; k < s->nr_cells; k++) {
      c = &s->cells_top[k];
      for (j = 0; j < 3; j++)
        ind[j] = c->loc[j] / s->dim[j] * initial_partition->grid[j];
      c->nodeID = ind[0] + initial_partition->grid[0] *
                               (ind[1] + initial_partition->grid[1] * ind[2]);
    }

    /* The grid technique can fail, so check for this before proceeding. */
    if (!check_complete(s, (nodeID == 0), nr_nodes)) {
      if (nodeID == 0)
        message("Grid initial partition failed, using a vectorised partition");
      initial_partition->type = INITPART_VECTORIZE;
      partition_initial_partition(initial_partition, nodeID, nr_nodes, s);
      return;
    }

  } else if (initial_partition->type == INITPART_METIS_WEIGHT ||
             initial_partition->type == INITPART_METIS_WEIGHT_EDGE ||
             initial_partition->type == INITPART_METIS_NOWEIGHT) {
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
    /* Simple k-way partition selected by METIS using cell particle
     * counts as weights or not. Should be best when starting with a
     * inhomogeneous dist.
     */
    double *weights_v = NULL;
    double *weights_e = NULL;
    if (initial_partition->type == INITPART_METIS_WEIGHT) {
      /* Particles sizes per cell, which will be used as weights. */
      if ((weights_v = (double *)malloc(sizeof(double) * s->nr_cells)) == NULL)
        error("Failed to allocate weights_v buffer.");

      /* Check each particle and accumulate the sizes per cell. */
      accumulate_sizes(s, s->e->verbose, weights_v);

    } else if (initial_partition->type == INITPART_METIS_WEIGHT_EDGE) {

      /* Particle sizes also counted towards the edges. */

      if ((weights_v = (double *)malloc(sizeof(double) * s->nr_cells)) == NULL)
        error("Failed to allocate weights_v buffer.");
      if ((weights_e = (double *)malloc(sizeof(double) * s->nr_cells * 26)) ==
          NULL)
        error("Failed to allocate weights_e buffer.");

      /* Check each particle and accumulate the sizes per cell. */
      accumulate_sizes(s, s->e->verbose, weights_v);

      /* Spread these into edge weights. */
      sizes_to_edges(s, weights_v, weights_e);
    }

    /* Do the calculation. */
    int *celllist = NULL;
    if ((celllist = (int *)malloc(sizeof(int) * s->nr_cells)) == NULL)
      error("Failed to allocate celllist");
#ifdef HAVE_PARMETIS
    if (initial_partition->usemetis) {
      pick_metis(nodeID, s, nr_nodes, weights_v, weights_e, celllist);
    } else {
      pick_parmetis(nodeID, s, nr_nodes, weights_v, weights_e, 0, 0, 0.0f,
                    celllist);
    }
#else
    pick_metis(nodeID, s, nr_nodes, weights_v, weights_e, celllist);
#endif

    /* And apply to our cells */
    split_metis(s, nr_nodes, celllist);

    /* It's not known if this can fail, but check for this before
     * proceeding. */
    if (!check_complete(s, (nodeID == 0), nr_nodes)) {
      if (nodeID == 0)
        message("METIS initial partition failed, using a vectorised partition");
      initial_partition->type = INITPART_VECTORIZE;
      partition_initial_partition(initial_partition, nodeID, nr_nodes, s);
    }

    if (weights_v != NULL) free(weights_v);
    if (weights_e != NULL) free(weights_e);
    free(celllist);
#else
    error("SWIFT was not compiled with METIS or ParMETIS support");
#endif

  } else if (initial_partition->type == INITPART_VECTORIZE) {

#if defined(WITH_MPI)
    /* Vectorised selection, guaranteed to work for samples less than the
     * number of cells, but not very clumpy in the selection of regions. */
    int *samplecells = NULL;
    if ((samplecells = (int *)malloc(sizeof(int) * nr_nodes * 3)) == NULL)
      error("Failed to allocate samplecells");

    if (nodeID == 0) {
      pick_vector(s, nr_nodes, samplecells);
    }

    /* Share the samplecells around all the nodes. */
    int res = MPI_Bcast(samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS)
      mpi_error(res, "Failed to bcast the partition sample cells.");

    /* And apply to our cells */
    split_vector(s, nr_nodes, samplecells);
    free(samplecells);
#else
    error("SWIFT was not compiled with MPI support");
#endif
  }

  if (s->e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Initialises the partition and re-partition scheme from the parameter
 *        file.
 *
 * @param partition The #partition scheme to initialise.
 * @param repartition The #repartition scheme to initialise.
 * @param params The parsed parameter file.
 * @param nr_nodes The number of MPI nodes we are running on.
 */
void partition_init(struct partition *partition,
                    struct repartition *repartition,
                    struct swift_params *params, int nr_nodes) {

#ifdef WITH_MPI

  /* Defaults make use of METIS if available */
#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)
  const char *default_repart = "fullcosts";
  const char *default_part = "edgememory";
#else
  const char *default_repart = "none";
  const char *default_part = "grid";
#endif

  /* Set a default grid so that grid[0]*grid[1]*grid[2] == nr_nodes. */
  factor(nr_nodes, &partition->grid[0], &partition->grid[1]);
  factor(nr_nodes / partition->grid[1], &partition->grid[0],
         &partition->grid[2]);
  factor(partition->grid[0] * partition->grid[1], &partition->grid[1],
         &partition->grid[0]);

  /* Initialise the repartition celllist. */
  repartition->ncelllist = 0;
  repartition->celllist = NULL;

  /* Now let's check what the user wants as an initial domain. */
  char part_type[20];
  parser_get_opt_param_string(params, "DomainDecomposition:initial_type",
                              part_type, default_part);
  switch (part_type[0]) {
    case 'g':
      partition->type = INITPART_GRID;
      break;
    case 'v':
      partition->type = INITPART_VECTORIZE;
      break;
#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)
    case 'r':
      partition->type = INITPART_METIS_NOWEIGHT;
      break;
    case 'm':
      partition->type = INITPART_METIS_WEIGHT;
      break;
    case 'e':
      partition->type = INITPART_METIS_WEIGHT_EDGE;
      break;
    default:
      message("Invalid choice of initial partition type '%s'.", part_type);
      error(
          "Permitted values are: 'grid', 'region', 'memory', 'edgememory' or "
          "'vectorized'");
#else
    default:
      message("Invalid choice of initial partition type '%s'.", part_type);
      error(
          "Permitted values are: 'grid' or 'vectorized' when compiled "
          "without METIS or ParMETIS.");
#endif
  }

  /* In case of grid, read more parameters */
  if (part_type[0] == 'g') {
    parser_get_opt_param_int_array(params, "DomainDecomposition:initial_grid",
                                   3, partition->grid);
  }

  /* Now let's check what the user wants as a repartition strategy */
  parser_get_opt_param_string(params, "DomainDecomposition:repartition_type",
                              part_type, default_repart);

  if (strcmp("none", part_type) == 0) {
    repartition->type = REPART_NONE;

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)
  } else if (strcmp("fullcosts", part_type) == 0) {
    repartition->type = REPART_METIS_VERTEX_EDGE_COSTS;

  } else if (strcmp("edgecosts", part_type) == 0) {
    repartition->type = REPART_METIS_EDGE_COSTS;

  } else if (strcmp("memory", part_type) == 0) {
    repartition->type = REPART_METIS_VERTEX_COUNTS;

  } else if (strcmp("timecosts", part_type) == 0) {
    repartition->type = REPART_METIS_VERTEX_COSTS_TIMEBINS;

  } else {
    message("Invalid choice of re-partition type '%s'.", part_type);
    error(
        "Permitted values are: 'none', 'fullcosts', 'edgecosts' "
        "'memory' or 'timecosts'");
#else
  } else {
    message("Invalid choice of re-partition type '%s'.", part_type);
    error(
        "Permitted values are: 'none' when compiled without "
        "METIS or ParMETIS.");
#endif
  }

  /* Get the fraction CPU time difference between nodes (<1) or the number
   * of steps between repartitions (>1). */
  repartition->trigger =
      parser_get_opt_param_float(params, "DomainDecomposition:trigger", 0.05f);
  if (repartition->trigger <= 0)
    error("Invalid DomainDecomposition:trigger, must be greater than zero");
  if (repartition->trigger < 2 && repartition->trigger >= 1)
    error(
        "Invalid DomainDecomposition:trigger, must be 2 or greater or less"
        " than 1");

  /* Fraction of particles that should be updated before a repartition
   * based on CPU time is considered, needs to be high. */
  repartition->minfrac =
      parser_get_opt_param_float(params, "DomainDecomposition:minfrac", 0.95f);
  if (repartition->minfrac <= 0.5 || repartition->minfrac > 1)
    error(
        "Invalid DomainDecomposition:minfrac, must be greater than 0.5 "
        "and less than equal to 1");

  /* Use METIS or ParMETIS when ParMETIS is also available. */
  repartition->usemetis =
      parser_get_opt_param_int(params, "DomainDecomposition:usemetis", 0);
  partition->usemetis = repartition->usemetis;

  /* Use adaptive or simple refinement when repartitioning. */
  repartition->adaptive =
      parser_get_opt_param_int(params, "DomainDecomposition:adaptive", 1);

  /* Ratio of interprocess communication time to data redistribution time. */
  repartition->itr =
      parser_get_opt_param_float(params, "DomainDecomposition:itr", 100.0f);

  /* Do we have fixed costs available? These can be used to force
   * repartitioning at any time. Not required if not repartitioning.*/
  repartition->use_fixed_costs = parser_get_opt_param_int(
      params, "DomainDecomposition:use_fixed_costs", 0);
  if (repartition->type == REPART_NONE) repartition->use_fixed_costs = 0;

  /* Check if this is true or required and initialise them. */
  if (repartition->use_fixed_costs || repartition->trigger > 1) {
    if (!repart_init_fixed_costs()) {
      if (repartition->trigger <= 1) {
        if (engine_rank == 0)
          message(
              "WARNING: fixed cost repartitioning was requested but is"
              " not available.");
        repartition->use_fixed_costs = 0;
      } else {
        error(
            "Forced fixed cost repartitioning was requested but is"
            " not available.");
      }
    }
  }

#else
  error("SWIFT was not compiled with MPI support");
#endif
}

/**
 * @brief Clean up any allocated resources.
 *
 * @param partition The #partition
 * @param repartition The #repartition
 */
void partition_clean(struct partition *partition,
                     struct repartition *repartition) {
#ifdef WITH_MPI
  /* Only the celllist is dynamic. */
  if (repartition->celllist != NULL) free(repartition->celllist);

  /* Zero structs for reuse. */
  bzero(partition, sizeof(struct partition));
  bzero(repartition, sizeof(struct repartition));
#endif
}

#ifdef WITH_MPI
/**
 * @brief Set the fixed costs for repartition using METIS.
 *
 *  These are determined using a run with the -y flag on which produces
 *  a statistical analysis that is condensed into a .h file for inclusion.
 *
 *  If the default include file is used then no fixed costs are set and this
 *  function will return 0.
 */
static int repart_init_fixed_costs(void) {

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
  /* Set the default fixed cost. */
  for (int j = 0; j < task_type_count; j++) {
    for (int k = 0; k < task_subtype_count; k++) {
      repartition_costs[j][k] = 1.0;
    }
  }

#include <partition_fixed_costs.h>
  return HAVE_FIXED_COSTS;
#endif

  return 0;
}
#endif /* WITH_MPI */

/*  General support */
/*  =============== */

/**
 * @brief Check if all regions have been assigned a node in the
 *        cells of a space.
 *
 * @param s the space containing the cells to check.
 * @param nregions number of regions expected.
 * @param verbose if true report the missing regions.
 * @return true if all regions have been found, false otherwise.
 */
static int check_complete(struct space *s, int verbose, int nregions) {

  int *present = NULL;
  if ((present = (int *)malloc(sizeof(int) * nregions)) == NULL)
    error("Failed to allocate present array");

  int failed = 0;
  for (int i = 0; i < nregions; i++) present[i] = 0;
  for (int i = 0; i < s->nr_cells; i++) {
    if (s->cells_top[i].nodeID <= nregions)
      present[s->cells_top[i].nodeID]++;
    else
      message("Bad nodeID: s->cells_top[%d].nodeID = %d", i,
              s->cells_top[i].nodeID);
  }
  for (int i = 0; i < nregions; i++) {
    if (!present[i]) {
      failed = 1;
      if (verbose) message("Region %d is not present in partition", i);
    }
  }
  free(present);
  return (!failed);
}

/**
 * @brief Partition a space of cells based on another space of cells.
 *
 * The two spaces are expected to be at different cell sizes, so what we'd
 * like to do is assign the second space to geometrically closest nodes
 * of the first, with the effect of minimizing particle movement when
 * rebuilding the second space from the first.
 *
 * Since two spaces cannot exist simultaneously the old space is actually
 * required in a decomposed state. These are the old cells sizes and counts
 * per dimension, along with a list of the old nodeIDs. The old nodeIDs are
 * indexed by the cellid (see cell_getid()), so should be stored that way.
 *
 * On exit the new space cells will have their nodeIDs assigned.
 *
 * @param oldh the cell dimensions of old space.
 * @param oldcdim number of cells per dimension in old space.
 * @param oldnodeIDs the nodeIDs of cells in the old space, indexed by old
 *cellid.
 * @param s the space to be partitioned.
 *
 * @return 1 if the new space contains nodeIDs from all nodes, 0 otherwise.
 */
int partition_space_to_space(double *oldh, double *oldcdim, int *oldnodeIDs,
                             struct space *s) {

  /* Loop over all the new cells. */
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {

        /* Scale indices to old cell space. */
        const int ii = rint(i * s->iwidth[0] * oldh[0]);
        const int jj = rint(j * s->iwidth[1] * oldh[1]);
        const int kk = rint(k * s->iwidth[2] * oldh[2]);

        const int cid = cell_getid(s->cdim, i, j, k);
        const int oldcid = cell_getid(oldcdim, ii, jj, kk);
        s->cells_top[cid].nodeID = oldnodeIDs[oldcid];
      }
    }
  }

  /* Check we have all nodeIDs present in the resample. */
  return check_complete(s, 1, s->e->nr_nodes);
}

/**
 * @brief save the nodeIDs of the current top-level cells by adding them to a
 *             repartition struct. Used when restarting application.
 *
 * @param s the space with the top-level cells.
 * @param reparttype struct to update with the a list of nodeIDs.
 *
 */
void partition_store_celllist(struct space *s, struct repartition *reparttype) {
  if (reparttype->ncelllist != s->nr_cells) {
    free(reparttype->celllist);
    if ((reparttype->celllist = (int *)malloc(sizeof(int) * s->nr_cells)) ==
        NULL)
      error("Failed to allocate celllist");
    reparttype->ncelllist = s->nr_cells;
  }

  for (int i = 0; i < s->nr_cells; i++) {
    reparttype->celllist[i] = s->cells_top[i].nodeID;
  }
}

/**
 * @brief restore the saved list of nodeIDs by applying them to the
 *        top-level cells of a space. Used when restarting application.
 *
 * @param s the space with the top-level cells.
 * @param reparttype struct with the list of nodeIDs saved,
 *
 */
void partition_restore_celllist(struct space *s,
                                struct repartition *reparttype) {
  if (reparttype->ncelllist > 0) {
    if (reparttype->ncelllist == s->nr_cells) {
      for (int i = 0; i < s->nr_cells; i++) {
        s->cells_top[i].nodeID = reparttype->celllist[i];
      }
      if (!check_complete(s, 1, s->e->nr_nodes)) {
        error("Not all ranks are present in the restored partition");
      }
    } else {
      error(
          "Cannot apply the saved partition celllist as the "
          "number of top-level cells (%d) is different to the "
          "saved number (%d)",
          s->nr_cells, reparttype->ncelllist);
    }
  }
}

/**
 * @brief Write a repartition struct to the given FILE as a stream of bytes.
 *
 * @param reparttype the struct
 * @param stream the file stream
 */
void partition_struct_dump(struct repartition *reparttype, FILE *stream) {
  restart_write_blocks(reparttype, sizeof(struct repartition), 1, stream,
                       "repartition", "repartition params");

  /* Also save the celllist, if we have one. */
  if (reparttype->ncelllist > 0)
    restart_write_blocks(reparttype->celllist,
                         sizeof(int) * reparttype->ncelllist, 1, stream,
                         "celllist", "repartition celllist");
}

/**
 * @brief Restore a repartition struct from the given FILE as a stream of
 * bytes.
 *
 * @param reparttype the struct
 * @param stream the file stream
 */
void partition_struct_restore(struct repartition *reparttype, FILE *stream) {
  restart_read_blocks(reparttype, sizeof(struct repartition), 1, stream, NULL,
                      "repartition params");

  /* Also restore the celllist, if we have one. */
  if (reparttype->ncelllist > 0) {
    if ((reparttype->celllist =
             (int *)malloc(sizeof(int) * reparttype->ncelllist)) == NULL)
      error("Failed to allocate celllist");
    restart_read_blocks(reparttype->celllist,
                        sizeof(int) * reparttype->ncelllist, 1, stream, NULL,
                        "repartition celllist");
  }
}
