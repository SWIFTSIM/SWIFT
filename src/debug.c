/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013- 2015:
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk),
 *                    Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Peter W. Draper (p.w.draper@durham.ac.uk).
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

/* This object's header. */
#include "debug.h"

/* Some standard headers. */
#include <float.h>
#include <stdio.h>
#include <unistd.h>

/* Local includes. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "hydro.h"
#include "inline.h"
#include "part.h"
#include "space.h"

/* Import the right hydro definition */
#if defined(MINIMAL_SPH)
#include "./hydro/Minimal/hydro_debug.h"
#elif defined(GADGET2_SPH)
#include "./hydro/Gadget2/hydro_debug.h"
#elif defined(HOPKINS_PE_SPH)
#include "./hydro/PressureEntropy/hydro_debug.h"
#elif defined(HOPKINS_PU_SPH)
#include "./hydro/PressureEnergy/hydro_debug.h"
#elif defined(HOPKINS_PU_SPH_MONAGHAN)
#include "./hydro/PressureEnergyMorrisMonaghanAV/hydro_debug.h"
#elif defined(DEFAULT_SPH)
#include "./hydro/Default/hydro_debug.h"
#elif defined(GIZMO_MFV_SPH)
#include "./hydro/GizmoMFV/hydro_debug.h"
#elif defined(GIZMO_MFM_SPH)
#include "./hydro/GizmoMFM/hydro_debug.h"
#elif defined(SHADOWFAX_SPH)
#include "./hydro/Shadowswift/hydro_debug.h"
#elif defined(PLANETARY_SPH)
#include "./hydro/Planetary/hydro_debug.h"
#elif defined(ANARCHY_PU_SPH)
#include "./hydro/AnarchyPU/hydro_debug.h"
#else
#error "Invalid choice of SPH variant"
#endif

/* Import the right gravity definition */
#if defined(DEFAULT_GRAVITY)
#include "./gravity/Default/gravity_debug.h"
#elif defined(POTENTIAL_GRAVITY)
#include "./gravity/Potential/gravity_debug.h"
#else
#error "Invalid choice of gravity variant"
#endif

/**
 * @brief Looks for the particle with the given id and prints its information to
 *the standard output.
 *
 * @param parts The array of particles.
 * @param xparts The array of particle extended data.
 * @param id The id too look for.
 * @param N The size of the array of particles.
 *
 * (Should be used for debugging only as it runs in O(N).)
 */
void printParticle(const struct part *parts, const struct xpart *xparts,
                   long long int id, size_t N) {

  int found = 0;

  /* Look for the particle. */
  for (size_t i = 0; i < N; i++)
    if (parts[i].id == id) {
      printf("## Particle[%zu]:\n id=%lld ", i, parts[i].id);
      hydro_debug_particle(&parts[i], &xparts[i]);
      found = 1;
      break;
    }

  if (!found) printf("## Particles[???] id=%lld not found\n", id);
}

/**
 * @brief Looks for the g-particle with the given id and prints its information
 * to
 * the standard output.
 *
 * @param gparts The array of g-particles.
 * @param parts The array of particles.
 * @param id The id too look for.
 * @param N The size of the array of g-particles.
 *
 * (Should be used for debugging only as it runs in O(N).)
 */
void printgParticle(const struct gpart *gparts, const struct part *parts,
                    long long int id, size_t N) {

  int found = 0;

  /* Look for the particle. */
  for (size_t i = 0; i < N; i++)
    if (gparts[i].id_or_neg_offset == id) {
      printf("## gParticle[%zu] (DM) :\n id=%lld", i, id);
      gravity_debug_particle(&gparts[i]);
      found = 1;
      break;
    } else if (gparts[i].id_or_neg_offset < 0 &&
               parts[-gparts[i].id_or_neg_offset].id == id) {
      printf("## gParticle[%zu] (hydro) :\n id=%lld", i, id);
      gravity_debug_particle(&gparts[i]);
      found = 1;
      break;
    }

  if (!found) printf("## Particles[???] id=%lld not found\n", id);
}

/**
 * @brief Prints the details of a given particle to stdout
 *
 * @param p The particle to print
 * @param xp The extended data ot the particle to print
 */
void printParticle_single(const struct part *p, const struct xpart *xp) {

  printf("## Particle: id=%lld ", p->id);
  hydro_debug_particle(p, xp);
  printf("\n");
}

/**
 * @brief Prints the details of a given particle to stdout
 *
 * @param gp The g-particle to print
 */
void printgParticle_single(struct gpart *gp) {

  printf("## g-Particle: id=%lld ", gp->id_or_neg_offset);
  gravity_debug_particle(gp);
  printf("\n");
}

/**
 * @brief Check that the cells and particles of a space have consistent h_max
 *        values.
 *
 * @param s the space.
 * @result 1 or 0
 */
int checkSpacehmax(struct space *s) {

  /* Loop over local cells. */
  float cell_h_max = 0.0f;
  for (int k = 0; k < s->nr_cells; k++) {
    if (s->cells_top[k].nodeID == s->e->nodeID &&
        s->cells_top[k].hydro.h_max > cell_h_max) {
      cell_h_max = s->cells_top[k].hydro.h_max;
    }
  }

  float cell_stars_h_max = 0.0f;
  for (int k = 0; k < s->nr_cells; k++) {
    if (s->cells_top[k].nodeID == s->e->nodeID &&
        s->cells_top[k].stars.h_max > cell_stars_h_max) {
      cell_stars_h_max = s->cells_top[k].stars.h_max;
    }
  }

  /* Now all particles. */
  float part_h_max = 0.0f;
  for (size_t k = 0; k < s->nr_parts; k++) {
    if (s->parts[k].h > part_h_max) {
      part_h_max = s->parts[k].h;
    }
  }

  /* Now all the sparticles. */
  float spart_h_max = 0.0f;
  for (size_t k = 0; k < s->nr_sparts; k++) {
    if (s->sparts[k].h > spart_h_max) {
      spart_h_max = s->sparts[k].h;
    }
  }

  /*  If within some epsilon we are OK. */
  if (fabsf(cell_h_max - part_h_max) <= FLT_EPSILON &&
      fabsf(cell_stars_h_max - spart_h_max) <= FLT_EPSILON)
    return 1;

  /* There is a problem. Hunt it down. */
  /* part */
  for (int k = 0; k < s->nr_cells; k++) {
    if (s->cells_top[k].nodeID == s->e->nodeID) {
      if (s->cells_top[k].hydro.h_max > part_h_max) {
        message("cell %d is inconsistent (%f > %f)", k,
                s->cells_top[k].hydro.h_max, part_h_max);
      }
    }
  }

  for (size_t k = 0; k < s->nr_parts; k++) {
    if (s->parts[k].h > cell_h_max) {
      message("part %lld is inconsistent (%f > %f)", s->parts[k].id,
              s->parts[k].h, cell_h_max);
    }
  }

  /* spart */
  for (int k = 0; k < s->nr_cells; k++) {
    if (s->cells_top[k].nodeID == s->e->nodeID) {
      if (s->cells_top[k].stars.h_max > spart_h_max) {
        message("cell %d is inconsistent (%f > %f)", k,
                s->cells_top[k].stars.h_max, spart_h_max);
      }
    }
  }

  for (size_t k = 0; k < s->nr_sparts; k++) {
    if (s->sparts[k].h > cell_stars_h_max) {
      message("spart %lld is inconsistent (%f > %f)", s->sparts[k].id,
              s->sparts[k].h, cell_stars_h_max);
    }
  }

  return 0;
}

/**
 * @brief Check if the h_max and dx_max values of a cell's hierarchy are
 * consistent with the particles. Also checks if particles are correctly
 * in a cell. Report verbosely if not.
 *
 * @param c the top cell of the hierarchy.
 * @param depth the recursion depth for use in messages. Set to 0 initially.
 * @result 1 or 0
 */
int checkCellhdxmax(const struct cell *c, int *depth) {

  *depth = *depth + 1;

  float h_max = 0.0f;
  float dx_max = 0.0f;
  float stars_h_max = 0.0f;
  float stars_dx_max = 0.0f;
  int result = 1;

  const double loc_min[3] = {c->loc[0], c->loc[1], c->loc[2]};
  const double loc_max[3] = {c->loc[0] + c->width[0], c->loc[1] + c->width[1],
                             c->loc[2] + c->width[2]};

  const size_t nr_parts = c->hydro.count;
  struct part *parts = c->hydro.parts;
  struct xpart *xparts = c->hydro.xparts;
  for (size_t k = 0; k < nr_parts; k++) {

    struct part *const p = &parts[k];
    struct xpart *const xp = &xparts[k];

    if (p->x[0] < loc_min[0] || p->x[0] >= loc_max[0] || p->x[1] < loc_min[1] ||
        p->x[1] >= loc_max[1] || p->x[2] < loc_min[2] ||
        p->x[2] >= loc_max[2]) {

      message(
          "Inconsistent part position p->x=[%e %e %e], c->loc=[%e %e %e] "
          "c->width=[%e %e %e]",
          p->x[0], p->x[1], p->x[2], c->loc[0], c->loc[1], c->loc[2],
          c->width[0], c->width[1], c->width[2]);

      result = 0;
    }

    const float dx2 = xp->x_diff[0] * xp->x_diff[0] +
                      xp->x_diff[1] * xp->x_diff[1] +
                      xp->x_diff[2] * xp->x_diff[2];

    h_max = max(h_max, p->h);
    dx_max = max(dx_max, sqrt(dx2));
  }

  const size_t nr_sparts = c->stars.count;
  struct spart *sparts = c->stars.parts;
  for (size_t k = 0; k < nr_sparts; k++) {

    struct spart *const sp = &sparts[k];

    if (sp->x[0] < loc_min[0] || sp->x[0] >= loc_max[0] ||
        sp->x[1] < loc_min[1] || sp->x[1] >= loc_max[1] ||
        sp->x[2] < loc_min[2] || sp->x[2] >= loc_max[2]) {

      message(
          "Inconsistent part position p->x=[%e %e %e], c->loc=[%e %e %e] "
          "c->width=[%e %e %e]",
          sp->x[0], sp->x[1], sp->x[2], c->loc[0], c->loc[1], c->loc[2],
          c->width[0], c->width[1], c->width[2]);

      result = 0;
    }

    const float dx2 = sp->x_diff[0] * sp->x_diff[0] +
                      sp->x_diff[1] * sp->x_diff[1] +
                      sp->x_diff[2] * sp->x_diff[2];

    stars_h_max = max(stars_h_max, sp->h);
    stars_dx_max = max(stars_dx_max, sqrt(dx2));
  }

  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        checkCellhdxmax(cp, depth);
      }
    }
  }

  /* Check. */
  if (c->hydro.h_max != h_max) {
    message("%d Inconsistent h_max: cell %f != parts %f", *depth,
            c->hydro.h_max, h_max);
    message("location: %f %f %f", c->loc[0], c->loc[1], c->loc[2]);
    result = 0;
  }
  if (c->hydro.dx_max_part != dx_max) {
    message("%d Inconsistent dx_max: %f != %f", *depth, c->hydro.dx_max_part,
            dx_max);
    message("location: %f %f %f", c->loc[0], c->loc[1], c->loc[2]);
    result = 0;
  }

  if (c->stars.h_max != stars_h_max) {
    message("%d Inconsistent stars_h_max: cell %f != parts %f", *depth,
            c->stars.h_max, stars_h_max);
    message("location: %f %f %f", c->loc[0], c->loc[1], c->loc[2]);
    result = 0;
  }
  if (c->stars.dx_max_part != stars_dx_max) {
    message("%d Inconsistent stars_dx_max: %f != %f", *depth,
            c->stars.dx_max_part, stars_dx_max);
    message("location: %f %f %f", c->loc[0], c->loc[1], c->loc[2]);
    result = 0;
  }

  return result;
}

/**
 * @brief map function for dumping cells.
 */
static void dumpCells_map(struct cell *c, void *data) {
  size_t *ldata = (size_t *)data;
  FILE *file = (FILE *)ldata[0];
  struct engine *e = (struct engine *)ldata[1];
  float ntasks = c->nr_tasks;
  int super = (int)ldata[2];
  int active = (int)ldata[3];
  int mpiactive = (int)ldata[4];
  int pactive = (int)ldata[5];

#if SWIFT_DEBUG_CHECKS
  /* The c->nr_tasks field does not include all the tasks. So let's check this
   * the hard way. Note pairs share the task 50/50 with the other cell. */
  ntasks = 0.0f;
  struct task *tasks = e->sched.tasks;
  int nr_tasks = e->sched.nr_tasks;
  for (int k = 0; k < nr_tasks; k++) {
    if (tasks[k].cj == NULL) {
      if (c == tasks[k].ci) {
        ntasks = ntasks + 1.0f;
      }
    } else {
      if (c == tasks[k].ci || c == tasks[k].cj) {
        ntasks = ntasks + 0.5f;
      }
    }
  }
#endif

  /* Only cells with particles are dumped. */
  if (c->hydro.count > 0 || c->grav.count > 0 || c->stars.count > 0) {

    /* In MPI mode we may only output cells with foreign partners.
     * These define the edges of the partitions. */
    int ismpiactive = 0;
#if WITH_MPI
    ismpiactive = (c->mpi.hydro.send_xv != NULL);
    if (mpiactive)
      mpiactive = ismpiactive;
    else
      mpiactive = 1;
#else
    mpiactive = 1;
#endif

    /* Active cells, otherwise all. */
    if (active)
      active = cell_is_active_hydro(c, e);
    else
      active = 1;

    /* So output local super cells or top-level cells that are active and have
     * MPI
     * tasks as requested. */
    if (c->nodeID == e->nodeID &&
        (!super || ((super && c->super == c) || (c->parent == NULL))) &&
        active && mpiactive) {

      /* If requested we work out how many particles are active in this cell. */
      int pactcount = 0;
      if (pactive) {
        const struct part *parts = c->hydro.parts;
        for (int k = 0; k < c->hydro.count; k++)
          if (part_is_active(&parts[k], e)) pactcount++;
        struct gpart *gparts = c->grav.parts;
        for (int k = 0; k < c->grav.count; k++)
          if (gpart_is_active(&gparts[k], e)) pactcount++;
        struct spart *sparts = c->stars.parts;
        for (int k = 0; k < c->stars.count; k++)
          if (spart_is_active(&sparts[k], e)) pactcount++;
      }

      fprintf(file,
              "  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6d %6d %6d %6d %6d %6d "
              "%6.1f %20lld %6d %6d %6d %6d %6d %6d %6d\n",
              c->loc[0], c->loc[1], c->loc[2], c->width[0], c->width[1],
              c->width[2], e->step, c->hydro.count, c->grav.count,
              c->stars.count, pactcount, c->depth, ntasks, c->hydro.ti_end_min,
              get_time_bin(c->hydro.ti_end_min), (c->super == c),
              (c->parent == NULL), cell_is_active_hydro(c, e), c->nodeID,
              c->nodeID == e->nodeID, ismpiactive);
    }
  }
}

/**
 * @brief Dump the location, depth, task counts and timebins and active state,
 * for all cells to a simple text file. A more costly count of the active
 * particles in a cell can also be output.
 *
 * @param prefix base output filename, result is written to
 *               %prefix%_%rank%_%step%.dat
 * @param super just output the super cells.
 * @param active just output active cells.
 * @param mpiactive just output MPI active cells, i.e. those with foreign cells.
 * @param pactive also output a count of active particles.
 * @param s the space holding the cells to dump.
 * @param rank node ID of MPI rank, or 0 if not relevant.
 * @param step the current engine step, or some unique integer.
 */
void dumpCells(const char *prefix, int super, int active, int mpiactive,
               int pactive, struct space *s, int rank, int step) {

  FILE *file = NULL;

  /* Name of output file. */
  char fname[200];
  sprintf(fname, "%s_%03d_%03d.dat", prefix, rank, step);
  file = fopen(fname, "w");

  /* Header. */
  fprintf(file,
          "# %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s "
          "%20s %6s %6s %6s %6s %6s %6s %6s\n",
          "x", "y", "z", "xw", "yw", "zw", "step", "count", "gcount", "scount",
          "actcount", "depth", "tasks", "ti_end_min", "timebin", "issuper",
          "istop", "active", "rank", "local", "mpiactive");

  size_t data[6];
  data[0] = (size_t)file;
  data[1] = (size_t)s->e;
  data[2] = (size_t)super;
  data[3] = (size_t)active;
  data[4] = (size_t)mpiactive;
  data[5] = (size_t)pactive;
  space_map_cells_pre(s, 1, dumpCells_map, &data);
  fclose(file);
}

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

/**
 * @brief Dump a graph in METIS standard format, simple format and weights
 * only, to a file.
 *
 * The standard format output can be read into the METIS and some ParMETIS
 * command-line tools. The simple format is just the cell connectivity (this
 * should not change between calls).  The weights format is the standard one,
 * minus the cell connectivity.
 *
 * The output filenames are generated from the prefix and the sequence number
 * of calls. So the first is called {prefix}_std_001.dat,
 *{prefix}_simple_001.dat,
 * {prefix}_weights_001.dat, etc.
 *
 * @param prefix base output filename
 * @param nvertices the number of vertices
 * @param nvertexweights the number vertex weights
 * @param cellconruns first part of cell connectivity info (CSR)
 * @param cellcon second part of cell connectivity info (CSR)
 * @param vertexweights weights of vertices
 * @param vertexsizes size of vertices
 * @param edgeweights weights of edges
 */
void dumpMETISGraph(const char *prefix, idx_t nvertices, idx_t nvertexweights,
                    idx_t *cellconruns, idx_t *cellcon, idx_t *vertexweights,
                    idx_t *vertexsizes, idx_t *edgeweights) {
  FILE *stdfile = NULL;
  FILE *simplefile = NULL;
  FILE *weightfile = NULL;
  char fname[200];
  int haveedgeweight = 0;
  int havevertexsize = 0;
  int havevertexweight = 0;
  static int nseq = 0;
  nseq++;

  if (vertexweights != NULL) {
    for (idx_t i = 0; i < nvertices * nvertexweights; i++) {
      if (vertexweights[i] != 1) {
        havevertexweight = 1;
        break;
      }
    }
  }

  if (vertexsizes != NULL) {
    for (idx_t i = 0; i < nvertices; i++) {
      if (vertexsizes[i] != 1) {
        havevertexsize = 1;
        break;
      }
    }
  }

  if (edgeweights != NULL) {
    for (idx_t i = 0; i < cellconruns[nvertices]; i++) {
      if (edgeweights[i] != 1) {
        haveedgeweight = 1;
        break;
      }
    }
  }

  /*  Open output files. */
  sprintf(fname, "%s_std_%03d.dat", prefix, nseq);
  stdfile = fopen(fname, "w");

  sprintf(fname, "%s_simple_%03d.dat", prefix, nseq);
  simplefile = fopen(fname, "w");

  if (havevertexweight || havevertexsize || haveedgeweight) {
    sprintf(fname, "%s_weights_%03d.dat", prefix, nseq);
    weightfile = fopen(fname, "w");
  }

  /*  Write the header lines. */
  fprintf(stdfile, "%" PRIDX " %" PRIDX, nvertices, cellconruns[nvertices] / 2);
  fprintf(simplefile, "%" PRIDX " %" PRIDX, nvertices,
          cellconruns[nvertices] / 2);
  if (havevertexweight || havevertexsize || haveedgeweight) {
    fprintf(weightfile, "%" PRIDX " %" PRIDX, nvertices,
            cellconruns[nvertices] / 2);

    fprintf(stdfile, " %d%d%d", havevertexsize, havevertexweight,
            haveedgeweight);
    fprintf(weightfile, " %d%d%d", havevertexsize, havevertexweight,
            haveedgeweight);

    if (havevertexweight) {
      fprintf(stdfile, " %d", (int)nvertexweights);
      fprintf(weightfile, " %d", (int)nvertexweights);
    }
  }

  /*  Write the rest of the graph. */
  for (idx_t i = 0; i < nvertices; i++) {
    fprintf(stdfile, "\n");
    fprintf(simplefile, "\n");
    if (weightfile != NULL) {
      fprintf(weightfile, "\n");
    }

    if (havevertexsize) {
      fprintf(stdfile, " %" PRIDX, vertexsizes[i]);
      fprintf(weightfile, " %" PRIDX, vertexsizes[i]);
    }

    if (havevertexweight) {
      for (idx_t j = 0; j < nvertexweights; j++) {
        fprintf(stdfile, " %" PRIDX, vertexweights[i * nvertexweights + j]);
        fprintf(weightfile, " %" PRIDX, vertexweights[i * nvertexweights + j]);
      }
    }

    for (idx_t j = cellconruns[i]; j < cellconruns[i + 1]; j++) {
      fprintf(stdfile, " %" PRIDX, cellcon[j] + 1);
      fprintf(simplefile, " %" PRIDX, cellcon[j] + 1);
      if (haveedgeweight) {
        fprintf(stdfile, " %" PRIDX, edgeweights[j]);
        fprintf(weightfile, " %" PRIDX, edgeweights[j]);
      }
    }
  }
  fprintf(stdfile, "\n");
  fprintf(simplefile, "\n");
  if (weightfile != NULL) {
    fprintf(weightfile, "\n");
  }

  fclose(stdfile);
  fclose(simplefile);
  if (weightfile != NULL) {
    fclose(weightfile);
  }
}

#endif /* HAVE_METIS || HAVE_PARMETIS */

#ifdef HAVE_MPI
/**
 * @brief Dump the positions and MPI ranks of the given top-level cells
 *        to a simple text file.
 *
 * Can be used to visualise the partitioning of an MPI run. Note should
 * be used immediately after repartitioning when the top-level cells
 * have been assigned their nodes. Each time this is called a new file
 * with the given prefix, a unique integer and type of .dat is created.
 *
 * @param prefix base output filename
 * @param cells_top the top-level cells.
 * @param nr_cells the number of cells.
 */
void dumpCellRanks(const char *prefix, struct cell *cells_top, int nr_cells) {

  FILE *file = NULL;

  /* Name of output file. */
  static int nseq = 0;
  char fname[200];
  sprintf(fname, "%s_%03d.dat", prefix, nseq);
  nseq++;

  file = fopen(fname, "w");

  /* Header. */
  fprintf(file, "# %6s %6s %6s %6s %6s %6s %6s\n", "x", "y", "z", "xw", "yw",
          "zw", "rank");

  /* Output */
  for (int i = 0; i < nr_cells; i++) {
    struct cell *c = &cells_top[i];
    fprintf(file, "  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6d\n", c->loc[0],
            c->loc[1], c->loc[2], c->width[0], c->width[1], c->width[2],
            c->nodeID);
  }

  fclose(file);
}

#endif /* HAVE_MPI */
