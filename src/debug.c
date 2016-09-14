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
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <stdio.h>

/* This object's header. */
#include "debug.h"

/* Local includes. */
#include "cell.h"
#include "const.h"
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
#elif defined(DEFAULT_SPH)
#include "./hydro/Default/hydro_debug.h"
#elif defined(GIZMO_SPH)
#include "./hydro/Gizmo/hydro_debug.h"
#else
#error "Invalid choice of SPH variant"
#endif

#include "./gravity/Default/gravity_debug.h"

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
        s->cells_top[k].h_max > cell_h_max) {
      cell_h_max = s->cells_top[k].h_max;
    }
  }

  /* Now all particles. */
  float part_h_max = 0.0f;
  for (size_t k = 0; k < s->nr_parts; k++) {
    if (s->parts[k].h > part_h_max) {
      part_h_max = s->parts[k].h;
    }
  }

  /*  If within some epsilon we are OK. */
  if (abs(cell_h_max - part_h_max) <= FLT_EPSILON) return 1;

  /* There is a problem. Hunt it down. */
  for (int k = 0; k < s->nr_cells; k++) {
    if (s->cells_top[k].nodeID == s->e->nodeID) {
      if (s->cells_top[k].h_max > part_h_max) {
        message("cell %d is inconsistent (%f > %f)", k, s->cells_top[k].h_max,
                part_h_max);
      }
    }
  }

  for (size_t k = 0; k < s->nr_parts; k++) {
    if (s->parts[k].h > cell_h_max) {
      message("part %lld is inconsistent (%f > %f)", s->parts[k].id,
              s->parts[k].h, cell_h_max);
    }
  }

  return 0;
}

#ifdef HAVE_METIS

/**
 * @brief Dump the METIS graph in standard format, simple format and weights
 * only, to a file.
 *
 * The standard format output can be read into the METIS
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

#endif
