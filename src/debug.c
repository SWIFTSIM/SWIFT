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

#include <stdio.h>

#include "config.h"
#include "const.h"
#include "part.h"
#include "debug.h"

/* Import the right hydro definition */
#if defined(MINIMAL_SPH)
#include "./hydro/Minimal/hydro_debug.h"
#elif defined(GADGET2_SPH)
#include "./hydro/Gadget2/hydro_debug.h"
#elif defined(DEFAULT_SPH)
#include "./hydro/Default/hydro_debug.h"
#elif defined(NO_SPH)
#include "./hydro/Gadget2/hydro_debug.h"
#else
#error "Invalid choice of SPH variant"
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

void printParticle(struct part *parts, struct xpart *xparts, long long int id,
                   int N) {

  int i, found = 0;

  /* Look for the particle. */
  for (i = 0; i < N; i++)
    if (parts[i].id == id) {
      printf("## Particle[%d]:\n id=%lld", i, parts[i].id);
      hydro_debug_particle(&parts[i], &xparts[i]);
      found = 1;
    }

  if (!found) printf("## Particles[???] id=%lld not found\n", id);
}

void printgParticle(struct gpart *parts, long long int id, int N) {

  int i, found = 0;

  /* Look for the particle. */
  for (i = 0; i < N; i++)
    if (parts[i].id == -id || (parts[i].id > 0 && parts[i].part->id == id)) {
      printf(
          "## gParticle[%d]: id=%lld, x=[%.16e,%.16e,%.16e], "
          "v=[%.3e,%.3e,%.3e], a=[%.3e,%.3e,%.3e], m=%.3e, t_begin=%d, "
          "t_end=%d\n",
          i, parts[i].part->id, parts[i].x[0], parts[i].x[1], parts[i].x[2],
          parts[i].v_full[0], parts[i].v_full[1], parts[i].v_full[2], parts[i].a[0],
          parts[i].a[1], parts[i].a[2], parts[i].mass, parts[i].ti_begin,
          parts[i].ti_end);
      found = 1;
    }

  if (!found) printf("## Particles[???] id=%lld not found\n", id);
}

/**
 * @brief Prints the details of a given particle to stdout
 *
 * @param p The particle to print
 * @param xp The extended data ot the particle to print
 *
 */

void printParticle_single(struct part *p, struct xpart *xp) {

  printf("## Particle: id=%lld", p->id);
  hydro_debug_particle(p, xp);
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
  idx_t i;
  idx_t j;
  int haveedgeweight = 0;
  int havevertexsize = 0;
  int havevertexweight = 0;
  static int nseq = 0;
  nseq++;

  if (vertexweights != NULL) {
    for (i = 0; i < nvertices * nvertexweights; i++) {
      if (vertexweights[i] != 1) {
        havevertexweight = 1;
        break;
      }
    }
  }

  if (vertexsizes != NULL) {
    for (i = 0; i < nvertices; i++) {
      if (vertexsizes[i] != 1) {
        havevertexsize = 1;
        break;
      }
    }
  }

  if (edgeweights != NULL) {
    for (i = 0; i < cellconruns[nvertices]; i++) {
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
  for (i = 0; i < nvertices; i++) {
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
      for (j = 0; j < nvertexweights; j++) {
        fprintf(stdfile, " %" PRIDX, vertexweights[i * nvertexweights + j]);
        fprintf(weightfile, " %" PRIDX, vertexweights[i * nvertexweights + j]);
      }
    }

    for (j = cellconruns[i]; j < cellconruns[i + 1]; j++) {
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
