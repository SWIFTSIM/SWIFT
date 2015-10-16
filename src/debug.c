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

/**
 * @brief Looks for the particle with the given id and prints its information to
 *the standard output.
 *
 * @param parts The array of particles.
 * @param id The id too look for.
 * @param N The size of the array of particles.
 *
 * (Should be used for debugging only as it runs in O(N).)
 */

void printParticle(struct part *parts, long long int id, int N) {

  int i, found = 0;

  /* Look for the particle. */
  for (i = 0; i < N; i++)
    if (parts[i].id == id) {
      printf(
          "## Particle[%d]: id=%lld, x=[%.16e,%.16e,%.16e], "
          "v=[%.3e,%.3e,%.3e], a=[%.3e,%.3e,%.3e], h=%.3e, h_dt=%.3e, "
          "wcount=%.3e, m=%.3e, rho=%.3e, rho_dh=%.3e, div_v=%.3e, u=%.3e, "
          "dudt=%.3e, bals=%.3e, POrho2=%.3e, v_sig=%.3e, dt=%.3e\n",
          i, parts[i].id, parts[i].x[0], parts[i].x[1], parts[i].x[2],
          parts[i].v[0], parts[i].v[1], parts[i].v[2], parts[i].a[0],
          parts[i].a[1], parts[i].a[2], parts[i].h, parts[i].force.h_dt,
          parts[i].density.wcount, parts[i].mass, parts[i].rho, parts[i].rho_dh,
          parts[i].density.div_v, parts[i].u, parts[i].force.u_dt,
          parts[i].force.balsara, parts[i].force.POrho2, parts[i].force.v_sig,
          parts[i].dt);
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
          "v=[%.3e,%.3e,%.3e], a=[%.3e,%.3e,%.3e], m=%.3e, dt=%.3e\n",
          i, (parts[i].id < 0) ? -parts[i].id : parts[i].part->id,
          parts[i].x[0], parts[i].x[1], parts[i].x[2], parts[i].v[0],
          parts[i].v[1], parts[i].v[2], parts[i].a[0], parts[i].a[1],
          parts[i].a[2], parts[i].mass, parts[i].dt);
      found = 1;
    }

  if (!found) printf("## Particles[???] id=%lld not found\n", id);
}

/**
 * @brief Prints the details of a given particle to stdout
 *
 * @param p The particle to print
 *
 */

void printParticle_single(struct part *p) {

  printf(
      "## Particle: id=%lld, x=[%e,%e,%e], v=[%.3e,%.3e,%.3e], "
      "a=[%.3e,%.3e,%.3e], h=%.3e, h_dt=%.3e, wcount=%.3e, m=%.3e, rho=%.3e, "
      "rho_dh=%.3e, div_v=%.3e, u=%.3e, dudt=%.3e, bals=%.3e, POrho2=%.3e, "
      "v_sig=%.3e, dt=%.3e\n",
      p->id, p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2], p->a[0],
      p->a[1], p->a[2], p->h, p->force.h_dt, p->density.wcount, p->mass, p->rho,
      p->rho_dh, p->density.div_v, p->u, p->force.u_dt, p->force.balsara,
      p->force.POrho2, p->force.v_sig, p->dt);
}

#ifdef HAVE_METIS

/**
 * @brief Dump the METIS graph in standard format to stdout
 *
 * @param nvtxs the number of vertices
 * @param ncon the number vertex weights
 * @param xadj first part of adjacency info
 * @param adjncy second part of adjacency info
 * @param vwgt weights of vertices
 * @param vsize size of vertices
 * @param adjwgt weights of edges
 */
void printMETISGraph(idx_t nvtxs, idx_t ncon, idx_t *xadj, idx_t *adjncy,
                     idx_t *vwgt, idx_t *vsize, idx_t *adjwgt) {
  idx_t i;
  idx_t j;
  int hasvwgt = 0;
  int hasewgt = 0;
  int hasvsize = 0;

  /* Check for vwgt, vsize and adjwgt values. */
  if (vwgt) {
    for (i = 0; i < nvtxs * ncon; i++) {
      if (vwgt[i] != 1) {
        hasvwgt = 1;
        break;
      }
    }
  }
  if (vsize) {
    for (i = 0; i < nvtxs; i++) {
      if (vsize[i] != 1) {
        hasvsize = 1;
        break;
      }
    }
  }
  if (adjwgt) {
    for (i = 0; i < xadj[nvtxs]; i++) {
      if (adjwgt[i] != 1) {
        hasewgt = 1;
        break;
      }
    }
  }

  /*  Write the header line. */
  printf("METIS: ");
  printf("%" PRIDX " %" PRIDX, nvtxs, xadj[nvtxs] / 2);
  if (hasvwgt || hasvsize || hasewgt) {
    printf(" %d%d%d", hasvsize, hasvwgt, hasewgt);
    if (hasvwgt) {
      printf(" %d", (int)ncon);
    }
  }

  /*  Write the rest of the graph. */
  for (i = 0; i < nvtxs; i++) {
    printf("\nMETIS: ");

    if (hasvsize) {
      printf(" %" PRIDX, vsize[i]);
    }

    if (hasvwgt) {
      for (j = 0; j < ncon; j++) {
        printf(" %" PRIDX, vwgt[i * ncon + j]);
      }
    }

    for (j = xadj[i]; j < xadj[i + 1]; j++) {
      printf(" %" PRIDX, adjncy[j] + 1);
      if (hasewgt) {
        printf(" %" PRIDX, adjwgt[j]);
      }
    }
  }
  printf("\n");
}

#endif
