/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* This object's header. */
#include "map.h"

/* Local headers. */
#include "atomic.h"
#include "error.h"

/**
 * @brief Mapping function to draw a specific cell (gnuplot).
 */
void map_cells_plot(struct cell *c, void *data) {

  int depth = *(int *)data;
  double *l = c->loc, *h = c->width;

  if (c->depth <= depth) {

    printf("%.16e %.16e %.16e\n", l[0], l[1], l[2]);
    printf("%.16e %.16e %.16e\n", l[0] + h[0], l[1], l[2]);
    printf("%.16e %.16e %.16e\n", l[0] + h[0], l[1] + h[1], l[2]);
    printf("%.16e %.16e %.16e\n\n\n", l[0], l[1] + h[1], l[2]);

    printf("%.16e %.16e %.16e\n", l[0], l[1], l[2] + h[2]);
    printf("%.16e %.16e %.16e\n", l[0] + h[0], l[1], l[2] + h[2]);
    printf("%.16e %.16e %.16e\n", l[0] + h[0], l[1] + h[1], l[2] + h[2]);
    printf("%.16e %.16e %.16e\n\n\n", l[0], l[1] + h[1], l[2] + h[2]);

    printf("%.16e %.16e %.16e\n", l[0], l[1], l[2]);
    printf("%.16e %.16e %.16e\n", l[0], l[1] + h[1], l[2]);
    printf("%.16e %.16e %.16e\n", l[0], l[1] + h[1], l[2] + h[2]);
    printf("%.16e %.16e %.16e\n\n\n", l[0], l[1], l[2] + h[2]);

    printf("%.16e %.16e %.16e\n", l[0] + h[0], l[1], l[2]);
    printf("%.16e %.16e %.16e\n", l[0] + h[0], l[1] + h[1], l[2]);
    printf("%.16e %.16e %.16e\n", l[0] + h[0], l[1] + h[1], l[2] + h[2]);
    printf("%.16e %.16e %.16e\n\n\n", l[0] + h[0], l[1], l[2] + h[2]);

    printf("%.16e %.16e %.16e\n", l[0], l[1], l[2]);
    printf("%.16e %.16e %.16e\n", l[0], l[1], l[2] + h[2]);
    printf("%.16e %.16e %.16e\n", l[0] + h[0], l[1], l[2] + h[2]);
    printf("%.16e %.16e %.16e\n\n\n", l[0] + h[0], l[1], l[2]);

    printf("%.16e %.16e %.16e\n", l[0], l[1] + h[1], l[2]);
    printf("%.16e %.16e %.16e\n", l[0], l[1] + h[1], l[2] + h[2]);
    printf("%.16e %.16e %.16e\n", l[0] + h[0], l[1] + h[1], l[2] + h[2]);
    printf("%.16e %.16e %.16e\n\n\n", l[0] + h[0], l[1] + h[1], l[2]);

    if (!c->split) {
      for (int k = 0; k < c->hydro.count; k++)
        printf("0 0 0 %.16e %.16e %.16e\n", c->hydro.parts[k].x[0],
               c->hydro.parts[k].x[1], c->hydro.parts[k].x[2]);
      printf("\n\n");
    }
    /* else
        for ( int k = 0 ; k < 8 ; k++ )
            if ( c->progeny[k] != NULL )
                map_cells_plot( c->progeny[k] , data ); */
  }
}

/**
 * @brief Mapping function for checking if each part is in its box.
 */
void map_check(struct part *p, struct cell *c, void *data) {

  if (p->x[0] < c->loc[0] || p->x[0] > c->loc[0] + c->width[0] ||
      p->x[0] < c->loc[0] || p->x[0] > c->loc[0] + c->width[0] ||
      p->x[0] < c->loc[0] || p->x[0] > c->loc[0] + c->width[0])
    printf("map_check: particle %lld is outside of its box.\n", p->id);
}

/**
 * @brief Mapping function for neighbour count.
 */
void map_cellcheck(struct cell *c, void *data) {

  int *count = (int *)data;
  atomic_add(count, c->hydro.count);

  /* Loop over all parts and check if they are in the cell. */
  for (int k = 0; k < c->hydro.count; k++) {
    struct part *p = &c->hydro.parts[k];
    if (p->x[0] < c->loc[0] || p->x[1] < c->loc[1] || p->x[2] < c->loc[2] ||
        p->x[0] > c->loc[0] + c->width[0] ||
        p->x[1] > c->loc[1] + c->width[1] ||
        p->x[2] > c->loc[2] + c->width[2]) {
      printf(
          "map_cellcheck: particle at [ %.16e %.16e %.16e ] outside of cell [ "
          "%.16e %.16e %.16e ] - [ %.16e %.16e %.16e ].\n",
          p->x[0], p->x[1], p->x[2], c->loc[0], c->loc[1], c->loc[2],
          c->loc[0] + c->width[0], c->loc[1] + c->width[1],
          c->loc[2] + c->width[2]);
      error("particle out of bounds!");
    }
  }

  /* Loop over all gparts and check if they are in the cell. */
  for (int k = 0; k < c->grav.count; k++) {
    struct gpart *p = &c->grav.parts[k];
    if (p->x[0] < c->loc[0] || p->x[1] < c->loc[1] || p->x[2] < c->loc[2] ||
        p->x[0] > c->loc[0] + c->width[0] ||
        p->x[1] > c->loc[1] + c->width[1] ||
        p->x[2] > c->loc[2] + c->width[2]) {
      printf(
          "map_cellcheck: g-particle at [ %.16e %.16e %.16e ] outside of cell "
          "[ %.16e %.16e %.16e ] - [ %.16e %.16e %.16e ].\n",
          p->x[0], p->x[1], p->x[2], c->loc[0], c->loc[1], c->loc[2],
          c->loc[0] + c->width[0], c->loc[1] + c->width[1],
          c->loc[2] + c->width[2]);
      error("g-particle out of bounds!");
    }
  }
}

/**
 * @brief Mapping function for maxdepth cell count.
 */
void map_maxdepth(struct cell *c, void *data) {

  int maxdepth = ((int *)data)[0];
  int *count = &((int *)data)[1];

  // printf( "%e\n" , p->count );

  if (c->depth == maxdepth) *count += 1;
}

/**
 * @brief Mapping function for neighbour count.
 */
void map_count(struct part *p, struct cell *c, void *data) {

  double *wcount = (double *)data;

  // printf( "%i %e %e\n" , p->id , p->count , p->count_dh );

  *wcount += p->density.wcount;
}
void map_wcount_min(struct part *p, struct cell *c, void *data) {

  struct part **p2 = (struct part **)data;

  if (p->density.wcount < (*p2)->density.wcount) *p2 = p;
}

void map_wcount_max(struct part *p, struct cell *c, void *data) {

  struct part **p2 = (struct part **)data;

  if (p->density.wcount > (*p2)->density.wcount) *p2 = p;
}

void map_h_min(struct part *p, struct cell *c, void *data) {

  struct part **p2 = (struct part **)data;

  if (p->h < (*p2)->h) *p2 = p;
}

void map_h_max(struct part *p, struct cell *c, void *data) {

  struct part **p2 = (struct part **)data;

  if (p->h > (*p2)->h) *p2 = p;
}

void map_stars_h_max(struct spart *p, struct cell *c, void *data) {

  struct spart **p2 = (struct spart **)data;

  if (p->h > (*p2)->h) *p2 = p;
}

/**
 * @brief Mapping function for neighbour count.
 */
void map_icount(struct part *p, struct cell *c, void *data) {

  // int *count = (int *)data;

  // printf( "%i\n" , p->icount );

  // *count += p->icount;
}

/**
 * @brief Mapping function to print the particle position.
 */
void map_dump(struct part *p, struct cell *c, void *data) {

  double *shift = (double *)data;

  printf("%g\t%g\t%g\n", p->x[0] - shift[0], p->x[1] - shift[1],
         p->x[2] - shift[2]);
}
