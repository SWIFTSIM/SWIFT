/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 SWIFT developers.
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
#include "cell_boundary.h"
#include "swift.h"

#include <assert.h>
#include <string.h>

/* Reflects a small two-level cell in a wall and checks the mirror cell. */
static void test(const unsigned int wall, const double loc_mirror) {

  struct space s;
  struct cell *cells_sub = NULL;
  bzero(&s, sizeof(s));
  s.cells_sub = &cells_sub;
  s.dim[0] = s.dim[1] = s.dim[2] = 1.;

  /* A unit cell filling the box, with two particles in its lower-x half and
   * one in its upper-x half. */
  struct part parts[3], mirrors[3];
  bzero(parts, sizeof(parts));
  const double x[3][3] = {{0.1, 0.2, 0.3}, {0.2, 0.7, 0.6}, {0.9, 0.4, 0.8}};
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      parts[i].x[k] = x[i][k];
      parts[i].v[k] = 1. + i + k;
    }
    parts[i].h = 0.1f * (i + 1);
#ifdef REMIX_SPH
    for (int k = 0; k < 3; k++) {
      parts[i].grad_m0[k] = 1.f;
      parts[i].force.B[k] = 1.f;
      for (int j = 0; j < 3; j++) parts[i].dv_norm_kernel[k][j] = 1.f;
    }
#endif
  }

  struct cell c, lower, upper, cm;
  bzero(&c, sizeof(c));
  bzero(&lower, sizeof(lower));
  bzero(&upper, sizeof(upper));
  bzero(&cm, sizeof(cm));
  c.width[0] = c.width[1] = c.width[2] = 1.;
  c.split = 1;
  c.hydro.parts = parts;
  c.hydro.count = 3;
  c.progeny[0] = &lower;
  c.progeny[4] = &upper;
  lower.hydro.parts = parts;
  lower.hydro.count = 2;
  upper.hydro.parts = parts + 2;
  upper.hydro.count = 1;
  lower.hydro.h_max = 0.2f;
  upper.hydro.h_max = 0.3f;
  lower.hydro.dx_max_part = 0.02f;
  upper.hydro.dx_max_part = 0.01f;
  lower.hydro.ti_old_part = 8;
  upper.hydro.ti_old_part = 4;

  /* The mirror cell on the other side of the wall */
  cm.loc[0] = loc_mirror;
  cm.width[0] = cm.width[1] = cm.width[2] = 1.;
  cm.dmin = 1.;
  cm.subtype = cell_subtype_mirror;
  cm.mirror_walls = wall;
  cm.top = &cm;

  space_mirror_clone(&s, &cm, &c, mirrors);
  cell_mirror_refresh(&cm, s.dim, /*moved=*/1);

  /* Same content, with the halves swapped along x only */
  assert(cm.split && cm.hydro.count == 3 && cm.hydro.parts == mirrors);
  for (int k = 0; k < 8; k++)
    assert((cm.progeny[k] != NULL) == (k == 0 || k == 4));
  const struct cell *ml = cm.progeny[0], *mu = cm.progeny[4];
  assert(ml->mirror_of == &upper && mu->mirror_of == &lower);
  assert(ml->hydro.count == 1 && ml->hydro.parts == mirrors);
  assert(mu->hydro.count == 2 && mu->hydro.parts == mirrors + 1);
  assert(ml->subtype == cell_subtype_mirror && ml->mirror_walls == wall);
  assert(ml->depth == 1 && ml->parent == &cm && ml->top == &cm);
  assert(ml->loc[0] == loc_mirror && mu->loc[0] == loc_mirror + 0.5);
  assert(ml->width[0] == 0.5 && !ml->split);

  /* The mirror cells are never active and carry no extended particle data */
  assert(cm.hydro.ti_end_min == max_nr_timesteps && cm.hydro.xparts == NULL);
  assert(ml->hydro.ti_end_min == max_nr_timesteps && ml->hydro.xparts == NULL);

  /* They take the state of the cells they reflect, summarised up the tree */
  assert(ml->hydro.h_max == 0.3f && mu->hydro.h_max == 0.2f);
  assert(cm.hydro.h_max == 0.3f && cm.hydro.dx_max_part == 0.02f);
  assert(ml->hydro.ti_old_part == 4 && cm.hydro.ti_old_part == 4);

  /* Only a refresh that moved the particles invalidates the sorts */
  cm.hydro.sorted = cm.progeny[0]->hydro.sorted = 1;
  cell_mirror_refresh(&cm, s.dim, /*moved=*/0);
  assert(cm.hydro.sorted == 1 && ml->hydro.sorted == 1);
  cell_mirror_refresh(&cm, s.dim, /*moved=*/1);
  assert(cm.hydro.sorted == 0 && ml->hydro.sorted == 0);

  /* Each mirror particle is the reflection of its counterpart */
  const double wall_x = (wall == boundary_wall_x_min) ? 0. : 1.;
  const int parent[3] = {2, 0, 1};
  for (int i = 0; i < 3; i++) {
    const struct part *p = &parts[parent[i]], *m = &mirrors[i];
    assert(fabs(m->x[0] - (2. * wall_x - p->x[0])) < 1e-15);
    assert(m->x[1] == p->x[1] && m->x[2] == p->x[2]);
    assert(m->v[0] == -p->v[0] && m->v[1] == p->v[1] && m->v[2] == p->v[2]);
    assert(m->h == p->h && m->gpart == NULL);
#ifdef REMIX_SPH
    /* Vectors flip along x; tensors where one index only is x */
    assert(m->grad_m0[0] == -1.f && m->grad_m0[1] == 1.f);
    assert(m->force.B[0] == -1.f && m->force.B[2] == 1.f);
    assert(m->dv_norm_kernel[0][0] == 1.f && m->dv_norm_kernel[1][2] == 1.f);
    assert(m->dv_norm_kernel[0][1] == -1.f && m->dv_norm_kernel[2][0] == -1.f);
#endif
  }

  /* And lies inside the mirror cell holding it */
  assert(mirrors[0].x[0] >= ml->loc[0] && mirrors[0].x[0] < mu->loc[0]);
  assert(mirrors[1].x[0] >= mu->loc[0] && mirrors[2].x[0] >= mu->loc[0]);
}

/* Sum of the kernel weights of all the particles around a position. */
static double kernel_sum(const double x[3], const float h,
                         const struct part *parts, const int count) {
  double sum = 0.;
  for (int i = 0; i < count; i++) {
    double r2 = 0.;
    for (int k = 0; k < 3; k++) {
      r2 += (x[k] - parts[i].x[k]) * (x[k] - parts[i].x[k]);
    }
    const float u = sqrt(r2) / h;
    float w = 0.f;
    if (u < kernel_gamma) kernel_eval(u, &w);
    sum += w;
  }
  return sum;
}

/* Checks that a particle on a lattice next to a mirror wall has the kernel
 * sum of a particle in the interior. */
#define N 8
static void test_kernel_sum(void) {

  struct space s;
  bzero(&s, sizeof(s));
  s.dim[0] = s.dim[1] = s.dim[2] = 1.;

  /* A lattice in a unit cell, half a spacing away from the walls */
  static struct part parts[2 * N * N * N];
  struct part *mirrors = parts + N * N * N;
  bzero(parts, sizeof(parts));
  for (int i = 0; i < N * N * N; i++) {
    parts[i].x[0] = (i / (N * N) + 0.5) / N;
    parts[i].x[1] = ((i / N) % N + 0.5) / N;
    parts[i].x[2] = (i % N + 0.5) / N;
  }

  struct cell c, cm;
  bzero(&c, sizeof(c));
  bzero(&cm, sizeof(cm));
  c.hydro.parts = parts;
  c.hydro.count = N * N * N;
  cm.subtype = cell_subtype_mirror;
  cm.mirror_walls = boundary_wall_x_min;
  space_mirror_clone(&s, &cm, &c, mirrors);
  cell_mirror_refresh(&cm, s.dim, /*moved=*/1);

  /* One particle in the middle and one against the wall, neither within
   * reach of another face */
  const float h = 1.5f / N;
  const double *x_in = parts[(N / 2) * (N * N + N + 1)].x;
  const double *x_wall = parts[(N / 2) * (N + 1)].x;
  const double sum_in = kernel_sum(x_in, h, parts, N * N * N);
  const double sum_cut = kernel_sum(x_wall, h, parts, N * N * N);
  const double sum_wall = kernel_sum(x_wall, h, parts, 2 * N * N * N);
  assert(sum_cut < 0.9 * sum_in);
  assert(fabs(sum_wall - sum_in) < 1e-5 * sum_in);
}

int main(void) {
  test(boundary_wall_x_min, -1.);
  test(boundary_wall_x_max, 1.);
  test_kernel_sum();
  return 0;
}
