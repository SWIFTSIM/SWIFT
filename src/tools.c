/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 * Copyright (c) 2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include <ctype.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>

/* This object's header. */
#include "tools.h"

/* Local includes. */
#include "active.h"
#include "cell.h"
#include "chemistry.h"
#include "cosmology.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "part.h"
#include "periodic.h"
#include "runner.h"
#include "star_formation_iact.h"
#include "stars.h"

/**
 *  Factorize a given integer, attempts to keep larger pair of factors.
 */
void factor(int value, int *f1, int *f2) {
  int j;
  int i;

  j = (int)sqrt(value);
  for (i = j; i > 0; i--) {
    if ((value % i) == 0) {
      *f1 = i;
      *f2 = value / i;
      break;
    }
  }
}

/**
 * @brief Compute the average number of pairs per particle using
 *      a brute-force O(N^2) computation.
 *
 * @param dim The space dimensions.
 * @param parts The #part array.
 * @param N The number of parts.
 * @param periodic Periodic boundary conditions flag.
 */
void pairs_n2(double *dim, struct part *restrict parts, int N, int periodic) {
  int i, j, k, count = 0;
  // int mj, mk;
  // double maxratio = 1.0;
  double r2, dx[3], rho = 0.0;
  double rho_max = 0.0, rho_min = 100;
  float a = 1.f, H = 0.f;

  /* Loop over all particle pairs. */
  for (j = 0; j < N; j++) {
    if (j % 1000 == 0) {
      printf("pairs_n2: j=%i.\n", j);
      fflush(stdout);
    }
    for (k = j + 1; k < N; k++) {
      for (i = 0; i < 3; i++) {
        dx[i] = parts[j].x[i] - parts[k].x[i];
        if (periodic) {
          if (dx[i] < -dim[i] / 2)
            dx[i] += dim[i];
          else if (dx[i] > dim[i] / 2)
            dx[i] -= dim[i];
        }
      }
      r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      if (r2 < parts[j].h * parts[j].h || r2 < parts[k].h * parts[k].h) {
        runner_iact_density(r2, NULL, parts[j].h, parts[k].h, &parts[j],
                            &parts[k], a, H);
        /* if ( parts[j].h / parts[k].h > maxratio )
            {
            maxratio = parts[j].h / parts[k].h;
            mj = j; mk = k;
            }
        else if ( parts[k].h / parts[j].h > maxratio )
            {
            maxratio = parts[k].h / parts[j].h;
            mj = j; mk = k;
            } */
      }
    }
  }

  /* Aggregate the results. */
  for (k = 0; k < N; k++) {
    // count += parts[k].icount;
    rho += parts[k].density.wcount;
    rho_min = fmin(parts[k].density.wcount, rho_min);
    rho_min = fmax(parts[k].density.wcount, rho_max);
  }

  /* Dump the result. */
  printf("pairs_n2: avg. density per part is %.3f (nr. pairs %.3f).\n",
         rho / N + 32.0 / 3, ((double)count) / N);
  printf("pairs_n2: densities are in [ %e , %e ].\n", rho_min / N + 32.0 / 3,
         rho_max / N + 32.0 / 3);
  /* printf( "pairs_n2: maximum ratio between parts %i [%e,%e,%e] and %i
     [%e,%e,%e] is %.3f/%.3f\n" ,
      mj , parts[mj].x[0] , parts[mj].x[1] , parts[mj].x[2] ,
      mk , parts[mk].x[0] , parts[mk].x[1] , parts[mk].x[2] ,
      parts[mj].h , parts[mk].h ); fflush(stdout); */
  fflush(stdout);
}

void pairs_single_density(double *dim, long long int pid,
                          struct part *restrict parts, int N, int periodic) {
  int i, k;
  double r2, dx[3];
  float fdx[3];
  struct part p;
  float a = 1.f, H = 0.f;

  /* Find "our" part. */
  for (k = 0; k < N && parts[k].id != pid; k++)
    ;
  if (k == N) error("Part not found.");
  p = parts[k];
  printf("pairs_single: part[%i].id == %lli.\n", k, pid);

  hydro_init_part(&p, NULL);

  /* Loop over all particle pairs. */
  for (k = 0; k < N; k++) {
    if (parts[k].id == p.id) continue;
    for (i = 0; i < 3; i++) {
      dx[i] = p.x[i] - parts[k].x[i];
      if (periodic) {
        if (dx[i] < -dim[i] / 2)
          dx[i] += dim[i];
        else if (dx[i] > dim[i] / 2)
          dx[i] -= dim[i];
      }
      fdx[i] = dx[i];
    }
    r2 = fdx[0] * fdx[0] + fdx[1] * fdx[1] + fdx[2] * fdx[2];
    if (r2 < p.h * p.h) {
      runner_iact_nonsym_density(r2, fdx, p.h, parts[k].h, &p, &parts[k], a, H);
      /* printf( "pairs_simple: interacting particles %lli [%i,%i,%i] and %lli
         [%i,%i,%i], r=%e.\n" ,
          pid , (int)(p.x[0]*ih) , (int)(p.x[1]*ih) , (int)(p.x[2]*ih) ,
          parts[k].id , (int)(parts[k].x[0]*ih) , (int)(parts[k].x[1]*ih) ,
         (int)(parts[k].x[2]*ih) ,
          sqrtf(r2) ); */
    }
  }

  /* Dump the result. */
  printf("pairs_single: wcount of part %lli (h=%e) is %f.\n", p.id, p.h,
         p.density.wcount + 32.0 / 3);
  fflush(stdout);
}

void pairs_all_density(struct runner *r, struct cell *ci, struct cell *cj) {

  float r2, hi, hj, hig2, hjg2, dx[3];
  struct part *pi, *pj;
  const double dim[3] = {r->e->s->dim[0], r->e->s->dim[1], r->e->s->dim[2]};
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->hydro.count; ++i) {

    pi = &ci->hydro.parts[i];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(pi, e)) continue;

    for (int j = 0; j < cj->hydro.count; ++j) {

      pj = &cj->hydro.parts[j];

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = ci->hydro.parts[i].x[k] - cj->hydro.parts[j].x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 && !part_is_inhibited(pj, e)) {

        /* Interact */
        runner_iact_nonsym_density(r2, dx, hi, pj->h, pi, pj, a, H);
        runner_iact_nonsym_chemistry(r2, dx, hi, pj->h, pi, pj, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hi, pj->h, pi, pj, a, H);
      }
    }
  }

  /* Reverse double-for loop and checks every interaction */
  for (int j = 0; j < cj->hydro.count; ++j) {

    pj = &cj->hydro.parts[j];
    hj = pj->h;
    hjg2 = hj * hj * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(pj, e)) continue;

    for (int i = 0; i < ci->hydro.count; ++i) {

      pi = &ci->hydro.parts[i];

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = cj->hydro.parts[j].x[k] - ci->hydro.parts[i].x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hjg2 && !part_is_inhibited(pi, e)) {

        /* Interact */
        runner_iact_nonsym_density(r2, dx, hj, pi->h, pj, pi, a, H);
        runner_iact_nonsym_chemistry(r2, dx, hj, pi->h, pj, pi, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hj, pi->h, pj, pi, a, H);
      }
    }
  }
}

#ifdef EXTRA_HYDRO_LOOP
void pairs_all_gradient(struct runner *r, struct cell *ci, struct cell *cj) {

  float r2, hi, hj, hig2, hjg2, dx[3];
  struct part *pi, *pj;
  const double dim[3] = {r->e->s->dim[0], r->e->s->dim[1], r->e->s->dim[2]};
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->hydro.count; ++i) {

    pi = &ci->hydro.parts[i];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(pi, e)) continue;

    for (int j = 0; j < cj->hydro.count; ++j) {

      pj = &cj->hydro.parts[j];
      hj = pj->h;
      hjg2 = hj * hj * kernel_gamma2;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = ci->hydro.parts[i].x[k] - cj->hydro.parts[j].x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 && !part_is_inhibited(pj, e)) {

        /* Interact */
        runner_iact_nonsym_gradient(r2, dx, hi, hj, pi, pj, a, H);
      }
    }
  }

  /* Reverse double-for loop and checks every interaction */
  for (int j = 0; j < cj->hydro.count; ++j) {

    pj = &cj->hydro.parts[j];
    hj = pj->h;
    hjg2 = hj * hj * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(pj, e)) continue;

    for (int i = 0; i < ci->hydro.count; ++i) {

      pi = &ci->hydro.parts[i];
      hi = pi->h;
      hig2 = hi * hi * kernel_gamma2;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = cj->hydro.parts[j].x[k] - ci->hydro.parts[i].x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hjg2 && !part_is_inhibited(pi, e)) {

        /* Interact */
        runner_iact_nonsym_gradient(r2, dx, hj, pi->h, pj, pi, a, H);
      }
    }
  }
}
#endif /* EXTRA_HDYRO_LOOP */

void pairs_all_force(struct runner *r, struct cell *ci, struct cell *cj) {

  float r2, hi, hj, hig2, hjg2, dx[3];
  struct part *pi, *pj;
  const double dim[3] = {r->e->s->dim[0], r->e->s->dim[1], r->e->s->dim[2]};
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->hydro.count; ++i) {

    pi = &ci->hydro.parts[i];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(pi, e)) continue;

    for (int j = 0; j < cj->hydro.count; ++j) {

      pj = &cj->hydro.parts[j];
      hj = pj->h;
      hjg2 = hj * hj * kernel_gamma2;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = ci->hydro.parts[i].x[k] - cj->hydro.parts[j].x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 || r2 < hjg2) {

        /* Interact */
        runner_iact_nonsym_force(r2, dx, hi, hj, pi, pj, a, H);
      }
    }
  }

  /* Reverse double-for loop and checks every interaction */
  for (int j = 0; j < cj->hydro.count; ++j) {

    pj = &cj->hydro.parts[j];
    hj = pj->h;
    hjg2 = hj * hj * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(pj, e)) continue;

    for (int i = 0; i < ci->hydro.count; ++i) {

      pi = &ci->hydro.parts[i];
      hi = pi->h;
      hig2 = hi * hi * kernel_gamma2;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = cj->hydro.parts[j].x[k] - ci->hydro.parts[i].x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hjg2 || r2 < hig2) {

        /* Interact */
        runner_iact_nonsym_force(r2, dx, hj, pi->h, pj, pi, a, H);
      }
    }
  }
}

void pairs_all_stars_density(struct runner *r, struct cell *ci,
                             struct cell *cj) {

  float r2, dx[3];
  const double dim[3] = {r->e->s->dim[0], r->e->s->dim[1], r->e->s->dim[2]};
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->stars.count; ++i) {
    struct spart *spi = &ci->stars.parts[i];

    float hi = spi->h;
    float hig2 = hi * hi * kernel_gamma2;

    /* Skip inactive particles. */
    if (!spart_is_active(spi, e)) continue;

    for (int j = 0; j < cj->hydro.count; ++j) {

      struct part *pj = &cj->hydro.parts[j];

      /* Early abort? */
      if (part_is_inhibited(pj, e)) continue;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = spi->x[k] - pj->x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2) {
        /* Interact */
        runner_iact_nonsym_stars_density(r2, dx, hi, pj->h, spi, pj, a, H);
      }
    }
  }

  /* Reverse double-for loop and checks every interaction */
  for (int j = 0; j < cj->stars.count; ++j) {

    struct spart *spj = &cj->stars.parts[j];
    float hj = spj->h;
    float hjg2 = hj * hj * kernel_gamma2;

    /* Skip inactive particles. */
    if (!spart_is_active(spj, e)) continue;

    for (int i = 0; i < ci->hydro.count; ++i) {

      struct part *pi = &ci->hydro.parts[i];

      /* Early abort? */
      if (part_is_inhibited(pi, e)) continue;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = spj->x[k] - pi->x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hjg2) {
        /* Interact */
        runner_iact_nonsym_stars_density(r2, dx, hj, pi->h, spj, pi, a, H);
      }
    }
  }
}

void self_all_density(struct runner *r, struct cell *ci) {
  float r2, hi, hj, hig2, hjg2, dxi[3];  //, dxj[3];
  struct part *pi, *pj;
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->hydro.count; ++i) {

    pi = &ci->hydro.parts[i];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    for (int j = i + 1; j < ci->hydro.count; ++j) {

      pj = &ci->hydro.parts[j];
      hj = pj->h;
      hjg2 = hj * hj * kernel_gamma2;

      if (pi == pj) continue;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dxi[k] = ci->hydro.parts[i].x[k] - ci->hydro.parts[j].x[k];
        r2 += dxi[k] * dxi[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 && part_is_active(pi, e) && !part_is_inhibited(pj, e)) {

        /* Interact */
        runner_iact_nonsym_density(r2, dxi, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_chemistry(r2, dxi, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_star_formation(r2, dxi, hi, hj, pi, pj, a, H);
      }

      /* Hit or miss? */
      if (r2 < hjg2 && part_is_active(pj, e) && !part_is_inhibited(pi, e)) {

        dxi[0] = -dxi[0];
        dxi[1] = -dxi[1];
        dxi[2] = -dxi[2];

        /* Interact */
        runner_iact_nonsym_density(r2, dxi, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_chemistry(r2, dxi, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_star_formation(r2, dxi, hj, hi, pj, pi, a, H);
      }
    }
  }
}

#ifdef EXTRA_HYDRO_LOOP
void self_all_gradient(struct runner *r, struct cell *ci) {
  float r2, hi, hj, hig2, hjg2, dxi[3];  //, dxj[3];
  struct part *pi, *pj;
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->hydro.count; ++i) {

    pi = &ci->hydro.parts[i];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    for (int j = i + 1; j < ci->hydro.count; ++j) {

      pj = &ci->hydro.parts[j];
      hj = pj->h;
      hjg2 = hj * hj * kernel_gamma2;

      if (pi == pj) continue;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dxi[k] = ci->hydro.parts[i].x[k] - ci->hydro.parts[j].x[k];
        r2 += dxi[k] * dxi[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 && part_is_active(pi, e) && !part_is_inhibited(pj, e)) {

        /* Interact */
        runner_iact_nonsym_gradient(r2, dxi, hi, hj, pi, pj, a, H);
      }

      /* Hit or miss? */
      if (r2 < hjg2 && part_is_active(pj, e) && !part_is_inhibited(pi, e)) {

        dxi[0] = -dxi[0];
        dxi[1] = -dxi[1];
        dxi[2] = -dxi[2];

        /* Interact */
        runner_iact_nonsym_gradient(r2, dxi, hj, hi, pj, pi, a, H);
      }
    }
  }
}
#endif /* EXTRA_HYDRO_LOOP */

void self_all_force(struct runner *r, struct cell *ci) {
  float r2, hi, hj, hig2, hjg2, dxi[3];  //, dxj[3];
  struct part *pi, *pj;
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->hydro.count; ++i) {

    pi = &ci->hydro.parts[i];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    for (int j = i + 1; j < ci->hydro.count; ++j) {

      pj = &ci->hydro.parts[j];
      hj = pj->h;
      hjg2 = hj * hj * kernel_gamma2;

      if (pi == pj) continue;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dxi[k] = ci->hydro.parts[i].x[k] - ci->hydro.parts[j].x[k];
        r2 += dxi[k] * dxi[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 || r2 < hjg2) {

        /* Interact */
        runner_iact_force(r2, dxi, hi, hj, pi, pj, a, H);
      }
    }
  }
}

void self_all_stars_density(struct runner *r, struct cell *ci) {
  float r2, hi, hj, hig2, dxi[3];
  struct spart *spi;
  struct part *pj;
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->stars.count; ++i) {

    spi = &ci->stars.parts[i];
    hi = spi->h;
    hig2 = hi * hi * kernel_gamma2;

    if (!spart_is_active(spi, e)) continue;

    for (int j = 0; j < ci->hydro.count; ++j) {

      pj = &ci->hydro.parts[j];
      hj = pj->h;

      /* Early abort? */
      if (part_is_inhibited(pj, e)) continue;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dxi[k] = spi->x[k] - pj->x[k];
        r2 += dxi[k] * dxi[k];
      }

      /* Hit or miss? */
      if (r2 < hig2) {
        /* Interact */
        runner_iact_nonsym_stars_density(r2, dxi, hi, hj, spi, pj, a, H);
      }
    }
  }
}

/**
 * @brief Compute the force on a single particle brute-force.
 */
void engine_single_density(double *dim, long long int pid,
                           struct part *restrict parts, int N, int periodic,
                           const struct cosmology *cosmo) {
  double r2, dx[3];
  float fdx[3];
  struct part p;
  float a = 1.f, H = 0.f;

  /* Find "our" part. */
  int k;
  for (k = 0; k < N && parts[k].id != pid; k++)
    ;
  if (k == N) error("Part not found.");
  p = parts[k];

  /* Clear accumulators. */
  hydro_init_part(&p, NULL);

  /* Loop over all particle pairs (force). */
  for (k = 0; k < N; k++) {
    if (parts[k].id == p.id) continue;
    for (int i = 0; i < 3; i++) {
      dx[i] = p.x[i] - parts[k].x[i];
      if (periodic) {
        if (dx[i] < -dim[i] / 2)
          dx[i] += dim[i];
        else if (dx[i] > dim[i] / 2)
          dx[i] -= dim[i];
      }
      fdx[i] = dx[i];
    }
    r2 = fdx[0] * fdx[0] + fdx[1] * fdx[1] + fdx[2] * fdx[2];
    if (r2 < p.h * p.h * kernel_gamma2) {
      runner_iact_nonsym_density(r2, fdx, p.h, parts[k].h, &p, &parts[k], a, H);
    }
  }

  /* Dump the result. */
  hydro_end_density(&p, cosmo);
  message("part %lli (h=%e) has wcount=%e, rho=%e.", p.id, p.h,
          p.density.wcount, hydro_get_comoving_density(&p));
  fflush(stdout);
}

void engine_single_force(double *dim, long long int pid,
                         struct part *restrict parts, int N, int periodic) {
  int i, k;
  double r2, dx[3];
  float fdx[3];
  struct part p;
  float a = 1.f, H = 0.f;

  /* Find "our" part. */
  for (k = 0; k < N && parts[k].id != pid; k++)
    ;
  if (k == N) error("Part not found.");
  p = parts[k];

  /* Clear accumulators. */
  hydro_reset_acceleration(&p);

  /* Loop over all particle pairs (force). */
  for (k = 0; k < N; k++) {
    // for ( k = N-1 ; k >= 0 ; k-- ) {
    if (parts[k].id == p.id) continue;
    for (i = 0; i < 3; i++) {
      dx[i] = p.x[i] - parts[k].x[i];
      if (periodic) {
        if (dx[i] < -dim[i] / 2)
          dx[i] += dim[i];
        else if (dx[i] > dim[i] / 2)
          dx[i] -= dim[i];
      }
      fdx[i] = dx[i];
    }
    r2 = fdx[0] * fdx[0] + fdx[1] * fdx[1] + fdx[2] * fdx[2];
    if (r2 < p.h * p.h * kernel_gamma2 ||
        r2 < parts[k].h * parts[k].h * kernel_gamma2) {
      hydro_reset_acceleration(&p);
      runner_iact_nonsym_force(r2, fdx, p.h, parts[k].h, &p, &parts[k], a, H);
    }
  }

  /* Dump the result. */
  message("part %lli (h=%e) has a=[%.3e,%.3e,%.3e]", p.id, p.h, p.a_hydro[0],
          p.a_hydro[1], p.a_hydro[2]);
  fflush(stdout);
}

/**
 * Returns a random number (uniformly distributed) in [a,b[
 */
double random_uniform(double a, double b) {
  return (rand() / (double)RAND_MAX) * (b - a) + a;
}

/**
 * @brief Randomly shuffle an array of particles.
 */
void shuffle_particles(struct part *parts, const int count) {
  if (count > 1) {
    for (int i = 0; i < count - 1; i++) {
      int j = i + random_uniform(0., (double)(count - 1 - i));

      struct part particle = parts[j];

      parts[j] = parts[i];

      parts[i] = particle;
    }
  }
}

/**
 * @brief Randomly shuffle an array of sparticles.
 */
void shuffle_sparticles(struct spart *sparts, const int scount) {
  if (scount > 1) {
    for (int i = 0; i < scount - 1; i++) {
      int j = i + random_uniform(0., (double)(scount - 1 - i));

      struct spart sparticle = sparts[j];

      sparts[j] = sparts[i];

      sparts[i] = sparticle;
    }
  }
}

/**
 * @brief Compares two values based on their relative difference: |a - b|/|a +
 * b|
 *
 * @param a Value a
 * @param b Value b
 * @param threshold The limit on the relative difference between the two values
 * @param absDiff Absolute difference: |a - b|
 * @param absSum Absolute sum: |a + b|
 * @param relDiff Relative difference: |a - b|/|a + b|
 *
 * @return 1 if difference found, 0 otherwise
 */
int compare_values(double a, double b, double threshold, double *absDiff,
                   double *absSum, double *relDiff) {

  int result = 0;
  *absDiff = 0.0, *absSum = 0.0, *relDiff = 0.0;

  *absDiff = fabs(a - b);
  *absSum = fabs(a + b);
  if (*absSum > 0.f) {
    *relDiff = *absDiff / *absSum;
  }

  if (*relDiff > threshold) {
    result = 1;
  }

  return result;
}

/**
 * @brief Compares two particles' properties using the relative difference and a
 * threshold.
 *
 * @param a Particle A
 * @param b Particle B
 * @param threshold The limit on the relative difference between the two values
 *
 * @return 1 if difference found, 0 otherwise
 */
int compare_particles(struct part *a, struct part *b, double threshold) {

#ifdef GADGET2_SPH

  int result = 0;
  double absDiff = 0.0, absSum = 0.0, relDiff = 0.0;

  for (int k = 0; k < 3; k++) {
    if (compare_values(a->x[k], b->x[k], threshold, &absDiff, &absSum,
                       &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for x[%d] of "
          "particle %lld.",
          relDiff, threshold, k, a->id);
      message("a = %e, b = %e", a->x[k], b->x[k]);
      result = 1;
    }
  }
  for (int k = 0; k < 3; k++) {
    if (compare_values(a->v[k], b->v[k], threshold, &absDiff, &absSum,
                       &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for v[%d] of "
          "particle %lld.",
          relDiff, threshold, k, a->id);
      message("a = %e, b = %e", a->v[k], b->v[k]);
      result = 1;
    }
  }
  for (int k = 0; k < 3; k++) {
    if (compare_values(a->a_hydro[k], b->a_hydro[k], threshold, &absDiff,
                       &absSum, &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for a_hydro[%d] "
          "of particle %lld.",
          relDiff, threshold, k, a->id);
      message("a = %e, b = %e", a->a_hydro[k], b->a_hydro[k]);
      result = 1;
    }
  }
  if (compare_values(a->rho, b->rho, threshold, &absDiff, &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for rho of "
        "particle %lld.",
        relDiff, threshold, a->id);
    message("a = %e, b = %e", a->rho, b->rho);
    result = 1;
  }
  if (compare_values(a->density.rho_dh, b->density.rho_dh, threshold, &absDiff,
                     &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for rho_dh of "
        "particle %lld.",
        relDiff, threshold, a->id);
    message("a = %e, b = %e", a->density.rho_dh, b->density.rho_dh);
    result = 1;
  }
  if (compare_values(a->density.wcount, b->density.wcount, threshold, &absDiff,
                     &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for wcount of "
        "particle %lld.",
        relDiff, threshold, a->id);
    message("a = %e, b = %e", a->density.wcount, b->density.wcount);
    result = 1;
  }
  if (compare_values(a->density.wcount_dh, b->density.wcount_dh, threshold,
                     &absDiff, &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for wcount_dh of "
        "particle %lld.",
        relDiff, threshold, a->id);
    message("a = %e, b = %e", a->density.wcount_dh, b->density.wcount_dh);
    result = 1;
  }
  if (compare_values(a->force.h_dt, b->force.h_dt, threshold, &absDiff, &absSum,
                     &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for h_dt of "
        "particle %lld.",
        relDiff, threshold, a->id);
    message("a = %e, b = %e", a->force.h_dt, b->force.h_dt);
    result = 1;
  }
  if (compare_values(a->force.v_sig, b->force.v_sig, threshold, &absDiff,
                     &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for v_sig of "
        "particle %lld.",
        relDiff, threshold, a->id);
    message("a = %e, b = %e", a->force.v_sig, b->force.v_sig);
    result = 1;
  }
  if (compare_values(a->entropy_dt, b->entropy_dt, threshold, &absDiff, &absSum,
                     &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for entropy_dt of "
        "particle %lld.",
        relDiff, threshold, a->id);
    message("a = %e, b = %e", a->entropy_dt, b->entropy_dt);
    result = 1;
  }
  if (compare_values(a->density.div_v, b->density.div_v, threshold, &absDiff,
                     &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for div_v of "
        "particle %lld.",
        relDiff, threshold, a->id);
    message("a = %e, b = %e", a->density.div_v, b->density.div_v);
    result = 1;
  }
  for (int k = 0; k < 3; k++) {
    if (compare_values(a->density.rot_v[k], b->density.rot_v[k], threshold,
                       &absDiff, &absSum, &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for rot_v[%d] "
          "of particle %lld.",
          relDiff, threshold, k, a->id);
      message("a = %e, b = %e", a->density.rot_v[k], b->density.rot_v[k]);
      result = 1;
    }
  }

  return result;

#else

  error("Function not supported for this flavour of SPH");
  return 0;

#endif
}

/**
 * @brief return the resident memory use of the process and its children.
 *
 * @result memory use in Kb.
 */
long get_maxrss() {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return usage.ru_maxrss;
}

/**
 * @brief trim leading white space from a string.
 *
 * Returns pointer to first character.
 *
 * @param s the string.
 * @result the result.
 */
char *trim_leading(char *s) {
  if (s == NULL || strlen(s) < 2) return s;
  while (isspace(*s)) s++;
  return s;
}

/**
 * @brief trim trailing white space from a string.
 *
 * Modifies the string by adding a NULL to the end.
 *
 * @param s the string.
 * @result the result.
 */
char *trim_trailing(char *s) {
  if (s == NULL || strlen(s) < 2) return s;
  char *end = s + strlen(s) - 1;
  while (isspace(*end)) end--;
  *(end + 1) = '\0';
  return s;
}

/**
 * @brief trim leading and trailing white space from a string.
 *
 * Can modify the string by adding a NULL to the end.
 *
 * @param s the string.
 * @result the result.
 */
char *trim_both(char *s) {
  if (s == NULL || strlen(s) < 2) return s;
  return trim_trailing(trim_leading(s));
}
