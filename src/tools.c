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
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

/* This object's header. */
#include "tools.h"

/* Local includes. */
#include "active.h"
#include "cell.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "part.h"
#include "periodic.h"
#include "runner.h"

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
                            &parts[k]);
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
  // int mj, mk;
  // double maxratio = 1.0;
  double r2, dx[3];
  float fdx[3];
  struct part p;
  // double ih = 12.0/6.25;

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
      runner_iact_nonsym_density(r2, fdx, p.h, parts[k].h, &p, &parts[k]);
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

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->count; ++i) {

    pi = &ci->parts[i];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(pi, e)) continue;

    for (int j = 0; j < cj->count; ++j) {

      pj = &cj->parts[j];

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = ci->parts[i].x[k] - cj->parts[j].x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2) {

        /* Interact */
        runner_iact_nonsym_density(r2, dx, hi, pj->h, pi, pj);
      }
    }
  }

  /* Reverse double-for loop and checks every interaction */
  for (int j = 0; j < cj->count; ++j) {

    pj = &cj->parts[j];
    hj = pj->h;
    hjg2 = hj * hj * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(pj, e)) continue;

    for (int i = 0; i < ci->count; ++i) {

      pi = &ci->parts[i];

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = cj->parts[j].x[k] - ci->parts[i].x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hjg2) {

        /* Interact */
        runner_iact_nonsym_density(r2, dx, hj, pi->h, pj, pi);
      }
    }
  }
}

void pairs_all_force(struct runner *r, struct cell *ci, struct cell *cj) {

  float r2, hi, hj, hig2, hjg2, dx[3];
  struct part *pi, *pj;
  const double dim[3] = {r->e->s->dim[0], r->e->s->dim[1], r->e->s->dim[2]};

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->count; ++i) {

    pi = &ci->parts[i];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    for (int j = 0; j < cj->count; ++j) {

      pj = &cj->parts[j];
      hj = pj->h;
      hjg2 = hj * hj * kernel_gamma2;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = ci->parts[i].x[k] - cj->parts[j].x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 || r2 < hjg2) {

        /* Interact */
        runner_iact_nonsym_force(r2, dx, hi, hj, pi, pj);
      }
    }
  }

  /* Reverse double-for loop and checks every interaction */
  for (int j = 0; j < cj->count; ++j) {

    pj = &cj->parts[j];
    hj = pj->h;
    hjg2 = hj * hj * kernel_gamma2;

    for (int i = 0; i < ci->count; ++i) {

      pi = &ci->parts[i];
      hi = pi->h;
      hig2 = hi * hi * kernel_gamma2;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = cj->parts[j].x[k] - ci->parts[i].x[k];
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hjg2 || r2 < hig2) {

        /* Interact */
        runner_iact_nonsym_force(r2, dx, hj, pi->h, pj, pi);
      }
    }
  }
}

void self_all_density(struct runner *r, struct cell *ci) {
  float r2, hi, hj, hig2, hjg2, dxi[3];  //, dxj[3];
  struct part *pi, *pj;
  const struct engine *e = r->e;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->count; ++i) {

    pi = &ci->parts[i];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    for (int j = i + 1; j < ci->count; ++j) {

      pj = &ci->parts[j];
      hj = pj->h;
      hjg2 = hj * hj * kernel_gamma2;

      if (pi == pj) continue;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dxi[k] = ci->parts[i].x[k] - ci->parts[j].x[k];
        r2 += dxi[k] * dxi[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 && part_is_active(pi, e)) {

        /* Interact */
        runner_iact_nonsym_density(r2, dxi, hi, hj, pi, pj);
      }

      /* Hit or miss? */
      if (r2 < hjg2 && part_is_active(pj, e)) {

        dxi[0] = -dxi[0];
        dxi[1] = -dxi[1];
        dxi[2] = -dxi[2];

        /* Interact */
        runner_iact_nonsym_density(r2, dxi, hj, hi, pj, pi);
      }
    }
  }
}

void self_all_force(struct runner *r, struct cell *ci) {
  float r2, hi, hj, hig2, hjg2, dxi[3];  //, dxj[3];
  struct part *pi, *pj;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->count; ++i) {

    pi = &ci->parts[i];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    for (int j = i + 1; j < ci->count; ++j) {

      pj = &ci->parts[j];
      hj = pj->h;
      hjg2 = hj * hj * kernel_gamma2;

      if (pi == pj) continue;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dxi[k] = ci->parts[i].x[k] - ci->parts[j].x[k];
        r2 += dxi[k] * dxi[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 || r2 < hjg2) {

        /* Interact */
        runner_iact_force(r2, dxi, hi, hj, pi, pj);
      }
    }
  }
}

/**
 * @brief Compute the force on a single particle brute-force.
 */
void engine_single_density(double *dim, long long int pid,
                           struct part *restrict parts, int N, int periodic) {
  double r2, dx[3];
  float fdx[3];
  struct part p;

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
      runner_iact_nonsym_density(r2, fdx, p.h, parts[k].h, &p, &parts[k]);
    }
  }

  /* Dump the result. */
  hydro_end_density(&p);
  message("part %lli (h=%e) has wcount=%e, rho=%e.", p.id, p.h,
          p.density.wcount, hydro_get_density(&p));
  fflush(stdout);
}

void engine_single_force(double *dim, long long int pid,
                         struct part *restrict parts, int N, int periodic) {
  int i, k;
  double r2, dx[3];
  float fdx[3];
  struct part p;

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
      runner_iact_nonsym_force(r2, fdx, p.h, parts[k].h, &p, &parts[k]);
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

  } else
    error("Array not big enough to shuffle!");
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
int compare_particles(struct part a, struct part b, double threshold) {

#ifdef GADGET2_SPH

  int result = 0;
  double absDiff = 0.0, absSum = 0.0, relDiff = 0.0;

  for (int k = 0; k < 3; k++) {
    if (compare_values(a.x[k], b.x[k], threshold, &absDiff, &absSum,
                       &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for x[%d] of "
          "particle %lld.",
          relDiff, threshold, k, a.id);
      message("a = %e, b = %e", a.x[k], b.x[k]);
      result = 1;
    }
  }
  for (int k = 0; k < 3; k++) {
    if (compare_values(a.v[k], b.v[k], threshold, &absDiff, &absSum,
                       &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for v[%d] of "
          "particle %lld.",
          relDiff, threshold, k, a.id);
      message("a = %e, b = %e", a.v[k], b.v[k]);
      result = 1;
    }
  }
  for (int k = 0; k < 3; k++) {
    if (compare_values(a.a_hydro[k], b.a_hydro[k], threshold, &absDiff, &absSum,
                       &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for a_hydro[%d] "
          "of particle %lld.",
          relDiff, threshold, k, a.id);
      message("a = %e, b = %e", a.a_hydro[k], b.a_hydro[k]);
      result = 1;
    }
  }
  if (compare_values(a.rho, b.rho, threshold, &absDiff, &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for rho of "
        "particle %lld.",
        relDiff, threshold, a.id);
    message("a = %e, b = %e", a.rho, b.rho);
    result = 1;
  }
  if (compare_values(a.density.rho_dh, b.density.rho_dh, threshold, &absDiff,
                     &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for rho_dh of "
        "particle %lld.",
        relDiff, threshold, a.id);
    message("a = %e, b = %e", a.density.rho_dh, b.density.rho_dh);
    result = 1;
  }
  if (compare_values(a.density.wcount, b.density.wcount, threshold, &absDiff,
                     &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for wcount of "
        "particle %lld.",
        relDiff, threshold, a.id);
    message("a = %e, b = %e", a.density.wcount, b.density.wcount);
    result = 1;
  }
  if (compare_values(a.density.wcount_dh, b.density.wcount_dh, threshold,
                     &absDiff, &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for wcount_dh of "
        "particle %lld.",
        relDiff, threshold, a.id);
    message("a = %e, b = %e", a.density.wcount_dh, b.density.wcount_dh);
    result = 1;
  }
  if (compare_values(a.force.h_dt, b.force.h_dt, threshold, &absDiff, &absSum,
                     &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for h_dt of "
        "particle %lld.",
        relDiff, threshold, a.id);
    message("a = %e, b = %e", a.force.h_dt, b.force.h_dt);
    result = 1;
  }
  if (compare_values(a.force.v_sig, b.force.v_sig, threshold, &absDiff, &absSum,
                     &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for v_sig of "
        "particle %lld.",
        relDiff, threshold, a.id);
    message("a = %e, b = %e", a.force.v_sig, b.force.v_sig);
    result = 1;
  }
  if (compare_values(a.entropy_dt, b.entropy_dt, threshold, &absDiff, &absSum,
                     &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for entropy_dt of "
        "particle %lld.",
        relDiff, threshold, a.id);
    message("a = %e, b = %e", a.entropy_dt, b.entropy_dt);
    result = 1;
  }
  if (compare_values(a.density.div_v, b.density.div_v, threshold, &absDiff,
                     &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for div_v of "
        "particle %lld.",
        relDiff, threshold, a.id);
    message("a = %e, b = %e", a.density.div_v, b.density.div_v);
    result = 1;
  }
  for (int k = 0; k < 3; k++) {
    if (compare_values(a.density.rot_v[k], b.density.rot_v[k], threshold,
                       &absDiff, &absSum, &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for rot_v[%d] "
          "of particle %lld.",
          relDiff, threshold, k, a.id);
      message("a = %e, b = %e", a.density.rot_v[k], b.density.rot_v[k]);
      result = 1;
    }
  }

  return result;

#else

  error("Function not supported for this flavour of SPH");
  return 0;

#endif
}
