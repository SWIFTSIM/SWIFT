/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#include <config.h>

/* Some standard headers. */
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "runner_doiact_grav.h"
#include "swift.h"

const int L = 8;
const float THETA = 0.5;
const float EPSILON = 0.001;
const int USE_MAC = 1;
 
struct cell *make_cell(const size_t n, const double offset[3],
                       const double size, long long *partId) {

  const size_t count = n * n * n;
  struct cell *cell = NULL;
  if (posix_memalign((void **)&cell, cell_align, sizeof(struct cell)) != 0) {
    error("Couldn't allocate the cell");
  }
  bzero(cell, sizeof(struct cell));

  if (posix_memalign((void **)&cell->grav.parts, gpart_align,
                     count * sizeof(struct gpart)) != 0) {
    error("couldn't allocate particles, no. of particles: %d", (int)count);
  }
  bzero(cell->grav.parts, count * sizeof(struct gpart));

  if (posix_memalign((void **)&cell->grav.multipole, multipole_align,
                     sizeof(struct gravity_tensors)) != 0) {
    error("Couldn't allocate the mulipole");
  }

  /* Construct the parts */
  struct gpart *gpart = cell->grav.parts;
  for (size_t x = 0; x < n; ++x) {
    for (size_t y = 0; y < n; ++y) {
      for (size_t z = 0; z < n; ++z) {
	
	gpart->x[0] = offset[0] + random_uniform(0., size);
        gpart->x[1] = offset[1] + random_uniform(0., size);
        gpart->x[2] = offset[2] + random_uniform(0., size);

        gpart->v_full[0] = 0.f;
        gpart->v_full[1] = 0.f;
        gpart->v_full[2] = 0.f;

        gpart->a_grav[0] = 0.f;
        gpart->a_grav[1] = 0.f;
        gpart->a_grav[2] = 0.f;

        gpart->id_or_neg_offset = -(++(*partId));
        gpart->mass = 1.f;
        gpart->time_bin = 1;
	gpart->type = swift_type_dark_matter;
	gpart->epsilon = 0.001;

#ifdef SWIFT_DEBUG_CHECKS
        gpart->ti_drift = 8;
        gpart->ti_kick = 8;
#endif

	gpart++;
      }
    }
  }

  /* Cell properties */
  cell->split = 0;
  cell->grav.count = count;
  cell->width[0] = size;
  cell->width[1] = size;
  cell->width[2] = size;
  cell->loc[0] = offset[0];
  cell->loc[1] = offset[1];
  cell->loc[2] = offset[2];

  cell->grav.super = cell;
  cell->grav.ti_old_part = 8;
  cell->grav.ti_end_min = 8;

  return cell;
}

void clean_up(struct cell *c) {
  free(c->grav.parts);
  free(c->grav.multipole);
  free(c);
}

void end_force(struct cell *c) {
  for (int i = 0; i < c->grav.count; ++i) {
    gravity_end_force(&c->grav.parts[i], /*G=*/1.f, 0.f, /*periodic=*/0, /*self_gravity=*/1);
  }
}

void init(struct cell *c) {
  for (int i = 0; i < c->grav.count; ++i) {
    gravity_init_gpart(&c->grav.parts[i]);
  }
}

void grav_down(struct cell *c) {
  const struct grav_tensor *pot = &c->grav.multipole->pot;
  const double CoM[3] = {c->grav.multipole->CoM[0], c->grav.multipole->CoM[1],
    c->grav.multipole->CoM[2]};  
  for (int i = 0; i < c->grav.count; ++i) {
    gravity_L2P(pot, CoM, &c->grav.parts[i]);
  }
}


void compute_exact_forces(struct cell *ci, const struct cell *cj) {

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  
  for (int i = 0; i < ci->grav.count; ++i) {

    struct gpart *gpi = &ci->grav.parts[i];
    const double pix[3] = {gpi->x[0], gpi->x[1], gpi->x[2]};

    double a_grav[3] = {0., 0., 0.};
    
    for (int j = 0; j < cj->grav.count; ++j) {

      const struct gpart *gpj = &cj->grav.parts[j];
      const double pjx[3] = {gpj->x[0], gpj->x[1], gpj->x[2]};

      /* Compute the pairwise distance. */
      const double dx = pjx[0] - pix[0];
      const double dy = pjx[1] - pix[1];
      const double dz = pjx[2] - pix[2];

      const double r2 = dx * dx + dy * dy + dz * dz;
      const double r_inv = 1. / sqrt(r2);
      const double mj = gpj->mass;

      /* Exact Newtonian gravity */
      const double f = mj * r_inv * r_inv * r_inv;
      a_grav[0] += f * dx;
      a_grav[1] += f * dy;
      a_grav[2] += f * dz;
    }

    gpi->a_grav_exact[0] = a_grav[0];
    gpi->a_grav_exact[1] = a_grav[1];
    gpi->a_grav_exact[2] = a_grav[2];
  }

#endif
}

int compare(const void *p1, const void *p2) {

  const double a = *(double *)p1;
  const double b = *(double *)p2;
  if ( a< b)
    return -1;
  else if ( a > b)
    return 1;
  else
    return 0;
    
}

void check_error(const struct cell *ci, double *ret_mean, double *ret_std, double *ret_median, double *ret_per99) {

  double *error = (double*)malloc(ci->grav.count * sizeof(double));
  bzero(error, ci->grav.count * sizeof(double));
  
  for (int i = 0; i < ci->grav.count; ++i) {

    const struct gpart *gpi = &ci->grav.parts[i];

    const double a_swift[3] = {gpi->a_grav[0], gpi->a_grav[1], gpi->a_grav[2]}; 
    const double a_exact[3] = {gpi->a_grav_exact[0], gpi->a_grav_exact[1], gpi->a_grav_exact[2]};   

    double a_swift_norm = a_swift[0] * a_swift[0] + a_swift[1] * a_swift[1] + a_swift[2] * a_swift[2];
    double a_exact_norm = a_exact[0] * a_exact[0] + a_exact[1] * a_exact[1] + a_exact[2] * a_exact[2];

    a_swift_norm = sqrt(a_swift_norm);
    a_exact_norm = sqrt(a_exact_norm);

    const double err_rel = fabs(a_swift_norm - a_exact_norm) / a_exact_norm;

    error[i] = err_rel;
  }

  qsort(error, ci->grav.count, sizeof(double), compare);

  double mean = 0., mean2 = 0.;
  for (int i = 0; i < ci->grav.count; ++i) {
    mean += error[i];
    mean2 += error[i] * error[i];
  }
  mean /= (double) ci->grav.count;
  mean2 /= (double) ci->grav.count;
  const double std = sqrt(mean2 - mean * mean);
  const double median = error[ci->grav.count / 2];
  const double per99 = error[(int)(ci->grav.count * 0.99)];

  *ret_mean = mean;
  *ret_std = std;
  *ret_median = median;
  *ret_per99 = per99;
  
  free(error);
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  srand(time(NULL));
  
/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Initialise a few things to get us going */
  struct engine e;
  bzero(&e, sizeof(struct engine));
  e.max_active_bin = num_time_bins;
  e.time = 0.1f;
  e.ti_current = 8;
  e.time_base = 1e-10;

  struct pm_mesh mesh;
  mesh.periodic = 0;
  mesh.dim[0] = 10.;
  mesh.dim[1] = 10.;
  mesh.dim[2] = 10.;
  mesh.r_s_inv = 0.;
  mesh.r_cut_min = 0.;
  e.mesh = &mesh;

  struct gravity_props props;
  props.a_smooth = 1.25;
  props.epsilon_DM_cur = 1e-3;
  props.theta_crit = THETA;
  props.adaptive_tolerance = EPSILON;
  props.use_advanced_MAC = USE_MAC;
  e.gravity_properties = &props;

  struct runner runner;
  bzero(&runner, sizeof(struct runner));
  runner.e = &e;

  /* Init the cache for gravity interaction */
  gravity_cache_init(&runner.ci_gravity_cache, L * L * L);
  gravity_cache_init(&runner.cj_gravity_cache, L * L * L);

  for (int i = 0; i < 100; i+=5) {
    
    long long partID = 0LL;
    
    const double r = 1. + 0.1 * (double) i;
  
    const double offset_i[3] = {0., 0., 0.};
    const double offset_j[3] = {r, 0., 0.};
    
    /* Build cells */
    struct cell *ci = make_cell(L, offset_i, 1.0, &partID);
    struct cell *cj = make_cell(L, offset_j, 1.0, &partID);

    init(ci);
    init(ci);
    
    /* Build multipoles */
    gravity_P2M(ci->grav.multipole, ci->grav.parts, ci->grav.count, &props);
    gravity_P2M(cj->grav.multipole, cj->grav.parts, cj->grav.count, &props);

    gravity_multipole_compute_power(&ci->grav.multipole->m_pole);
    gravity_multipole_compute_power(&cj->grav.multipole->m_pole);
    
    /* Init the field tensors */
    gravity_field_tensors_init(&ci->grav.multipole->pot, e.ti_current);
    gravity_field_tensors_init(&cj->grav.multipole->pot, e.ti_current);
    
    /* Interact the two cells the SWIFT way */
    runner_dopair_recursive_grav(&runner, ci, cj, 0);

    /* Propagate the field tensor to the particles */
    grav_down(ci);
    grav_down(cj);

    /* Finish the calculation */
    end_force(ci);
    end_force(cj);

    if (USE_MAC) {
      
      init(ci);
      init(ci);
      
      /* Build multipoles */
      gravity_P2M(ci->grav.multipole, ci->grav.parts, ci->grav.count, &props);
      gravity_P2M(cj->grav.multipole, cj->grav.parts, cj->grav.count, &props);

      gravity_multipole_compute_power(&ci->grav.multipole->m_pole);
      gravity_multipole_compute_power(&cj->grav.multipole->m_pole);

      /* Init the field tensors */
      gravity_field_tensors_init(&ci->grav.multipole->pot, e.ti_current);
      gravity_field_tensors_init(&cj->grav.multipole->pot, e.ti_current);
      
      /* Interact the two cells the SWIFT way */
      runner_dopair_recursive_grav(&runner, ci, cj, 0);

      /* Propagate the field tensor to the particles */
      grav_down(ci);
      grav_down(cj);

      /* Finish the calculation */
      end_force(ci);
      end_force(cj);
    }
      
    /* Direct summation comparison */
    compute_exact_forces(ci, cj);
    
    /* Compute differences */
    double mean, std, median, per99;
    check_error(ci, &mean, &std, &median, &per99);
    
    message("%d %e %e %e %e %e %e %e",
	    USE_MAC, THETA, EPSILON,
	    r, mean, std, median, per99);
  
    clean_up(ci);
    clean_up(cj);
  }

  return 0;
}
