/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

/* This object's header. */
#include "projected_kernel.h"

/* Local headers */
#include "error.h"
#include "kernel_hydro.h"
#include "math.h"

#ifdef HAVE_LIBGSL
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#endif

/**
 * @brief Integrand used in evaluating the projected kernel
 *
 * See section 4.3.1 in Price et al. 2007:
 * https://ui.adsabs.harvard.edu/abs/2007PASA...24..159P/abstract
 *
 * This function is used to carry out the integral in equation 30.
 *
 * @param qz z coordinate at which to evaluate the kernel, in units of h
 * @param param Ratio of distance in the xy plane in units of h
 */
static double projected_kernel_integrand(double qz, void *param) {

  const double qxy = *((double *)param);
  const double q = sqrt(pow(qxy, 2.0) + pow(qz, 2.0));
  double W;
  kernel_eval_double(q, &W);
  return W;
}

/**
 * @brief Computes 2D projection of the 3D kernel function.
 *
 * Given a distance in the xy plane, we integrate along the
 * z axis to evaluate the projected kernel.
 *
 * @param u The ratio of the (2D) distance to the smoothing length
 */
double projected_kernel_integrate(double u) {

#ifdef HAVE_LIBGSL

  /* Swift's hydro kernel can be evaluated with kernel_eval(u, W)
     where u = r / h and W returns the result. The kernel goes to
     zero at u=kernel_gamma. Projection is only implemented in 3D.*/
#ifndef HYDRO_DIMENSION_3D
  error("projected_kernel_eval() is only defined for the 3D case.");
#endif

  /* Initalise the GSL workspace */
  const size_t workspace_size = 100000;
  gsl_integration_workspace *space =
      gsl_integration_workspace_alloc(workspace_size);

  /* Compute the integral */
  double result;
  double abserr;
  double qxy = u;
  const double qz_max = sqrt(pow(kernel_gamma, 2.0) - pow(qxy, 2.0));
  const double qz_min = -qz_max;
  gsl_function F = {&projected_kernel_integrand, &qxy};
  gsl_integration_qag(&F, qz_min, qz_max, 1.0e-10, 1.0e-10, workspace_size,
                      GSL_INTEG_GAUSS61, space, &result, &abserr);

  /* Free the workspace */
  gsl_integration_workspace_free(space);

  return result;

#else
  error("Need GSL library to evaluate the projected kernel");
  return 0.0;
#endif
}

/**
 * @brief Tabulate the projected kernel
 *
 * @param tab The projected_kernel_table struct
 */
void projected_kernel_init(struct projected_kernel_table *tab) {

  /* Allocate storage */
  tab->n = PROJECTED_KERNEL_NTAB;
  tab->value = malloc(sizeof(double) * tab->n);

  /* Determine range to tabulate */
  tab->u_max = kernel_gamma;
  tab->du = tab->u_max / (tab->n - 1);
  tab->inv_du = 1.0 / tab->du;

  /* Evaluate the kernel at points in the table */
  for (int i = 0; i < tab->n - 1; i += 1)
    tab->value[i] = projected_kernel_integrate(i * tab->du);
  tab->value[tab->n - 1] = 0.0;
}

/**
 * @brief Deallocate the projected kernel table
 */
void projected_kernel_clean(struct projected_kernel_table *tab) {
  free(tab->value);
}

void projected_kernel_dump(void) {

  struct projected_kernel_table tab;
  projected_kernel_init(&tab);

  const int N = 5000;
  const double du = kernel_gamma / (N - 1);
  FILE *fd;

  fd = fopen("projected_kernel.txt", "w");
  fprintf(fd, "u, kernel, projected kernel\n");
  for (int i = 0; i < N; i += 1) {
    double u = i * du;
    float kernel;
    kernel_eval(u, &kernel);
    double kernel_proj = projected_kernel_eval(&tab, u);
    double kernel_proj_int = projected_kernel_integrate(u);

    fprintf(fd, "%e, %e, %e, %e\n", u, (double)kernel, kernel_proj,
            kernel_proj_int);
  }

  fclose(fd);
  projected_kernel_clean(&tab);
}
