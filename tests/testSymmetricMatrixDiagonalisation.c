/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2025 Darwin Roduit (darwin.roduit@epfl.ch).
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
#include <time.h>

/* Local headers */
#include "src/chemistry/GEAR_MF_DIFFUSION/chemistry_utils.h"
#include "const.h"
#include "clocks.h"
#include "dimension.h"
#include "error.h"
#include "tools.h"

/**
 * @brief Verification logic for a diagonalized matrix.
 * Checks Orthonormality (V^T V = I) and Eigenpair validity (Av = lambda v).
 */
void verify_diagonalization(const double A[3][3], const double eigenvalues[3],
			    const double ev0[3], const double ev1[3],
			    const double ev2[3], const char* test_name) {

  const double tolerance = 1e-9;
  const double *vecs[3] = {ev0, ev1, ev2};

  for (int i = 0; i < 3; i++) {
    /* 1. Check Orthonormality */
    for (int j = 0; j < 3; j++) {
      double dot = chemistry_utils_dot_product(vecs[i], vecs[j]);
      double target = (i == j) ? 1.0 : 0.0;
      if (fabs(dot - target) > tolerance) {
	error("[%s] Orthonormality failed: dot(ev%d, ev%d) = %.10e", test_name, i, j, dot);
      }
    }

    /* 2. Check Reconstruction: A * v_i = lambda_i * v_i */
    double Avi[3];
    for (int row = 0; row < 3; row++) {
      Avi[row] = A[row][0] * vecs[i][0] + A[row][1] * vecs[i][1] + A[row][2] * vecs[i][2];
    }

    for (int row = 0; row < 3; row++) {
      double lambda_vi = eigenvalues[i] * vecs[i][row];
      
      /* Get the maximum scale of the matrix/eigenvalues */
      double max_val = 0.0;
      for(int k=0; k<3; k++) {
	for(int l=0; l<3; l++) max_val = max(max_val, fabs(A[k][l]));
      }

      /* Use a relative tolerance based on the matrix scale. This allows for
       * precision loss proportional to the dynamic range */      
      double allowed_diff = tolerance * max(1.0, max_val) * 100.0; 

      if (fabs(Avi[row] - lambda_vi) > allowed_diff) {
        error("[%s] Eigenpair failed: i=%d, row=%d, diff=%e (allowed=%e)", 
              test_name, i, row, Avi[row] - lambda_vi, allowed_diff);
      }
    }    
  }
}

/**
 * Generate a random symmetric 3x3 matrix.
 */
void setup_symmetric_matrix(double A[3][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      double val = ((double)rand() / (double)RAND_MAX) * 2.0 - 1.0;
      A[i][j] = val;
      A[j][i] = val;
    }
  }
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency */
  clocks_set_cpufreq(0);

#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  double A[3][3], evals[3], ev0[3], ev1[3], ev2[3];

  /* --- Part 1: Selected Simmple Stress Tests --- */

  /* Identity Matrix (Triple Multiplicity) */
  double I[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  chemistry_utils_diagonalize_3x3(I, evals, ev0, ev1, ev2);
  verify_diagonalization(I, evals, ev0, ev1, ev2, "Identity");

  message("%e %e %e | %e %e %e | %e %e %e | %e %e %e", evals[0], evals[1],
	  evals[2], ev0[0], ev0[1], ev0[2],
	  ev1[0], ev1[1], ev1[2], ev2[0], ev2[1], ev2[2]);

  /* Double Multiplicity (Eigenvalues 2, 4, 5) */
  double M_mult[3][3] = {{3, -1, 0}, {-1, 3, 0}, {0, 0, 5}};
  chemistry_utils_diagonalize_3x3(M_mult, evals, ev0, ev1, ev2);
  verify_diagonalization(M_mult, evals, ev0, ev1, ev2, "Multiplicity");

  message("%e %e %e | %e %e %e | %e %e %e | %e %e %e", evals[0], evals[1],
	  evals[2], ev0[0], ev0[1], ev0[2],
	  ev1[0], ev1[1], ev1[2], ev2[0], ev2[1], ev2[2]);

  /* Preconditioning Test (Large Numbers) */
  double large = 1e12;
  double M_large[3][3] = {{large, 0, 0}, {0, large, 0}, {0, 0, large}};
  chemistry_utils_diagonalize_3x3(M_large, evals, ev0, ev1, ev2);
  verify_diagonalization(M_large, evals, ev0, ev1, ev2, "Large Numbers 1");

  message("%e %e %e | %e %e %e | %e %e %e | %e %e %e", evals[0], evals[1],
	  evals[2], ev0[0], ev0[1], ev0[2], ev1[0], ev1[1], ev1[2], ev2[0],
	  ev2[1], ev2[2]);

  double M_large2[3][3] = {{3*large, -1*large, 0}, {-1*large, 3*large, 0}, {0, 0, 5*large}};
  chemistry_utils_diagonalize_3x3(M_large2, evals, ev0, ev1, ev2);
  verify_diagonalization(M_large2, evals, ev0, ev1, ev2, "Large Numbers 2");

  message("%e %e %e | %e %e %e | %e %e %e | %e %e %e", evals[0], evals[1],
	  evals[2], ev0[0], ev0[1], ev0[2],
	  ev1[0], ev1[1], ev1[2], ev2[0], ev2[1], ev2[2]);

  double M_large3[3][3] = {{large, large*0.2, 0}, {large*0.2, large*0.8, 0}, {0, 0, large*0.5}};
  chemistry_utils_diagonalize_3x3(M_large3, evals, ev0, ev1, ev2);
  verify_diagonalization(M_large3, evals, ev0, ev1, ev2, "Large Numbers 3");

  message("%e %e %e | %e %e %e | %e %e %e | %e %e %e", evals[0], evals[1],
	  evals[2], ev0[0], ev0[1], ev0[2],
	  ev1[0], ev1[1], ev1[2], ev2[0], ev2[1], ev2[2]);

  /* Zero Matrix */
  double Z[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  chemistry_utils_diagonalize_3x3(Z, evals, ev0, ev1, ev2);
  verify_diagonalization(Z, evals, ev0, ev1, ev2, "Zero Matrix");

  message("%e %e %e | %e %e %e | %e %e %e | %e %e %e", evals[0], evals[1],
          evals[2], ev0[0], ev0[1], ev0[2], ev1[0], ev1[1], ev1[2], ev2[0],
          ev2[1], ev2[2]);

  /* --- Part 1.5: Advanced Edge Cases --- */

  /* Nearly Diagonal (tests the threshold of the 'norm > 0' check) */
  double M_near_diag[3][3] = {{1.0, 1e-15, 0}, {1e-15, 2.0, 0}, {0, 0, 3.0}};
  chemistry_utils_diagonalize_3x3(M_near_diag, evals, ev0, ev1, ev2);
  verify_diagonalization(M_near_diag, evals, ev0, ev1, ev2, "Nearly Diagonal");

  message("%e %e %e | %e %e %e | %e %e %e | %e %e %e", evals[0], evals[1],
          evals[2], ev0[0], ev0[1], ev0[2], ev1[0], ev1[1], ev1[2], ev2[0],
          ev2[1], ev2[2]);  

  /* Rotated High-Contrast Matrix */
  /* This matrix has eigenvalues 1, 1e12, 1e12 but is not axis-aligned. Given
     its high dynamic range, floating points errors are expected to slithly
     alter the solution. However, the verification takes this into account. */
  double L = 1e12;
  double M_rotated[3][3] = {
      { (L+1.0)*0.5, (L-1.0)*0.5, 0.0 },
      { (L-1.0)*0.5, (L+1.0)*0.5, 0.0 },
      { 0.0,         0.0,         L   }
  };
  chemistry_utils_diagonalize_3x3(M_rotated, evals, ev0, ev1, ev2);
  verify_diagonalization(M_rotated, evals, ev0, ev1, ev2,
                         "Rotated High-Contrast");

  message("%e %e %e | %e %e %e | %e %e %e | %e %e %e", evals[0], evals[1],
          evals[2], ev0[0], ev0[1], ev0[2], ev1[0], ev1[1], ev1[2], ev2[0],
          ev2[1], ev2[2]);  

  /* The "Negative Eigenvalue" case */
  /* Ensure the acos logic handles cases where trace/3 is negative */
  double M_neg[3][3] = {{-2.0, 1.0, 0.0}, {1.0, -2.0, 1.0}, {0.0, 1.0, -2.0}};
  chemistry_utils_diagonalize_3x3(M_neg, evals, ev0, ev1, ev2);
  verify_diagonalization(M_neg, evals, ev0, ev1, ev2, "Negative Eigenvalues");

  message("%e %e %e | %e %e %e | %e %e %e | %e %e %e", evals[0], evals[1],
          evals[2], ev0[0], ev0[1], ev0[2], ev1[0], ev1[1], ev1[2], ev2[0],
          ev2[1], ev2[2]);

  /* Filamentary Case: Eigenvalues (1e10, 1e10, 1e-5) */
  /* This tests if the solver can find a very small eigenvalue in the presence 
   * of very large ones without losing the eigenvector direction. */
  double M_planar[3][3] = {{1e10, 0, 0}, {0, 1e10, 0}, {0, 0, 1e-5}};
  chemistry_utils_diagonalize_3x3(M_planar, evals, ev0, ev1, ev2);
  verify_diagonalization(M_planar, evals, ev0, ev1, ev2, "Planar/Filament");

  message("%e %e %e | %e %e %e | %e %e %e | %e %e %e", evals[0], evals[1],
          evals[2], ev0[0], ev0[1], ev0[2], ev1[0], ev1[1], ev1[2], ev2[0],
          ev2[1], ev2[2]);  

  /* --- Part 2: Randomized Monte Carlo Tests --- */

  for (int test = 0; test < 10000; ++test) {
    setup_symmetric_matrix(A);
    chemistry_utils_diagonalize_3x3(A, evals, ev0, ev1, ev2);
    verify_diagonalization(A, evals, ev0, ev1, ev2, "Random Matrix");
  }

  /* --- Part 3: Extreme Dynamic Range Monte Carlo --- */
  for (int test = 0; test < 1000; ++test) {
    double scale = pow(10.0, ((double)rand() / (double)RAND_MAX) * 20.0 - 10.0); // 1e-10 to 1e10
    for (int i = 0; i < 3; i++) {
      for (int j = i; j < 3; j++) {
	double val = (((double)rand() / (double)RAND_MAX) * 2.0 - 1.0) * scale;
	A[i][j] = val;
	A[j][i] = val;
      }
    }
    chemistry_utils_diagonalize_3x3(A, evals, ev0, ev1, ev2);
    verify_diagonalization(A, evals, ev0, ev1, ev2, "Extreme Random Matrix");
  }

  /* --- Part 4: Symmetry Sensitivity --- */
  /* Test with a matrix that is *almost symmetric but has 1 epsilon
   * difference. */  
  double M_nosym[3][3] = {{1.0, 0.2, 0.3}, {0.200000000000001, 1.0, 0.1}, {0.3, 0.1, 1.0}};
  chemistry_utils_diagonalize_3x3(M_nosym, evals, ev0, ev1, ev2);
  verify_diagonalization(M_nosym, evals, ev0, ev1, ev2, "Near-Symmetry");

  message("%e %e %e | %e %e %e | %e %e %e | %e %e %e", evals[0], evals[1],
          evals[2], ev0[0], ev0[1], ev0[2], ev1[0], ev1[1], ev1[2], ev2[0],
          ev2[1], ev2[2]);

  /* --- Part 5: Stability over time --- */
  /* Run a very tight loop to ensure there are no memory leaks or 
   * instruction-level cache issues (mostly for profiling). */
  ticks tic = getticks();
  for (int i = 0; i < 100000; i++) {
    chemistry_utils_diagonalize_3x3(M_rotated, evals, ev0, ev1, ev2);
  }
  ticks toc = getticks();
  message("Performance: 100k diagonalizations in %.3f ms", 
          clocks_from_ticks(toc - tic) * 1000.0);  

  message("All 11,000 random tests and specialized cases passed!");
  return 0;
}
