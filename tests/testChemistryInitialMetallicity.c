/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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

/* Local headers. */
#include "swift.h"

/* System headers. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * @file testChemistryInitialMetallicity.c
 * @brief Check the metal mass of a first-initialised FVPM particle for a
 * given GEARChemistry:initial_metallicity.
 *
 * The last metal entry of the FVPM schemes holds the untracked metals and the
 * total metal mass is the sum of all entries. A value z >= 0 must give every
 * tracked element the mass fraction z and the untracked entry 0, so the total
 * is (N - 1) * z. A negative value must keep the entries read from the ICs.
 */

#if defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) || \
    defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)

#define N GEAR_CHEMISTRY_ELEMENT_COUNT

/**
 * @brief Set up a particle of mass 2 and first-initialise its chemistry.
 */
static void first_init(struct part *p, struct xpart *xp,
                       const struct chemistry_global_data *cd,
                       const double *ic_fractions) {
  bzero(p, sizeof(struct part));
  bzero(xp, sizeof(struct xpart));
  p->id = 1;
  p->h = 1.f;
#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)
  p->conserved.mass = 2.f;
  p->conserved.energy = 1.f;
  hydro_first_init_part(p, xp);
#else
  p->mass = 2.f;
  p->u = 1.f;
#endif
  p->rho = 1.f;
  for (int i = 0; i < N; i++) p->chemistry_data.metal_mass[i] = ic_fractions[i];
  chemistry_first_init_part(NULL, NULL, NULL, cd, p, xp);
}

int main(int argc, char *argv[]) {

  struct chemistry_global_data cd;
  bzero(&cd, sizeof(struct chemistry_global_data));
  struct part p;
  struct xpart xp;
  int pass = 1;

  /* z >= 0: tracked entries z * mass, untracked entry 0. */
  const float z = 0.01f;
  const double ic[N] = {0.};
  chemistry_set_initial_metallicities(&cd, z);
  first_init(&p, &xp, &cd, ic);
  const double mass = hydro_get_mass(&p);
  double total = 0.;
  for (int i = 0; i < N; i++) total += p.chemistry_data.metal_mass[i];
  const double expected_total = (N - 1) * z * mass;
  for (int i = 0; i < N - 1; i++)
    pass &= (fabs(p.chemistry_data.metal_mass[i] - z * mass) < 1e-6 * mass);
  pass &= (p.chemistry_data.metal_mass[N - 1] == 0.);
  pass &= (fabs(total - expected_total) < 1e-6 * mass);
  message(
      "z=%g: total Z=%.9g (expected %.9g), untracked entry %.9g (must be 0)", z,
      total / mass, expected_total / mass,
      p.chemistry_data.metal_mass[N - 1] / mass);

  /* z < 0: the mass fractions of the ICs are kept in every entry. */
  double ic_fractions[N];
  for (int i = 0; i < N; i++) ic_fractions[i] = 0.001 * (i + 1);
  chemistry_set_initial_metallicities(&cd, -1.f);
  first_init(&p, &xp, &cd, ic_fractions);
  int ic_kept = 1;
  for (int i = 0; i < N; i++)
    ic_kept &= (fabs(p.chemistry_data.metal_mass[i] -
                     ic_fractions[i] * hydro_get_mass(&p)) <
                1e-6 * hydro_get_mass(&p));
  message("z<0: IC entries kept: %s", ic_kept ? "yes" : "NO");
  pass &= ic_kept;

  if (!pass) error("testChemistryInitialMetallicity: FAILED");
  message("testChemistryInitialMetallicity: PASSED");
  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif
