/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Willem Elbers (willem.h.elbers@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_NEUTRINO_H
#define SWIFT_DEFAULT_NEUTRINO_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "../../engine.h"
#include "fermi_dirac.h"
#include "neutrino_properties.h"
#include "neutrino_response.h"
#include "relativity.h"

/* Riemann function zeta(3) */
#define M_ZETA_3 1.2020569031595942853997

/**
 * @brief Shared information for delta-f neutrino weighting of a cell.
 */
struct neutrino_model {
  char use_delta_f_mesh_only;
  double *M_nu_eV;
  double *deg_nu;
  int N_nu;
  double fac;
  double inv_mass_factor;
  long long neutrino_seed;
};

void gather_neutrino_consts(const struct space *s, struct neutrino_model *nm);
void gpart_neutrino_weight_mesh_only(const struct gpart *gp,
                                     const struct neutrino_model *nm,
                                     double *weight);
void gpart_neutrino_mass_weight(const struct gpart *gp,
                                const struct neutrino_model *nm, double *mass,
                                double *weight);

/* Compute the ratio of macro particle mass in internal mass units to
 * the mass of one microscopic neutrino in eV.
 *
 * @param cosmo The #cosmology used for this run.
 * @param internal_units The system of units used internally.
 * @param physical_constants The #phys_const used for this run.
 * @param volume The volume occupied by neutrino particles
 * @param nr_nuparts The number of macro neutrino particles
 */
INLINE static double neutrino_mass_factor(
    const struct cosmology *cosmo, const struct unit_system *internal_units,
    const struct phys_const *physical_constants, double volume,
    double nr_nuparts) {
  /* Some constants */
  const double k_b = physical_constants->const_boltzmann_k;
  const double hbar = physical_constants->const_planck_hbar;
  const double c = physical_constants->const_speed_light_c;
  const double eV = physical_constants->const_electron_volt;
  const double eV_mass = eV / (c * c);  // 1 eV/c^2 in internal mass units
  const double prefactor = (1.5 * M_ZETA_3) / (M_PI * M_PI);
  const double T_nu = cosmo->T_nu_0;
  const int N_nu = cosmo->N_nu;

  /* Compute the comoving number density per flavour */
  const double kThc = k_b * T_nu / (hbar * c);
  const double n = prefactor * kThc * kThc * kThc;

  /* Compute the conversion factor per flavour */
  const double mass_factor = nr_nuparts / (n * volume * N_nu);

  /* Convert to eV */
  const double mass_factor_eV = mass_factor / eV_mass;

  return mass_factor_eV;
}

/**
 * @brief Initialises the neutrino g-particles for the first time. This is done
 * in addition to gravity_first_init_gpart().
 *
 * This function is called only once just after the ICs have been read in
 * and after IDs have been remapped (if used) by space_remap_ids().
 *
 * @param gp The particle to act upon
 * @param engine The engine of the run
 */
__attribute__((always_inline)) INLINE static void gravity_first_init_neutrino(
    struct gpart *gp, const struct engine *e) {

  /* Do we need to do anything? */
  if (!e->neutrino_properties->generate_ics) return;

  /* Retrieve physical and cosmological constants */
  const double c_vel = e->physical_constants->const_speed_light_c;
  const double *m_eV_array = e->cosmology->M_nu_eV;
  const double *deg_array = e->cosmology->deg_nu;
  const int N_nu = e->cosmology->N_nu;
  const double T_eV = e->cosmology->T_nu_0_eV;
  const double inv_fac = c_vel * T_eV;
  const double inv_mass_factor = 1. / e->neutrino_mass_conversion_factor;
  const long long neutrino_seed = e->neutrino_properties->neutrino_seed;

  /* Use a particle id dependent seed */
  const long long seed = gp->id_or_neg_offset + neutrino_seed;

  /* Compute the initial dimensionless momentum from the seed */
  const double pi = neutrino_seed_to_fermi_dirac(seed);

  /* The neutrino mass and degeneracy (we cycle based on the neutrino seed) */
  const double m_eV = neutrino_seed_to_mass(N_nu, m_eV_array, seed);
  const double deg = neutrino_seed_to_degeneracy(N_nu, deg_array, seed);

  /* Compute the initial direction of the momentum vector from the seed */
  double n[3];
  neutrino_seed_to_direction(seed, n);

  /* Set the initial velocity */
  const double vi = pi * inv_fac / m_eV;
  gp->v_full[0] = n[0] * vi;
  gp->v_full[1] = n[1] * vi;
  gp->v_full[2] = n[2] * vi;

  /* If running with the delta-f method, set the weight to (almost) zero */
  if (e->neutrino_properties->use_delta_f) {
    gp->mass = FLT_MIN;
  }
  /* Otherwise, set the mass based on the mass factor */
  else {
    gp->mass = deg * m_eV * inv_mass_factor;
  }
}

void compute_neutrino_diagnostics(
    const struct space *s, const struct cosmology *cosmo,
    const struct phys_const *physical_constants,
    const struct neutrino_props *neutrino_properties, const int rank, double *r,
    double *I_df, double *mass_tot);
void neutrino_check_cosmology(const struct space *s,
                              const struct cosmology *cosmo,
                              const struct phys_const *physical_constants,
                              struct swift_params *params,
                              const struct neutrino_props *neutrino_props,
                              const int rank, const int verbose);
double lightcone_map_neutrino_baseline_value(
    const struct cosmology *c, const struct lightcone_props *lightcone_props,
    const struct lightcone_map *map);
#endif /* SWIFT_DEFAULT_NEUTRINO_H */
