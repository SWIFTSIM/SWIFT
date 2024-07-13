/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
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

/* This object's header. */
#include "neutrino.h"

/* Standard headers */
#include <math.h>

/* Local includes */
#include "lightcone/lightcone.h"
#include "lightcone/lightcone_map_types.h"

/* Compute the dimensionless neutrino momentum (units of kb*T).
 *
 * @param v The internal 3-velocity
 * @param m_eV The neutrino mass in electron-volts
 * @param fac Conversion factor = 1. / (speed_of_light * T_nu_eV)
 */
INLINE static double neutrino_momentum(const float v[3], const double m_eV,
                                       const double fac) {

  float v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  float vmag = sqrtf(v2);
  double p = vmag * fac * m_eV;
  return p;
}

/**
 * @brief Gather neutrino constants
 *
 * @param s The #space for this run.
 * @param nm Struct with neutrino constants
 */
void gather_neutrino_consts(const struct space *s, struct neutrino_model *nm) {
  nm->use_delta_f_mesh_only = s->e->neutrino_properties->use_delta_f_mesh_only;
  nm->M_nu_eV = s->e->cosmology->M_nu_eV;
  nm->deg_nu = s->e->cosmology->deg_nu;
  nm->N_nu = s->e->cosmology->N_nu;
  nm->fac = 1.0 / (s->e->physical_constants->const_speed_light_c *
                   s->e->cosmology->T_nu_0_eV);
  nm->inv_mass_factor = 1. / s->e->neutrino_mass_conversion_factor;
  nm->neutrino_seed = s->e->neutrino_properties->neutrino_seed;
}

/**
 * @brief Compute delta-f weight of a neutrino particle, but *only* when using
 * the delta-f method exclusively on the mesh (otherwise the mass is already
 * weighted).
 *
 * @param gp The #gpart.
 * @param nm Properties of the neutrino model
 * @param weight The resulting weight (output)
 */
void gpart_neutrino_weight_mesh_only(const struct gpart *gp,
                                     const struct neutrino_model *nm,
                                     double *weight) {
  /* Anything to do? */
  if (!nm->use_delta_f_mesh_only) return;

  /* Use a particle id dependent seed */
  const long long seed = gp->id_or_neg_offset + nm->neutrino_seed;

  /* Compute the initial dimensionless momentum from the seed */
  const double pi = neutrino_seed_to_fermi_dirac(seed);

  /* The neutrino mass and degeneracy (we cycle based on the seed) */
  const double m_eV = neutrino_seed_to_mass(nm->N_nu, nm->M_nu_eV, seed);

  /* Compute the current dimensionless momentum */
  double p = neutrino_momentum(gp->v_full, m_eV, nm->fac);

  /* Compute the initial and current background phase-space density */
  double fi = fermi_dirac_density(pi);
  double f = fermi_dirac_density(p);
  *weight = 1.0 - f / fi;
}

/**
 * @brief Compute the mass and delta-f weight of a neutrino particle
 *
 * @param gp The #gpart.
 * @param nm Properties of the neutrino model
 * @param mass The mass (output)
 * @param weight The resulting weight (output)
 */
void gpart_neutrino_mass_weight(const struct gpart *gp,
                                const struct neutrino_model *nm, double *mass,
                                double *weight) {

  /* Use a particle id dependent seed */
  const long long seed = gp->id_or_neg_offset + nm->neutrino_seed;

  /* Compute the initial dimensionless momentum from the seed */
  const double pi = neutrino_seed_to_fermi_dirac(seed);

  /* The neutrino mass and degeneracy (we cycle based on the seed) */
  const double m_eV = neutrino_seed_to_mass(nm->N_nu, nm->M_nu_eV, seed);
  const double deg = neutrino_seed_to_degeneracy(nm->N_nu, nm->deg_nu, seed);
  *mass = deg * m_eV * nm->inv_mass_factor;

  /* Compute the current dimensionless momentum */
  const double p = neutrino_momentum(gp->v_full, m_eV, nm->fac);

  /* Compute the initial and current background phase-space density */
  const double fi = fermi_dirac_density(pi);
  const double f = fermi_dirac_density(p);
  *weight = 1.0 - f / fi;
}

/**
 * @brief Compute diagnostics for the neutrino delta-f method, including
 * the mean squared weight.
 *
 * @param s The #space.
 * @param cosmo The current cosmology model.
 * @param physical_constants The #phys_const used for this run.
 * @param neutrino_props The #neutrino_props used for this run.
 * @param rank The MPI rank of this #space.
 * @param r Output: correlation coefficient between current and sampled momenta
 * @param I_df Output: half the mean squared weight
 * @param mass_tot Output: the total mass in neutrino particles
 */
void compute_neutrino_diagnostics(
    const struct space *s, const struct cosmology *cosmo,
    const struct phys_const *physical_constants,
    const struct neutrino_props *neutrino_properties, const int rank, double *r,
    double *I_df, double *mass_tot) {

  int use_df = neutrino_properties->use_delta_f;
  int use_df_mesh = neutrino_properties->use_delta_f_mesh_only;

  if (!use_df && !use_df_mesh) {
    error("Neutrino diagnostics only defined when using the delta-f method.");
  }

  struct gpart *gparts = s->gparts;
  const size_t nr_gparts = s->nr_gparts;

  /* Retrieve physical and cosmological constants */
  const double c_vel = physical_constants->const_speed_light_c;
  const double *m_eV_array = cosmo->M_nu_eV;
  const double *deg_array = cosmo->deg_nu;
  const int N_nu = cosmo->N_nu;
  const double T_eV = cosmo->T_nu_0_eV;
  const double fac = 1.0 / (c_vel * T_eV);
  const double inv_mass_factor = 1. / s->e->neutrino_mass_conversion_factor;
  const long long total_nr_neutrinos = s->e->total_nr_neutrino_gparts;
  const long long neutrino_seed = neutrino_properties->neutrino_seed;

  /* Sum up the masses, weights, and momenta for the neutrinos in this space */
  double mass_sum = 0;
  double weight2_sum = 0;
  double p_sum = 0;  // current momenta
  double p2_sum = 0;
  double pi_sum = 0;  // sampled momenta
  double pi2_sum = 0;
  double ppi_sum = 0;
  for (size_t i = 0; i < nr_gparts; ++i) {

    /* Skip extra and non-neutrino particles */
    if (gparts[i].time_bin == time_bin_not_created) continue;
    if (gparts[i].type != swift_type_neutrino) continue;

    /* Use a particle id dependent seed */
    const long long seed = gparts[i].id_or_neg_offset + neutrino_seed;

    /* Compute the initial dimensionless momentum from the seed */
    const double pi = neutrino_seed_to_fermi_dirac(seed);

    /* The neutrino mass and degeneracy (we cycle based on the seed) */
    const double m_eV = neutrino_seed_to_mass(N_nu, m_eV_array, seed);
    const double deg = neutrino_seed_to_degeneracy(N_nu, deg_array, seed);
    const double mass = deg * m_eV * inv_mass_factor;

    /* Compute the current dimensionless momentum */
    double v2 = gparts[i].v_full[0] * gparts[i].v_full[0] +
                gparts[i].v_full[1] * gparts[i].v_full[1] +
                gparts[i].v_full[2] * gparts[i].v_full[2];
    double p = sqrt(v2) * m_eV * fac;

    /* Compute the initial and current background phase-space density */
    double fi = fermi_dirac_density(pi);
    double f = fermi_dirac_density(p);
    double weight = 1.0 - f / fi;

    p_sum += p;
    p2_sum += p * p;
    pi_sum += pi;
    pi2_sum += pi * pi;
    ppi_sum += p * pi;
    mass_sum += mass;
    weight2_sum += weight * weight;
  }

/* Reduce the total mass, weights, and momenta */
#ifdef WITH_MPI
  double sums[7] = {p_sum,   p2_sum,   pi_sum,     pi2_sum,
                    ppi_sum, mass_sum, weight2_sum};
  double total_sums[7];

  MPI_Reduce(sums, total_sums, 7, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  double total_p = total_sums[0];
  double total_p2 = total_sums[1];
  double total_pi = total_sums[2];
  double total_pi2 = total_sums[3];
  double total_ppi = total_sums[4];
  double total_mass = total_sums[5];
  double total_weight2 = total_sums[6];
#else
  double total_p = p_sum;
  double total_p2 = p2_sum;
  double total_pi = pi_sum;
  double total_pi2 = pi2_sum;
  double total_ppi = ppi_sum;
  double total_mass = mass_sum;
  double total_weight2 = weight2_sum;
#endif

  if (rank == 0) {
    /* Compute the correlation coefficient */
    double Sp = total_nr_neutrinos * total_p2 - total_p * total_p;
    double Spi = total_nr_neutrinos * total_pi2 - total_pi * total_pi;
    double Sppi = sqrt(Sp * Spi);
    *r = (total_nr_neutrinos * total_ppi - total_p * total_pi) / Sppi;

    /* Half the mean squared weight */
    *I_df = 0.5 * total_weight2 / total_nr_neutrinos;

    /* The total mass in neutrino particles */
    *mass_tot = total_mass;
  }
}

/**
 * @brief Verify that the neutrino content matches the cosmological model
 * and has been correctly loaded, when using the delta-f method.
 *
 * @param s The #space.
 * @param cosmo The current cosmology model.
 * @param physical_constants The #phys_const used for this run.
 * @param params The parsed parameters.
 * @param neutrino_props The #neutrino_props used for this run.
 * @param rank The MPI rank of this #space.
 * @param with_neutrinos Are we running with neutrino particles?
 * @param verbose Are we verbose?
 */
void neutrino_check_cosmology(const struct space *s,
                              const struct cosmology *cosmo,
                              const struct phys_const *physical_constants,
                              struct swift_params *params,
                              const struct neutrino_props *neutrino_props,
                              const int rank, const int verbose) {

  /* Check that we have neutrino particles if and only if we need them */
  int use_df = neutrino_props->use_delta_f;
  int use_df_mesh = neutrino_props->use_delta_f_mesh_only;
  int use_linres = neutrino_props->use_linear_response;
  int use_none = neutrino_props->use_model_none;
  int genics = neutrino_props->generate_ics;
  int with_neutrinos = s->with_neutrinos;

  if ((use_df || use_df_mesh || genics) && !with_neutrinos) {
    error(
        "Running without neutrino particles, but specified a neutrino "
        "model that requires them.");
  } else if ((use_linres || use_none) && with_neutrinos) {
    error(
        "Running with neutrino particles, but specified a neutrino "
        "model that is incompatible with particles.");
  } else if (cosmo->Omega_nu_0 > 0. && !(use_linres || use_none) &&
             !with_neutrinos) {
    error(
        "Running without neutrino particles, but specified neutrinos "
        "in the background cosmology and not using a neutrino model that does "
        "not use particles.");
  }

  /* We are done if the delta-f method is not used, since the total mass
   * has otherwise already been checked in space_check_cosmology. */
  if (!use_df && !use_df_mesh) return;

  /* Compute neutrino diagnostics, including the total mass */
  double r, I_df, total_mass;
  compute_neutrino_diagnostics(s, cosmo, physical_constants, neutrino_props,
                               rank, &r, &I_df, &total_mass);

  if (rank == 0) {
    /* Check the correlation coefficient */
    if (r < 0.1)
      error(
          "There is no correlation between current and sampled neutrino "
          "momenta (r = %e, I = %e). Most likely, the neutrino seed is "
          "incorrect or the particle IDs have been scrambled.",
          r, I_df);

    /* Check the mean squared weight */
    else if (I_df > 0.1)
      error(
          "The neutrino particle weights are very large (r = %e, I = %e). "
          "Most likely, the particle velocities are incorrectly normalised.",
          r, I_df);

    if (verbose) message("Neutrino delta-f diagnostic: I = %e", I_df);

    /* Check the neutrino mass */
    const double volume = s->dim[0] * s->dim[1] * s->dim[2];

    /* Current Hubble constant */
    const double H = cosmo->H;

    /* z=0 Hubble parameter */
    const double H0 = cosmo->H0;

    /* Critical density at z=0 */
    const double rho_crit0 = cosmo->critical_density * H0 * H0 / (H * H);

    /* Density in neutrino particles */
    const double Omega_particles_nu = (total_mass / volume) / rho_crit0;

    if (fabs(Omega_particles_nu - cosmo->Omega_nu_0) > 1e-4)
      error(
          "The massive neutrino content of the simulation does not match the "
          "cosmology in the parameter file: cosmo.Omega_nu = %e particles "
          "Omega_nu = %e",
          cosmo->Omega_nu_0, Omega_particles_nu);
  }
}

/*
  Lightcone map of neutrino mass perturbation
*/

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_neutrino_mass_type_contributes(int ptype) {

  switch (ptype) {
    case swift_type_neutrino:
      return 1;
    default:
      return 0;
  }
}

/**
 * @brief Make a healpix map of the neutrino mass perturbation
 *
 * When a neutrino particle crosses the lightcone this function
 * should return the value to accumulate to the corresponding
 * pixel in the healpix map.
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_neutrino_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {

  switch (gp->type) {
    case swift_type_neutrino: {
      struct neutrino_model nu_model;
      bzero(&nu_model, sizeof(struct neutrino_model));
      if (e->neutrino_properties->use_delta_f_mesh_only)
        gather_neutrino_consts(e->s, &nu_model);
      double weight = 1.0;
      gpart_neutrino_weight_mesh_only(gp, &nu_model, &weight);
      return gp->mass * weight;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0; /* Prevent 'missing return' error */
  }
}

/**
 * @brief Return baseline value for neutrino mass lightcone maps.
 *
 * This is the mean neutrino density integrated over the volume of the pixel.
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param map The lightcone map
 */
double lightcone_map_neutrino_baseline_value(
    const struct cosmology *c, const struct lightcone_props *lightcone_props,
    const struct lightcone_map *map) {

  /* Fetch the area of healpix pixels */
  const double area = lightcone_props->pixel_area_steradians;

  /* Fetch the inner and outer radii */
  const double r_inner = map->r_min;
  const double r_outer = map->r_max;
  const double r_inner_3 = r_inner * r_inner * r_inner;
  const double r_outer_3 = r_outer * r_outer * r_outer;

  /* The volume mapped into a healpix pixel */
  const double volume = area * (r_outer_3 - r_inner_3) / 3.0;

  /* The mean comoving neutrino density at z = 0 */
  const double rho_nu_0 = c->critical_density_0 * c->Omega_nu_0;

  return rho_nu_0 * volume;
}
