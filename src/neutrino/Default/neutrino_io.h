/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Willem Elbers (willem.h.elbers@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_NEUTRINO_IO_H
#define SWIFT_DEFAULT_NEUTRINO_IO_H

#include "../../gravity.h"

/* Local includes */
#include "fermi_dirac.h"
#include "neutrino.h"
#include "neutrino_properties.h"

/**
 * @brief Recover and store the initial Fermi-Dirac speed, vi, for a neutrino
 * particle. These are used in the weight calculation of the delta-f method.
 *
 * @param e The engine of the run
 * @param gp The neutrino gpart in question
 * @param ret Output
 */
INLINE static void convert_gpart_vi(const struct engine* e,
                                    const struct gpart* gp, float* ret) {

  /* When we are running with the delta-f method, resample the momentum */
  if (e->neutrino_properties->use_delta_f) {
    /* Retrieve physical constants, including the neutrino mass array */
    const double a_scale = e->cosmology->a;
    const double N_nu = e->cosmology->N_nu;
    const double* m_eV_array = e->cosmology->M_nu_eV;
    const double c_vel = e->physical_constants->const_speed_light_c;
    const double T_eV = e->cosmology->T_nu_0_eV;
    const double a_fac = c_vel * T_eV / a_scale;

    /* Use a particle id dependent seed (sum of global seed and ID) */
    const long long neutrino_seed = e->neutrino_properties->neutrino_seed;
    const long long seed = gp->id_or_neg_offset + neutrino_seed;

    /* Convert momentum in electronvolts to speed in internal units */
    double pi_eV = neutrino_seed_to_fermi_dirac(seed);
    double m_eV = neutrino_seed_to_mass(N_nu, m_eV_array, seed);
    double vi = pi_eV / m_eV * a_fac;  // scales like peculiar velocity

    ret[0] = vi;
  } else {
    /* We don't know what the initial momentum was and we don't need it */
    ret[0] = 0.f;
  }
}

/**
 * @brief Recover and store the microscopic mass of a neutrino particle
 *
 * @param e The engine of the run
 * @param gp The neutrino gpart in question
 * @param ret Output
 */
INLINE static void convert_gpart_mnu(const struct engine* e,
                                     const struct gpart* gp, double* ret) {

  /* Physical constants */
  const double c = e->physical_constants->const_speed_light_c;
  const double eV = e->physical_constants->const_electron_volt;
  const double eV_mass = eV / (c * c);

  double micro_mass;

  /* Resample if running with the delta-f method or neutrino ic generation */
  if (e->neutrino_properties->use_delta_f ||
      e->neutrino_properties->generate_ics) {

    /* Use a particle id dependent seed (sum of global seed and ID) */
    const long long neutrino_seed = e->neutrino_properties->neutrino_seed;
    const long long seed = gp->id_or_neg_offset + neutrino_seed;

    /* Fetch neutrino masses defined in the cosmology */
    const int N_nu = e->cosmology->N_nu;
    const double* m_eV_array = e->cosmology->M_nu_eV;

    micro_mass = neutrino_seed_to_mass(N_nu, m_eV_array, seed);  // eV
  } else {
    /* Otherwise, simply use the mass implied by the conversion factor and
     * total degeneracy */
    const double deg_nu_tot = e->cosmology->deg_nu_tot;
    const double mass_factor = e->neutrino_mass_conversion_factor;

    micro_mass = gp->mass * mass_factor / deg_nu_tot;  // eV
  }

  /* Convert units and store the answer */
  ret[0] = micro_mass * eV_mass;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param gparts The particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int neutrino_write_particles(
    const struct gpart* gparts, struct io_props* list) {

  list[0] = io_make_output_field_convert_gpart(
      "SampledSpeeds", FLOAT, 1, UNIT_CONV_SPEED, 0.f, gparts, convert_gpart_vi,
      "Initial Fermi-Dirac speed sampled at infinity. This is a * |dx/dt| "
      "where x is the co-moving position of the particles.");

  list[1] = io_make_output_field_convert_gpart(
      "MicroscopicMasses", DOUBLE, 1, UNIT_CONV_MASS, 0.f, gparts,
      convert_gpart_mnu, "Microscopic masses of individual neutrino particles");

  return 2;
}

#endif /* SWIFT_DEFAULT_NEUTRINO_IO_H */
