/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (schaller@strw.leidenuniv.nl).
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

#if defined(CHEMISTRY_EAGLE) && defined(COOLING_EAGLE) && defined(GADGET2_SPH)

#include "../src/cooling/EAGLE/cooling_rates.h"

/*
 * @brief Assign particle density and entropy corresponding to the
 * hydrogen number density and internal energy specified.
 *
 * @param p Particle data structure
 * @param xp extra particle structure
 * @param us unit system struct
 * @param cooling Cooling function data structure
 * @param cosmo Cosmology data structure
 * @param phys_const Physical constants data structure
 * @param nh_cgs Hydrogen number density (cgs units)
 * @param u Internal energy (cgs units)
 * @param ti_current integertime to set cosmo quantities
 */
void set_quantities(struct part *restrict p, struct xpart *restrict xp,
                    const struct unit_system *restrict us,
                    const struct cooling_function_data *restrict cooling,
                    struct cosmology *restrict cosmo,
                    const struct phys_const *restrict phys_const, float nh_cgs,
                    double u_cgs, integertime_t ti_current) {
  /* calculate density */
  double hydrogen_number_density = nh_cgs / cooling->number_density_to_cgs;
  p->rho = hydrogen_number_density * phys_const->const_proton_mass /
           p->chemistry_data.metal_mass_fraction[chemistry_element_H];

  /* update entropy based on internal energy */
  float pressure = (u_cgs)*cooling->internal_energy_from_cgs * p->rho *
                   (hydro_gamma_minus_one);
  p->entropy = pressure * (pow(p->rho, -hydro_gamma));
  xp->entropy_full = p->entropy;

  p->entropy_dt = 0.f;
}

/*
 * @brief Tests cooling integration scheme by comparing EAGLE
 * integration to subcycled explicit equation.
 */
int main(int argc, char **argv) {
  // Declare relevant structs
  struct swift_params *params = malloc(sizeof(struct swift_params));
  struct unit_system us;
  struct chemistry_global_data chem_data;
  struct part p;
  struct xpart xp;
  struct phys_const phys_const;
  struct hydro_props hydro_properties;
  struct entropy_floor_properties floor_props;
  struct pressure_floor_props pressure_floor;
  struct cooling_function_data cooling;
  struct cosmology cosmo;
  char *parametersFileName = "./testCooling.yml";

  float nh_cgs;  // hydrogen number density
  double u_cgs;  // internal energy

  const float seconds_per_year = 3.154e7;

  /* Number of values to test for in redshift,
   * hydrogen number density and internal energy */
  const int n_z = 10;
  const int n_nh = 10;
  const int n_u = 10;

  /* Number of subcycles and tolerance used to compare
   * subcycled and implicit solution. Note, high value
   * of tolerance due to mismatch between explicit and
   * implicit solution for large timesteps */
  const float integration_tolerance = 0.2;

  /* Read the parameter file */
  if (params == NULL) error("Error allocating memory for the parameter file.");
  message("Reading runtime parameters from file '%s'", parametersFileName);
  parser_read_file(parametersFileName, params);

  /* Init units */
  units_init_from_params(&us, params, "InternalUnitSystem");
  phys_const_init(&us, params, &phys_const);

  /* Init chemistry */
  chemistry_init(params, &us, &phys_const, &chem_data);
  chemistry_first_init_part(&phys_const, &us, &cosmo, &chem_data, &p, &xp);
  chemistry_part_has_no_neighbours(&p, &xp, &chem_data, &cosmo);
  chemistry_print(&chem_data);

  /* Init cosmology */
  cosmology_init(params, &us, &phys_const, &cosmo);
  cosmology_print(&cosmo);

  /* Init hydro_props */
  hydro_props_init(&hydro_properties, &phys_const, &us, params);

  /* Init entropy floor */
  entropy_floor_init(&floor_props, &phys_const, &us, &hydro_properties, params);

  /* Init the pressure floor */
  pressure_floor_init(&pressure_floor, &phys_const, &us, &hydro_properties,
                      params);

  /* Set dt */
  const int timebin = 38;
  float dt_cool, dt_therm;

  /* Init cooling */
  cooling_init(params, &us, &phys_const, &hydro_properties, &cooling);
  cooling_print(&cooling);
  cooling_update(&cosmo, &pressure_floor, &cooling, 0);

  /* Cooling function needs to know the minimal energy. Set it to the lowest
   * internal energy in the cooling table. */
  hydro_properties.minimal_internal_energy =
      exp(M_LN10 * cooling.Therm[0]) * cooling.internal_energy_from_cgs;

  /* Calculate abundance ratios */
  float *abundance_ratio;
  abundance_ratio = malloc((chemistry_element_count + 2) * sizeof(float));
  abundance_ratio_to_solar(&p, &cooling, abundance_ratio);

  /* extract mass fractions, calculate table indices and offsets */
  float XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac =
      p.chemistry_data.metal_mass_fraction[chemistry_element_He] /
      (XH + p.chemistry_data.metal_mass_fraction[chemistry_element_He]);
  int He_i;
  float d_He;
  get_index_1d(cooling.HeFrac, eagle_cooling_N_He_frac, HeFrac, &He_i, &d_He);

  /* calculate spacing in nh and u */
  const float log_u_min_cgs = 11, log_u_max_cgs = 17;
  const float log_nh_min_cgs = -6, log_nh_max_cgs = 3;
  const float delta_log_nh_cgs = (log_nh_max_cgs - log_nh_min_cgs) / n_nh;
  const float delta_log_u_cgs = (log_u_max_cgs - log_u_min_cgs) / n_u;

  /* Declare variables we will be checking */
  double du_dt_implicit, du_dt_check;
  integertime_t ti_current;

  /* Loop over values of nh and u */
  for (int nh_i = 0; nh_i < n_nh; nh_i++) {
    nh_cgs = exp(M_LN10 * log_nh_min_cgs + delta_log_nh_cgs * nh_i);
    for (int u_i = 0; u_i < n_u; u_i++) {
      u_cgs = exp(M_LN10 * log_u_min_cgs + delta_log_u_cgs * u_i);

      /* Calculate cooling solution at redshift zero if we're doing the comoving
       * check */
      /* reset quantities to nh, u, and z that we want to test */
      ti_current = max_nr_timesteps;
      cosmology_update(&cosmo, &phys_const, ti_current);
      set_quantities(&p, &xp, &us, &cooling, &cosmo, &phys_const, nh_cgs, u_cgs,
                     ti_current);

      /* Set dt */
      const integertime_t ti_step = get_integer_timestep(timebin);
      const integertime_t ti_begin =
          get_integer_time_begin(ti_current - 1, timebin);
      dt_cool = cosmology_get_delta_time(&cosmo, ti_begin, ti_begin + ti_step);
      dt_therm =
          cosmology_get_therm_kick_factor(&cosmo, ti_begin, ti_begin + ti_step);

      cooling_init(params, &us, &phys_const, &hydro_properties, &cooling);
      cooling_update(&cosmo, &pressure_floor, &cooling, 0);

      /* compute implicit solution */
      cooling_cool_part(&phys_const, &us, &cosmo, &hydro_properties,
                        &floor_props, &pressure_floor, &cooling, &p, &xp,
                        dt_cool, dt_therm, 0);
      du_dt_check = hydro_get_physical_internal_energy_dt(&p, &cosmo);

      /* Now we can test the cooling at various redshifts and compare the result
       * to the redshift zero solution */
      for (int z_i = 0; z_i <= n_z; z_i++) {
        ti_current = max_nr_timesteps / n_z * z_i + 1;

        /* reset to get the comoving density */
        cosmology_update(&cosmo, &phys_const, ti_current);
        cosmo.z = 0.f;
        set_quantities(&p, &xp, &us, &cooling, &cosmo, &phys_const,
                       nh_cgs * cosmo.a * cosmo.a * cosmo.a,
                       u_cgs / cosmo.a2_inv, ti_current);

        /* Load the appropriate tables */
        cooling_init(params, &us, &phys_const, &hydro_properties, &cooling);
        cooling_update(&cosmo, &pressure_floor, &cooling, 0);

        /* compute implicit solution */
        cooling_cool_part(&phys_const, &us, &cosmo, &hydro_properties,
                          &floor_props, &pressure_floor, &cooling, &p, &xp,
                          dt_cool, dt_therm, 0.);
        du_dt_implicit = hydro_get_physical_internal_energy_dt(&p, &cosmo);

        /* check if the two solutions are consistent */
        if (fabs((du_dt_implicit - du_dt_check) / du_dt_check) >
                integration_tolerance ||
            (du_dt_check == 0.0 && du_dt_implicit != 0.0))
          error(
              "Solutions do not match. scale factor %.5e z %.5e nh_cgs %.5e "
              "u_cgs %.5e dt (years) %.5e du cgs implicit %.5e reference %.5e "
              "error %.5e",
              cosmo.a, cosmo.z, nh_cgs, u_cgs,
              dt_cool * units_cgs_conversion_factor(&us, UNIT_CONV_TIME) /
                  seconds_per_year,
              du_dt_implicit *
                  units_cgs_conversion_factor(&us,
                                              UNIT_CONV_ENERGY_PER_UNIT_MASS) *
                  dt_therm,
              du_dt_check *
                  units_cgs_conversion_factor(&us,
                                              UNIT_CONV_ENERGY_PER_UNIT_MASS) *
                  dt_therm,
              fabs((du_dt_implicit - du_dt_check) / du_dt_check));
      }
    }
  }
  message("done comoving cooling test");

  free(params);
  return 0;
}

#else

int main(int argc, char **argv) { return 0; }

#endif
