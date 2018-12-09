/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#include "../config.h"

/* Local headers. */
#include "swift.h"

#if 0

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
                    double u, integertime_t ti_current) {

  /* Update cosmology quantities */
  cosmology_update(cosmo, phys_const, ti_current);

  /* calculate density */
  double hydrogen_number_density = nh_cgs / cooling->number_density_scale;
  p->rho = hydrogen_number_density * phys_const->const_proton_mass /
           p->chemistry_data.metal_mass_fraction[chemistry_element_H] *
           (cosmo->a * cosmo->a * cosmo->a);

  /* update entropy based on internal energy */
  float pressure = (u * cosmo->a * cosmo->a) / cooling->internal_energy_scale *
                   p->rho * (hydro_gamma_minus_one);
  p->entropy = pressure * (pow(p->rho, -hydro_gamma));
  xp->entropy_full = p->entropy;
}

/*
 * @brief Produces contributions to cooling rates for different
 * hydrogen number densities, from different metals,
 * tests 1d and 4d table interpolations produce
 * same results for cooling rate, dlambda/du and temperature.
 */
int main(int argc, char **argv) {
  // Declare relevant structs
  struct swift_params *params = malloc(sizeof(struct swift_params));
  struct unit_system us;
  struct chemistry_global_data chem_data;
  struct part p;
  struct xpart xp;
  struct phys_const phys_const;
  struct cooling_function_data cooling;
  struct cosmology cosmo;
  char *parametersFileName = "./testCooling.yml";

  float nh;  // hydrogen number density
  double u;  // internal energy

  /* Number of values to test for in redshift,
   * hydrogen number density and internal energy */
  const int n_z = 50;
  const int n_nh = 50;
  const int n_u = 50;

  /* Number of subcycles and tolerance used to compare
   * subcycled and implicit solution. Note, high value
   * of tolerance due to mismatch between explicit and
   * implicit solution for large timesteps */
  const int n_subcycle = 1000;
  const float integration_tolerance = 0.2;

  /* Set dt */
  const float dt_cool = 1.0e-5;
  const float dt_therm = 1.0e-5;

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
  chemistry_print(&chem_data);

  /* Init cosmology */
  cosmology_init(params, &us, &phys_const, &cosmo);
  cosmology_print(&cosmo);

  /* Init cooling */
  cooling_init(params, &us, &phys_const, &cooling);
  cooling_print(&cooling);
  cooling_update(&cosmo, &cooling, 0);

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
  get_index_1d(cooling.HeFrac, cooling.N_He, HeFrac, &He_i, &d_He);

  /* Cooling function needs to know the minimal energy. Set it to the lowest
   * internal energy in the cooling table. */
  struct hydro_props hydro_properties;
  hydro_properties.minimal_internal_energy =
      exp(M_LN10 * cooling.Therm[0]) / cooling.internal_energy_scale;

  /* calculate spacing in nh and u */
  const float delta_nh = (cooling.nH[cooling.N_nH - 1] - cooling.nH[0]) / n_nh;
  const float delta_u =
      (cooling.Therm[cooling.N_Temp - 1] - cooling.Therm[0]) / n_u;

  for (int z_i = 0; z_i < n_z; z_i++) {
    integertime_t ti_current = max_nr_timesteps / n_z * z_i;
    for (int nh_i = 0; nh_i < n_nh; nh_i++) {
      nh = exp(M_LN10 * cooling.nH[0] + delta_nh * nh_i);
      for (int u_i = 0; u_i < n_u; u_i++) {
        u = exp(M_LN10 * cooling.Therm[0] + delta_u * u_i);

        /* update nh, u, z */
        set_quantities(&p, &xp, &us, &cooling, &cosmo, &phys_const, nh, u,
                       ti_current);

        /* calculate subcycled solution */
        for (int t_subcycle = 0; t_subcycle < n_subcycle; t_subcycle++) {
          p.entropy_dt = 0;
          cooling_cool_part(&phys_const, &us, &cosmo, &hydro_properties,
                            &cooling, &p, &xp, dt_cool / n_subcycle,
                            dt_therm / n_subcycle);
          xp.entropy_full += p.entropy_dt * dt_therm / n_subcycle;
        }
        double u_subcycled =
            hydro_get_physical_internal_energy(&p, &xp, &cosmo) *
            cooling.internal_energy_scale;

        /* reset quantities to nh, u, and z that we want to test */
        set_quantities(&p, &xp, &us, &cooling, &cosmo, &phys_const, nh, u,
                       ti_current);

        /* compute implicit solution */
        cooling_cool_part(&phys_const, &us, &cosmo, &hydro_properties, &cooling,
                          &p, &xp, dt_cool, dt_therm);
        double u_implicit =
            hydro_get_physical_internal_energy(&p, &xp, &cosmo) *
            cooling.internal_energy_scale;

        /* check if the two solutions are consistent */
        if (fabs((u_implicit - u_subcycled) / u_subcycled) >
            integration_tolerance)
          message(
              "implicit and subcycled solutions do not match. z_i %d nh_i %d "
              "u_i %d implicit %.5e subcycled %.5e error %.5e",
              z_i, nh_i, u_i, u_implicit, u_subcycled,
              fabs((u_implicit - u_subcycled) / u_subcycled));
      }
    }
  }
  message("done test");

  free(params);
  return 0;
}

#else

int main(int argc, char **argv) { return 0; }

#endif
