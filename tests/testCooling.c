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
#include "../src/cooling/EAGLE/cooling_rates.h"
#include "../src/cooling/EAGLE/interpolate.h"
#include "../src/cooling/EAGLE/cooling_tables.h"

#if 1

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
  double hydrogen_number_density = nh_cgs / cooling->number_density_to_cgs;
  p->rho = hydrogen_number_density * phys_const->const_proton_mass /
           p->chemistry_data.metal_mass_fraction[chemistry_element_H] *
           (cosmo->a * cosmo->a * cosmo->a);

  /* update entropy based on internal energy */
  float pressure = (u * cosmo->a * cosmo->a) * cooling->internal_energy_from_cgs *
                   p->rho * (hydro_gamma_minus_one);
  p->entropy = pressure * (pow(p->rho, -hydro_gamma));
  xp->entropy_full = p->entropy;
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
  struct cooling_function_data cooling;
  struct cosmology cosmo;
  char *parametersFileName = "./testCooling.yml";

  float nh;  // hydrogen number density
  double u;  // internal energy

  /* switch between checking comoving cooling and subcycling */
  const int comoving_check = 1;

  /* Number of values to test for in redshift,
   * hydrogen number density and internal energy */
  const int n_z = 50;
  const int n_nh = 50;
  const int n_u = 50;
  //const int n_z = 2;
  //const int n_nh = 2;
  //const int n_u = 2;

  /* Number of subcycles and tolerance used to compare
   * subcycled and implicit solution. Note, high value
   * of tolerance due to mismatch between explicit and
   * implicit solution for large timesteps */
  const int n_subcycle = 10000;
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
  chemistry_part_has_no_neighbours(&p,&xp,&chem_data,&cosmo);
  chemistry_print(&chem_data);

  /* Init cosmology */
  cosmology_init(params, &us, &phys_const, &cosmo);
  cosmology_print(&cosmo);
  
  /* Init hydro_props */
  struct hydro_props hydro_properties;
  hydro_props_init(&hydro_properties, &phys_const, &us, params);

  /* Init cooling */
  cooling_init(params, &us, &phys_const, &hydro_properties, &cooling);
  cooling_print(&cooling);
  cooling_update(&cosmo, &cooling, 0);

  /* Init entropy floor */
  struct entropy_floor_properties floor_props;
  entropy_floor_init(&floor_props, &phys_const, &us, &hydro_properties, params);
  
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
  const float delta_nh = (cooling.nH[eagle_cooling_N_density - 1] - cooling.nH[0]) / n_nh;
  const float delta_u =
      (cooling.Therm[eagle_cooling_N_temperature - 1] - cooling.Therm[0]) / n_u;

  /* Declare variables we will be checking */
  double u_implicit_cgs, u_check_cgs;
  integertime_t ti_current;

  for (int nh_i = 0; nh_i < n_nh; nh_i++) {
    nh = exp(M_LN10 * cooling.nH[0] + delta_nh * nh_i);
    for (int u_i = 0; u_i < n_u; u_i++) {
      u = exp(M_LN10 * cooling.Therm[0] + delta_u * u_i);
      //u = exp(M_LN10 * cooling.Therm[0] + delta_u * u_i + 3.f);

      if (comoving_check) {
          /* reset quantities to nh, u, and z that we want to test */
          ti_current = 0;
          set_quantities(&p, &xp, &us, &cooling, &cosmo, &phys_const, nh, u,
                         ti_current);

          /* compute implicit solution */
          cooling_cool_part(&phys_const, &us, &cosmo, &hydro_properties, &floor_props, 
                            &cooling, &p, &xp, dt_cool, dt_therm);
          u_check_cgs =
              hydro_get_physical_internal_energy(&p, &xp, &cosmo) *
              cooling.internal_energy_to_cgs;
      }
      for (int z_i = 0; z_i < n_z; z_i++) {
        ti_current = max_nr_timesteps / n_z * z_i;

        if (!comoving_check) {
          /* update nh, u, z */
          set_quantities(&p, &xp, &us, &cooling, &cosmo, &phys_const, nh, u,
                         ti_current);

          /* calculate subcycled solution */
          for (int t_subcycle = 0; t_subcycle < n_subcycle; t_subcycle++) {
            p.entropy_dt = 0;
            cooling_cool_part(&phys_const, &us, &cosmo, &hydro_properties,
                              &floor_props, &cooling, &p, &xp, dt_cool / n_subcycle,
                              dt_therm / n_subcycle);
            xp.entropy_full += p.entropy_dt * dt_therm / n_subcycle;
          }
          u_check_cgs =
              hydro_get_physical_internal_energy(&p, &xp, &cosmo) *
              cooling.internal_energy_to_cgs;
          
          /* reset quantities to nh, u, and z that we want to test */
          set_quantities(&p, &xp, &us, &cooling, &cosmo, &phys_const, nh, u,
                         ti_current);

          /* compute implicit solution */
          cooling_cool_part(&phys_const, &us, &cosmo, &hydro_properties, &floor_props, 
                            &cooling, &p, &xp, dt_cool, dt_therm);
          u_implicit_cgs =
              hydro_get_physical_internal_energy(&p, &xp, &cosmo) *
              cooling.internal_energy_to_cgs;

        } else {
          /* use set_quantities to set the redshift (and scalefactor) */
          set_quantities(&p, &xp, &us, &cooling, &cosmo, &phys_const, nh, u,
                         ti_current);
          float density_1_physical = hydro_get_physical_density(&p, &cosmo);
          float density_1_comoving = hydro_get_comoving_density(&p);
          float internal_energy_1_physical = hydro_get_physical_internal_energy(&p, &xp, &cosmo);
          float internal_energy_1_comoving = hydro_get_comoving_internal_energy(&p, &xp);

          /* reset to get the comoving density */
          set_quantities(&p, &xp, &us, &cooling, &cosmo, &phys_const, nh * cosmo.a*cosmo.a*cosmo.a, u / cosmo.a2_inv,
                         ti_current);
          float density_2_physical = hydro_get_physical_density(&p, &cosmo);
          float density_2_comoving = hydro_get_comoving_density(&p);
          float internal_energy_2_physical = hydro_get_physical_internal_energy(&p, &xp, &cosmo);
          float internal_energy_2_comoving = hydro_get_comoving_internal_energy(&p, &xp);

          message("scale factor %.5e density 1 2 physical comoving %.5e %.5e %.5e %.5e", cosmo.a, density_1_physical, density_1_comoving, density_2_physical, density_2_comoving);
          message("scale factor %.5e internal_energy 1 2 physical comoving %.5e %.5e %.5e %.5e", cosmo.a, internal_energy_1_physical, internal_energy_1_comoving, internal_energy_2_physical, internal_energy_2_comoving);

          /* compute implicit solution */
          cooling_cool_part(&phys_const, &us, &cosmo, &hydro_properties, &floor_props, 
                            &cooling, &p, &xp, dt_cool, dt_therm);
          u_implicit_cgs =
              hydro_get_physical_internal_energy(&p, &xp, &cosmo) *
              cooling.internal_energy_to_cgs;
        }

        /* check if the two solutions are consistent */
        if (fabs((u_implicit_cgs - u_check_cgs) / u_check_cgs) >
            integration_tolerance) {
          message(
              "implicit and reference solutions do not match. z_i %d nh_i %d "
              "u_i %d implicit %.5e reference %.5e error %.5e",
              z_i, nh_i, u_i, u_implicit_cgs, u_check_cgs,
              fabs((u_implicit_cgs - u_check_cgs) / u_check_cgs));
        } else {
          //message(
          //    "implicit and reference solutions match. z_i %d nh_i %d "
          //    "u_i %d implicit %.5e reference %.5e error %.5e",
          //    z_i, nh_i, u_i, u_implicit_cgs, u_check_cgs,
          //    fabs((u_implicit_cgs - u_check_cgs) / u_check_cgs));
        }
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
