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

#include <unistd.h>
#include "cooling.h"
#include "cooling_struct.h"
#include "hydro.h"
#include "physical_constants.h"
#include "swift.h"
#include "units.h"

/*
 * @brief Assign particle density and entropy corresponding to the
 * hydrogen number density and internal energy specified.
 *
 * @param p Particle data structure
 * @param cooling Cooling function data structure
 * @param cosmo Cosmology data structure
 * @param phys_const Physical constants data structure
 * @param nh Hydrogen number density (cgs units)
 * @param u Internal energy (cgs units)
 */
void set_quantities(struct part *restrict p,
		    struct xpart *restrict xp,
                    const struct unit_system *restrict us,
                    const struct cooling_function_data *restrict cooling,
                    struct cosmology *restrict cosmo,
                    const struct phys_const *restrict phys_const, float nh_cgs,
                    double u,
		    float z) {

  cosmo->z = z;
  float scale_factor = 1.0 / (1.0 + cosmo->z);
  double hydrogen_number_density =
      nh_cgs * units_cgs_conversion_factor(us, UNIT_CONV_LENGTH) * 
               units_cgs_conversion_factor(us, UNIT_CONV_LENGTH) * 
               units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  p->rho = hydrogen_number_density * phys_const->const_proton_mass /
           p->chemistry_data.metal_mass_fraction[chemistry_element_H] * (scale_factor * scale_factor * scale_factor);

  float pressure = (u * scale_factor * scale_factor) / cooling->internal_energy_scale *
                   p->rho * (hydro_gamma_minus_one);
  p->entropy = pressure * (pow(p->rho, -hydro_gamma));
  xp->entropy_full = p->entropy;

  // Using hydro_set_init_internal_energy seems to work better for higher z for
  // setting the internal energy correctly However, with Gadget2 this just sets
  // the entropy to the internal energy, which needs to be converted somehow
  if (cosmo->z >= 1)
    hydro_set_init_internal_energy(
        p, (u * scale_factor * scale_factor) / cooling->internal_energy_scale);
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
  char *parametersFileName = "./coolingBox.yml";

  float nh;              // hydrogen number density
  double u;              // internal energy

  // Read options
  int param;
  float redshift = -1.0,
        log_10_nh =
            100;  // unreasonable values will be changed if not set in options
  while ((param = getopt(argc, argv, "z:d:m:t")) != -1) switch (param) {
      case 'z':
        // read redshift
        redshift = atof(optarg);
        break;
      case 'd':
        // read log10 of hydrogen number density
        log_10_nh = atof(optarg);
        break;
      case 'm':
        // read which yml file we need to use
        parametersFileName = optarg;
        break;
      case '?':
        if (optopt == 'z')
          printf("Option -%c requires an argument.\n", optopt);
        else
          printf("Unknown option character `\\x%x'.\n", optopt);
        error("invalid option(s) to testCooling");
    }

  // Read the parameter file
  if (params == NULL) error("Error allocating memory for the parameter file.");
  message("Reading runtime parameters from file '%s'", parametersFileName);
  parser_read_file(parametersFileName, params);

  // Init units
  units_init_from_params(&us, params, "InternalUnitSystem");
  phys_const_init(&us, params, &phys_const);

  // Init chemistry
  chemistry_init(params, &us, &phys_const, &chem_data);
  chemistry_first_init_part(&phys_const, &us, &cosmo, &chem_data, &p, &xp);
  chemistry_print(&chem_data);

  // Init cosmology
  cosmology_init(params, &us, &phys_const, &cosmo);
  cosmology_print(&cosmo);
  if (redshift == -1.0) {
    cosmo.z = 7.0;
  } else {
    cosmo.z = redshift;
  }
  message("redshift %.5e", cosmo.z);

  // Init cooling
  cooling_init(params, &us, &phys_const, &cooling);
  cooling_print(&cooling);
  cooling_update(&cosmo, &cooling, 0);

  // Calculate abundance ratios
  float *abundance_ratio;
  abundance_ratio = malloc((chemistry_element_count + 2) * sizeof(float));
  abundance_ratio_to_solar(&p, &cooling, abundance_ratio);

  // extract mass fractions, calculate table indices and offsets
  float XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac =
      p.chemistry_data.metal_mass_fraction[chemistry_element_He] /
      (XH + p.chemistry_data.metal_mass_fraction[chemistry_element_He]);
  int He_i;//, n_h_i;
  float d_He;//, d_n_h;
  get_index_1d(cooling.HeFrac, cooling.N_He, HeFrac, &He_i, &d_He);

  int n_z = 100;
  int n_nh = 100;
  int n_u = 100; 
  
  struct hydro_props hydro_properties;
  hydro_properties.minimal_internal_energy = exp(M_LN10 * cooling.Therm[0]) / cooling.internal_energy_scale;
  float dt_cool = 0.1;
  float dt_therm = 0.1;
  
  float delta_z = (cooling.Redshifts[cooling.N_Redshifts - 1] - cooling.Redshifts[0])/n_z;
  float delta_nh = (cooling.nH[cooling.N_nH - 1] - cooling.nH[0])/n_nh;
  float delta_u = (cooling.Therm[cooling.N_Temp - 1] - cooling.Therm[0])/n_u;
  for(int i = 0; i < n_z; i++) {
    for(int j = 0; j < n_nh; j++) {
      for (int k = 0; k < n_u; k++) {
        float z = cooling.Redshifts[0] + delta_z*i;
	nh = exp(M_LN10 * cooling.nH[0] + delta_nh*j);
	u = exp(M_LN10 * cooling.Therm[0] + delta_u*k);
        set_quantities(&p,&xp,&us,&cooling,&cosmo,&phys_const, nh, u, z);

	cooling_cool_part(&phys_const,&us,&cosmo,&hydro_properties,&cooling,&p,&xp,dt_cool,dt_therm);
      }
    }
  }
  
  message("done test");

  free(params);
  return 0;
}

