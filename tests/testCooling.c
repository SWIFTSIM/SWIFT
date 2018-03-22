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

#include "swift.h"
#include "cooling.h"
#include "physical_constants.h"
#include "units.h"
#include "hydro.h"

int main() {

  /* Read the parameter file */
  struct swift_params *params = malloc(sizeof(struct swift_params));
  struct unit_system us;
  struct chemistry_data chemistry_data;
  struct part p;
  struct phys_const internal_const;
  struct cooling_function_data cooling;
  struct cosmology cosmo;
  char *parametersFileName = "../examples/CoolingBox/coolingBox.yml";
  enum table_index {EAGLE_Hydrogen=0,EAGLE_Helium,EAGLE_Carbon,EAGLE_Nitrogen,
			EAGLE_Oxygen,EAGLE_Neon,EAGLE_Magnesium,EAGLE_Silicon,
			EAGLE_Iron};

  if (params == NULL) error("Error allocating memory for the parameter file.");
  message("Reading runtime parameters from file '%s'", parametersFileName);
  parser_read_file(parametersFileName, params);

  /* And dump the parameters as used. */
  // parser_print_params(&params);
  parser_write_params_to_file(params, "used_parameters.yml");

  //units_init_cgs(&us);
  units_init(&us, params, "InternalUnitSystem");
  phys_const_init(&us, params, &internal_const);

  double number_density_cgs = 0.1;
  double temperature_cgs = 1.0e6;
  //double power_scale = units_cgs_conversion_factor(&us,UNIT_CONV_POWER)/units_cgs_conversion_factor(&us,UNIT_CONV_MASS);
  double power_per_num_density_factor = units_cgs_conversion_factor(&us,UNIT_CONV_POWER)*pow(units_cgs_conversion_factor(&us,UNIT_CONV_LENGTH),3)/number_density_cgs;

  double gamma = 5.0/3.0;

  double number_density = number_density_cgs*pow(units_cgs_conversion_factor(&us,UNIT_CONV_LENGTH),3);
  p.rho = number_density*(1.0/0.6*internal_const.const_proton_mass);
  
  for (int i = 0; i < 9; i++) p.chemistry_data.metal_mass_fraction[i] = 0.0;
  p.chemistry_data.metal_mass_fraction[EAGLE_Hydrogen] = 0.752;
  p.chemistry_data.metal_mass_fraction[EAGLE_Helium] = 0.248;

  double temperature = temperature_cgs/units_cgs_conversion_factor(&us,UNIT_CONV_TEMPERATURE);
  double pressure = number_density*internal_const.const_boltzmann_k*temperature;
  printf("non-dim number density, code number density, number density scale %.5e, %.5e, %.5e\n",number_density, chemistry_get_number_density(&p,chemistry_element_H,&internal_const), pow(units_cgs_conversion_factor(&us,UNIT_CONV_LENGTH),3));
  printf("number density, boltzmann constant, temperature: %.5e, %.5e, %.5e\n",number_density_cgs, internal_const.const_boltzmann_k, temperature);
  printf("proton mass, boltzmann constant, %.5e, %.5e\n", internal_const.const_proton_mass, internal_const.const_boltzmann_k);
  double internal_energy = pressure/(p.rho*(gamma - 1.0));
  //p.entropy = internal_energy/((gamma - 1.0)*pow(p.rho,gamma - 1.0));
  p.entropy = pressure/(pow(p.rho,gamma));
  printf("double check pressure, actual pressure: %.5e, %.5e\n", hydro_get_comoving_pressure(&p), pressure);
  printf("temperature, pressure, internal energy, entropy: %.5e, %.5e, %.5e, %.5e, %.5e\n", temperature,pressure,internal_energy,p.entropy,internal_const.const_boltzmann_k);

  cosmology_init(params, &us, &internal_const, &cosmo);
  cosmology_print(&cosmo);
  printf("testCooling.c redshift %.5e\n", cosmo.z);

  cooling_init(params, &us, &internal_const, &cooling);
  cooling_print(&cooling);

  chemistry_init(params, &us, &internal_const, &chemistry_data);
  chemistry_print(&chemistry_data);

  double cooling_du_dt = eagle_cooling_rate(&p,&cooling,&cosmo,&internal_const)*power_per_num_density_factor;
  printf("cooling rate: %.5e\n",cooling_du_dt);

  free(params);

  return 0;
}
