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
  struct xpart xp;
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
  parser_write_params_to_file(params, "used_parameters.yml");

  units_init(&us, params, "InternalUnitSystem");
  phys_const_init(&us, params, &internal_const);

  double hydrogen_number_density_cgs = 1e-4;
  double u, hydrogen_number_density, pressure, gamma, cooling_du_dt, temperature_cgs;
  int n_t_i = 2000;

  char *output_filename = "cooling_output.dat";
  FILE *output_file = fopen(output_filename, "w");
  if (output_file == NULL)
  {
      printf("Error opening file!\n");
      exit(1);
  }
  char *output_filename2 = "temperature_output.dat";
  FILE *output_file2 = fopen(output_filename2, "w");
  if (output_file2 == NULL)
  {
      printf("Error opening file!\n");
      exit(1);
  }

  gamma = 5.0/3.0;

  chemistry_init(params, &us, &internal_const, &chemistry_data);
  chemistry_first_init_part(&p,&xp,&chemistry_data);
  chemistry_print(&chemistry_data);
    
  u = 1.0*pow(10.0,11)/(units_cgs_conversion_factor(&us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(&us,UNIT_CONV_MASS));
  pressure = u*p.rho*(gamma -1.0);
  hydrogen_number_density = hydrogen_number_density_cgs*pow(units_cgs_conversion_factor(&us,UNIT_CONV_LENGTH),3);
  p.rho = hydrogen_number_density*internal_const.const_proton_mass*(1.0+p.chemistry_data.metal_mass_fraction[EAGLE_Helium]/p.chemistry_data.metal_mass_fraction[EAGLE_Hydrogen]);
  p.entropy = pressure/(pow(p.rho,gamma));

  cosmology_init(params, &us, &internal_const, &cosmo);
  cosmology_print(&cosmo);

  cooling_init(params, &us, &internal_const, &cooling);
  cooling_print(&cooling);


  for(int t_i = 0; t_i < n_t_i; t_i++){
    
    u = 1.0*pow(10.0,11 + t_i*6.0/n_t_i)/(units_cgs_conversion_factor(&us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(&us,UNIT_CONV_MASS));
    pressure = u*p.rho*(gamma -1.0);
    hydrogen_number_density = hydrogen_number_density_cgs*pow(units_cgs_conversion_factor(&us,UNIT_CONV_LENGTH),3);
    p.rho = hydrogen_number_density*internal_const.const_proton_mass*(1.0+p.chemistry_data.metal_mass_fraction[EAGLE_Helium]/p.chemistry_data.metal_mass_fraction[EAGLE_Hydrogen]);
    p.entropy = pressure/(pow(p.rho,gamma));

    cooling_du_dt = eagle_print_metal_cooling_rate(&p,&cooling,&cosmo,&internal_const);
    temperature_cgs = eagle_convert_u_to_temp(&p,&cooling,&cosmo,&internal_const);
    fprintf(output_file,"%.5e %.5e\n",temperature_cgs,cooling_du_dt);
    fprintf(output_file2,"%.5e %.5e\n",u*(units_cgs_conversion_factor(&us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(&us,UNIT_CONV_MASS)), temperature_cgs);
  }
  fclose(output_file);
  fclose(output_file2);

  free(params);

  return 0;
}
