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

  double hydrogen_number_density_cgs = 1e-1;
  double u,u_cgs, hydrogen_number_density, pressure, gamma, cooling_du_dt, temperature_cgs, newton_func, u_ini_cgs, ratefact, dt_cgs;
  u_cgs = 0; cooling_du_dt = 0; temperature_cgs = 0; newton_func = 0;
  //int n_t_i = 2000;
  dt_cgs = 4.0e-8*units_cgs_conversion_factor(&us,UNIT_CONV_TIME);

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
  char *output_filename3 = "newton_output.dat";
  FILE *output_file3 = fopen(output_filename3, "w");
  if (output_file3 == NULL)
  {
      printf("Error opening file!\n");
      exit(1);
  }

  gamma = 5.0/3.0;

  chemistry_init(params, &us, &internal_const, &chemistry_data);
  chemistry_first_init_part(&p,&xp,&chemistry_data);
  chemistry_print(&chemistry_data);
  
  cosmology_init(params, &us, &internal_const, &cosmo);
  cosmology_print(&cosmo);

    
  u = 1.0*pow(10.0,11)/(units_cgs_conversion_factor(&us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(&us,UNIT_CONV_MASS));
  pressure = u*p.rho*(gamma -1.0);
  hydrogen_number_density = hydrogen_number_density_cgs*pow(units_cgs_conversion_factor(&us,UNIT_CONV_LENGTH),3);
  p.rho = hydrogen_number_density*internal_const.const_proton_mass*(1.0+p.chemistry_data.metal_mass_fraction[EAGLE_Helium]/p.chemistry_data.metal_mass_fraction[EAGLE_Hydrogen]);
  p.entropy = pressure/(pow(p.rho,gamma));


  cooling_init(params, &us, &internal_const, &cooling);
  cooling_print(&cooling);

  // construct 1d table of cooling rates wrt temperature
  float H_plus_He_heat_table[176];                      // WARNING sort out how it is declared/allocated
  float H_plus_He_electron_abundance_table[176];        // WARNING sort out how it is declared/allocated
  float temp_table[176];        			// WARNING sort out how it is declared/allocated
  float element_cooling_table[9*176];                   // WARNING sort out how it is declared/allocated
  float element_electron_abundance_table[176];          // WARNING sort out how it is declared/allocated
  construct_1d_table_from_4d(&p,&cooling,&cosmo,&internal_const,cooling.table.element_cooling.H_plus_He_heating,H_plus_He_heat_table);
  construct_1d_table_from_4d(&p,&cooling,&cosmo,&internal_const,cooling.table.element_cooling.H_plus_He_electron_abundance,H_plus_He_electron_abundance_table);
  construct_1d_table_from_4d(&p,&cooling,&cosmo,&internal_const,cooling.table.element_cooling.temperature,temp_table);
  construct_1d_table_from_4d_elements(&p,&cooling,&cosmo,&internal_const,cooling.table.element_cooling.metal_heating,element_cooling_table);
  construct_1d_table_from_3d(&p,&cooling,&cosmo,&internal_const,cooling.table.element_cooling.electron_abundance,element_electron_abundance_table);

  float XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];
  float inn_h = chemistry_get_number_density(&p,&cosmo,chemistry_element_H,&internal_const)*cooling.number_density_scale;
  ratefact = inn_h * (XH / eagle_proton_mass_cgs);
  u_ini_cgs = pow(10.0,14);

  // Compute contributions to cooling rate from different metals
  //for(int t_i = 0; t_i < n_t_i; t_i++){
  //  
  //  u_cgs = pow(10.0,11 + t_i*6.0/n_t_i);
  //  u = u_cgs/(units_cgs_conversion_factor(&us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(&us,UNIT_CONV_MASS));
  //  pressure = u*p.rho*(gamma -1.0);
  //  hydrogen_number_density = hydrogen_number_density_cgs*pow(units_cgs_conversion_factor(&us,UNIT_CONV_LENGTH),3);
  //  p.rho = hydrogen_number_density*internal_const.const_proton_mass*(1.0+p.chemistry_data.metal_mass_fraction[EAGLE_Helium]/p.chemistry_data.metal_mass_fraction[EAGLE_Hydrogen]);
  //  p.entropy = pressure/(pow(p.rho,gamma));

  //  //cooling_du_dt = eagle_print_metal_cooling_rate(&p,&cooling,&cosmo,&internal_const);
  //  cooling_du_dt = eagle_print_metal_cooling_rate_1d_table(H_plus_He_heat_table,H_plus_He_electron_abundance_table,element_cooling_table,element_electron_abundance_table,temp_table,&p,&cooling,&cosmo,&internal_const);
  //  temperature_cgs = eagle_convert_u_to_temp(u_cgs,&p,&cooling,&cosmo,&internal_const);
  //  fprintf(output_file,"%.5e %.5e\n",temperature_cgs,cooling_du_dt);
  //  fprintf(output_file2,"%.5e %.5e\n",u*(units_cgs_conversion_factor(&us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(&us,UNIT_CONV_MASS)), temperature_cgs);

  //  newton_func = u_cgs - u_ini_cgs - cooling_du_dt*ratefact*dt_cgs;
  //  fprintf(output_file3,"%.5e %.5e\n",u_cgs,newton_func);
  //  //printf("u,u_ini,cooling_du_dt,ratefact,dt %.5e %.5e %.5e %.5e %.5e \n",u_cgs,u_ini_cgs,cooling_du_dt,ratefact,dt_cgs);
  //}
  //fclose(output_file);
  //fclose(output_file2);
  //fclose(output_file3);

  double dLambdaNet_du, LambdaNet, LambdaNext;
  for(int j = 0; j < 5; j++){
    float u_ini = eagle_convert_temp_to_u_1d_table(pow(10.0,j+4),temp_table,&p,&cooling,&cosmo,&internal_const),x,du;
    float dt = 2.0e-2*units_cgs_conversion_factor(&us,UNIT_CONV_TIME);
    LambdaNet = eagle_cooling_rate_1d_table(u_ini, &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, &p, &cooling, &cosmo, &internal_const);
    float u_temp = u_ini + LambdaNet*ratefact*dt;
    if (u_temp > 0) LambdaNext = eagle_cooling_rate_1d_table(u_temp, &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table, element_cooling_table, element_electron_abundance_table, temp_table, &p, &cooling, &cosmo, &internal_const);
    if (fabs(LambdaNet - LambdaNext)/LambdaNet < 0.5) {
      u_temp = u_ini;
    } else {
      u_temp = eagle_convert_temp_to_u_1d_table(1.0e4,temp_table,&p,&cooling,&cosmo,&internal_const);
    }
    float x_init = log(u_temp);

    x = newton_iter(x_init,u_ini,H_plus_He_heat_table,H_plus_He_electron_abundance_table,element_cooling_table,element_electron_abundance_table,temp_table,&p,&cosmo,&cooling,&internal_const,dt);
    printf("testing newton integration, u_ini, u %.5e %.5e, temperature initial, final %.5e %.5e\n", u_ini, exp(x), eagle_convert_u_to_temp_1d_table(u_ini,&du,temp_table,&p,&cooling,&cosmo,&internal_const), eagle_convert_u_to_temp_1d_table(exp(x),&du,temp_table,&p,&cooling,&cosmo,&internal_const));
  }

  free(params);

  return 0;
}
