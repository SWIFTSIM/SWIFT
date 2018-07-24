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

#include "cooling.h"
#include "hydro.h"
#include "physical_constants.h"
#include "swift.h"
#include "units.h"

int main(int argc, char **argv) {

  /* Read the parameter file */
  struct swift_params *params = malloc(sizeof(struct swift_params));
  struct unit_system us;
  struct chemistry_global_data chem_data;
  struct part p;
  struct xpart xp;
  struct phys_const internal_const;
  struct cooling_function_data cooling;
  struct cosmology cosmo;
  char *parametersFileName = "./testCooling.yml";
  enum table_index {
    EAGLE_Hydrogen = 0,
    EAGLE_Helium,
    EAGLE_Carbon,
    EAGLE_Nitrogen,
    EAGLE_Oxygen,
    EAGLE_Neon,
    EAGLE_Magnesium,
    EAGLE_Silicon,
    EAGLE_Iron
  };

  if (params == NULL) error("Error allocating memory for the parameter file.");
  message("Reading runtime parameters from file '%s'", parametersFileName);
  parser_read_file(parametersFileName, params);

  /* And dump the parameters as used. */
  parser_write_params_to_file(params, "used_parameters.yml");

  double UnitMass_in_cgs = 1.989e43;
  double UnitLength_in_cgs = 3.085678e24;
  double UnitVelocity_in_cgs = 1e5;
  double UnitCurrent_in_cgs = 1;
  double UnitTemp_in_cgs = 1;
  units_init(&us, UnitMass_in_cgs, UnitLength_in_cgs, UnitVelocity_in_cgs, UnitCurrent_in_cgs, UnitTemp_in_cgs);
  phys_const_init(&us, params, &internal_const);

  double hydrogen_number_density_cgs = 1e-4;
  double u, u_cgs, hydrogen_number_density, pressure, gamma, cooling_du_dt,
      temperature_cgs, newton_func, u_ini_cgs, ratefact, dt_cgs;
  u_cgs = 0;
  cooling_du_dt = 0;
  temperature_cgs = 0;
  newton_func = 0;
  // int n_t_i = 2000;
  dt_cgs = 4.0e-8 * units_cgs_conversion_factor(&us, UNIT_CONV_TIME);

  char *output_filename = "cooling_output.dat";
  FILE *output_file = fopen(output_filename, "w");
  if (output_file == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }
  char *output_filename2 = "temperature_output.dat";
  FILE *output_file2 = fopen(output_filename2, "w");
  if (output_file2 == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }
  char *output_filename3 = "newton_output.dat";
  FILE *output_file3 = fopen(output_filename3, "w");
  if (output_file3 == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }

  gamma = 5.0 / 3.0;

  chemistry_init(params, &us, &internal_const, &chem_data);
  chemistry_first_init_part(&internal_const, &us, &cosmo, &chem_data, &p, &xp);
  chemistry_print(&chem_data);

  cosmology_init(params, &us, &internal_const, &cosmo);
  cosmology_print(&cosmo);

  u = 1.0 * pow(10.0, 11) /
      (units_cgs_conversion_factor(&us, UNIT_CONV_ENERGY) /
       units_cgs_conversion_factor(&us, UNIT_CONV_MASS));
  pressure = u * p.rho * (gamma - 1.0);
  hydrogen_number_density =
      hydrogen_number_density_cgs *
      pow(units_cgs_conversion_factor(&us, UNIT_CONV_LENGTH), 3);
  p.rho = hydrogen_number_density * internal_const.const_proton_mass *
          (1.0 +
           p.chemistry_data.metal_mass_fraction[EAGLE_Helium] /
               p.chemistry_data.metal_mass_fraction[EAGLE_Hydrogen]);
  p.entropy = pressure / (pow(p.rho, gamma));

  cooling_init(params, &us, &internal_const, &cooling);
  cooling_print(&cooling);

  float XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac = p.chemistry_data.metal_mass_fraction[chemistry_element_He] /
         (XH + p.chemistry_data.metal_mass_fraction[chemistry_element_He]);
  float inn_h = chemistry_get_number_density(&p, &cosmo, chemistry_element_H,
                                             &internal_const) *
                cooling.number_density_scale;
  ratefact = inn_h * (XH / eagle_proton_mass_cgs);
  
  float *abundance_ratio;
  abundance_ratio = malloc((chemistry_element_count + 2)*sizeof(float));
  abundance_ratio_to_solar(&p, &cooling, abundance_ratio);

  // construct 1d table of cooling rates wrt temperature
  double H_plus_He_heat_table[176];              
  double H_plus_He_electron_abundance_table[176];
  double temp_table[176];  
  double element_cooling_table[176];         
  double element_print_cooling_table[9 * 176];  
  double element_electron_abundance_table[176];

  int z_index,He_i,n_h_i;
  float dz,d_He,d_n_h;
  get_redshift_index(cosmo.z, &z_index, &dz, &cooling);
  get_index_1d(cooling.HeFrac, cooling.N_He, HeFrac, &He_i, &d_He);
  get_index_1d(cooling.nH, cooling.N_nH, log10(inn_h), &n_h_i, &d_n_h);

  construct_1d_table_from_4d(
      &p, &cooling, &cosmo, &internal_const,
      cooling.table.element_cooling.temperature,
      z_index, dz, cooling.N_Redshifts, n_h_i, d_n_h,
      cooling.N_nH, He_i, d_He, cooling.N_He, cooling.N_Temp, temp_table);
  construct_1d_table_from_4d(
      &p, &cooling, &cosmo, &internal_const,
      cooling.table.element_cooling.H_plus_He_heating, z_index, dz, cooling.N_Redshifts, n_h_i, d_n_h,
      cooling.N_nH, He_i, d_He, cooling.N_He, cooling.N_Temp, H_plus_He_heat_table);
  construct_1d_table_from_4d(
      &p, &cooling, &cosmo, &internal_const,
      cooling.table.element_cooling.H_plus_He_electron_abundance, z_index,
      dz, cooling.N_Redshifts, n_h_i, d_n_h, cooling.N_nH, He_i, d_He,
      cooling.N_He, cooling.N_Temp, H_plus_He_electron_abundance_table);
  construct_1d_table_from_4d_elements(
      &p, &cooling, &cosmo, &internal_const,
      cooling.table.element_cooling.metal_heating, z_index, dz, cooling.N_Redshifts,
      n_h_i, d_n_h, cooling.N_nH, cooling.N_Temp, element_cooling_table,abundance_ratio);
  construct_1d_print_table_from_4d_elements(
      &p, &cooling, &cosmo, &internal_const,
      cooling.table.element_cooling.metal_heating, z_index, dz, cooling.N_Redshifts,
      n_h_i, d_n_h, cooling.N_nH, cooling.N_Temp, element_print_cooling_table,abundance_ratio);
  construct_1d_table_from_3d(
      &p, &cooling, &cosmo, &internal_const,
      cooling.table.element_cooling.electron_abundance, z_index, dz, cooling.N_Redshifts,
      n_h_i, d_n_h, cooling.N_nH, cooling.N_Temp, element_electron_abundance_table);


  //
  // printf("cooling table values \n");
  // for( int j=0; j < 176; j++) {
  //    printf("   %.5e",  H_plus_He_heat_table[j]+element_cooling_table[j]
  u_ini_cgs = pow(10.0, 14);
  float delta_u;


   //Compute contributions to cooling rate from different metals
  // int n_t_i = 100;
  // for(int t_i = 0; t_i < n_t_i; t_i++){
  //   u_cgs = pow(10.0,11 + t_i*6.0/n_t_i);
  //   //u = u_cgs/(units_cgs_conversion_factor(&us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(&us,UNIT_CONV_MASS));
  //   u = u_cgs/cooling.internal_energy_scale;
  //   hydrogen_number_density =
  //   hydrogen_number_density_cgs*pow(units_cgs_conversion_factor(&us,UNIT_CONV_LENGTH),3);
  //   p.rho =
  //   hydrogen_number_density*internal_const.const_proton_mass*(1.0+p.chemistry_data.metal_mass_fraction[EAGLE_Helium]/p.chemistry_data.metal_mass_fraction[EAGLE_Hydrogen]);
  //   pressure = u*p.rho*(gamma -1.0);
  //   p.entropy = pressure/(pow(p.rho,gamma));
  //   double u_swift = hydro_get_physical_internal_energy(&p, &cosmo) *
  //                cooling.internal_energy_scale;

  //   cooling_du_dt = eagle_print_metal_cooling_rate_1d_table(H_plus_He_heat_table,
  //                     H_plus_He_electron_abundance_table,
  //      	       element_print_cooling_table,
  //      	       element_electron_abundance_table,
  //      	       temp_table,
  //      	       z_index, dz, n_h_i, d_n_h, He_i, d_He,
  //      	       &p,&cooling,&cosmo,&internal_const,
  //      	       abundance_ratio);
  //   float logT = eagle_convert_u_to_temp_1d_table(log10(u_swift), &delta_u, temp_table,
  //                  &p,&cooling, &cosmo, &internal_const);
  //   float temperature_swift = pow(10.0,logT);

  //   fprintf(output_file,"%.5e %.5e\n", temperature_swift,cooling_du_dt);
  //   fprintf(output_file2,"%.5e %.5e\n",u_swift*(units_cgs_conversion_factor(&us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(&us,UNIT_CONV_MASS)),
  //   temperature_cgs);

  //   newton_func = u_cgs - u_ini_cgs - cooling_du_dt*ratefact*dt_cgs;
  //   fprintf(output_file3,"%.5e %.5e\n",u_cgs,newton_func);
  //   //printf("u,u_ini,cooling_du_dt,ratefact,dt %.5e %.5e %.5e %.5e %.5e
  //   // \n",u_cgs,u_ini_cgs,cooling_du_dt,ratefact,dt_cgs);
  //}
   
   int n_t_i = 100;
   for(int nh_i = 1; nh_i < 7; nh_i++){
     char output_filename4[21];
     sprintf(output_filename4, "%s%d%s", "cooling_output_", nh_i, ".dat");
     FILE *output_file4 = fopen(output_filename4, "w");
     if (output_file4 == NULL) {
       printf("Error opening file!\n");
       exit(1);
     }
     hydrogen_number_density = pow(10.0,-nh_i)*pow(units_cgs_conversion_factor(&us,UNIT_CONV_LENGTH),3)*(1.0/(pow(1.0+cosmo.z,3)));
     p.rho =
     hydrogen_number_density*internal_const.const_proton_mass*(1.0+p.chemistry_data.metal_mass_fraction[EAGLE_Helium]/p.chemistry_data.metal_mass_fraction[EAGLE_Hydrogen]);

     XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];
     HeFrac = p.chemistry_data.metal_mass_fraction[chemistry_element_He] /
      (XH + p.chemistry_data.metal_mass_fraction[chemistry_element_He]);
     inn_h = chemistry_get_number_density(&p, &cosmo, chemistry_element_H,
                                                &internal_const) *
                   cooling.number_density_scale;
     ratefact = inn_h * (XH / eagle_proton_mass_cgs);
     printf("test cooling hydrogen num dens inn_h %.5e %.5e\n", hydrogen_number_density*cooling.number_density_scale, inn_h);
     
     abundance_ratio = malloc((chemistry_element_count + 2)*sizeof(float));
     abundance_ratio_to_solar(&p, &cooling, abundance_ratio);

     get_redshift_index(cosmo.z, &z_index, &dz, &cooling);
     get_index_1d(cooling.HeFrac, cooling.N_He, HeFrac, &He_i, &d_He);
     get_index_1d(cooling.nH, cooling.N_nH, log10(inn_h), &n_h_i, &d_n_h);

     construct_1d_table_from_4d(
         &p, &cooling, &cosmo, &internal_const,
         cooling.table.element_cooling.temperature,
         z_index, dz, cooling.N_Redshifts, n_h_i, d_n_h,
         cooling.N_nH, He_i, d_He, cooling.N_He, cooling.N_Temp, temp_table);
     construct_1d_table_from_4d(
         &p, &cooling, &cosmo, &internal_const,
         cooling.table.element_cooling.H_plus_He_heating, z_index, dz, cooling.N_Redshifts, n_h_i, d_n_h,
         cooling.N_nH, He_i, d_He, cooling.N_He, cooling.N_Temp, H_plus_He_heat_table);
     construct_1d_table_from_4d(
         &p, &cooling, &cosmo, &internal_const,
         cooling.table.element_cooling.H_plus_He_electron_abundance, z_index,
         dz, cooling.N_Redshifts, n_h_i, d_n_h, cooling.N_nH, He_i, d_He,
         cooling.N_He, cooling.N_Temp, H_plus_He_electron_abundance_table);
     construct_1d_table_from_4d_elements(
         &p, &cooling, &cosmo, &internal_const,
         cooling.table.element_cooling.metal_heating, z_index, dz, cooling.N_Redshifts,
         n_h_i, d_n_h, cooling.N_nH, cooling.N_Temp, element_cooling_table,abundance_ratio);
     construct_1d_print_table_from_4d_elements(
         &p, &cooling, &cosmo, &internal_const,
         cooling.table.element_cooling.metal_heating, z_index, dz, cooling.N_Redshifts,
         n_h_i, d_n_h, cooling.N_nH, cooling.N_Temp, element_print_cooling_table,abundance_ratio);
     construct_1d_table_from_3d(
         &p, &cooling, &cosmo, &internal_const,
         cooling.table.element_cooling.electron_abundance, z_index, dz, cooling.N_Redshifts,
         n_h_i, d_n_h, cooling.N_nH, cooling.N_Temp, element_electron_abundance_table);
     for(int t_i = 0; t_i < n_t_i; t_i++){
       u_cgs = pow(10.0,11 + t_i*6.0/n_t_i);
       u = u_cgs/cooling.internal_energy_scale;
       pressure = u*p.rho*(gamma -1.0);
       p.entropy = pressure/(pow(p.rho,gamma));
       double u_swift = hydro_get_physical_internal_energy(&p, &cosmo) *
                    cooling.internal_energy_scale;

       double dlambda_du;
       cooling_du_dt = eagle_metal_cooling_rate_1d_table(log10(u_swift),&dlambda_du,H_plus_He_heat_table,
                         H_plus_He_electron_abundance_table,
          	       element_print_cooling_table,
          	       element_electron_abundance_table,
          	       temp_table,
          	       z_index, dz, n_h_i, d_n_h, He_i, d_He,
          	       &p,&cooling,&cosmo,&internal_const,NULL,
          	       abundance_ratio);
       float logT = eagle_convert_u_to_temp_1d_table(log10(u_swift), &delta_u, temp_table,
                      &p,&cooling, &cosmo, &internal_const);
       float temperature_swift = pow(10.0,logT);

       fprintf(output_file4,"%.5e %.5e\n", temperature_swift,cooling_du_dt);
       fprintf(output_file2,"%.5e %.5e\n",u_swift*(units_cgs_conversion_factor(&us,UNIT_CONV_ENERGY)/units_cgs_conversion_factor(&us,UNIT_CONV_MASS)),
       temperature_cgs);

       newton_func = u_cgs - u_ini_cgs - cooling_du_dt*ratefact*dt_cgs;
       fprintf(output_file3,"%.5e %.5e\n",u_cgs,newton_func);
     }
     fclose(output_file4);
   }
   fclose(output_file);
   fclose(output_file2);
   fclose(output_file3);

  //double dLambdaNet_du, LambdaNet;  //, LambdaNext;
  //float x_init, u_eq = 2.0e12;
  //for (int j = 5; j < 6; j++) {
  //  float u_ini = eagle_convert_temp_to_u_1d_table(pow(10.0, 0.5 * (j + 5)),
  //                                                 temp_table, &p, &cooling,
  //                                                 &cosmo, &internal_const),
  //        x, du;
  //  float dt = 2.0e-4 * units_cgs_conversion_factor(&us, UNIT_CONV_TIME);
  //  LambdaNet = eagle_cooling_rate_1d_table(
  //      u_ini, &dLambdaNet_du, H_plus_He_heat_table,
  //      H_plus_He_electron_abundance_table, element_cooling_table,
  //      element_electron_abundance_table, temp_table, &p, &cooling, &cosmo,
  //      &internal_const);
  //  float u_temp = u_ini + LambdaNet * ratefact * dt;
  //  /* RGB removed this **
  //  if (u_temp > 0) LambdaNext = eagle_cooling_rate_1d_table(u_temp,
  //  &dLambdaNet_du, H_plus_He_heat_table, H_plus_He_electron_abundance_table,
  //  element_cooling_table, element_electron_abundance_table, temp_table, &p,
  //  &cooling, &cosmo, &internal_const);
  //  if (fabs(LambdaNet - LambdaNext)/LambdaNet < 0.5) {
  //    u_temp = u_ini;
  //  } else {
  //    u_temp =
  //  eagle_convert_temp_to_u_1d_table(1.0e4,temp_table,&p,&cooling,&cosmo,&internal_const);
  //              } */
  //  if (u_temp > u_eq) {
  //    x_init = log(u_temp);
  //  } else {
  //    x_init = log(u_eq);
  //  }
  //  x = newton_iter(x_init, u_ini, H_plus_He_heat_table,
  //                  H_plus_He_electron_abundance_table, element_cooling_table,
  //                  element_electron_abundance_table, temp_table, &p, &cosmo,
  //                  &cooling, &internal_const, dt);
  //  printf(
  //      "testing newton integration, u_ini, u %.5e %.5e, temperature initial, "
  //      "final %.5e %.5e\n",
  //      u_ini, exp(x),
  //      eagle_convert_u_to_temp_1d_table(u_ini, &du, temp_table, &p, &cooling,
  //                                       &cosmo, &internal_const),
  //      eagle_convert_u_to_temp_1d_table(exp(x), &du, temp_table, &p, &cooling,
  //                                       &cosmo, &internal_const));
  //}

  free(params);

  return 0;
}
