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

void set_quantities(struct part *restrict p, 
                    const struct unit_system *restrict us,
		    const struct cooling_function_data *restrict cooling,
		    const struct cosmology *restrict cosmo,
		    const struct phys_const *restrict internal_const,
		    float nh,
		    double u){
     
  const float gamma = 5.0/3.0;
  double hydrogen_number_density = nh*pow(units_cgs_conversion_factor(us,UNIT_CONV_LENGTH),3)*(1.0/(pow(1.0+cosmo->z,3)));
  p->rho = hydrogen_number_density*internal_const->const_proton_mass*
          (1.0+p->chemistry_data.metal_mass_fraction[EAGLE_Helium]/p->chemistry_data.metal_mass_fraction[EAGLE_Hydrogen]);

  float scale_factor = 1.0/(1.0+cosmo->z);
  float pressure = (u*pow(scale_factor,3))/cooling->internal_energy_scale*p->rho*(gamma -1.0);
  p->entropy = pressure/(pow(p->rho,gamma));
}

void construct_1d_tables_test(const struct part *restrict p,
			 const struct cooling_function_data *restrict cooling,
			 const struct cosmology *restrict cosmo,
			 const struct phys_const *restrict internal_const,
			 double *temp_table,
			 double *H_plus_He_heat_table,
			 double *H_plus_He_electron_abundance_table,
			 double *element_electron_abundance_table,
			 double *element_print_cooling_table,
			 double *element_cooling_table,
			 float *abundance_ratio,
			 float *ub, float *lb){

  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] /
                 (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  float inn_h = chemistry_get_number_density(p, cosmo, chemistry_element_H,
                                             internal_const) *
                cooling->number_density_scale;

  int z_index,He_i,n_h_i;
  float dz,d_He,d_n_h;
  
  get_redshift_index(cosmo->z, &z_index, &dz, cooling);
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);
  printf("testCooling 1d z_index dz nh_i d_nh He_i d_He %d %.5e %d %.5e %d %.5e\n", z_index,dz, n_h_i, d_n_h, He_i, d_He);

  construct_1d_table_from_4d(
      p, cooling, cosmo, internal_const,
      cooling->table.element_cooling.temperature,
      z_index, dz, cooling->N_Redshifts, n_h_i, d_n_h,
      cooling->N_nH, He_i, d_He, cooling->N_He, cooling->N_Temp, temp_table, ub, lb);
  construct_1d_table_from_4d(
      p, cooling, cosmo, internal_const,
      cooling->table.element_cooling.H_plus_He_heating, z_index, dz, cooling->N_Redshifts, n_h_i, d_n_h,
      cooling->N_nH, He_i, d_He, cooling->N_He, cooling->N_Temp, H_plus_He_heat_table, ub, lb);
  construct_1d_table_from_4d(
      p, cooling, cosmo, internal_const,
      cooling->table.element_cooling.H_plus_He_electron_abundance, z_index,
      dz, cooling->N_Redshifts, n_h_i, d_n_h, cooling->N_nH, He_i, d_He,
      cooling->N_He, cooling->N_Temp, H_plus_He_electron_abundance_table, ub, lb);
  construct_1d_table_from_4d_elements(
      p, cooling, cosmo, internal_const,
      cooling->table.element_cooling.metal_heating, z_index, dz, cooling->N_Redshifts,
      n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_cooling_table,abundance_ratio, ub, lb);
  construct_1d_print_table_from_4d_elements(
      p, cooling, cosmo, internal_const,
      cooling->table.element_cooling.metal_heating, z_index, dz, cooling->N_Redshifts,
      n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_print_cooling_table,abundance_ratio, ub, lb);
  construct_1d_table_from_3d(
      p, cooling, cosmo, internal_const,
      cooling->table.element_cooling.electron_abundance, z_index, dz, cooling->N_Redshifts,
      n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_electron_abundance_table, ub, lb);
    
}

void compare_dlambda_du(
  const struct unit_system *restrict us,
  struct part *restrict p,
  const struct xpart *restrict xp,
  const struct phys_const *restrict internal_const,
  const struct cooling_function_data *restrict cooling,
  const struct cosmology *restrict cosmo){

  double H_plus_He_heat_table[176];              
  double H_plus_He_electron_abundance_table[176];
  double temp_table[176];  
  double element_cooling_table[176];         
  double element_print_cooling_table[9 * 176];  
  double element_electron_abundance_table[176];
  double rate_element_table[11];
  float *abundance_ratio, cooling_du_dt1, cooling_du_dt2;
  double dlambda_du1, dlambda_du2;

  abundance_ratio = malloc((chemistry_element_count + 2)*sizeof(float));
  abundance_ratio_to_solar(p, cooling, abundance_ratio);
  
  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] /
                 (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  
  int z_index,He_i,n_h_i;
  float dz,d_He,d_n_h;
  get_redshift_index(cosmo->z, &z_index, &dz, cooling);
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);

  float nh = 1.0e-1, u = 1.0e9;
  set_quantities(p, us, cooling, cosmo, internal_const, nh, u);

  const float log_10_e = 0.43429448190325182765;
  float upper_bound = cooling->Temp[cooling->N_Temp-1]/log_10_e;
  float lower_bound = cooling->Temp[0]/log_10_e;
  construct_1d_tables_test(p, cooling, cosmo, internal_const, temp_table, 
                           H_plus_He_heat_table, H_plus_He_electron_abundance_table, 
			   element_electron_abundance_table, element_print_cooling_table, 
			   element_cooling_table, abundance_ratio, &upper_bound, &lower_bound);
  
  float inn_h = chemistry_get_number_density(p, cosmo, chemistry_element_H,
                                           internal_const) *
              cooling->number_density_scale;
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);
  printf("testCooling.c 4d z_i dz nh_i d_nh He_i d_He %d %.5e %d %.5e %d %.5e \n", z_index, dz, n_h_i, d_n_h, He_i, d_He);
  int nt = 10;
  for(int i = 0; i < nt; i++){
    u = pow(10.0,9 + i);
    set_quantities(p, us, cooling, cosmo, internal_const, nh, u);
    cooling_du_dt1 = eagle_metal_cooling_rate_1d_table(log10(u),&dlambda_du1,H_plus_He_heat_table,
                      H_plus_He_electron_abundance_table,
       	              element_print_cooling_table,
       	              element_electron_abundance_table,
       	              temp_table,
       	              z_index, dz, n_h_i, d_n_h, He_i, d_He,
       	              p,cooling,cosmo,internal_const,rate_element_table,
       	              abundance_ratio);
    cooling_du_dt2 = eagle_metal_cooling_rate(log10(u),
       	              &dlambda_du2,z_index, dz, n_h_i, d_n_h, He_i, d_He,
       	              p,cooling,cosmo,internal_const,NULL,
       	              abundance_ratio);
    printf("u du_dt_1d du_dt_4d dlambda_du_1d dlambda_du_4d %.5e %.5e %.5e %.5e %.5e\n",u,cooling_du_dt1,cooling_du_dt2,dlambda_du1,dlambda_du2);
  }
}

void compare_temp(
  const struct unit_system *restrict us,
  struct part *restrict p,
  const struct xpart *restrict xp,
  const struct phys_const *restrict internal_const,
  const struct cooling_function_data *restrict cooling,
  const struct cosmology *restrict cosmo){

  double H_plus_He_heat_table[176];              
  double H_plus_He_electron_abundance_table[176];
  double temp_table[176];  
  double element_cooling_table[176];         
  double element_print_cooling_table[9 * 176];  
  double element_electron_abundance_table[176];
  float *abundance_ratio;

  abundance_ratio = malloc((chemistry_element_count + 2)*sizeof(float));
  abundance_ratio_to_solar(p, cooling, abundance_ratio);
  
  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] /
                 (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  
  int z_index,He_i,n_h_i;
  float dz,d_He,d_n_h;
  get_redshift_index(cosmo->z, &z_index, &dz, cooling);
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);

  float nh = 1.0e-1, u = 1.0e9;
  set_quantities(p, us, cooling, cosmo, internal_const, nh, u);
  const float log_10_e = 0.43429448190325182765;
  float upper_bound = cooling->Temp[cooling->N_Temp-1]/log_10_e;
  float lower_bound = cooling->Temp[0]/log_10_e;

  construct_1d_tables_test(p, cooling, cosmo, internal_const, temp_table, 
                           H_plus_He_heat_table, H_plus_He_electron_abundance_table, 
			   element_electron_abundance_table, element_print_cooling_table, 
			   element_cooling_table, abundance_ratio, &upper_bound, &lower_bound);
  float T1d, T4d, delta_u;
  float inn_h = chemistry_get_number_density(p, cosmo, chemistry_element_H,
                                           internal_const) *
              cooling->number_density_scale;
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);
  printf("testCooling.c z_i dz nh_i d_nh He_i d_He %d %.5e %d %.5e %d %.5e \n", z_index, dz, n_h_i, d_n_h, He_i, d_He);
  int nt = 10;
  for(int i = 0; i < nt; i++){
    u = pow(10.0,9 + i);
    set_quantities(p, us, cooling, cosmo, internal_const, nh, u);
    T1d = eagle_convert_u_to_temp_1d_table(log10(u), &delta_u, temp_table,
                   p,cooling, cosmo, internal_const);
    T4d = eagle_convert_u_to_temp(log10(u), &delta_u, z_index, n_h_i, He_i,
        		 dz, d_n_h, d_He, p,cooling, cosmo, internal_const);
    printf("u T1d T4d %.5e %.5e %.5e\n",u,pow(10.0,T1d),pow(10.0,T4d));
  }
}

int main(int argc, char **argv) {
  /* Declare relevant structs */
  struct swift_params *params = malloc(sizeof(struct swift_params));
  struct unit_system us;
  struct chemistry_global_data chem_data;
  struct part p;
  struct xpart xp;
  struct phys_const internal_const;
  struct cooling_function_data cooling;
  struct cosmology cosmo;
  char *parametersFileName = "./testCooling.yml";

  int tables = 0;

  /* Read the parameter file */
  if (params == NULL) error("Error allocating memory for the parameter file.");
  message("Reading runtime parameters from file '%s'", parametersFileName);
  parser_read_file(parametersFileName, params);

  /* And dump the parameters as used. */
  parser_write_params_to_file(params, "used_parameters.yml");

  /* Init units */
  units_init_from_params(&us, params, "InternalUnitSystem");
  phys_const_init(&us, params, &internal_const);

  /* Init chemistry */
  chemistry_init(params, &us, &internal_const, &chem_data);
  chemistry_first_init_part(&internal_const, &us, &cosmo, &chem_data, &p, &xp);
  chemistry_print(&chem_data);

  /* Init cosmology */
  cosmology_init(params, &us, &internal_const, &cosmo);
  cosmology_print(&cosmo);

  /* Init cooling */
  cooling_init(params, &us, &internal_const, &cooling);
  cooling_print(&cooling);

  /* Calculate abundance ratios */
  float *abundance_ratio;
  abundance_ratio = malloc((chemistry_element_count + 2)*sizeof(float));
  abundance_ratio_to_solar(&p, &cooling, abundance_ratio);
  
  /* Declare 1D tables */
  double H_plus_He_heat_table[176];              
  double H_plus_He_electron_abundance_table[176];
  double temp_table[176];  
  double element_cooling_table[176];         
  double element_print_cooling_table[9 * 176];  
  double element_electron_abundance_table[176];
  
  float nh;

  float XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac = p.chemistry_data.metal_mass_fraction[chemistry_element_He] /
                 (XH + p.chemistry_data.metal_mass_fraction[chemistry_element_He]);
  int z_index,He_i,n_h_i;
  float dz,d_He,d_n_h;
  get_redshift_index(cosmo.z, &z_index, &dz, &cooling);
  get_index_1d(cooling.HeFrac, cooling.N_He, HeFrac, &He_i, &d_He);

  const float log_10_e = 0.43429448190325182765;
  float upper_bound = cooling.Temp[cooling.N_Temp-1]/log_10_e;
  float lower_bound = cooling.Temp[0]/log_10_e;

  /* Loop over densities */
  int nt = 250, n_nh = 6;
  double u = pow(10.0,10);
  if (argc == 1 || strcmp(argv[1], "nh") == 0){
    for(int i = 0; i < n_nh; i++){
      /* Open files */
      char output_filename[21];
      sprintf(output_filename, "%s%d%s", "cooling_output_", i, ".dat");
      FILE *output_file = fopen(output_filename, "w");
      if (output_file == NULL) {
        printf("Error opening file!\n");
        exit(1);
      }

      nh = pow(10.0,-i);
      set_quantities(&p, &us, &cooling, &cosmo, &internal_const, nh, u);
      construct_1d_tables_test(&p, &cooling, &cosmo, &internal_const, 
        			temp_table, H_plus_He_heat_table, 
          		H_plus_He_electron_abundance_table, 
          		element_electron_abundance_table, 
          		element_print_cooling_table, 
          		element_cooling_table, 
          		abundance_ratio, &upper_bound, &lower_bound);
      get_index_1d(cooling.nH, cooling.N_nH, log10(nh), &n_h_i, &d_n_h);

      for(int j = 0; j < nt; j++){
        set_quantities(&p, &us, &cooling, &cosmo, &internal_const, nh, pow(10.0,11.0 + j*8.0/nt));
	u = hydro_get_physical_internal_energy(&p,&cosmo)*cooling.internal_energy_scale;
        double dlambda_du;
        float delta_u, cooling_du_dt, logT;
        if (tables == 1) {
	  cooling_du_dt = eagle_metal_cooling_rate_1d_table(log10(u),&dlambda_du,H_plus_He_heat_table,
                                  H_plus_He_electron_abundance_table,
             	              element_print_cooling_table,
             	              element_electron_abundance_table,
             	              temp_table,
             	              z_index, dz, n_h_i, d_n_h, He_i, d_He,
             	              &p,&cooling,&cosmo,&internal_const,NULL,
             	              abundance_ratio);
          logT = eagle_convert_u_to_temp_1d_table(log10(u), &delta_u, temp_table,
                         &p,&cooling, &cosmo, &internal_const);
	} else {
          cooling_du_dt = eagle_metal_cooling_rate(log10(u),
             	              &dlambda_du,z_index, dz, n_h_i, d_n_h, He_i, d_He,
             	              &p,&cooling,&cosmo,&internal_const,NULL,
             	              abundance_ratio);
          logT = eagle_convert_u_to_temp(log10(u), &delta_u, z_index, n_h_i, He_i,
	      		 dz, d_n_h, d_He, &p,&cooling, &cosmo, &internal_const);
	}
        float temperature_swift = pow(10.0,logT);

        fprintf(output_file,"%.5e %.5e\n", temperature_swift,cooling_du_dt);
      }
      fclose(output_file);
    }
  }
  if (argc == 1 || strcmp(argv[1],"metals") == 0) {
    printf("here \n");
    char output_filename[21];
    sprintf(output_filename, "%s", "cooling_output.dat");
    FILE *output_file = fopen(output_filename, "w");
    if (output_file == NULL) {
      printf("Error opening file!\n");
      exit(1);
    }
    nh = pow(10.0,0);
    u = pow(10.0,14.0);
    set_quantities(&p, &us, &cooling, &cosmo, &internal_const, nh, u);
    float nh_swift = chemistry_get_number_density(&p, &cosmo, chemistry_element_H, &internal_const)*cooling.number_density_scale;
    printf("hydrogen number density %.5e\n", nh_swift);
    construct_1d_tables_test(&p, &cooling, &cosmo, &internal_const, 
      			temp_table, H_plus_He_heat_table, 
        		H_plus_He_electron_abundance_table, 
        		element_electron_abundance_table, 
        		element_print_cooling_table, 
        		element_cooling_table, 
        		abundance_ratio, &upper_bound, &lower_bound);
    float inn_h = chemistry_get_number_density(&p, &cosmo, chemistry_element_H,
                                             &internal_const) *
                cooling.number_density_scale;
    get_index_1d(cooling.nH, cooling.N_nH, log10(inn_h), &n_h_i, &d_n_h);
    printf("testcooling.c z_i dz nh_i d_nh He_i d_He %d %.5e %d %.5e %d %.5e\n",z_index, dz, n_h_i, d_n_h, He_i, d_He);
    for(int j = 0; j < nt; j++){
      set_quantities(&p, &us, &cooling, &cosmo, &internal_const, nh, pow(10.0,10.0 + j*8.0/nt));
      u = hydro_get_physical_internal_energy(&p,&cosmo)*cooling.internal_energy_scale;
      float delta_u, cooling_du_dt, logT;
      if (tables == 1) {
        cooling_du_dt = eagle_print_metal_cooling_rate_1d_table(H_plus_He_heat_table,
                              H_plus_He_electron_abundance_table,
           	              element_print_cooling_table,
           	              element_electron_abundance_table,
           	              temp_table,
           	              z_index, dz, n_h_i, d_n_h, He_i, d_He,
           	              &p,&cooling,&cosmo,&internal_const,
           	              abundance_ratio);
        logT = eagle_convert_u_to_temp_1d_table(log10(u), &delta_u, temp_table,
                       &p,&cooling, &cosmo, &internal_const);
      } else {
        cooling_du_dt = eagle_print_metal_cooling_rate(
			      z_index, dz, n_h_i, d_n_h, He_i, d_He,
           	              &p,&cooling,&cosmo,&internal_const,
           	              abundance_ratio);
      }
      fprintf(output_file,"%.5e %.5e\n", u,cooling_du_dt);
    }
    fclose(output_file);
  }

  if (argc == 1 || strcmp(argv[1],"compare_temp") == 0) compare_temp(&us, &p, &xp, &internal_const, &cooling, &cosmo);
  if (argc == 1 || strcmp(argv[1],"compare_dlambda_du") == 0) compare_dlambda_du(&us, &p, &xp, &internal_const, &cooling, &cosmo);

  free(params);
  return 0;
}


