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

/*
 * @brief Assign particle density and entropy corresponding to the 
 * hydrogen number density and internal energy specified.
 *
 * @param p Particle data structure
 * @param cooling Cooling function data structure
 * @param cosmo Cosmology data structure
 * @param internal_const Physical constants data structure
 * @param nh Hydrogen number density (cgs units)
 * @param u Internal energy (cgs units)
 */
void set_quantities(struct part *restrict p, 
                    const struct unit_system *restrict us,
		    const struct cooling_function_data *restrict cooling,
		    const struct cosmology *restrict cosmo,
		    const struct phys_const *restrict internal_const,
		    float nh,
		    double u){
     
  const float gamma = 5.0/3.0;
  float scale_factor = 1.0/(1.0+cosmo->z);
  double hydrogen_number_density = nh*pow(units_cgs_conversion_factor(us,UNIT_CONV_LENGTH),3);
  p->rho = hydrogen_number_density*internal_const->const_proton_mass/p->chemistry_data.metal_mass_fraction[EAGLE_Hydrogen];

  float pressure = (u*pow(scale_factor,2))/cooling->internal_energy_scale*p->rho*(gamma -1.0); 
  p->entropy = pressure/(pow(p->rho,gamma));

  // Using hydro_set_init_internal_energy seems to work better for higher z for setting the internal energy correctly
  // However, with Gadget2 this just sets the entropy to the internal energy, which needs to be converted somehow
  if(cosmo->z >= 1) hydro_set_init_internal_energy(p,(u*pow(scale_factor,2))/cooling->internal_energy_scale);
}

/*
 * @brief Construct 1d tables from 4d EAGLE tables by 
 * interpolating over redshift, hydrogen number density
 * and helium fraction. 
 *
 * @param p Particle data structure
 * @param cooling Cooling function data structure
 * @param cosmo Cosmology data structure
 * @param internal_const Physical constants data structure
 * @param temp_table Pointer to 1d interpolated table of temperature values
 * @param H_plus_He_heat_table Pointer to 1d interpolated table of cooling rates
 * due to hydrogen and helium
 * @param H_plus_He_electron_abundance_table Pointer to 1d interpolated table
 * of electron abundances due to hydrogen and helium
 * @param element_electron_abundance_table Pointer to 1d interpolated table 
 * of electron abundances due to metals
 * @param element_print_cooling_table Pointer to 1d interpolated table of
 * cooling rates due to each of the metals
 * @param element_cooling_table Pointer to 1d interpolated table of cooling
 * rates due to the contribution of all the metals
 * @param abundance_ratio Pointer to array of ratios of metal abundances to solar
 * @param ub Upper bound in temperature on table construction 
 * @param lb Lower bound in temperature on table construction
 */
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

  // obtain mass fractions and number density for particle
  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] /
                 (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  float inn_h = chemistry_get_number_density(p, cosmo, chemistry_element_H,
                                             internal_const) *
                cooling->number_density_scale;

  // find redshift, helium fraction, hydrogen number density indices and offsets
  int z_index,He_i,n_h_i;
  float dz,d_He,d_n_h;
  get_redshift_index(cosmo->z, &z_index, &dz, cooling);
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);

  if (cosmo->z > cooling->reionisation_redshift) { 
    // Photodissociation table 
    printf("Eagle testCooling.c photodissociation table redshift reionisation redshift %.5e %.5e\n", cosmo->z, cooling->reionisation_redshift);
    construct_1d_table_from_3d(p, cooling, cosmo, internal_const, 
                cooling->table.photodissociation_cooling.temperature, 
                n_h_i, d_n_h, cooling->N_nH, He_i, d_He, cooling->N_He, cooling->N_Temp, temp_table, ub, lb); 
    construct_1d_table_from_3d( 
                p, cooling, cosmo, internal_const, 
                cooling->table.photodissociation_cooling.H_plus_He_heating, n_h_i, d_n_h, cooling->N_nH,   
                He_i, d_He, cooling->N_He, cooling->N_Temp, H_plus_He_heat_table, ub, lb); 
    construct_1d_table_from_3d( 
                p, cooling, cosmo, internal_const, 
                cooling->table.photodissociation_cooling.H_plus_He_electron_abundance, 
                n_h_i, d_n_h, cooling->N_nH, He_i, d_He, cooling->N_He, 
                cooling->N_Temp, H_plus_He_electron_abundance_table, ub, lb);
    construct_1d_table_from_3d_elements(
                p, cooling, cosmo, internal_const,
                cooling->table.photodissociation_cooling.metal_heating, 
                n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_cooling_table, abundance_ratio, ub, lb);
    construct_1d_print_table_from_3d_elements(
                p, cooling, cosmo, internal_const,
                cooling->table.photodissociation_cooling.metal_heating, 
                n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_print_cooling_table, abundance_ratio, ub, lb);
    construct_1d_table_from_2d(
                p, cooling, cosmo, internal_const,
                cooling->table.photodissociation_cooling.electron_abundance,
                n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_electron_abundance_table, ub, lb);
  } else if (cosmo->z > cooling->Redshifts[cooling->N_Redshifts - 1]) {
    printf("Eagle testCooling.c no compton table redshift max redshift %.5e %.5e\n", cosmo->z, cooling->Redshifts[cooling->N_Redshifts - 1]);
    // High redshift table
    construct_1d_table_from_3d(p, cooling, cosmo, internal_const,
                cooling->table.no_compton_cooling.temperature,
                n_h_i, d_n_h, cooling->N_nH, He_i, d_He, cooling->N_He, cooling->N_Temp, temp_table, ub, lb);
    construct_1d_table_from_3d(
                p, cooling, cosmo, internal_const,
                cooling->table.no_compton_cooling.H_plus_He_heating, n_h_i, d_n_h, cooling->N_nH,  
                He_i, d_He, cooling->N_He, cooling->N_Temp, H_plus_He_heat_table, ub, lb);
    construct_1d_table_from_3d(
                p, cooling, cosmo, internal_const,
                cooling->table.no_compton_cooling.H_plus_He_electron_abundance,
                n_h_i, d_n_h, cooling->N_nH, He_i, d_He, cooling->N_He,
                cooling->N_Temp, H_plus_He_electron_abundance_table, ub, lb);
    construct_1d_table_from_3d_elements(
                p, cooling, cosmo, internal_const,
                cooling->table.no_compton_cooling.metal_heating, 
                n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_cooling_table, abundance_ratio, ub, lb);
    construct_1d_print_table_from_3d_elements(
                p, cooling, cosmo, internal_const,
                cooling->table.no_compton_cooling.metal_heating, 
                n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_print_cooling_table, abundance_ratio, ub, lb);
    construct_1d_table_from_2d(
                p, cooling, cosmo, internal_const,
                cooling->table.no_compton_cooling.electron_abundance,
                n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_electron_abundance_table, ub, lb);
  } else {
    // Normal tables 
    printf("Eagle testCooling.c normal table redshift %.5e\n", cosmo->z);
    construct_1d_table_from_4d(p, cooling, cosmo, internal_const,
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
                n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_cooling_table, abundance_ratio, ub, lb);
    construct_1d_print_table_from_4d_elements(
                p, cooling, cosmo, internal_const,
                cooling->table.element_cooling.metal_heating, z_index, dz, cooling->N_Redshifts, 
                n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_print_cooling_table, abundance_ratio, ub, lb);
    construct_1d_table_from_3d(
                p, cooling, cosmo, internal_const,
                cooling->table.element_cooling.electron_abundance, z_index, dz, cooling->N_Redshifts,
                n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_electron_abundance_table, ub, lb);
  } 
    
}

/*
 * @brief Compare calculation of dlambda/du (gradient of cooling rate, 
 * needed for Newton's method) between interpolating 1d tables and 
 * 4d tables
 *
 * @param us Units data structure
 * @param p Particle data structure
 * @param xp Extended particle data structure
 * @param internal_const Physical constants data structure
 * @param cooling Cooling function data structure
 * @param cosmo Cosmology data structure
 */
void compare_dlambda_du(
  const struct unit_system *restrict us,
  struct part *restrict p,
  const struct xpart *restrict xp,
  const struct phys_const *restrict internal_const,
  const struct cooling_function_data *restrict cooling,
  const struct cosmology *restrict cosmo){

  // allocate tables
  double H_plus_He_heat_table[cooling->N_Temp];              
  double H_plus_He_electron_abundance_table[cooling->N_Temp];
  double temp_table[cooling->N_Temp];  
  double element_cooling_table[cooling->N_Temp];         
  double element_print_cooling_table[cooling->N_Elements * cooling->N_Temp];  
  double element_electron_abundance_table[cooling->N_Temp];
  double rate_element_table[cooling->N_Elements+2];
  float *abundance_ratio, cooling_du_dt1, cooling_du_dt2;
  double dlambda_du1, dlambda_du2;

  // calculate ratio of particle metal abundances to solar 
  abundance_ratio = malloc((chemistry_element_count + 2)*sizeof(float));
  abundance_ratio_to_solar(p, cooling, abundance_ratio);
  
  // set hydrogen number density and internal energy
  float nh = 1.0e-1, u = 1.0e9;
  set_quantities(p, us, cooling, cosmo, internal_const, nh, u);
  
  // extract hydrogen and helium mass fractions
  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] /
                 (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  float inn_h = chemistry_get_number_density(p, cosmo, chemistry_element_H,
                                           internal_const) *
              cooling->number_density_scale;
  
  // find redshift, hydrogen number density and helium fraction indices and offsets
  int z_index,He_i,n_h_i;
  float dz,d_He,d_n_h;
  get_redshift_index(cosmo->z, &z_index, &dz, cooling);
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);


  // construct tables
  float upper_bound = cooling->Temp[cooling->N_Temp-1]/eagle_log_10_e;
  float lower_bound = cooling->Temp[0]/eagle_log_10_e;
  construct_1d_tables_test(p, cooling, cosmo, internal_const, temp_table, 
                           H_plus_He_heat_table, H_plus_He_electron_abundance_table, 
			   element_electron_abundance_table, element_print_cooling_table, 
			   element_cooling_table, abundance_ratio, &upper_bound, &lower_bound);
  
  // calculate dlambda/du for different values of internal energy
  int nt = 10;
  for(int i = 0; i < nt; i++){
    u = pow(10.0,9 + i);
    set_quantities(p, us, cooling, cosmo, internal_const, nh, u);
    cooling_du_dt1 = eagle_metal_cooling_rate_1d_table(log10(u),&dlambda_du1,H_plus_He_heat_table,
                      H_plus_He_electron_abundance_table,
       	              element_print_cooling_table,
       	              element_electron_abundance_table,
       	              temp_table,
       	              p,cooling,cosmo,internal_const,rate_element_table);
    cooling_du_dt2 = eagle_metal_cooling_rate(log10(u),
       	              &dlambda_du2,z_index, dz, n_h_i, d_n_h, He_i, d_He,
       	              p,cooling,cosmo,internal_const,NULL,
       	              abundance_ratio);
    printf("u du_dt_1d du_dt_4d dlambda_du_1d dlambda_du_4d %.5e %.5e %.5e %.5e %.5e\n",u,cooling_du_dt1,cooling_du_dt2,dlambda_du1,dlambda_du2);
  }
}

/*
 * @brief Compare calculating temperature between interpolating 1d tables and 
 * 4d tables
 *
 * @param us Units data structure
 * @param p Particle data structure
 * @param xp Extended particle data structure
 * @param internal_const Physical constants data structure
 * @param cooling Cooling function data structure
 * @param cosmo Cosmology data structure
 */
void compare_temp(
  const struct unit_system *restrict us,
  struct part *restrict p,
  const struct xpart *restrict xp,
  const struct phys_const *restrict internal_const,
  const struct cooling_function_data *restrict cooling,
  const struct cosmology *restrict cosmo){

  // allocate tables
  double H_plus_He_heat_table[cooling->N_Temp];              
  double H_plus_He_electron_abundance_table[cooling->N_Temp];
  double temp_table[cooling->N_Temp];  
  double element_cooling_table[cooling->N_Temp];         
  double element_print_cooling_table[cooling->N_Elements * cooling->N_Temp];  
  double element_electron_abundance_table[cooling->N_Temp];
  float *abundance_ratio;

  // calculate ratio of particle metal abundances to solar 
  abundance_ratio = malloc((chemistry_element_count + 2)*sizeof(float));
  abundance_ratio_to_solar(p, cooling, abundance_ratio);
  
  // set hydrogen number density and internal energy
  float nh = 1.0e-1, u = 1.0e9;
  set_quantities(p, us, cooling, cosmo, internal_const, nh, u);

  // extract hydrogen and helium mass fractions
  float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac = p->chemistry_data.metal_mass_fraction[chemistry_element_He] /
                 (XH + p->chemistry_data.metal_mass_fraction[chemistry_element_He]);
  float inn_h = chemistry_get_number_density(p, cosmo, chemistry_element_H,
                                           internal_const) *
              cooling->number_density_scale;
  
  // find redshift, hydrogen number density and helium fraction indices and offsets
  int z_index,He_i,n_h_i;
  float dz,d_He,d_n_h;
  get_redshift_index(cosmo->z, &z_index, &dz, cooling);
  get_index_1d(cooling->HeFrac, cooling->N_He, HeFrac, &He_i, &d_He);
  get_index_1d(cooling->nH, cooling->N_nH, log10(inn_h), &n_h_i, &d_n_h);

  // construct tables
  float upper_bound = cooling->Temp[cooling->N_Temp-1]/eagle_log_10_e;
  float lower_bound = cooling->Temp[0]/eagle_log_10_e;
  construct_1d_tables_test(p, cooling, cosmo, internal_const, temp_table, 
                           H_plus_He_heat_table, H_plus_He_electron_abundance_table, 
			   element_electron_abundance_table, element_print_cooling_table, 
			   element_cooling_table, abundance_ratio, &upper_bound, &lower_bound);
  
  // calculate temperature for different values of internal energy
  float T1d, T4d, delta_u;
  int nt = 10;
  for(int i = 0; i < nt; i++){
    u = pow(10.0,9 + i);
    set_quantities(p, us, cooling, cosmo, internal_const, nh, u);
    T1d = eagle_convert_u_to_temp_1d_table(log10(u), &delta_u, temp_table,
                   cooling);
    T4d = eagle_convert_u_to_temp(log10(u), &delta_u, z_index, n_h_i, He_i,
        		 dz, d_n_h, d_He,cooling, cosmo);
    printf("u T1d T4d %.5e %.5e %.5e\n",u,pow(10.0,T1d),pow(10.0,T4d));
  }
}

/*
 * @brief Produces contributions to cooling rates for different 
 * hydrogen number densities, from different metals, 
 * tests 1d and 4d table interpolations produce 
 * same results for cooling rate, dlambda/du and temperature.
 *
 * @param pass string "nh" to produce cooling rates for one internal
 * energy, different redshifts
 * @param pass string "metals" to produce contribution to cooling
 * rates from different metals
 * @param pass string "compare_temp" to compare temperature interpolation
 * from 1d and 4d tables
 * @param pass string "compare_dlambda_du" to compare dlambda/du calculation
 * from 1d and 4d tables
 * @param if no string passed, all of the above are performed.
 */
int main(int argc, char **argv) {
  // Declare relevant structs
  struct swift_params *params = malloc(sizeof(struct swift_params));
  struct unit_system us;
  struct chemistry_global_data chem_data;
  struct part p;
  struct xpart xp;
  struct phys_const internal_const;
  struct cooling_function_data cooling;
  struct cosmology cosmo;
  char *parametersFileName = "./testCooling.yml";

  float nh;

  // Read options
  int param, tables = 0;
  float redshift = -1.0, log_10_nh = 100;
  while ((param = getopt(argc, argv, "z:d:t")) != -1)
  switch(param){
    case 'z':
      redshift = atof(optarg);
      break;
    case 'd':
      log_10_nh = atof(optarg);
      break;
    case 't':
      tables = 1;
      break;
    case '?':
      if (optopt == 'z')
        printf ("Option -%c requires an argument.\n", optopt);
      else
        printf ("Unknown option character `\\x%x'.\n",
                 optopt);
      error("invalid option(s) to testCooling");
  }

  // Read the parameter file
  if (params == NULL) error("Error allocating memory for the parameter file.");
  message("Reading runtime parameters from file '%s'", parametersFileName);
  parser_read_file(parametersFileName, params);

  // Init units 
  units_init_from_params(&us, params, "InternalUnitSystem");
  phys_const_init(&us, params, &internal_const);

  // Init chemistry 
  chemistry_init(params, &us, &internal_const, &chem_data);
  chemistry_first_init_part(&internal_const, &us, &cosmo, &chem_data, &p, &xp);
  chemistry_print(&chem_data);

  // Init cosmology 
  cosmology_init(params, &us, &internal_const, &cosmo);
  cosmology_print(&cosmo);
  if (redshift == -1.0) {
    cosmo.z = 3.0;
  } else {
    cosmo.z = redshift;
  }

  // Init cooling 
  cooling_init(params, &us, &internal_const, &cooling);
  cooling_print(&cooling);

  // Calculate abundance ratios 
  float *abundance_ratio;
  abundance_ratio = malloc((chemistry_element_count + 2)*sizeof(float));
  abundance_ratio_to_solar(&p, &cooling, abundance_ratio);
  
  // Declare 1D tables 
  double H_plus_He_heat_table[cooling.N_Temp];              
  double H_plus_He_electron_abundance_table[cooling.N_Temp];
  double temp_table[cooling.N_Temp];  
  double element_cooling_table[cooling.N_Temp];         
  double element_print_cooling_table[cooling.N_Elements * cooling.N_Temp];  
  double element_electron_abundance_table[cooling.N_Temp];
  
  // extract mass fractions, calculate table indices and offsets 
  float XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac = p.chemistry_data.metal_mass_fraction[chemistry_element_He] /
                 (XH + p.chemistry_data.metal_mass_fraction[chemistry_element_He]);
  int z_index,He_i,n_h_i;
  float dz,d_He,d_n_h;
  get_redshift_index(cosmo.z, &z_index, &dz, &cooling);
  get_index_1d(cooling.HeFrac, cooling.N_He, HeFrac, &He_i, &d_He);

  float upper_bound = cooling.Temp[cooling.N_Temp-1]/eagle_log_10_e;
  float lower_bound = cooling.Temp[0]/eagle_log_10_e;

  int nt = 250;//, n_nh = 6;
  double u = pow(10.0,10);
  //if (argc == 1 || strcmp(argv[1], "nh") == 0){
  // Calculate cooling rates at different densities 
    //for(int i = 0; i < n_nh; i++){
    //  // Open files 
    //  char output_filename[21];
    //  sprintf(output_filename, "%s%d%s", "cooling_output_", i, ".dat");
    //  FILE *output_file = fopen(output_filename, "w");
    //  if (output_file == NULL) {
    //    printf("Error opening file!\n");
    //    exit(1);
    //  }

    //  // set hydrogen number density, construct tables
    //  nh = pow(10.0,-i);
    //  set_quantities(&p, &us, &cooling, &cosmo, &internal_const, nh, u);
    //  construct_1d_tables_test(&p, &cooling, &cosmo, &internal_const, 
    //    			temp_table, H_plus_He_heat_table, 
    //      		H_plus_He_electron_abundance_table, 
    //      		element_electron_abundance_table, 
    //      		element_print_cooling_table, 
    //      		element_cooling_table, 
    //      		abundance_ratio, &upper_bound, &lower_bound);
    //  get_index_1d(cooling.nH, cooling.N_nH, log10(nh), &n_h_i, &d_n_h);

    //  for(int j = 0; j < nt; j++){
    //    // set internal energy
    //    set_quantities(&p, &us, &cooling, &cosmo, &internal_const, nh, pow(10.0,11.0 + j*8.0/nt));
    //    u = hydro_get_physical_internal_energy(&p,&cosmo)*cooling.internal_energy_scale;
    //    double dlambda_du;
    //    float delta_u, cooling_du_dt, logT;

    //    // calculate cooling rates using 1d tables
    //    if (tables == 1) {
    //      cooling_du_dt = eagle_metal_cooling_rate_1d_table(log10(u),&dlambda_du,H_plus_He_heat_table,
    //                              H_plus_He_electron_abundance_table,
    //         	              element_print_cooling_table,
    //         	              element_electron_abundance_table,
    //         	              temp_table,
    //         	              z_index, dz, n_h_i, d_n_h, He_i, d_He,
    //         	              &p,&cooling,&cosmo,&internal_const,NULL,
    //         	              abundance_ratio);
    //      logT = eagle_convert_u_to_temp_1d_table(log10(u), &delta_u, temp_table,
    //                     &p,&cooling, &cosmo, &internal_const);
    //    } else {
    //    // calculate cooling rates using 4d tables
    //      cooling_du_dt = eagle_metal_cooling_rate(log10(u),
    //         	              &dlambda_du,z_index, dz, n_h_i, d_n_h, He_i, d_He,
    //         	              &p,&cooling,&cosmo,&internal_const,NULL,
    //         	              abundance_ratio);
    //      logT = eagle_convert_u_to_temp(log10(u), &delta_u, z_index, n_h_i, He_i,
    //          		 dz, d_n_h, d_He, &p,&cooling, &cosmo, &internal_const);
    //    }
    //    float temperature_swift = pow(10.0,logT);

    //    fprintf(output_file,"%.5e %.5e\n", temperature_swift,cooling_du_dt);
    //  }
    //  fclose(output_file);
    //}
  //}
  // Calculate contributions from metals to cooling rate
  //if (argc >= 1 || strcmp(argv[1],"metals") == 0) {
    // open file
    char output_filename[21];
    sprintf(output_filename, "%s", "cooling_output.dat");
    FILE *output_file = fopen(output_filename, "w");
    if (output_file == NULL) {
      printf("Error opening file!\n");
      exit(1);
    }

    // set hydrogen number density, construct 1d tables
    if (log_10_nh == 100) {
      nh = 1.0e-1;
    } else {
      nh = pow(10.0,log_10_nh);
    }
    //if (argc > 2) nh = pow(10.0,strtod(argv[2],NULL));
    printf("Eagle testcooling.c nh %.5e\n", nh);
    u = pow(10.0,14.0);
    set_quantities(&p, &us, &cooling, &cosmo, &internal_const, nh, u);
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
    printf("Eagle testcooling.c inn_h %.5e\n", inn_h);
    get_index_1d(cooling.nH, cooling.N_nH, log10(inn_h), &n_h_i, &d_n_h);

    // Loop over internal energy
    for(int j = 0; j < nt; j++){
      set_quantities(&p, &us, &cooling, &cosmo, &internal_const, nh, pow(10.0,10.0 + j*8.0/nt));
      u = hydro_get_physical_internal_energy(&p,&cosmo)*cooling.internal_energy_scale;
      float delta_u, cooling_du_dt, logT;

      // calculate cooling rates using 1d tables
      if (tables == 1) {
        cooling_du_dt = eagle_print_metal_cooling_rate_1d_table(H_plus_He_heat_table,
                              H_plus_He_electron_abundance_table,
           	              element_print_cooling_table,
           	              element_electron_abundance_table,
           	              temp_table,
           	              &p,&cooling,&cosmo,&internal_const);
        logT = eagle_convert_u_to_temp_1d_table(log10(u), &delta_u, temp_table,
                       &cooling);
      } else {
      // calculate cooling rates using 4d tables
        cooling_du_dt = eagle_print_metal_cooling_rate(
			      z_index, dz, n_h_i, d_n_h, He_i, d_He,
           	              &p,&cooling,&cosmo,&internal_const,
           	              abundance_ratio);
          logT = eagle_convert_u_to_temp(log10(u), &delta_u, z_index, n_h_i, He_i,
	      		 dz, d_n_h, d_He,&cooling, &cosmo);
      }
      //float temperature_swift = pow(10.0,logT);
      //fprintf(output_file,"%.5e %.5e\n", temperature_swift,cooling_du_dt);
      fprintf(output_file,"%.5e %.5e\n", u,cooling_du_dt);
    }
    fclose(output_file);
  //}

  // compare temperatures and dlambda/du calculated from 1d and 4d tables
  //compare_temp(&us, &p, &xp, &internal_const, &cooling, &cosmo);
  //compare_dlambda_du(&us, &p, &xp, &internal_const, &cooling, &cosmo);

  free(params);
  return 0;
}


