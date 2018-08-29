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
  
  // extract mass fractions, calculate table indices and offsets 
  float XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac = p.chemistry_data.metal_mass_fraction[chemistry_element_He] /
                 (XH + p.chemistry_data.metal_mass_fraction[chemistry_element_He]);
  int z_index,He_i,n_h_i;
  float dz,d_He,d_n_h;
  get_redshift_index(cosmo.z, &z_index, &dz, &cooling);
  get_index_1d(cooling.HeFrac, cooling.N_He, HeFrac, &He_i, &d_He);

  int nt = 250;
  double u = pow(10.0,10);

  // Calculate contributions from metals to cooling rate
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
  u = pow(10.0,14.0);
  set_quantities(&p, &us, &cooling, &cosmo, &internal_const, nh, u);
  float inn_h = chemistry_get_number_density(&p, &cosmo, chemistry_element_H,
                                           &internal_const) *
              cooling.number_density_scale;
  get_index_1d(cooling.nH, cooling.N_nH, log10(inn_h), &n_h_i, &d_n_h);

  // Loop over internal energy
  for(int j = 0; j < nt; j++){
    set_quantities(&p, &us, &cooling, &cosmo, &internal_const, nh, pow(10.0,10.0 + j*8.0/nt));
    u = hydro_get_physical_internal_energy(&p,&cosmo)*cooling.internal_energy_scale;
    float cooling_du_dt;

    // calculate cooling rates using 4d tables
    cooling_du_dt = eagle_print_metal_cooling_rate(
    		      z_index, dz, n_h_i, d_n_h, He_i, d_He,
       	              &p,&cooling,&cosmo,&internal_const,
       	              abundance_ratio);
    fprintf(output_file,"%.5e %.5e\n", u,cooling_du_dt);
  }
  fclose(output_file);
  printf("done\n");

  free(params);
  return 0;
}


