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

/** 
 * @brief compute the relative error between two floats
 *
 * @param a, b numbers between which to compute the relative difference
 */
float relative_error(float a, float b){
  return fabs((a - b)/b);
}

/**
 * @brief Test function to check feedback is working correctly. Produces a table
 * of values for the cumulative amount of mass enrichment up to a given time from
 * the birth of the star (normalized by the initial stellar mass). Compares these
 * results to those produced by an analogous test in EAGLE.
 */
int main(int argc, char *argv[]) {
  /* Declare relevant structs */
  struct swift_params *params = malloc(sizeof(struct swift_params));
  struct unit_system us;
  struct chemistry_global_data chem_data;
  struct part p;
  struct xpart xp;
  struct spart sp;
  struct phys_const phys_const;
  struct cosmology cosmo;
  struct hydro_props hydro_properties;
  struct stars_props stars_properties;
  char *parametersFileName = "./testFeedback.yml";

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

  /* Init hydro properties */
  hydro_props_init(&hydro_properties, &phys_const, &us, params);

  /* Init star properties */
  stars_props_init(&stars_properties, &phys_const, &us, params,
                   &hydro_properties, &cosmo);

  /* Read yield tables */
  stars_evolve_init(params, &stars_properties);

  /* Init spart */
  stars_first_init_spart(&sp);

  sp.mass_init = 4.706273e-5;

  /* Set metal mass fractions */
  for (int i = 0; i < chemistry_element_count; i++)
    sp.chemistry_data.metal_mass_fraction[i] = 0.f;
  sp.chemistry_data.metal_mass_fraction[0] = 0.752;
  sp.chemistry_data.metal_mass_fraction[1] = 0.248;

  sp.chemistry_data.metal_mass_fraction_total = 0.01;

  float Gyr_to_s = 3.154e16;
  float dt = 0.1 * Gyr_to_s / units_cgs_conversion_factor(&us, UNIT_CONV_TIME);
  float max_age =
      13.f * Gyr_to_s / units_cgs_conversion_factor(&us, UNIT_CONV_TIME);

  //for (int i = 0; i < chemistry_element_count; i++) sp.to_distribute.metal_mass[i] = 0.f;
  //sp.to_distribute.metal_mass_from_AGB = 0.f;
  //sp.to_distribute.total_metal_mass = 0.f;
  //sp.to_distribute.mass = 0.f;
  //
  //FILE *AGB_output;
  //const char AGB_fname[25] = "test_feedback_AGB.txt";
  //if (!(AGB_output = fopen(AGB_fname, "w"))) {
  //  error("error in opening file '%s'\n", AGB_fname);
  //}
  //fprintf(AGB_output,
  //        "# time[Gyr] | total mass | metal mass: total | H | He | C | N  | O  "
  //        "| Ne | Mg | Si | Fe | per solar mass\n");

  //for (float age = 0; age <= max_age; age += dt) {
  //  compute_stellar_evolution(&stars_properties, &sp, &us, age, dt);
  //  float age_Gyr =
  //      age * units_cgs_conversion_factor(&us, UNIT_CONV_TIME) / Gyr_to_s;
  //  fprintf(AGB_output, "%f %e %e ", age_Gyr,
  //          sp.to_distribute.mass_from_AGB / sp.mass_init,
  //          sp.to_distribute.metal_mass_from_AGB / sp.mass_init);
  //  for (int i = 0; i < chemistry_element_count; i++)
  //    fprintf(AGB_output, "%e ", sp.to_distribute.metal_mass[i] / sp.mass_init);
  //  fprintf(AGB_output, "\n");
  //}

  //for (int i = 0; i < chemistry_element_count; i++) sp.to_distribute.metal_mass[i] = 0.f;
  //sp.to_distribute.metal_mass_from_SNII = 0.f;
  //sp.to_distribute.mass = 0.f;
  //
  //FILE *SNII_output;
  //const char SNII_fname[25] = "test_feedback_SNII.txt";
  //if (!(SNII_output = fopen(SNII_fname, "w"))) {
  //  error("error in opening file '%s'\n", SNII_fname);
  //}
  //fprintf(SNII_output,
  //        "# time[Gyr] | total mass | metal mass: total | H | He | C | N  | O  "
  //        "| Ne | Mg | Si | Fe | Number of SNII per solar mass\n");

  //for (float age = 0; age <= max_age; age += dt) {
  //  sp.to_distribute.num_SNII = 0.f;
  //  compute_stellar_evolution(&stars_properties, &sp, &us, age, dt);
  //  float age_Gyr =
  //      age * units_cgs_conversion_factor(&us, UNIT_CONV_TIME) / Gyr_to_s;
  //  fprintf(SNII_output, "%f %e %e ", age_Gyr,
  //          sp.to_distribute.mass_from_SNII / sp.mass_init,
  //          sp.to_distribute.metal_mass_from_SNII / sp.mass_init);
  //  for (int i = 0; i < chemistry_element_count; i++)
  //    fprintf(SNII_output, "%e ", sp.to_distribute.metal_mass[i] / sp.mass_init);
  //  fprintf(SNII_output, "%e ", sp.to_distribute.num_SNII );
  //  fprintf(SNII_output, "\n");
  //}
  
  //for (int i = 0; i < chemistry_element_count; i++) sp.to_distribute.metal_mass[i] = 0.f;
  //sp.to_distribute.metal_mass_from_SNIa = 0.f;
  //sp.to_distribute.mass = 0.f;

  //FILE *SNIa_output;
  //const char SNIa_fname[25] = "test_feedback_SNIa.txt";
  //if (!(SNIa_output = fopen(SNIa_fname, "w"))) {
  //  error("error in opening file '%s'\n", SNIa_fname);
  //}
  //fprintf(SNIa_output,
  //        "# time[Gyr] | total mass | metal mass: total | H | He | C | N  | O  "
  //        "| Ne | Mg | Si | Fe | Number of SNIa per solar mass\n");

  //for (float age = 0; age <= max_age; age += dt) {
  //  sp.to_distribute.mass_from_SNIa = 0.f;
  //  sp.to_distribute.num_SNIa = 0.f;
  //  compute_stellar_evolution(&stars_properties, &sp, &us, age, dt);
  //  float age_Gyr =
  //      age * units_cgs_conversion_factor(&us, UNIT_CONV_TIME) / Gyr_to_s;
  //  fprintf(SNIa_output, "%f %e %e ", age_Gyr,
  //          sp.to_distribute.mass / sp.mass_init,
  //          sp.to_distribute.metal_mass_from_SNIa / sp.mass_init);
  //  for (int i = 0; i < chemistry_element_count; i++)
  //    fprintf(SNIa_output, "%e ", sp.to_distribute.metal_mass[i] / sp.mass_init);
  //  fprintf(SNIa_output, "%e ", sp.to_distribute.num_SNIa / (sp.mass_init * stars_properties.feedback.const_solar_mass));
  //  fprintf(SNIa_output, "\n");
  //}
  
  for (int i = 0; i < chemistry_element_count; i++) sp.to_distribute.metal_mass[i] = 0.f;
  sp.to_distribute.metal_mass_from_SNIa = 0.f;
  sp.to_distribute.metal_mass_from_SNII = 0.f;
  sp.to_distribute.metal_mass_from_AGB = 0.f;
  sp.to_distribute.mass_from_AGB = 0.f;
  sp.to_distribute.mass_from_SNII = 0.f;
  sp.to_distribute.mass_from_SNIa = 0.f;
  sp.to_distribute.Fe_mass_from_SNIa = 0.f;
  sp.to_distribute.total_metal_mass = 0.f;
  sp.to_distribute.mass = 0.f;
  
  /* Open EAGLE test file for reading  */
  FILE *EAGLE_test;
  const char EAGLE_fname[75] = "/cosma/home/dp004/dc-bori1/Eagle/data1/z_0.01/StellarEvolutionTotal.txt";
  if (!(EAGLE_test = fopen(EAGLE_fname, "r"))) {
    error("error in opening file '%s'\n", EAGLE_fname);
  }

  /* Declare constants necessary for reading EAGLE data */
  char *line = NULL;
  size_t len = 0;
  const int n_fields = 19;
  float tol = 1e-5;

  /* Declare array to store one line of data from EAGLE test */
  float eagle_data[n_fields];

  /* read first line */
  if (getline(&line, &len, EAGLE_test) == -1) error("failed to read first line of EAGLE test file");

  /* Open file for writing SWIFT feedback test output */
  FILE *Total_output;
  const char Total_fname[25] = "test_feedback_total.txt";
  if (!(Total_output = fopen(Total_fname, "w"))) {
    error("error in opening file '%s'\n", Total_fname);
  }
  fprintf(Total_output,
          "# time[Gyr] | total mass | metal mass: total | H | He | C | N  | O  "
          "| Ne | Mg | Si | Fe | per solar mass (m,z)_AGB (m,z)_SNII (m,z,M_fe)_SNIa \n");

  for (float age = 0; age <= max_age; age += dt) {
    compute_stellar_evolution(&stars_properties, &sp, &us, age, dt);

    float age_Gyr =
        age * units_cgs_conversion_factor(&us, UNIT_CONV_TIME) / Gyr_to_s;
    fprintf(Total_output, "%f %e %e ", age_Gyr,
            sp.to_distribute.mass / sp.mass_init,
            sp.to_distribute.total_metal_mass / sp.mass_init);
    for (int i = 0; i < chemistry_element_count; i++)
      fprintf(Total_output, "%e ", sp.to_distribute.metal_mass[i] / sp.mass_init);
    fprintf(Total_output, " %e %e %e %e %e %e %e", sp.to_distribute.mass_from_AGB / sp.mass_init, 
            sp.to_distribute.metal_mass_from_AGB / sp.mass_init,
            sp.to_distribute.mass_from_SNII / sp.mass_init, 
            sp.to_distribute.metal_mass_from_SNII / sp.mass_init, 
            sp.to_distribute.mass_from_SNIa / sp.mass_init, 
            sp.to_distribute.metal_mass_from_SNIa / sp.mass_init,
            sp.to_distribute.Fe_mass_from_SNIa / sp.mass_init);
    fprintf(Total_output, "\n");

    /* Read data from EAGLE test and compare it to what we have */
    if(!feof(EAGLE_test)){
      fscanf(EAGLE_test, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ", 
        &eagle_data[0], &eagle_data[1], &eagle_data[2], &eagle_data[3], &eagle_data[4], &eagle_data[5], 
        &eagle_data[6], &eagle_data[7], &eagle_data[8], &eagle_data[9], &eagle_data[10], &eagle_data[11], 
        &eagle_data[12], &eagle_data[13], &eagle_data[14], &eagle_data[15], &eagle_data[16], &eagle_data[17], 
        &eagle_data[18]);
      if (relative_error(age_Gyr,eagle_data[0]) > tol) 
        error("relative error in age greater than tolerance. Swift: %e Eagle %e", age_Gyr, eagle_data[0]);
      if (relative_error(sp.to_distribute.mass/sp.mass_init,eagle_data[1]) > tol) 
        error("relative error in total mass greater than tolerance. Swift: %e Eagle %e", sp.to_distribute.mass/sp.mass_init, eagle_data[1]);
      if (relative_error(sp.to_distribute.total_metal_mass/sp.mass_init,eagle_data[2]) > tol) 
        error("relative error in total metal mass greater than tolerance. Swift: %e Eagle %e", sp.to_distribute.total_metal_mass/sp.mass_init, eagle_data[2]);
      for (int i = 0; i < chemistry_element_count; i++)
        if (relative_error(sp.to_distribute.metal_mass[i]/sp.mass_init,eagle_data[3+i]) > tol) 
          error("relative error in mass released for element %d greater than tolerance. Swift: %e Eagle %e", i, age_Gyr, eagle_data[3+i]);
      if (relative_error(sp.to_distribute.mass_from_AGB/sp.mass_init,eagle_data[12]) > tol) 
        error("relative error in mass from AGB greater than tolerance. Swift: %e Eagle %e", sp.to_distribute.mass_from_AGB/sp.mass_init, eagle_data[12]);
      if (relative_error(sp.to_distribute.metal_mass_from_AGB/sp.mass_init,eagle_data[13]) > tol) 
        error("relative error in metal mass from AGB greater than tolerance. Swift: %e Eagle %e", sp.to_distribute.metal_mass_from_AGB/sp.mass_init, eagle_data[13]);
      if (relative_error(sp.to_distribute.mass_from_SNII/sp.mass_init,eagle_data[14]) > tol) 
        error("relative error in mass from SNII greater than tolerance. Swift: %e Eagle %e", sp.to_distribute.mass_from_SNII/sp.mass_init, eagle_data[14]);
      if (relative_error(sp.to_distribute.metal_mass_from_SNII/sp.mass_init,eagle_data[15]) > tol) 
        error("relative error in metal mass from SNII greater than tolerance. Swift: %e Eagle %e", sp.to_distribute.metal_mass_from_SNII/sp.mass_init, eagle_data[15]);
      if (relative_error(sp.to_distribute.mass_from_SNIa/sp.mass_init,eagle_data[16]) > tol) 
        error("relative error in mass from SNIa greater than tolerance. Swift: %e Eagle %e", sp.to_distribute.mass_from_SNIa/sp.mass_init, eagle_data[16]);
      if (relative_error(sp.to_distribute.metal_mass_from_SNIa/sp.mass_init,eagle_data[17]) > tol) 
        error("relative error in metal mass from SNIa greater than tolerance. Swift: %e Eagle %e", sp.to_distribute.metal_mass_from_SNIa/sp.mass_init, eagle_data[17]);
      if (relative_error(sp.to_distribute.Fe_mass_from_SNIa/sp.mass_init,eagle_data[18]) > tol) 
        error("relative error in iron mass from SNIa greater than tolerance. Swift: %e Eagle %e", sp.to_distribute.Fe_mass_from_SNIa/sp.mass_init, eagle_data[18]);
    } else {
      error("Failed to read line of EAGLE test file data");
    }
  }


  return 0;
}
