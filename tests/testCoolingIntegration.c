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
  // char *parametersFileName = "../examples/CoolingBox/coolingBox.yml";
  char *parametersFileName = "../examples/EAGLE_12/eagle_12.yml";
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

  char output_filename[40];
  FILE **output_file = malloc(10 * sizeof(FILE *));

  /* And dump the parameters as used. */
  parser_write_params_to_file(params, "used_parameters.yml");

  units_init(&us, params, "InternalUnitSystem");
  phys_const_init(&us, params, &internal_const);

  double hydrogen_number_density_cgs = 1.582e-3;
  // double hydrogen_number_density_cgs = 1.0e-1;
  double u, u_cgs, hydrogen_number_density, pressure, gamma, cooling_du_dt,
      temperature_cgs, newton_func, ratefact;
  u_cgs = 0;
  cooling_du_dt = 0;
  temperature_cgs = 0;
  newton_func = 0;
  // int n_t_i = 2000;

  gamma = 5.0 / 3.0;

  chemistry_init(params, &us, &internal_const, &chemistry_data);
  chemistry_first_init_part(&p, &xp, &chemistry_data);
  chemistry_print(&chemistry_data);

  cosmology_init(params, &us, &internal_const, &cosmo);
  cosmology_print(&cosmo);

  float u_ini = 3.357e15;
  u = u_ini / (units_cgs_conversion_factor(&us, UNIT_CONV_ENERGY) /
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
  cosmo.z = 0.0999744;

  cooling_init(params, &us, &internal_const, &cooling);
  cooling_print(&cooling);

  // construct 1d table of cooling rates wrt temperature
  float H_plus_He_heat_table[176];                // WARNING sort out how it is
                                                  // declared/allocated
  float H_plus_He_electron_abundance_table[176];  // WARNING sort out how it is
                                                  // declared/allocated
  float temp_table[176];  // WARNING sort out how it is declared/allocated
  float element_cooling_table[176];             // WARNING sort out how it is
                                                // declared/allocated
  float element_electron_abundance_table[176];  // WARNING sort out how it is
                                                // declared/allocated
  for (int k = 0; k < cooling.N_Temp; k++) H_plus_He_heat_table[k] = 0.0;
  for (int k = 0; k < cooling.N_Temp; k++)
    H_plus_He_electron_abundance_table[k] = 0.0;
  for (int k = 0; k < cooling.N_Temp; k++) temp_table[k] = 0.0;
  for (int k = 0; k < cooling.N_Temp; k++) element_cooling_table[k] = 0.0;
  for (int k = 0; k < cooling.N_Temp; k++)
    element_electron_abundance_table[k] = 0.0;
  // construct_1d_table_from_4d(&p,&cooling,&cosmo,&internal_const,cooling.table.element_cooling.H_plus_He_heating,H_plus_He_heat_table);

  float XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];
  // float inn_h =
  // chemistry_get_number_density(&p,&cosmo,chemistry_element_H,&internal_const)*cooling.number_density_scale;
  float inn_h = hydro_get_physical_density(&p, &cosmo) *
                units_cgs_conversion_factor(&us, UNIT_CONV_DENSITY) * XH /
                eagle_proton_mass_cgs;
  ratefact = inn_h * (XH / eagle_proton_mass_cgs);
  printf("XH inn_h ratefact %.5e %.5e %.5e\n", XH, inn_h, ratefact);

  double dLambdaNet_du, LambdaNet;
  float x_init;
  for (int j = 0; j < 1; j++) {
    // float u_ini =
    // eagle_convert_temp_to_u_1d_table(pow(10.0,0.5*(j+5)),temp_table,&p,&cooling,&cosmo,&internal_const);
    // float u_ini =  3.357e15;
    u = u_ini / (units_cgs_conversion_factor(&us, UNIT_CONV_ENERGY) /
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

    float x, du;
    double u_temp;
    // float dt = 2.0e-4*units_cgs_conversion_factor(&us,UNIT_CONV_TIME);
    float dt = 1.73798e-3 * units_cgs_conversion_factor(&us, UNIT_CONV_TIME);
    construct_1d_table_from_4d_H_He(
        &p, &cooling, &cosmo, &internal_const,
        cooling.table.element_cooling.H_plus_He_heating, H_plus_He_heat_table,
        &u_temp, u_ini, dt);
    construct_1d_table_from_4d(
        &p, &cooling, &cosmo, &internal_const,
        cooling.table.element_cooling.H_plus_He_electron_abundance,
        H_plus_He_electron_abundance_table);
    construct_1d_table_from_4d(&p, &cooling, &cosmo, &internal_const,
                               cooling.table.element_cooling.temperature,
                               temp_table);
    construct_1d_table_from_4d_elements(
        &p, &cooling, &cosmo, &internal_const,
        cooling.table.element_cooling.metal_heating, element_cooling_table);
    construct_1d_table_from_3d(&p, &cooling, &cosmo, &internal_const,
                               cooling.table.element_cooling.electron_abundance,
                               element_electron_abundance_table);

    sprintf(output_filename, "%s%d%s", "cooling_integration_output_", j,
            ".dat");
    output_file[j] = fopen(output_filename, "w");
    if (output_file[j] == NULL) {
      printf("Error opening file!\n");
      exit(1);
    }
    for (int k = 0; k < cooling.N_Temp - 1; k++) {
      float lambda1 = H_plus_He_heat_table[k] +
                      element_cooling_table[k] *
                          H_plus_He_electron_abundance_table[k] /
                          element_electron_abundance_table[k];
      float lambda2 = H_plus_He_heat_table[k + 1] +
                      element_cooling_table[k + 1] *
                          H_plus_He_electron_abundance_table[k + 1] /
                          element_electron_abundance_table[k + 1];
      // printf("temperature %.5e, internal energy
      // %.5e\n",pow(10.0,temp_table[k]),
      // eagle_convert_temp_to_u_1d_table(pow(10.0,temp_table[k]),temp_table,&p,&cooling,&cosmo,&internal_const));
      float u2 = eagle_convert_temp_to_u_1d_table(pow(10.0, temp_table[k + 1]),
                                                  temp_table, &p, &cooling,
                                                  &cosmo, &internal_const);
      float u1 = eagle_convert_temp_to_u_1d_table(pow(10.0, temp_table[k]),
                                                  temp_table, &p, &cooling,
                                                  &cosmo, &internal_const);
      float delta_u = u2 - u1;
      // fprintf(output_file[j], "%.5e %.5e %.5e\n", pow(10.0,temp_table[k]), u1
      // - u_ini - lambda1*ratefact*dt, (lambda2 - lambda1)/delta_u);
      fprintf(output_file[j], "%.5e %.5e %.5e\n", u1,
              u1 - u_ini - lambda1 * ratefact * dt,
              (lambda2 - lambda1) / delta_u);
    }
    fclose(output_file[j]);

    LambdaNet = eagle_cooling_rate_1d_table(
        u_ini, &dLambdaNet_du, H_plus_He_heat_table,
        H_plus_He_electron_abundance_table, element_cooling_table,
        element_electron_abundance_table, temp_table, &p, &cooling, &cosmo,
        &internal_const);
    double u_eq = 5.0e12;
    double u_temp_guess = u_ini + LambdaNet * ratefact * dt;
    printf(
        "u_guess, u_temp_guess, u_ini, LambdaNet, dLambdaNet_du %.5e %.5e %.5e "
        "%.5e %.5e %.5e %.5e \n",
        u_temp, u_temp_guess, u_ini, LambdaNet, dLambdaNet_du, ratefact, dt);

    // if ((LambdaNet < 0 && u_temp < u_temp_guess) ||
    //    (LambdaNet >= 0 && u_temp > u_temp_guess))
    //  u_temp = u_temp_guess;
    u_temp = u_temp_guess;
    if ((LambdaNet < 0 && u_temp < u_eq) || (LambdaNet >= 0 && u_temp > u_eq))
      u_temp = u_eq;

    x_init = log(u_temp);
    // int printflag = 1;
    // x =
    // newton_iter(x_init,u_ini,H_plus_He_heat_table,H_plus_He_electron_abundance_table,element_cooling_table,element_electron_abundance_table,temp_table,&p,&cosmo,&cooling,&internal_const,dt,&printflag);
    x = bisection_iter(x_init, u_ini, H_plus_He_heat_table,
                       H_plus_He_electron_abundance_table,
                       element_cooling_table, element_electron_abundance_table,
                       temp_table, &p, &cosmo, &cooling, &internal_const, dt);
    // printf("testing newton integration, u_ini, u %.5e %.5e, temperature
    // initial, final %.5e %.5e\n", u_ini, exp(x),
    // eagle_convert_u_to_temp_1d_table(u_ini,&du,temp_table,&p,&cooling,&cosmo,&internal_const),
    // eagle_convert_u_to_temp_1d_table(exp(x),&du,temp_table,&p,&cooling,&cosmo,&internal_const));

    u = u_ini;
    double u_next;
    int nt = 10000;
    float dt_sub = dt / nt;
    for (int t = 0; t < nt; t++) {
      LambdaNet = eagle_cooling_rate_1d_table(
          u, &dLambdaNet_du, H_plus_He_heat_table,
          H_plus_He_electron_abundance_table, element_cooling_table,
          element_electron_abundance_table, temp_table, &p, &cooling, &cosmo,
          &internal_const);
      u_next = u + LambdaNet * ratefact * dt_sub;
      printf(
          "here u_next u lambda_net ratefact dt_sub, du t %.5e %.5e %.5e %.5e "
          "%.5e %.5e %d\n",
          u_next, u, LambdaNet, ratefact, dt_sub, LambdaNet * ratefact * dt_sub,
          t);
      u = u_next;
    }
    printf(
        "testing newton integration, u_ini, u, u subcycle %.5e %.5e %.5e, "
        "temperature initial, final, subcycled %.5e %.5e %.5e\n",
        u_ini, exp(x), u,
        eagle_convert_u_to_temp_1d_table(u_ini, &du, temp_table, &p, &cooling,
                                         &cosmo, &internal_const),
        eagle_convert_u_to_temp_1d_table(exp(x), &du, temp_table, &p, &cooling,
                                         &cosmo, &internal_const),
        eagle_convert_u_to_temp_1d_table(u, &du, temp_table, &p, &cooling,
                                         &cosmo, &internal_const));
  }

  free(params);

  return 0;
}
