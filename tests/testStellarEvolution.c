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
#include "../config.h"
#include "hydro.h"
#include "physical_constants.h"
#include "stars.h"
#include "swift.h"
#include "units.h"

int main(int argc, char **argv) {
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
  char *parametersFileName = "./testStellarEvolution.yml";

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

  /* Init spart */
  stars_first_init_spart(&sp);

  /* Evolve spart */
  double dt = 1.0e-6;
  float current_time = 0.f;
  stars_evolve_spart(&sp, &stars_properties, &cosmo, &us, current_time, dt);

  for (int i = 0; i < 9; i++) {
    message("element %d to distribute fraction %.5e", i,
            sp.to_distribute.chemistry_data.metal_mass_fraction[i]);
  }
  message("to distribute mass %.5e", sp.to_distribute.mass);
  message("to distribute num_SNIa %.5e", sp.to_distribute.num_SNIa);
  message("to distribute metal_mass_fraction_total %.5e",
          sp.to_distribute.chemistry_data.metal_mass_fraction_total);
  message("to distribute mass_from_AGB %.5e",
          sp.to_distribute.chemistry_data.mass_from_AGB);
  message("to distribute metal_mass_fraction_from_AGB %.5e",
          sp.to_distribute.chemistry_data.metal_mass_fraction_from_AGB);
  message("to distribute mass_from_SNII %.5e",
          sp.to_distribute.chemistry_data.mass_from_SNII);
  message("to distribute metal_mass_fraction_from_SNII %.5e",
          sp.to_distribute.chemistry_data.metal_mass_fraction_from_SNII);
  message("to distribute mass_from_SNIa %.5e",
          sp.to_distribute.chemistry_data.mass_from_SNIa);
  message("to distribute metal_mass_fraction_from_SNIa %.5e",
          sp.to_distribute.chemistry_data.metal_mass_fraction_from_SNIa);
  message("to distribute iron_mass_fraction_from_SNIa %.5e",
          sp.to_distribute.chemistry_data.iron_mass_fraction_from_SNIa);

  message("done test");

  free(params);
  return 0;
}
