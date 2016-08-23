#error : this file is no longer in use!
/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Richard Bower (r.g.bower@durham.ac.uk)
 *                    Stefan Arridge  (stefan.arridge@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"
#include "../sourceterms.h"

void supernova_feedback_init(const struct swift_params* parameter_file,
                             struct UnitSystem* us, struct supernova* sn) {
  sn.time = parser_get_param_double(parameter_file, "SN:time");
  sn.energy = parser_get_param_double(parameter_file, "SN:energy");
  sn.x = parser_get_param_double(parameter_file, "SN:x");
  sn.y = parser_get_param_double(parameter_file, "SN:y");
  sn.z = parser_get_param_double(parameter_file, "SN:z");
}

void supernova_feedback_print(const struct supernova* sn) {
  message(
      " Single SNe of energy= %e will explode at time= %e at location "
      "(%e,%e,%e)",
      sn.energy, sn.time, sn.x, sn.y, sn.z);
};

__attribute__((always_inline)) INLINE static void do_supernova_feedback(
    const struct sourceterms* sourceterms, struct part* p){};

__attribute__((always_inline)) INLINE void update_entropy(
    const sourceterms* sourceterms, struct part* p) {

  /*updates the entropy of a particle due to feedback */
  float u_old;
  float u_new;
  float new_entropy;
  float old_entropy = p->entropy;
  float rho = p->rho;

  //  u_old = old_entropy/(GAMMA_MINUS1) * pow(rho,GAMMA_MINUS1);
  const float u_old =
      hydro_get_internal_energy(p, 0);  // dt = 0 because using current entropy
  const float u_new = u_old + sourceterms->supernova.energy;
  const float new_entropy =
      u_new * pow_minus_gamma_minus_one(-p > rho) * hydro_gamma_minus_one;
  p->entropy = new_entropy;
}
