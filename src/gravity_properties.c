/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* This object's header. */
#include "gravity_properties.h"

/* Standard headers */
#include <float.h>
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "dimension.h"
#include "error.h"
#include "gravity.h"
#include "kernel_gravity.h"

#define gravity_props_default_a_smooth 1.25f
#define gravity_props_default_r_cut_max 4.5f
#define gravity_props_default_r_cut_min 0.1f

void gravity_props_init(struct gravity_props *p,
                        const struct swift_params *params) {

  /* Tree-PM parameters */
  p->a_smooth = parser_get_opt_param_float(params, "Gravity:a_smooth",
                                           gravity_props_default_a_smooth);
  p->r_cut_max = parser_get_opt_param_float(params, "Gravity:r_cut_max",
                                            gravity_props_default_r_cut_max);
  p->r_cut_min = parser_get_opt_param_float(params, "Gravity:r_cut_min",
                                            gravity_props_default_r_cut_min);

  /* Time integration */
  p->eta = parser_get_param_float(params, "Gravity:eta");

  /* Opening angle */
  p->theta_crit = parser_get_param_double(params, "Gravity:theta");
  if (p->theta_crit >= 1.) error("Theta too large. FMM won't converge.");
  p->theta_crit2 = p->theta_crit * p->theta_crit;
  p->theta_crit_inv = 1. / p->theta_crit;

  /* Softening lengths */
  p->epsilon = 3. * parser_get_param_double(params, "Gravity:epsilon");
  p->epsilon2 = p->epsilon * p->epsilon;
  p->epsilon_inv = 1.f / p->epsilon;
  p->epsilon_inv3 = p->epsilon_inv * p->epsilon_inv * p->epsilon_inv;
}

void gravity_props_print(const struct gravity_props *p) {

  message("Self-gravity scheme: FMM-MM with m-poles of order %d",
          SELF_GRAVITY_MULTIPOLE_ORDER);

  message("Self-gravity time integration: eta=%.4f", p->eta);

  message("Self-gravity opening angle:  theta=%.4f", p->theta_crit);

  message("Self-gravity softening:    epsilon=%.4f (Plummer equivalent: %.4f)",
          p->epsilon, p->epsilon / 3.);

  message("Self-gravity mesh smoothing-scale: a_smooth=%f", p->a_smooth);

  message("Self-gravity tree cut-off: r_cut_max=%f", p->r_cut_max);
  message("Self-gravity truncation cut-off: r_cut_min=%f", p->r_cut_min);
}

#if defined(HAVE_HDF5)
void gravity_props_print_snapshot(hid_t h_grpgrav,
                                  const struct gravity_props *p) {

  io_write_attribute_f(h_grpgrav, "Time integration eta", p->eta);
  io_write_attribute_f(h_grpgrav, "Softening length", p->epsilon);
  io_write_attribute_f(h_grpgrav, "Softening length (Plummer equivalent)",
                       p->epsilon / 3.);
  io_write_attribute_f(h_grpgrav, "Opening angle", p->theta_crit);
  io_write_attribute_d(h_grpgrav, "MM order", SELF_GRAVITY_MULTIPOLE_ORDER);
  io_write_attribute_f(h_grpgrav, "Mesh a_smooth", p->a_smooth);
  io_write_attribute_f(h_grpgrav, "Mesh r_cut_max", p->r_cut_max);
  io_write_attribute_f(h_grpgrav, "Mesh r_cut_min", p->r_cut_min);
}
#endif
