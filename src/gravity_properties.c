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
#include "kernel_long_gravity.h"

#define gravity_props_default_a_smooth 1.25f
#define gravity_props_default_r_cut_max 4.5f
#define gravity_props_default_r_cut_min 0.1f
#define gravity_props_default_rebuild_frequency 0.01f

void gravity_props_init(struct gravity_props *p, struct swift_params *params,
                        const struct phys_const *phys_const,
                        const struct cosmology *cosmo, const int with_cosmology,
                        const int is_zoom_simulation, const int periodic) {

  /* Tree updates */
  p->rebuild_frequency =
      parser_get_opt_param_float(params, "Gravity:rebuild_frequency",
                                 gravity_props_default_rebuild_frequency);

  if (p->rebuild_frequency < 0.f || p->rebuild_frequency > 1.f)
    error("Invalid tree rebuild frequency. Must be in [0., 1.]");

  /* Tree-PM parameters */
  if (periodic) {
    p->mesh_size = parser_get_param_int(params, "Gravity:mesh_side_length");
    p->a_smooth = parser_get_opt_param_float(params, "Gravity:a_smooth",
                                             gravity_props_default_a_smooth);
    p->r_cut_max_ratio = parser_get_opt_param_float(
        params, "Gravity:r_cut_max", gravity_props_default_r_cut_max);
    p->r_cut_min_ratio = parser_get_opt_param_float(
        params, "Gravity:r_cut_min", gravity_props_default_r_cut_min);

    /* Some basic checks of what we read */
    if (p->mesh_size % 2 != 0)
      error("The mesh side-length must be an even number.");

    if (p->a_smooth <= 0.)
      error("The mesh smoothing scale 'a_smooth' must be > 0.");

    if (2. * p->a_smooth * p->r_cut_max_ratio > p->mesh_size)
      error("Mesh too small given r_cut_max. Should be at least %d cells wide.",
            (int)(2. * p->a_smooth * p->r_cut_max_ratio) + 1);
  } else {
    p->mesh_size = 0;
    p->a_smooth = 0.f;
    p->r_cut_min_ratio = 0.f;
    p->r_cut_max_ratio = 0.f;
  }

  /* Time integration */
  p->eta = parser_get_param_float(params, "Gravity:eta");

  /* Opening angle */
  p->theta_crit = parser_get_param_double(params, "Gravity:theta");
  if (p->theta_crit >= 1.) error("Theta too large. FMM won't converge.");
  p->theta_crit2 = p->theta_crit * p->theta_crit;
  p->theta_crit_inv = 1. / p->theta_crit;

  /* Softening parameters */
  if (with_cosmology) {

    /* Maximal physical softening taken straight from the parameter file */
    p->epsilon_DM_max_physical =
        parser_get_param_double(params, "Gravity:max_physical_DM_softening");
    p->epsilon_baryon_max_physical = parser_get_param_double(
        params, "Gravity:max_physical_baryon_softening");

    /* Co-moving softenings taken straight from the parameter file */
    p->epsilon_DM_comoving =
        parser_get_param_double(params, "Gravity:comoving_DM_softening");
    p->epsilon_baryon_comoving =
        parser_get_param_double(params, "Gravity:comoving_baryon_softening");

    if (is_zoom_simulation) {

      /* Compute the comoving softening length for background particles as
       * a fraction of the mean inter-particle density of the background DM
       * particles Since they have variable masses the mass factor will be
       * multiplied in later on. Note that we already multiply in the conversion
       * from Plummer -> real softening length */
      const double ratio_background =
          parser_get_param_double(params, "Gravity:softening_ratio_background");

      const double mean_matter_density =
          cosmo->Omega_m * cosmo->critical_density_0;

      p->epsilon_background_fac = kernel_gravity_softening_plummer_equivalent *
                                  ratio_background *
                                  cbrt(1. / mean_matter_density);
    }

  } else {

    p->epsilon_DM_max_physical =
        parser_get_param_double(params, "Gravity:max_physical_DM_softening");
    p->epsilon_baryon_max_physical = parser_get_param_double(
        params, "Gravity:max_physical_baryon_softening");

    p->epsilon_DM_comoving = p->epsilon_DM_max_physical;
    p->epsilon_baryon_comoving = p->epsilon_baryon_max_physical;
  }

  /* Copy over the gravitational constant */
  p->G_Newton = phys_const->const_newton_G;

  /* Set the softening to the current time */
  gravity_props_update(p, cosmo);
}

void gravity_props_update(struct gravity_props *p,
                          const struct cosmology *cosmo) {

  /* Current softening length for the high-res. DM particles. */
  double DM_softening, baryon_softening;
  if (p->epsilon_DM_comoving * cosmo->a > p->epsilon_DM_max_physical)
    DM_softening = p->epsilon_DM_max_physical / cosmo->a;
  else
    DM_softening = p->epsilon_DM_comoving;

  /* Current softening length for the high-res. baryon particles. */
  if (p->epsilon_baryon_comoving * cosmo->a > p->epsilon_baryon_max_physical)
    baryon_softening = p->epsilon_baryon_max_physical / cosmo->a;
  else
    baryon_softening = p->epsilon_baryon_comoving;

  /* Plummer equivalent -> internal */
  DM_softening *= kernel_gravity_softening_plummer_equivalent;
  baryon_softening *= kernel_gravity_softening_plummer_equivalent;

  /* Store things */
  p->epsilon_DM_cur = DM_softening;
  p->epsilon_baryon_cur = baryon_softening;
}

void gravity_props_print(const struct gravity_props *p) {

  message("Self-gravity scheme: %s", GRAVITY_IMPLEMENTATION);

  message("Self-gravity scheme: FMM-MM with m-poles of order %d",
          SELF_GRAVITY_MULTIPOLE_ORDER);

  message("Self-gravity time integration: eta=%.4f", p->eta);

  message("Self-gravity opening angle:  theta=%.4f", p->theta_crit);

  message("Self-gravity softening functional form: %s",
          kernel_gravity_softening_name);

  message(
      "Self-gravity DM comoving softening: epsilon=%.6f (Plummer equivalent: "
      "%.6f)",
      p->epsilon_DM_comoving * kernel_gravity_softening_plummer_equivalent,
      p->epsilon_DM_comoving);

  message(
      "Self-gravity DM maximal physical softening:    epsilon=%.6f (Plummer "
      "equivalent: %.6f)",
      p->epsilon_DM_max_physical * kernel_gravity_softening_plummer_equivalent,
      p->epsilon_DM_max_physical);

  message(
      "Self-gravity baryon comoving softening: epsilon=%.6f (Plummer "
      "equivalent: "
      "%.6f)",
      p->epsilon_baryon_comoving * kernel_gravity_softening_plummer_equivalent,
      p->epsilon_baryon_comoving);

  message(
      "Self-gravity baryon maximal physical softening:    epsilon=%.6f "
      "(Plummer "
      "equivalent: %.6f)",
      p->epsilon_baryon_max_physical *
          kernel_gravity_softening_plummer_equivalent,
      p->epsilon_baryon_max_physical);

  message("Self-gravity mesh side-length: N=%d", p->mesh_size);
  message("Self-gravity mesh smoothing-scale: a_smooth=%f", p->a_smooth);

  message("Self-gravity tree cut-off ratio: r_cut_max=%f", p->r_cut_max_ratio);
  message("Self-gravity truncation cut-off ratio: r_cut_min=%f",
          p->r_cut_min_ratio);

  message("Self-gravity mesh truncation function: %s",
          kernel_long_gravity_truncation_name);

  message("Self-gravity tree update frequency: f=%f", p->rebuild_frequency);
}

#if defined(HAVE_HDF5)
void gravity_props_print_snapshot(hid_t h_grpgrav,
                                  const struct gravity_props *p) {

  io_write_attribute_f(h_grpgrav, "Time integration eta", p->eta);
  io_write_attribute_s(h_grpgrav, "Softening style",
                       kernel_gravity_softening_name);

  io_write_attribute_f(
      h_grpgrav, "Comoving DM softening length [internal units]",
      p->epsilon_DM_comoving * kernel_gravity_softening_plummer_equivalent);
  io_write_attribute_f(
      h_grpgrav,
      "Comoving DM softening length (Plummer equivalent)  [internal units]",
      p->epsilon_DM_comoving);

  io_write_attribute_f(
      h_grpgrav, "Maximal physical DM softening length  [internal units]",
      p->epsilon_DM_max_physical * kernel_gravity_softening_plummer_equivalent);
  io_write_attribute_f(h_grpgrav,
                       "Maximal physical DM softening length (Plummer "
                       "equivalent) [internal units]",
                       p->epsilon_DM_max_physical);

  io_write_attribute_f(
      h_grpgrav, "Comoving baryon softening length [internal units]",
      p->epsilon_baryon_comoving * kernel_gravity_softening_plummer_equivalent);
  io_write_attribute_f(
      h_grpgrav,
      "Comoving baryon softening length (Plummer equivalent)  [internal units]",
      p->epsilon_baryon_comoving);

  io_write_attribute_f(
      h_grpgrav, "Maximal physical baryon softening length  [internal units]",
      p->epsilon_baryon_max_physical *
          kernel_gravity_softening_plummer_equivalent);
  io_write_attribute_f(h_grpgrav,
                       "Maximal physical baryon softening length (Plummer "
                       "equivalent) [internal units]",
                       p->epsilon_baryon_max_physical);

  io_write_attribute_f(h_grpgrav, "Opening angle", p->theta_crit);
  io_write_attribute_s(h_grpgrav, "Scheme", GRAVITY_IMPLEMENTATION);
  io_write_attribute_i(h_grpgrav, "MM order", SELF_GRAVITY_MULTIPOLE_ORDER);
  io_write_attribute_f(h_grpgrav, "Mesh a_smooth", p->a_smooth);
  io_write_attribute_f(h_grpgrav, "Mesh r_cut_max ratio", p->r_cut_max_ratio);
  io_write_attribute_f(h_grpgrav, "Mesh r_cut_min ratio", p->r_cut_min_ratio);
  io_write_attribute_f(h_grpgrav, "Tree update frequency",
                       p->rebuild_frequency);
  io_write_attribute_s(h_grpgrav, "Mesh truncation function",
                       kernel_long_gravity_truncation_name);
}
#endif

/**
 * @brief Write a gravity_props struct to the given FILE as a stream of bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void gravity_props_struct_dump(const struct gravity_props *p, FILE *stream) {
  restart_write_blocks((void *)p, sizeof(struct gravity_props), 1, stream,
                       "gravity", "gravity props");
}

/**
 * @brief Restore a gravity_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void gravity_props_struct_restore(struct gravity_props *p, FILE *stream) {
  restart_read_blocks((void *)p, sizeof(struct gravity_props), 1, stream, NULL,
                      "gravity props");
}
