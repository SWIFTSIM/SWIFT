/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <string.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "dimension.h"
#include "error.h"
#include "gravity.h"
#include "kernel_gravity.h"
#include "kernel_long_gravity.h"
#include "restart.h"

#define gravity_props_default_a_smooth 1.25f
#define gravity_props_default_r_cut_max 4.5f
#define gravity_props_default_r_cut_min 0.0f
#define gravity_props_default_rebuild_frequency 0.01f
#define gravity_props_default_rebuild_active_fraction 1.01f  // > 1 means never
#define gravity_props_default_distributed_mesh 0
#define gravity_props_default_max_adaptive_softening FLT_MAX
#define gravity_props_default_min_adaptive_softening 0.f

void gravity_props_init(struct gravity_props *p, struct swift_params *params,
                        const struct phys_const *phys_const,
                        const struct cosmology *cosmo, const int with_cosmology,
                        const int with_external_potential,
                        const int has_baryons, const int has_DM,
                        const int has_neutrinos, const int is_zoom_simulation,
                        const int periodic, const double dim[3],
                        const int cdim[3]) {

  /* Tree updates */
  p->rebuild_frequency =
      parser_get_opt_param_float(params, "Gravity:rebuild_frequency",
                                 gravity_props_default_rebuild_frequency);

  p->rebuild_active_fraction =
      parser_get_opt_param_float(params, "Gravity:rebuild_active_fraction",
                                 gravity_props_default_rebuild_active_fraction);

  if (p->rebuild_frequency < 0.f || p->rebuild_frequency > 1.f)
    error("Invalid tree rebuild frequency. Must be in [0., 1.]");

  /* Tree-PM parameters */
  if (periodic) {
    p->mesh_size = parser_get_param_int(params, "Gravity:mesh_side_length");
    p->distributed_mesh =
        parser_get_opt_param_int(params, "Gravity:distributed_mesh",
                                 gravity_props_default_distributed_mesh);
    p->mesh_uses_local_patches =
        parser_get_opt_param_int(params, "Gravity:mesh_uses_local_patches", 1);
    p->a_smooth = parser_get_opt_param_float(params, "Gravity:a_smooth",
                                             gravity_props_default_a_smooth);
    p->r_cut_max_ratio = parser_get_opt_param_float(
        params, "Gravity:r_cut_max", gravity_props_default_r_cut_max);
    p->r_cut_min_ratio = parser_get_opt_param_float(
        params, "Gravity:r_cut_min", gravity_props_default_r_cut_min);

    p->r_s = p->a_smooth * dim[0] / p->mesh_size;
    p->r_s_inv = 1. / p->r_s;

    /* Some basic checks of what we read */
    if (p->mesh_size % 2 != 0)
      error("The mesh side-length must be an even number.");

    if (p->a_smooth <= 0.)
      error("The mesh smoothing scale 'a_smooth' must be > 0.");

#if !defined(WITH_MPI) || !defined(HAVE_MPI_FFTW)
    if (p->distributed_mesh)
      error(
          "Need to use MPI and FFTW MPI library (i.e. compile with "
          "--enable-mpi-mesh-gravity) to run with distributed mesh.");
#endif

    if (2. * p->a_smooth * p->r_cut_max_ratio > p->mesh_size)
      error("Mesh too small given r_cut_max. Should be at least %d cells wide.",
            (int)(2. * p->a_smooth * p->r_cut_max_ratio) + 1);

    if (p->mesh_size < max3(cdim[0], cdim[1], cdim[2]))
      error(
          "Mesh too small given the number of top-level cells. Should be at "
          "least %d cells wide.",
          max3(cdim[0], cdim[1], cdim[2]));

  } else {
    p->mesh_size = 0;
    p->distributed_mesh = 0;
    p->a_smooth = 0.f;
    p->r_s = FLT_MAX;
    p->r_s_inv = 0.f;
    p->r_cut_min_ratio = 0.f;
    p->r_cut_max_ratio = 0.f;
  }

  /* Time integration */
  p->eta = parser_get_param_float(params, "Gravity:eta");

  /* Read the choice of multipole acceptance criterion */
  char buffer[32] = {0};
  parser_get_param_string(params, "Gravity:MAC", buffer);

  if (strcmp(buffer, "adaptive") == 0) {
    p->use_adaptive_tolerance = 1;
    p->use_gadget_tolerance = 0;
  } else if (strcmp(buffer, "gadget") == 0) {
    p->use_adaptive_tolerance = 1;
    p->use_gadget_tolerance = 1;
  } else if (strcmp(buffer, "geometric") == 0) {
    p->use_adaptive_tolerance = 0;
  } else {
    error(
        "Invalid choice of multipole acceptance criterion: '%s'. Should be "
        "'adaptive', 'gadget', or 'geometric'",
        buffer);
  }

  /* We always start with the geometric MAC */
  p->use_advanced_MAC = 0;

  /* Geometric opening angle */
  p->theta_crit = parser_get_param_double(params, "Gravity:theta_cr");
  if (p->theta_crit >= 1.) error("Theta too large. FMM won't converge.");

  /* Adaptive opening angle tolerance */
  if (p->use_adaptive_tolerance)
    p->adaptive_tolerance =
        parser_get_param_float(params, "Gravity:epsilon_fmm");

  /* Consider truncated forces in the MAC? */
  if (p->use_adaptive_tolerance)
    p->consider_truncation_in_MAC =
        parser_get_opt_param_int(params, "Gravity:allow_truncation_in_MAC", 0);

  /* Are we allowing tree use below softening? */
  p->use_tree_below_softening =
      parser_get_opt_param_int(params, "Gravity:use_tree_below_softening", 0);

#ifdef GADGET2_SOFTENING_CORRECTION
  if (p->use_tree_below_softening)
    error(
        "Cannot solve gravity via the tree below softening with the "
        "Gadget2-type softening kernel");
#endif

  /* Softening parameters */
  if (with_cosmology) {

    if (has_DM) {
      /* Maximal physical softening taken straight from the parameter file */
      p->epsilon_DM_max_physical =
          parser_get_param_double(params, "Gravity:max_physical_DM_softening");

      /* Co-moving softenings taken straight from the parameter file */
      p->epsilon_DM_comoving =
          parser_get_param_double(params, "Gravity:comoving_DM_softening");
    }

    if (has_baryons) {
      /* Maximal physical softening taken straight from the parameter file */
      p->epsilon_baryon_max_physical = parser_get_param_double(
          params, "Gravity:max_physical_baryon_softening");

      /* Co-moving softenings taken straight from the parameter file */
      p->epsilon_baryon_comoving =
          parser_get_param_double(params, "Gravity:comoving_baryon_softening");
    }

    if (has_neutrinos) {
      /* Maximal physical softening taken straight from the parameter file */
      p->epsilon_nu_max_physical =
          parser_get_param_double(params, "Gravity:max_physical_nu_softening");

      /* Co-moving softenings taken straight from the parameter file */
      p->epsilon_nu_comoving =
          parser_get_param_double(params, "Gravity:comoving_nu_softening");
    }

    if (is_zoom_simulation) {

      /* Compute the comoving softening length for background particles as
       * a fraction of the mean inter-particle density of the background DM
       * particles Since they have variable masses the mass factor will be
       * multiplied in later on. Note that we already multiply in the conversion
       * from Plummer -> real softening length */
      const double ratio_background =
          parser_get_param_double(params, "Gravity:softening_ratio_background");

      const double Omega_m = cosmo->Omega_cdm + cosmo->Omega_b;

      const double mean_matter_density = Omega_m * cosmo->critical_density_0;

      p->epsilon_background_fac = kernel_gravity_softening_plummer_equivalent *
                                  ratio_background *
                                  cbrt(1. / mean_matter_density);
    }

  } else {

    if (has_DM) {
      p->epsilon_DM_max_physical =
          parser_get_param_double(params, "Gravity:max_physical_DM_softening");
    }
    if (has_baryons) {
      p->epsilon_baryon_max_physical = parser_get_param_double(
          params, "Gravity:max_physical_baryon_softening");
    }
    if (has_neutrinos) {
      p->epsilon_nu_max_physical =
          parser_get_param_double(params, "Gravity:max_physical_nu_softening");
    }

    /* Some gravity models use the DM softening as the one and only softening
       length that exists. So, if we don't have DM (e.g. hydro test or planetary
       physics), we must have a non-zero epsilon. */
    if (!has_DM && has_baryons)
      p->epsilon_DM_max_physical = p->epsilon_baryon_max_physical;

    p->epsilon_DM_comoving = p->epsilon_DM_max_physical;
    p->epsilon_baryon_comoving = p->epsilon_baryon_max_physical;
  }

  /* Adaptive softening properties */
  p->max_adaptive_softening = parser_get_opt_param_float(
      params, "Gravity:max_adaptive_softening",
      gravity_props_default_max_adaptive_softening /
          kernel_gravity_softening_plummer_equivalent);
  p->min_adaptive_softening =
      parser_get_opt_param_float(params, "Gravity:min_adaptive_softening",
                                 gravity_props_default_min_adaptive_softening);

  p->max_adaptive_softening *= kernel_gravity_softening_plummer_equivalent;
  p->min_adaptive_softening *= kernel_gravity_softening_plummer_equivalent;

  /* Copy over the gravitational constant */
  p->G_Newton = phys_const->const_newton_G;

  /* Set the softening to the current time */
  gravity_props_update(p, cosmo);
}

void gravity_props_update_MAC_choice(struct gravity_props *p) {

  /* Now that we have run initial accelerations,
   * switch to the better MAC */
  if (p->use_adaptive_tolerance) p->use_advanced_MAC = 1;
}

void gravity_props_update(struct gravity_props *p,
                          const struct cosmology *cosmo) {

  /* Current softening length for the high-res. DM particles. */
  double DM_softening, baryon_softening, neutrino_softening;
  if (p->epsilon_DM_comoving * cosmo->a > p->epsilon_DM_max_physical)
    DM_softening = p->epsilon_DM_max_physical / cosmo->a;
  else
    DM_softening = p->epsilon_DM_comoving;

  /* Current softening length for the high-res. baryon particles. */
  if (p->epsilon_baryon_comoving * cosmo->a > p->epsilon_baryon_max_physical)
    baryon_softening = p->epsilon_baryon_max_physical / cosmo->a;
  else
    baryon_softening = p->epsilon_baryon_comoving;

  /* Current softening length for the neutrino DM particles. */
  if (p->epsilon_nu_comoving * cosmo->a > p->epsilon_nu_max_physical)
    neutrino_softening = p->epsilon_nu_max_physical / cosmo->a;
  else
    neutrino_softening = p->epsilon_nu_comoving;

  /* Plummer equivalent -> internal */
  DM_softening *= kernel_gravity_softening_plummer_equivalent;
  baryon_softening *= kernel_gravity_softening_plummer_equivalent;
  neutrino_softening *= kernel_gravity_softening_plummer_equivalent;

  /* Store things */
  p->epsilon_DM_cur = DM_softening;
  p->epsilon_baryon_cur = baryon_softening;
  p->epsilon_nu_cur = neutrino_softening;
}

void gravity_props_print(const struct gravity_props *p) {

  message("Self-gravity scheme: %s", GRAVITY_IMPLEMENTATION);

  message("Self-gravity scheme: FMM-MM with m-poles of order %d",
          SELF_GRAVITY_MULTIPOLE_ORDER);

  message("Self-gravity time integration: eta=%.4f", p->eta);

  if (p->use_adaptive_tolerance) {
    if (p->use_gadget_tolerance) {
      message("Self-gravity opening angle scheme:  Gadget");
      message("Self-gravity opening angle:  epsilon_fmm=%.6f",
              p->adaptive_tolerance);
    } else {
      message("Self-gravity opening angle scheme:  adaptive");
      message("Self-gravity opening angle:  epsilon_fmm=%.6f",
              p->adaptive_tolerance);
    }
  } else {
    message("Self-gravity opening angle scheme:  fixed");
    message("Self-gravity opening angle:  theta_cr=%.4f", p->theta_crit);
  }

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

  message(
      "Self-gravity neutrino DM comoving softening: epsilon=%.6f (Plummer "
      "equivalent: %.6f)",
      p->epsilon_nu_comoving * kernel_gravity_softening_plummer_equivalent,
      p->epsilon_nu_comoving);

  message(
      "Self-gravity neutrino DM maximal physical softening:    epsilon=%.6f "
      "(Plummer equivalent: %.6f)",
      p->epsilon_nu_max_physical * kernel_gravity_softening_plummer_equivalent,
      p->epsilon_nu_max_physical);

  message("Self-gravity mesh side-length: N=%d", p->mesh_size);
  message("Self-gravity mesh smoothing-scale: a_smooth=%f", p->a_smooth);
  message("Self-gravity distributed mesh enabled: %d", p->distributed_mesh);

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

  io_write_attribute_f(
      h_grpgrav, "Comoving neutrino softening length [internal units]",
      p->epsilon_nu_comoving * kernel_gravity_softening_plummer_equivalent);
  io_write_attribute_f(h_grpgrav,
                       "Comoving neutrino softening length (Plummer "
                       "equivalent)  [internal units]",
                       p->epsilon_nu_comoving);

  io_write_attribute_f(
      h_grpgrav, "Maximal physical neutrino softening length  [internal units]",
      p->epsilon_nu_max_physical * kernel_gravity_softening_plummer_equivalent);
  io_write_attribute_f(h_grpgrav,
                       "Maximal physical neutrino softening length (Plummer "
                       "equivalent) [internal units]",
                       p->epsilon_nu_max_physical);

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
