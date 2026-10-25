/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#ifndef SWIFT_MAGMA_HYDRO_PARAMETERS_H
#define SWIFT_MAGMA_HYDRO_PARAMETERS_H

/* Configuration file */
#include <config.h>

/* Global headers */
#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif
#include <math.h>

/* Local headers */
#include "common_io.h"
#include "error.h"
#include "inline.h"
#include "kernel_hydro.h"
#include "parser.h"

// #define TRADITIONAL_SPH_ACCELERATION_TERM /* MI1 (set by configure) */
// #define MAGMA_MI3_FORMULATION /* MI3 (set by configure): MI2's 1/(rho_a
// rho_b)
//                                  weights with P_a on G_a and P_b on G_b,
//                                  so that a vanishing W_ab(h_b) does not halve
//                                  the operator of particle a. */
// #define USE_ZEROTH_ORDER_VELOCITIES
// #define USE_STANDARD_KERNEL_GRADIENTS
#define GRAVITY_DIFF_VELOCITY

/**
 * @file MAGMA/hydro_parameters.h
 * @brief MAGMA-2 implementation of SPH (Rosswog. 2020)
 *
 *        This file defines a number of things that are used in
 *        hydro_properties.c as defaults for run-time parameters
 *        as well as a number of compile-time parameters.
 *
 *        All the constants of the scheme can be changed at run time through
 *        the (optional) parameters listed in viscosity_init() and
 *        diffusion_init(). The interaction functions have no access to the
 *        #hydro_props, so the values are also kept in the global variables
 *        magma_viscosity and magma_diffusion (defined in hydro.c), following
 *        the pattern of the equation-of-state parameters.
 */

/* Default values of the run-time parameters */

/*! Linear viscosity coefficient alpha (Rosswog 2020, eq. 14) */
#define hydro_props_default_viscosity_alpha 1.0f

/*! Quadratic viscosity coefficient beta (Rosswog 2020, eq. 14) */
#define hydro_props_default_viscosity_beta 2.0f

/*! Softening of the viscosity velocity jump, in units of h (eq. 15). Also
 * softens the approach velocity of the Courant condition (eq. 36). */
#define hydro_props_default_viscosity_epsilon 0.1f

/*! Distance (in units of h) below which the slope limiter switches the
 * reconstruction off (eq. 21-23). A negative value means the mean
 * inter-particle separation, 1 / eta_neighbours, as in the paper (whose
 * eta_crit = (32 pi / 3 N_ngb)^(1/3) is that separation in units of half the
 * kernel support). */
#define hydro_props_default_limiter_eta_crit -1.f

/*! Width (in units of h) of the exponential cut-off of the limiter (eq. 21).
 * A negative value means the paper's 0.2 in units of half the kernel support,
 * i.e. 0.1 * kernel_gamma. */
#define hydro_props_default_limiter_width -1.f

/*! Shear limiter of the viscous pressure: each particle's Q is scaled by its
 * Cullen & Dehnen (2010) limiter, the convergence fraction with the
 * shock-presence boost (their eq. 17-18). 0 off, 1 both terms of Q, 2 the
 * quadratic (beta mu^2) term only. */
#define hydro_props_default_viscosity_shear_limiter 1

/*! Floor of the shear limiter, xi = max(xi, floor): the fraction of the
 * viscosity kept in pure shear. */
#define hydro_props_default_viscosity_limiter_floor 0.f

/*! Mach number of the pair's velocity jump, |mu| / c, above which the shear
 * limiter also acts on the linear (alpha c mu) term of Q. Below it the
 * linear term keeps its full strength, so subsonic shear layers (which need
 * the viscosity to keep their particle order) are untouched and only cold,
 * supersonic shear (differential rotation) is relieved of it. The 2026-10
 * optimisation study of the scheme-comparison suite (one-at-a-time sweep,
 * factorial, replicated finalists on 27 tests, MI1/MI2/MI3) settled on 0.5:
 * the Keplerian rings and discs survive 6 and 20 orbits with nothing lost
 * on shocks, shear layers or interfaces once the pressure-residual
 * conduction is on (with the paper's conduction alone, 0.5 cost 8 per cent
 * on the Sedov pressure error, which is why 0.7 was the earlier choice). */
#define hydro_props_default_viscosity_limiter_mach 0.5f

/*! Time-dependent viscosity coefficient (Cullen & Dehnen 2010): 0 the
 * paper's constant alpha, 1 a per-particle alpha raised to
 * alpha_max A / (A + c_s^2) with A = xi H^2 max(0, -d(div v)/dt) where the
 * flow is converging (xi the shock-boosted limiter above, H the kernel
 * support) and decaying back otherwise, floored at
 * alpha_min; the decay time is Cullen & Dehnen's h / (2 length v_sig) with
 * their h the kernel support radius H: H / (2 length c_s). The slope-limited
 * reconstruction does not switch the
 * viscosity off in smooth flow to better than a few per cent of the
 * velocity jump (gradient noise and the first-order mid-point
 * reconstruction), which with a constant alpha = 1 damps a linear sound
 * wave by 4-5 per cent per period at 32 particles per wavelength. */
#define hydro_props_default_viscosity_switch 0

/*! Floor of the time-dependent viscosity coefficient. */
#define hydro_props_default_viscosity_alpha_min 0.1f

/*! Inverse decay time of the time-dependent viscosity coefficient in units
 * of the sound-crossing time of the kernel support (Cullen & Dehnen's l). */
#define hydro_props_default_viscosity_length 0.05f

/*! Maximal condition number of the inverse C-matrix. Directions in which the
 * neighbourhood is more degenerate than this (planar or filamentary
 * configurations) are regularised by clipping the spectrum of the inverse. */
#define hydro_props_default_gradient_max_condition_number 60.f

/*! Maximal angle (in radians) between the gradient function G and the
 * separation vector of a pair before falling back to kernel gradients for the
 * pair */
#define hydro_props_default_gradient_angle_limit 0.5f

/*! Thermal conduction coefficient alpha_u (Rosswog 2020, eq. 24) */
#define hydro_props_default_diffusion_alpha 0.05f

/*! Coefficient of the reconstructed-pressure-residual conduction term,
 * added to the paper's: the conduction speed is alpha_u |v~_ab| + alpha_P
 * sqrt(2 |dP_rec| / (rho_a + rho_b)); 0 switches it off. The residual
 * of the slope-limited mid-point reconstruction of P is ~0 across a
 * resolved pressure gradient (hydrostatic equilibrium under gravity) and the
 * full jump across the kernel-scale pressure blip of the density estimate
 * at an entropy jump, or a shock; the conduction then removes the blip and
 * switches itself off. With the threshold below, 0.5 takes the
 * static-contact pressure blip from 20 to 5 per cent, halves the Sod
 * contact blip, cures the Sedov central overheating and keeps the KH rolls
 * coherent, for MI2 and MI3 alike, at no cost on the sound wave, Gresho
 * vortex or hydrostatic halo; 1 (default, from the 2026-10 optimisation
 * study) further improves the Evrard entropy and pressure profiles, the
 * contact blips, the blob stripping and the square tests, at the price of
 * more entropy smearing at strong contacts (weak shock tube, 1D blast
 * waves). */
#define hydro_props_default_diffusion_alpha_pressure 1.0f

/*! Threshold of the pressure-residual conduction term: only the part of the
 * residual above this fraction of the pair's mean pressure counts, so that
 * the term switches off once a blip is reduced to the pressure noise of the
 * particle distribution (a few per cent) instead of diffusing the entropy
 * jump for ever at the rate the noise sets. */
#define hydro_props_default_diffusion_pressure_threshold 0.05f

/*! Growth rate of the per-particle conduction coefficient (SPHENIX's
 * beta_D): the coefficient grows at beta H laplace(u) / sqrt(u) where the
 * internal energy has a local minimum */
#define hydro_props_default_diffusion_beta 1.0f

/*! Floor of the per-particle conduction coefficient. The coefficient decays
 * towards this value over the signal crossing time of the kernel support and
 * is capped where the neighbours' viscosity coefficient is up (shocks, and
 * the surroundings of a particle a feedback event has just heated). With the
 * floor at 1 (default, 2026-10 feedback study) the conduction is the paper's
 * full term wherever the gas is not viscous, which keeps the suite's contact
 * results, and the cap alone protects feedback energy: the cooling Sedov
 * blast, the Balsara-Kim patch and the isolated galaxy then follow SPHENIX.
 * With the floor at 0 the coefficient is SPHENIX's: only the Laplacian
 * source switches it on, which costs 30 to 70 per cent larger contact blips
 * on the suite for no gain on the feedback tests. */
#define hydro_props_default_diffusion_alpha_min 1.0f

/* Structs that store the relevant variables */

/*! Artificial viscosity (and gradient-function) parameters */
struct viscosity_global_data {

  /*! Linear viscosity coefficient (eq. 14) */
  float alpha;

  /*! Quadratic viscosity coefficient (eq. 14) */
  float beta;

  /*! Softening of the velocity jump, in units of h (eq. 15) */
  float epsilon;

  /*! Limiter cut-off distance, in units of h (eq. 21-23) */
  float eta_crit;

  /*! Limiter cut-off width, in units of h (eq. 21) */
  float limiter_width;

  /*! Maximal condition number of the C-matrix */
  float max_condition_number;

  /*! Maximal angle between G and the pair separation (radians) */
  float angle_limit;

  /*! cos(angle_limit), pre-computed */
  float cos_angle_limit;

  /*! Binary switch (0 or 1, default 1): apply the distance cut-off of the
   * slope limiter (eq. 21, eta_crit) in proportion to the convergence
   * fraction of the flow.
   * The cut-off guards against particles converging past each other below
   * the resolution scale in a smooth compression; distance alone is a poor
   * proxy for that in shear flows, where close pairs approach along their
   * separation without ever colliding and the viscosity acting on their
   * un-reconstructed velocity difference dissipates the shear (differential
   * rotation). With the switch on, the cut-off keeps its full strength in
   * compressions and is inactive in pure shear. */
  int limiter_convergence_switch;

  /*! Scale each particle's viscous pressure by its own Cullen & Dehnen
   * (2010) limiter (convergence fraction with the shock-presence boost):
   * 0 off, 1 both terms of Q, 2 the quadratic (beta mu^2) term only
   * (default 1) */
  int shear_limiter;

  /*! Floor of the viscosity limiter, xi = max(xi, floor): keeps a fraction
   * of the viscosity in pure shear (default 0: none) */
  float limiter_floor;

  /*! The limiter acts on the linear (alpha c mu) term of Q only for pairs
   * with |mu| > limiter_mach * c, i.e. only in supersonic shear (default
   * 0.7; 0: always) */
  float limiter_mach;

  /*! Binary switch (0 or 1, default 0): time-dependent alpha (Cullen &
   * Dehnen 2010): 0 the constant alpha, 1 the per-particle alpha of
   * part::alpha_visc on the linear term of Q */
  int alpha_switch;

  /*! Floor of the time-dependent alpha (default 0.1) */
  float alpha_min;

  /*! Decay rate of the time-dependent alpha, Cullen & Dehnen's l: the decay
   * time is H / (2 l c_s) (default 0.05, i.e. 10 sound crossings of H) */
  float length;
};

/*! Thermal diffusion parameters */
struct diffusion_global_data {

  /*! Conduction coefficient (eq. 24) */
  float alpha;

  /*! Coefficient of the reconstructed-pressure-residual conduction term
   * (default 1; 0 switches it off) */
  float alpha_pressure;

  /*! Fraction of the mean pressure below which the residual is ignored
   * (default 0.05) */
  float pressure_threshold;

  /*! Growth rate of the per-particle conduction coefficient (default 1) */
  float beta;

  /*! Floor of the per-particle conduction coefficient (default 1) */
  float alpha_min;
};

/*! Copies of the parameters accessible to the interaction functions */
extern struct viscosity_global_data magma_viscosity;
extern struct diffusion_global_data magma_diffusion;

/* Functions for reading from parameter file */

/* Forward declarations */
struct swift_params;
struct phys_const;
struct unit_system;

/* Viscosity */

/**
 * @brief Sets the viscosity parameters to their default values.
 *
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 * @param eta_neighbours: resolution parameter (h in units of the mean
 * inter-particle separation) used to resolve the automatic limiter settings.
 **/
static INLINE void viscosity_set_defaults(
    struct viscosity_global_data *viscosity, const float eta_neighbours) {

  viscosity->alpha = hydro_props_default_viscosity_alpha;
  viscosity->beta = hydro_props_default_viscosity_beta;
  viscosity->epsilon = hydro_props_default_viscosity_epsilon;
  viscosity->eta_crit = hydro_props_default_limiter_eta_crit;
  viscosity->limiter_width = hydro_props_default_limiter_width;
  viscosity->max_condition_number =
      hydro_props_default_gradient_max_condition_number;
  viscosity->angle_limit = hydro_props_default_gradient_angle_limit;
  viscosity->limiter_convergence_switch = 1;
  viscosity->shear_limiter = hydro_props_default_viscosity_shear_limiter;
  viscosity->limiter_floor = hydro_props_default_viscosity_limiter_floor;
  viscosity->limiter_mach = hydro_props_default_viscosity_limiter_mach;
  viscosity->alpha_switch = hydro_props_default_viscosity_switch;
  viscosity->alpha_min = hydro_props_default_viscosity_alpha_min;
  viscosity->length = hydro_props_default_viscosity_length;

  if (viscosity->eta_crit < 0.f) viscosity->eta_crit = 1.f / eta_neighbours;
  if (viscosity->limiter_width < 0.f)
    viscosity->limiter_width = 0.1f * kernel_gamma;
  viscosity->cos_angle_limit = cosf(viscosity->angle_limit);
}

/**
 * @brief Initialises the viscosity parameters in the struct from
 *        the parameter file, or sets them to defaults.
 *
 * @param params: the pointer to the swift_params file
 * @param us: pointer to the internal unit system
 * @param phys_const: pointer to the physical constants system
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 **/
static INLINE void viscosity_init(struct swift_params *params,
                                  const struct unit_system *us,
                                  const struct phys_const *phys_const,
                                  struct viscosity_global_data *viscosity) {

  viscosity->alpha = parser_get_opt_param_float(
      params, "SPH:viscosity_alpha", hydro_props_default_viscosity_alpha);
  viscosity->beta = parser_get_opt_param_float(
      params, "SPH:viscosity_beta", hydro_props_default_viscosity_beta);
  viscosity->epsilon = parser_get_opt_param_float(
      params, "SPH:viscosity_epsilon", hydro_props_default_viscosity_epsilon);
  viscosity->eta_crit = parser_get_opt_param_float(
      params, "SPH:limiter_eta_crit", hydro_props_default_limiter_eta_crit);
  viscosity->limiter_width = parser_get_opt_param_float(
      params, "SPH:limiter_width", hydro_props_default_limiter_width);
  viscosity->max_condition_number = parser_get_opt_param_float(
      params, "SPH:gradient_max_condition_number",
      hydro_props_default_gradient_max_condition_number);
  viscosity->angle_limit =
      parser_get_opt_param_float(params, "SPH:gradient_angle_limit",
                                 hydro_props_default_gradient_angle_limit);
  viscosity->limiter_convergence_switch =
      parser_get_opt_param_int(params, "SPH:limiter_convergence_switch", 1);
  viscosity->shear_limiter =
      parser_get_opt_param_int(params, "SPH:viscosity_shear_limiter",
                               hydro_props_default_viscosity_shear_limiter);
  viscosity->limiter_floor =
      parser_get_opt_param_float(params, "SPH:viscosity_limiter_floor",
                                 hydro_props_default_viscosity_limiter_floor);
  viscosity->limiter_mach =
      parser_get_opt_param_float(params, "SPH:viscosity_limiter_mach",
                                 hydro_props_default_viscosity_limiter_mach);
  viscosity->alpha_switch = parser_get_opt_param_int(
      params, "SPH:viscosity_switch", hydro_props_default_viscosity_switch);
  viscosity->alpha_min =
      parser_get_opt_param_float(params, "SPH:viscosity_alpha_min",
                                 hydro_props_default_viscosity_alpha_min);
  viscosity->length = parser_get_opt_param_float(
      params, "SPH:viscosity_length", hydro_props_default_viscosity_length);

  /* Automatic limiter settings: the paper's prescription in SWIFT's h
   * convention (kernel support = kernel_gamma * h, the paper uses 2 h) */
  if (viscosity->eta_crit < 0.f) {
    const float eta_neighbours =
        parser_get_param_float(params, "SPH:resolution_eta");
    viscosity->eta_crit = 1.f / eta_neighbours;
  }
  if (viscosity->limiter_width < 0.f)
    viscosity->limiter_width = 0.1f * kernel_gamma;

  if (viscosity->alpha < 0.f || viscosity->beta < 0.f ||
      viscosity->epsilon <= 0.f || viscosity->limiter_width <= 0.f ||
      viscosity->max_condition_number < 1.f || viscosity->angle_limit < 0.f ||
      viscosity->angle_limit > M_PI_2 || viscosity->shear_limiter < 0 ||
      viscosity->shear_limiter > 2 || viscosity->limiter_floor < 0.f ||
      viscosity->limiter_floor > 1.f || viscosity->limiter_mach < 0.f ||
      viscosity->alpha_switch < 0 || viscosity->alpha_switch > 1 ||
      viscosity->limiter_convergence_switch < 0 ||
      viscosity->limiter_convergence_switch > 1 || viscosity->alpha_min < 0.f ||
      viscosity->alpha_min > viscosity->alpha || viscosity->length <= 0.f)
    error("Invalid MAGMA viscosity / gradient parameters.");

  viscosity->cos_angle_limit = cosf(viscosity->angle_limit);

  /* Make the values available to the interaction functions */
  magma_viscosity = *viscosity;
}

/**
 * @brief Initialises a viscosity struct to sensible numbers for mocking
 *        purposes.
 *
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 **/
static INLINE void viscosity_init_no_hydro(
    struct viscosity_global_data *viscosity) {

  viscosity_set_defaults(viscosity, /*eta_neighbours=*/1.2348f);
  magma_viscosity = *viscosity;
}

/**
 * @brief Prints out the viscosity parameters at the start of a run.
 *
 * @param viscosity: pointer to the viscosity_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void viscosity_print(
    const struct viscosity_global_data *viscosity) {
  message(
      "Artificial viscosity parameters set to alpha: %.3f, beta: %.3f, "
      "epsilon: %.3f.",
      viscosity->alpha, viscosity->beta, viscosity->epsilon);
  message("Slope limiter parameters set to eta_crit: %.4f h, width: %.4f h.",
          viscosity->eta_crit, viscosity->limiter_width);
  message(
      "Gradient-function regularisation: C-matrix inverse clipped to a "
      "condition number of %.1f (particle); fallback: angle > %.3f rad "
      "(pair).",
      viscosity->max_condition_number, viscosity->angle_limit);
  message(
      "Slope-limiter distance cut-off scaled by the convergence fraction "
      "(binary switch): %s.",
      viscosity->limiter_convergence_switch ? "on" : "off");
  if (viscosity->shear_limiter)
    message(
        "Viscous pressure scaled by the Cullen & Dehnen limiter (%s, floor "
        "%.2f, linear term limited for |mu| > %.2f c).",
        viscosity->shear_limiter == 2 ? "quadratic term only" : "both terms",
        viscosity->limiter_floor, viscosity->limiter_mach);
  else
    message("Viscous pressure not scaled by the Cullen & Dehnen limiter.");
  if (viscosity->alpha_switch)
    message(
        "Time-dependent alpha (Cullen & Dehnen 2010): alpha_max %.2f, "
        "alpha_min %.2f, decay length %.3f.",
        viscosity->alpha, viscosity->alpha_min, viscosity->length);
  else
    message("Constant viscosity alpha %.2f.", viscosity->alpha);
}

#if defined(HAVE_HDF5)
/**
 * @brief Prints the viscosity information to the snapshot when writing.
 *
 * @param h_grpsph: the SPH group in the ICs to write attributes to.
 * @param viscosity: pointer to the viscosity_global_data struct.
 **/
static INLINE void viscosity_print_snapshot(
    hid_t h_grpsph, const struct viscosity_global_data *viscosity) {

  io_write_attribute_f(h_grpsph, "Alpha viscosity", viscosity->alpha);
  io_write_attribute_f(h_grpsph, "Beta viscosity", viscosity->beta);
  io_write_attribute_f(h_grpsph, "Epsilon viscosity", viscosity->epsilon);
  io_write_attribute_f(h_grpsph, "Limiter eta_crit", viscosity->eta_crit);
  io_write_attribute_f(h_grpsph, "Limiter width", viscosity->limiter_width);
  io_write_attribute_i(h_grpsph, "Viscosity switch", viscosity->alpha_switch);
  io_write_attribute_f(h_grpsph, "Alpha viscosity (min)", viscosity->alpha_min);
  io_write_attribute_f(h_grpsph, "Viscosity decay length", viscosity->length);
  io_write_attribute_i(h_grpsph, "Limiter convergence switch",
                       viscosity->limiter_convergence_switch);
  io_write_attribute_f(h_grpsph, "Gradient max condition number",
                       viscosity->max_condition_number);
  io_write_attribute_f(h_grpsph, "Gradient angle limit",
                       viscosity->angle_limit);
}
#endif

/* Diffusion */

/**
 * @brief Initialises the diffusion parameters in the struct from
 *        the parameter file, or sets them to defaults.
 *
 * @param params: the pointer to the swift_params file
 * @param us: pointer to the internal unit system
 * @param phys_const: pointer to the physical constants system
 * @param diffusion: pointer to the diffusion struct to be filled.
 **/
static INLINE void diffusion_init(struct swift_params *params,
                                  const struct unit_system *us,
                                  const struct phys_const *phys_const,
                                  struct diffusion_global_data *diffusion) {

  diffusion->alpha = parser_get_opt_param_float(
      params, "SPH:diffusion_alpha", hydro_props_default_diffusion_alpha);
  diffusion->alpha_pressure =
      parser_get_opt_param_float(params, "SPH:diffusion_alpha_pressure",
                                 hydro_props_default_diffusion_alpha_pressure);

  diffusion->pressure_threshold = parser_get_opt_param_float(
      params, "SPH:diffusion_pressure_threshold",
      hydro_props_default_diffusion_pressure_threshold);
  diffusion->beta = parser_get_opt_param_float(
      params, "SPH:diffusion_beta", hydro_props_default_diffusion_beta);
  diffusion->alpha_min =
      parser_get_opt_param_float(params, "SPH:diffusion_alpha_min",
                                 hydro_props_default_diffusion_alpha_min);

  if (diffusion->alpha < 0.f || diffusion->alpha_pressure < 0.f ||
      diffusion->pressure_threshold < 0.f || diffusion->beta < 0.f ||
      diffusion->alpha_min < 0.f || diffusion->alpha_min > 1.f)
    error("Invalid MAGMA diffusion parameters.");

  /* Make the values available to the interaction functions */
  magma_diffusion = *diffusion;
}

/**
 * @brief Initialises a diffusion struct to sensible numbers for mocking
 *        purposes.
 *
 * @param diffusion: pointer to the diffusion_global_data struct to be filled.
 **/
static INLINE void diffusion_init_no_hydro(
    struct diffusion_global_data *diffusion) {

  diffusion->alpha = hydro_props_default_diffusion_alpha;
  diffusion->alpha_pressure = hydro_props_default_diffusion_alpha_pressure;
  diffusion->pressure_threshold =
      hydro_props_default_diffusion_pressure_threshold;
  diffusion->beta = hydro_props_default_diffusion_beta;
  diffusion->alpha_min = hydro_props_default_diffusion_alpha_min;
  magma_diffusion = *diffusion;
}

/**
 * @brief Prints out the diffusion parameters at the start of a run.
 *
 * @param diffusion: pointer to the diffusion_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void diffusion_print(
    const struct diffusion_global_data *diffusion) {
  message(
      "Artificial conduction parameters set to alpha_u: %.3f, "
      "pressure-residual term: %.2f (threshold %.2f of the pressure), "
      "per-particle coefficient with growth rate %.2f and floor %.2f.",
      diffusion->alpha, diffusion->alpha_pressure,
      diffusion->pressure_threshold, diffusion->beta, diffusion->alpha_min);
}

#ifdef HAVE_HDF5
/**
 * @brief Prints the diffusion information to the snapshot when writing.
 *
 * @param h_grpsph: the SPH group in the ICs to write attributes to.
 * @param diffusion: pointer to the diffusion_global_data struct.
 **/
static INLINE void diffusion_print_snapshot(
    hid_t h_grpsph, const struct diffusion_global_data *diffusion) {
  io_write_attribute_f(h_grpsph, "Alpha diffusion", diffusion->alpha);
  io_write_attribute_f(h_grpsph, "Alpha diffusion (pressure residual)",
                       diffusion->alpha_pressure);
  io_write_attribute_f(h_grpsph, "Beta diffusion", diffusion->beta);
  io_write_attribute_f(h_grpsph, "Alpha diffusion (min)", diffusion->alpha_min);
  io_write_attribute_f(h_grpsph, "Pressure residual threshold",
                       diffusion->pressure_threshold);
}
#endif

#endif /* SWIFT_MAGMA_HYDRO_PARAMETERS_H */
