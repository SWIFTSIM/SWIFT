/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Camila Correa (camila.correa@uva.nl)
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
#include "dark_matter_properties.h"

/* Standard headers */
#include <float.h>
#include <math.h>

/* Local headers. */
#include "common_io.h"
#include "dimension.h"
#include "parser.h"
#include "units.h"
#include "error.h"
#include "restart.h"
#include "kernel_dark_matter.h"
#include "gravity_properties.h"



#define sidm_props_default_max_iterations 30
#define sidm_props_default_h_max FLT_MAX
#define sidm_props_default_h_min_ratio 0.f
#define sidm_props_default_h_sidm FLT_MAX
#define sidm_props_default_h_tolerance 1e-4
#define sidm_props_default_volume_change 1.4f
#define sidm_props_default_sigma 1.f
#define sidm_props_default_mx 0.f
#define sidm_props_default_mphi 0.f

/**
 * @brief Initialize the global properties of the self-interacting dark matter scheme.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param sidm_props The already read-in properties of the SIDM scheme.
 * @param cosmo The cosmological model.
 */
void sidm_props_init(struct sidm_props* sidm_props,
                     const struct phys_const* phys_const,
                     const struct unit_system* us,
                     struct swift_params* params,
                     const struct cosmology* cosmo) {
    
    /* ------ SIDM scattering parameters ---------- */
    
    sidm_props->with_constant_sigma = parser_get_param_int(params, "SIDM:with_constant_cross_section");

    sidm_props->with_velocity_dependent_sigma = parser_get_param_int(params, "SIDM:with_velocity_dependent_cross_section");
    
    sidm_props->mx = parser_get_opt_param_double(params, "SIDM:mx_MeV", sidm_props_default_mx);
    
    sidm_props->mphi = parser_get_opt_param_double(params, "SIDM:mphi_MeV", sidm_props_default_mphi);

    /* Scattering cross section in physical units */
    sidm_props->sigma_cgs = parser_get_opt_param_double(params, "SIDM:sigma_cm2_g", sidm_props_default_sigma);

    /* Scattering cross section in internal units */
    sidm_props->sigma = sidm_props->sigma_cgs * units_cgs_conversion_factor(us, UNIT_CONV_MASS);
    sidm_props->sigma /= units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
    sidm_props->sigma /= units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
    
    sidm_props->with_isotropic_scattering = parser_get_param_int(params, "SIDM:with_isotropic_scattering");

    sidm_props->with_anisotropic_scattering = parser_get_param_int(params, "SIDM:with_anisotropic_scattering");


    /* ------ Smoothing lengths parameters ---------- */
    
    /* Kernel properties */
    sidm_props->eta_neighbours = parser_get_param_float(params, "SIDM:resolution_eta");
    
    /* Tolerance for the smoothing length Newton-Raphson scheme */
    sidm_props->h_tolerance = parser_get_opt_param_float(params, "SIDM:h_tolerance",
                                                sidm_props_default_h_tolerance);
    
    /* Get derived properties */
    sidm_props->target_neighbours = pow_dimension(sidm_props->eta_neighbours);
    
    const float delta_eta = sidm_props->eta_neighbours * (1.f + sidm_props->h_tolerance);

    sidm_props->delta_neighbours = (pow_dimension(delta_eta) - pow_dimension(sidm_props->eta_neighbours));
    
    /* Maximal smoothing length */
    sidm_props->h_max = parser_get_opt_param_float(params, "SIDM:h_max",
                                          sidm_props_default_h_max);
    
    /* Minimal smoothing length ratio to softening */
    sidm_props->h_min_ratio = parser_get_opt_param_float(params, "SIDM:h_min_ratio",
                                                sidm_props_default_h_min_ratio);
    
    /* Temporarily set the minimal softening to 0. */
    sidm_props->h_min = 0.f;
    
    /* Number of iterations to converge h */
    sidm_props->max_smoothing_iterations = parser_get_opt_param_int(
                                                           params, "SIDM:max_ghost_iterations", sidm_props_default_max_iterations);
    
    if (sidm_props->max_smoothing_iterations <= 10)
        error("The number of smoothing length iterations for DM density should be > 10");
    
    /* ------ Neighbour number definition ------------ */
    
    /* Non-conventional neighbour number definition */
    sidm_props->use_mass_weighted_num_ngb = parser_get_opt_param_int(params, "SIDM:use_mass_weighted_num_ngb", 0);
    
    /* ------ Time integration parameters ------------ */
    
    const float max_volume_change = parser_get_opt_param_float(params, "SPH:max_volume_change", sidm_props_default_volume_change);
    
    sidm_props->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));
}

/**
 * @brief Print the global properties of the SIDM scheme.
 *
 * @param sidm_props The #sidm_props.
 */
void sidm_props_print(struct sidm_props *sidm_props) {
    
    /* Now describe SIDM model */
    if (sidm_props->with_constant_sigma) message("Running SIDM scheme with constant cross section");
    if (sidm_props->with_velocity_dependent_sigma)  message("Running SIDM scheme with velocity-dependent cross section");
    if (sidm_props->with_isotropic_scattering) message("and isotropic scattering.");
    if (sidm_props->with_anisotropic_scattering) message("and anisotropic scattering.");

    message("SIDM kernel: %s with eta=%f (%.2f neighbours).", dm_kernel_name,
            sidm_props->eta_neighbours, sidm_props->target_neighbours);
    
    if (sidm_props->use_mass_weighted_num_ngb)
    message("Neighbour number definition: Mass-weighted.");
    else
    message("Neighbour number definition: Unweighted.");
    
    if (sidm_props->h_max != sidm_props_default_h_max)
    message("Maximal smoothing length allowed: %.4f", sidm_props->h_max);
}

/**
 * @brief Update the global properties of the hydro scheme for that time-step.
 *
 * @param p The properties to update.
 * @param gp The properties of the gravity scheme.
 * @param cosmo The cosmological model.
 */
void sidm_props_update(struct sidm_props *sidm_props, const struct gravity_props *gp,
                        const struct cosmology *cosmo) {
    
    /* Update the minimal allowed smoothing length
     *
     * We follow Gadget here and demand that the kernel support (h * gamma)
     * is a fixed fraction of the radius at which the softened forces
     * recover a Newtonian behaviour (i.e. 2.8 * Plummer equivalent softening
     * in the case of a cubic spline kernel). */
    sidm_props->h_min = sidm_props->h_min_ratio * gp->epsilon_DM_cur / dm_kernel_gamma;
}


#if defined(HAVE_HDF5)
void sidm_props_print_snapshot(hid_t h_grpsph, const struct sidm_props *p) {
    
    io_write_attribute_f(h_grpsph, "SIDM cross section [cgs units]", p->sigma_cgs);
    io_write_attribute_f(h_grpsph, "SIDM cross section [internal units]", p->sigma);
    
}
#endif


/**
 * @brief Write a sidm_props struct to the given FILE as a stream of bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void sidm_props_struct_dump(const struct sidm_props *p, FILE *stream) {
    restart_write_blocks((void *)p, sizeof(struct sidm_props), 1, stream,
                         "sidmprops", "sidm props");
}

/**
 * @brief Restore a sidm_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void sidm_props_struct_restore(const struct sidm_props *p, FILE *stream) {
    restart_read_blocks((void *)p, sizeof(struct sidm_props), 1, stream, NULL,
                        "sidm props");
}
