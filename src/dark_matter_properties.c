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

#define sidm_props_default_max_iterations 30
#define sidm_props_default_h_max FLT_MAX
#define sidm_props_default_h_min_ratio 0.f
#define sidm_props_default_h_sidm FLT_MAX
#define sidm_props_default_h_tolerance 1e-4
#define sidm_props_default_volume_change 1.4f



/**
 * @brief Initialize the global properties of the self-interacting dark matter scheme.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
void sidm_props_init(struct sidm_props* sidm_props,
                     const struct phys_const* phys_const,
                     const struct unit_system* us,
                     struct swift_params* params,
                     const struct cosmology* cosmo) {
    
    /* Scattering cross section in physical units */
    sidm_props->sigma_cgs = parser_get_param_double(params, "SIDM:sigma_cm2_g");
    
    sidm_props->sigma = sidm_props->sigma_cgs / units_cgs_conversion_factor(us, UNIT_CONV_MASS);
    
    sidm_props->sigma *= units_cgs_conversion_factor(us, UNIT_CONV_LENGTH) * units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
    
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
    
    sidm_props->h_search_radius = parser_get_opt_param_float(params, "SIDM:h_sidm",
                                                   sidm_props_default_h_sidm);
    
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
    sidm_props->use_mass_weighted_num_ngb =
    parser_get_opt_param_int(params, "SIDM:use_mass_weighted_num_ngb", 0);
    
    /* ------ Time integration parameters ------------ */
    
    /* Time integration properties */
    sidm_props->CFL_condition = parser_get_param_float(params, "SIDM:CFL_condition");
    
    const float max_volume_change = parser_get_opt_param_float(params, "SPH:max_volume_change", sidm_props_default_volume_change);
    
    sidm_props->log_max_h_change = logf(powf(max_volume_change, hydro_dimension_inv));

    
}

#if defined(HAVE_HDF5)
void sidm_props_print_snapshot(hid_t h_grpsph, const struct sidm_props *p) {
    
    io_write_attribute_f(h_grpsph, "SIDM cross section [cgs units]", p->sigma_cgs);
    io_write_attribute_f(h_grpsph, "SIDM cross section [internal units]", p->sigma);
    io_write_attribute_f(h_grpsph, "SIDM search radius [internal units]", p->h_search_radius);
    
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
