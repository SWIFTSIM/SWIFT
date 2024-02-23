/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2024 Camila Correa (camila.correa@cea.fr)
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
#define sidm_props_default_sigma 0
#define sidm_props_default_mx 0.f
#define sidm_props_default_mphi 0.f
#define sidm_props_default_alphax 0.01

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
                     const struct cosmology* cosmo) {}

/**
 * @brief Print the global properties of the SIDM scheme.
 *
 * @param sidm_props The #sidm_props.
 */
void sidm_props_print(struct sidm_props *sidm_props) {}

/**
 * @brief Update the global properties of the hydro scheme for that time-step.
 *
 * @param p The properties to update.
 * @param gp The properties of the gravity scheme.
 * @param cosmo The cosmological model.
 */
void sidm_props_update(struct sidm_props *sidm_props, const struct gravity_props *gp,
                        const struct cosmology *cosmo) {}


#if defined(HAVE_HDF5)
void sidm_props_print_snapshot(hid_t h_grpsph, const struct sidm_props *p) {}
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
