/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
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

#ifndef SWIFT_NONE_NEUTRINO_RESPONSE_H
#define SWIFT_NONE_NEUTRINO_RESPONSE_H

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

#include "cosmology.h"
#include "neutrino_properties.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Structure for handling the linear neutrino response on the mesh
 */
struct neutrino_response {};

INLINE static void neutrino_response_init(
    struct neutrino_response *numesh, struct swift_params *params,
    const struct unit_system *us, const double dim[3],
    const struct cosmology *c, const struct neutrino_props *np,
    const struct gravity_props *gp, int rank, int verbose) {}

INLINE static void neutrino_response_clean(struct neutrino_response *numesh) {}

INLINE static void neutrino_response_compute(
    const struct space *s, struct pm_mesh *mesh, struct threadpool *tp,
    fftw_complex *frho, const int slice_offset, const int slice_width,
    int verbose) {}

INLINE static void neutrino_response_struct_dump(
    const struct neutrino_response *nt, FILE *stream) {}

INLINE static void neutrino_response_struct_restore(
    const struct neutrino_response *nt, FILE *stream) {}

#endif /* SWIFT_NONE_NEUTRINO_RESPONSE_H */
