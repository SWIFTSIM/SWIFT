
/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EXTRA_IO_H
#define SWIFT_EXTRA_IO_H

/* Config parameters. */
#include "../config.h"

/* Import the i/o routines the user asked for */
#if defined(EXTRA_IO_EAGLE)
#include "extra_io/EAGLE/extra_io.h"
#elif defined(EXTRA_IO_NONE)

struct extra_io_properties {};

INLINE static int extra_io_write_particles(const struct part* parts,
                                           const struct xpart* xparts,
                                           struct io_props* list,
                                           const int with_cosmology) {
  return 0;
}

INLINE static int extra_io_write_sparticles(const struct spart* sparts,
                                            struct io_props* list,
                                            const int with_cosmology) {
  return 0;
}

INLINE static int extra_io_write_bparticles(const struct bpart* bparts,
                                            struct io_props* list,
                                            const int with_cosmology) {
  return 0;
}

#ifdef HAVE_HDF5
INLINE static void extra_io_write_flavour(hid_t h_grp, hid_t h_grp_columns) {}
#endif

INLINE static void extra_io_init(struct swift_params* parameter_file,
                                 const struct unit_system* us,
                                 const struct phys_const* phys_const,
                                 const struct cosmology* cosmo,
                                 struct extra_io_properties* props) {}

INLINE static void extra_io_clean(struct extra_io_properties* props) {}

INLINE static void extra_io_struct_dump(const struct extra_io_properties* props,
                                        FILE* stream) {}

INLINE static void extra_io_struct_restore(struct extra_io_properties* props,
                                           FILE* stream) {}

/* In this case there are no extra lightcone map types */
static const struct lightcone_map_type extra_lightcone_map_types[] = {
    {
        /* .name = */ "",
        /* .update_map = */ NULL,
        /* .ptype_contributes = */ NULL,
        /* .baseline_func = */ NULL,
        /* .units = */ UNIT_CONV_NO_UNITS,
        /* .smoothing = */ map_unsmoothed,
        /* .compression = */ compression_write_lossless,
        /* .buffer_scale_factor = */ 1.0,
    },
};

#else
#error "Invalid choice of extra-i/o."
#endif

#endif /* SWIFT_EXTRA_IO_H */
