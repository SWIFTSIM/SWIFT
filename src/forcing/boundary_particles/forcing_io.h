#ifndef SWIFT_FORCING_BOUNDARY_PARTICLES_IO_H
#define SWIFT_FORCING_BOUNDARY_PARTICLES_IO_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "forcing.h"
#include "engine.h"
#include "io_properties.h"

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int forcing_write_particles(
    const struct part *parts, const struct xpart *xparts,
    struct io_props *list) {
  return 0;
}

#endif /* SWIFT_FORCING_BOUNDARY_PARTICLES_IO_H */