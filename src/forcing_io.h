#ifndef SWIFT_FORCING_IO_H
#define SWIFT_FORCING_IO_H

/**
 * @file src/forcing_struct.h
 * @brief Branches between the different forcing functions 
 */

/* Config parameters. */
#include <config.h>

/* Import the right external potential definition */
#if defined(FORCING_NONE)
#include "./forcing/none/forcing_io.h"
#elif defined(FORCING_ROBERTS_FLOW)
#include "./forcing/roberts_flow/forcing_io.h"
#elif defined(FORCING_ROBERTS_FLOW_ACCELERATION)
#include "./forcing/roberts_flow_acceleration/forcing_io.h"
#elif defined(FORCING_ABC_FLOW)
#include "./forcing/ABC_flow/forcing_io.h"
#elif defined(FORCING_BALSARAKIM)
#include "./forcing/BalsaraKim/forcing_io.h"
#elif defined(FORCING_BOUNDARY_PARTICLES)
#include "./forcing/boundary_particles/forcing_io.h"
#else
#error "Invalid choice of forcing terms"
#endif

#endif /* SWIFT_FORCING_IO_H */