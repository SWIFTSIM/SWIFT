#ifndef SWIFT_FORCING_STRUCT_H
#define SWIFT_FORCING_STRUCT_H

/**
 * @file src/forcing_struct.h
 * @brief Branches between the different forcing functions 
 */

/* Config parameters. */
#include <config.h>

/* Import the right external potential definition */
#if defined(FORCING_NONE)
#include "./forcing/none/forcing_struct.h"
#elif defined(FORCING_ROBERTS_FLOW)
#include "./forcing/none/forcing_struct.h"
#elif defined(FORCING_ROBERTS_FLOW_ACCELERATION)
#include "./forcing/none/forcing_struct.h"
#elif defined(FORCING_ABC_FLOW)
#include "./forcing/none/forcing_struct.h"
#elif defined(FORCING_BALSARAKIM)
#include "./forcing/BalsaraKim/forcing_struct.h"
#else
#error "Invalid choice of forcing terms"
#endif

#endif /* SWIFT_FORCING_STRUCT_H */