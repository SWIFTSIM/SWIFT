#ifndef SWIFT_BOUNDARY_MOVEMENT_H
#define SWIFT_BOUNDARY_MOVEMENT_H

#include "../config.h"

#if defined(NO_BOUNDARY_MOVEMENT)
#include "boundary_motion/no_movement/boundary_movement.h"
#elif defined(SINUSOIDAL_BOUNDARY_MOVEMENT)
#include "boundary_motion/sinusoidal/boundary_movement.h"
#else
#error "Invalid choice of boundary movement"
#endif


#endif /*SWIFT_BOUNDARY_H*/
