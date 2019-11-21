#ifndef SWIFT_BOUNDARY_H
#define SWIFT_BOUNDARY_H

#include "../config.h"

#if defined(NO_BOUNDARY)
#include "boundary/no_boundary/boundary.h"
#elif defined(FIXED_FLUID_BOUNDARY)
#include "boundary/fixed_fluid_boundary/boundary.h"
#elif defined(MONAGHAN94_BOUNDARY)
#include "boundary/Monaghan94/boundary.h"
#else
#error "Invalid choice of boundary particle"
#endif


#endif /*SWIFT_BOUNDARY_H*/
