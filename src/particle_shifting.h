#ifndef SWIFT_PARTICLE_SHIFTING_H
#define SWIFT_PARTICLE_SHIFTING_H

#include "../config.h"

#if defined(NO_PARTICLE_SHIFTING)
#include "particle_shifting/no_particle_shifting/particle_shifting.h"
#elif defined(MOKOS13_SHIFTING)
#include "particle_shifting/Mokos13/particle_shifting.h"
#else
#error "Invalid choice of particle shifting"
#endif


#endif /*SWIFT_BOUNDARY_H*/
