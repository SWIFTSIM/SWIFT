#ifndef SWIFT_DENSITY_DIFFUSION_H
#define SWIFT_DENSITY_DIFFUSION_H

#include "../config.h"

#if defined(NO_DENSITY_DIFFUSION)
#include "density_diffusion/no_diffusion/density_diffusion.h"
#elif defined(MOLTENI09_DIFFUSION)
#include "density_diffusion/Molteni09/density_diffusion.h"
#else
#error "Invalid choice of boundary movement"
#endif


#endif /*SWIFT_BOUNDARY_H*/
