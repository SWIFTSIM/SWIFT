#ifndef SWIFT_DENSITY_DIFFUSION_PART_DATA_H
#define SWIFT_DENSITY_DIFFUSION_PART_DATA_H

#include "../config.h"

#if defined(NO_DENSITY_DIFFUSION)
#include "density_diffusion/no_diffusion/density_diffusion_part_data.h"
#elif defined(MOLTENI09_DIFFUSION)
#include "density_diffusion/Molteni09/density_diffusion_part_data.h"
#else
#error "Invalid choice of density diffusion"
#endif


#endif /*SWIFT_BOUNDARY_H*/

