#ifndef NO_BOUNDARY_H
#define NO_BOUNDARY_H

#include "part.h"

extern struct boundary_parameters boundary_params;

struct boundary_parameters{};

INLINE void init_boundary(struct boundary_parameters *bp, struct swift_params *params){}

INLINE void boundary_fluid_interaction(struct part *restrict pi, struct part * restrict pj, float r, float r2, const float *dx){}

INLINE void boundary_fluid_interaction_nonsym(struct part *restrict pi, struct part * restrict pj, float r, float r2, const float *dx){}

INLINE void boundary_end_force( struct part *p ){} 

#endif
