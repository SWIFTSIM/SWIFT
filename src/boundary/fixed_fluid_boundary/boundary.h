#ifndef FIXED_FLUID_BOUNDARY_H
#define FIXED_FLUID_BOUNDARY_H

#include "part.h"

extern struct boundary_parameters boundary_params;

struct boundary_parameters{};

INLINE void init_boundary(struct boundary_parameters *bp, struct swift_params *params){}

INLINE void boundary_fluid_interaction(struct part *restrict pi, struct part * restrict pj, float r, float r2, const float *dx){}

INLINE void boundary_fluid_interaction_nonsym(struct part *restrict pi, const struct part *restrict pj, float r, float r2, const float *dx){}

INLINE void boundary_end_force( struct part *p ){
  if(p->is_boundary){
    p->a_hydro[0] = 0.0;
    p->a_hydro[1] = 0.0;
    p->a_hydro[2] = 0.0;
  }
} 

#endif
