#ifndef NO_SHIFTING_H
#define NO_SHIFTING_H

#include "part.h"
#include "particle_shifting_part_data.h"

extern struct particle_shifting_parameters psp;

struct particle_shifting_parameters {
};


INLINE void reset_shifting_data(struct part *p){
}

//Must be called after engine_init
INLINE void init_shifting_parameters(struct particle_shifting_parameters *shift,
                            struct swift_params *params){

}

INLINE void compute_particle_shifting(struct part *p, float dt){
}

INLINE void compute_shifting_term(struct part *pi, struct part *pj,  float r2,  float r,  float wi_dx,  float wj_dx, const float *dx){


}


INLINE void compute_shifting_term_nonsym( struct part *pi, const struct part *pj,  float r2,  float r,  float wi_dx, const float *dx){

}
#endif
