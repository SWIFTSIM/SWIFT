#ifndef NO_DIFFUSIVE_H
#define NO_DIFFUSIVE_H

#include "part.h"

extern struct diffusive_term_parameters dtp;

struct diffusive_term_parameters {
};

//Must be called after engine_init
INLINE void init_diffusive_term(struct diffusive_term_parameters *diff,
                            struct swift_params *params){

}

INLINE void reset_diffusive_term(struct part *restrict p){
}

INLINE void diffusive_end_force(struct part *restrict p){
}


INLINE void compute_density_diffusive_term(struct part *restrict pi, struct part *restrict pj, float r2, float r, float wi_dx, float wj_dx, const float *dx){


}


INLINE void compute_density_diffusive_term_asym(struct part *restrict pi, const struct part *restrict pj, float r2, float r, float wi_dx, const float *dx){

}
#endif
