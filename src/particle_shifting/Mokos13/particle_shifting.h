#ifndef MOKOS13_SHIFTING_H
#define MOKOS13_SHIFTING_H

#include "part.h"

extern struct particle_shifting_parameters psp;

struct particle_shifting_parameters {
  float a_fsm;
  float a_fst;
  float A;
};


//Must be called after engine_init
INLINE void init_shifting_parameters(struct particle_shifting_parameters *shift,
                            struct swift_params *params){
#ifdef HYDRO_DIMENSION_2D
  shift->a_fsm = 2.0;
  shift->a_fst = 1.5;
#elif defined(HYDRO_DIMENSION_3D)
  shift->a_fsm = 3.0;
  shift->a_fst = 2.75;
#else
  shift->a_fsm = 0.0;
  shift->a_fst = 0.0;
#endif
  shift->A = parser_get_opt_param_float(params, "Particle Shifting: A", 2.0);
}

INLINE void reset_shifting_data(struct part *p){
  p->particle_shifting_data.particle_concentration[0] = 0.0f;
  p->particle_shifting_data.particle_concentration[1] = 0.0f;
  p->particle_shifting_data.particle_concentration[2] = 0.0f;
  p->particle_shifting_data.particle_divergence = 0.0f;
}

INLINE void compute_particle_shifting(struct part *p, float dt){
  if(p->particle_shifting_data.particle_concentration[0] != p->particle_shifting_data.particle_concentration[0]){
    return;
  }
  const float speed = sqrt(p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2]);
  const float D = psp.A * p->h * speed * dt;

  const float fsm_minus_fst = psp.a_fsm - psp.a_fst;
  const float inv_fsm_minus_fst = 1.0 / fsm_minus_fst;
  const float div_r_minus_fst = p->particle_shifting_data.particle_divergence - psp.a_fst;
  /*if( (p->x[1] > 13.0 && p->x[0] < 2.5) || (p->x[1] < 2.25 && p->x[0] < 2.5)){
    printf("id = %llu, x[1] = %f, div_r_minus_fst=%f,div=%f\n",p->id, p->x[1],div_r_minus_fst, p->particle_shifting_data.particle_divergence);
  } */
  if(div_r_minus_fst < 0.0){

    const float a_fsc = ( div_r_minus_fst ) * inv_fsm_minus_fst;
    p->x[0] -= a_fsc * D * p->particle_shifting_data.particle_concentration[0];
    p->x[1] -= a_fsc * D * p->particle_shifting_data.particle_concentration[1];
    p->x[2] -= a_fsc * D * p->particle_shifting_data.particle_concentration[2];

  }else if (div_r_minus_fst >= 0.0){

    p->x[0] -= D * p->particle_shifting_data.particle_concentration[0];
    p->x[1] -= D * p->particle_shifting_data.particle_concentration[1];
    p->x[2] -= D * p->particle_shifting_data.particle_concentration[2];

  }
}

INLINE void compute_shifting_term(struct part *restrict pi, struct part *restrict pj,  float r2,  float r,  float wi_dx,  float wj_dx, const float *dx){
  if(pi->is_boundary || pj->is_boundary){
    pi->particle_shifting_data.particle_concentration[0] = NAN;
    pj->particle_shifting_data.particle_concentration[0] = NAN;
    return;
  }
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rho_i = pi->rho;
  const float rho_j = pj->rho;


  const float mj_over_rho_j = mj/rho_j;
  //const float dx_dot_wi_dx = fabsf(dx[0] * wi_dx*dx[0]) + fabsf(dx[1] * wi_dx*dx[1]) + fabsf(dx[2]*wi_dx*dx[2]);
  pi->particle_shifting_data.particle_concentration[0] += mj_over_rho_j * wi_dx * dx[0] ;
  pi->particle_shifting_data.particle_concentration[1] += mj_over_rho_j * wi_dx * dx[1] ;
  pi->particle_shifting_data.particle_concentration[2] += mj_over_rho_j * wi_dx * dx[2] ;
/*  if( 0&& pi->id == 1218){
    printf("%e\n", -mj_over_rho_j*dx_dot_wi_dx);
  }*/ 
  pi->particle_shifting_data.particle_divergence -= mj_over_rho_j*r*wi_dx;

  const float mi_over_rho_i = mi/rho_i;
//  const float dx_dot_wj_dx = fabsf(dx[0] * wj_dx*dx[0]) + fabsf(dx[1] * wj_dx*dx[1]) + fabsf(dx[2]*wj_dx*dx[2]);
  pj->particle_shifting_data.particle_concentration[0] += mi_over_rho_i * wj_dx * -dx[0];
  pj->particle_shifting_data.particle_concentration[1] += mi_over_rho_i * wj_dx * -dx[1];
  pj->particle_shifting_data.particle_concentration[2] += mi_over_rho_i * wj_dx * -dx[2];

/*  if(0 && pj->id == 1218){
    printf("%e\n", -mi_over_rho_i*dx_dot_wj_dx);
  }*/ 
  pj->particle_shifting_data.particle_divergence -= mi_over_rho_i*r*wj_dx;


}


INLINE void compute_shifting_term_nonsym(struct part *restrict pi, const struct part *restrict pj,  float r2,  float r,  float wi_dx, const float *dx){

  if(pi->is_boundary || pj->is_boundary){
    pi->particle_shifting_data.particle_concentration[0] = NAN;
    return;
  }
  const float mj = pj->mass;
  const float rho_j = pj->rho;

  const float mj_over_rho_j = mj/rho_j;
  const float dx_dot_wi_dx = dx[0] * wi_dx*dx[0] + dx[1] * wi_dx*dx[1] + dx[2]*wi_dx*dx[2];
  pi->particle_shifting_data.particle_concentration[0] += mj_over_rho_j * wi_dx * dx[0] ;
  pi->particle_shifting_data.particle_concentration[1] += mj_over_rho_j * wi_dx * dx[1] ;
  pi->particle_shifting_data.particle_concentration[2] += mj_over_rho_j * wi_dx * dx[2] ;
  
  pi->particle_shifting_data.particle_divergence += mj_over_rho_j*dx_dot_wi_dx;
}





#endif
