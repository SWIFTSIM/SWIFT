#ifndef MOLTENI_DIFFUSIVE_H
#define MOLTENI_DIFFUSIVE_H

#include "part.h"
#include <math.h>

extern struct diffusive_term_parameters dtp;

struct diffusive_term_parameters {
  float delta;
  float soundspeed;
};

//Must be called after engine_init
INLINE void init_diffusive_term(struct diffusive_term_parameters *diff,
                            struct swift_params *params){

  diff->delta = parser_get_param_float(params, "Diffusive_term:delta"); 
  diff->soundspeed = parser_get_param_float(params, "EoS:soundspeed");

}



INLINE void reset_diffusive_term(struct part *restrict p){
  p->density_diffusive_data.delta = 0.0f;
}

INLINE void diffusive_end_force(struct part *restrict p){
  if(p->density_diffusive_data.delta == p->density_diffusive_data.delta){
    p->drho_dt += p->density_diffusive_data.delta;
  }
}

INLINE void compute_density_diffusive_term(struct part *restrict pi, struct part *restrict pj, float r2, float r, float wi_dx, float wj_dx, const float *dx){
  if(pi->is_boundary || pj->is_boundary){
    pi->density_diffusive_data.delta = NAN;
    pj->density_diffusive_data.delta = NAN; 
    return;
  }
#ifdef EOS_MULTIFLUID_TAIT
  if(pi->rho_base != pj->rho_base){
    return;
  }
#endif
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rho_i = pi->rho;
  const float rho_j = pj->rho;
  const float hi = pi->h;
  const float hj = pj->h;

  const float rhoi_over_rhoj_minus_1 = (rho_i / rho_j) - 1.f;
  const float r2_plus_eta2 = r2 + (0.01f*0.01f);
  const float inv_r2_plus_eta2 = 1.0f / r2_plus_eta2; 
  const float rr = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r_dot_grad_wi = wi_dx * rr;

  pi->density_diffusive_data.delta += mj * r_dot_grad_wi * 2.0 * dtp.delta * hi * dtp.soundspeed * (rhoi_over_rhoj_minus_1 * inv_r2_plus_eta2); 

/*  const float r_w_over_r2 = r * wi_dx / r2;
  const float mj_over_rhoj = mj / rho_j;
  pi->drho_dt += 2.0 * dtp.delta * hi * dtp.soundspeed * (rho_j - rho_i) * r_w_over_r2 * mj_over_rhoj;*/
//  if(pi->id == 4281){
//     printf("diffusive_term = %e %e %e\n",r2,  2.0 * dtp.delta * hi * dtp.soundspeed * (rho_j - rho_i) * r_w_over_r2 * mj_over_rhoj, rho_j-rho_i);
//  }

  const float rhoj_over_rhoi_minus_1 = (rho_j / rho_i) - 1.f;
  const float r_dot_grad_wj = wj_dx * rr;
  pj->density_diffusive_data.delta += mi * r_dot_grad_wj * 2.0 * dtp.delta * hj * dtp.soundspeed * (rhoj_over_rhoi_minus_1 * inv_r2_plus_eta2);

/*  const float mr_w_over_r2 = -r * wj_dx / r2;
  const float mi_over_rhoi = mi / rho_i;
  pj->drho_dt += 2.0 * dtp.delta * hj * dtp.soundspeed * (rho_i - rho_j) * mr_w_over_r2 * mi_over_rhoi;*/
//  if(pj->id == 4281){
//     printf("diffusive_term = %e %e %e\n",r2,  2.0 * dtp.delta * hj * dtp.soundspeed * (rho_i - rho_j) * mr_w_over_r2 * mi_over_rhoi, rho_j-rho_i);
//  }

}


INLINE void compute_density_diffusive_term_asym(struct part *restrict pi, const struct part *restrict pj, float r2, float r, float wi_dx, const float *dx){
  if(pi->is_boundary || pj->is_boundary){
    pi->density_diffusive_data.delta = NAN;
    return;
  }
#ifdef EOS_MULTIFLUID_TAIT
  if(pi->rho_base != pj->rho_base){
    return;
  }
#endif
  const float mj = pj->mass;
  const float rho_i = pi->rho;
  const float rho_j = pj->rho;
  const float hi = pi->h;


  const float rhoi_over_rhoj_minus_1 =  (rho_i / rho_j) - 1.f;
  const float r2_plus_eta2 = r2 + (0.01f*0.01f);
  const float inv_r2_plus_eta2 = 1.0f / r2_plus_eta2; 
  const float rr = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r_dot_grad_wi = wi_dx * rr;

  pi->density_diffusive_data.delta += mj * r_dot_grad_wi * 2.0 * dtp.delta * hi * dtp.soundspeed * (rhoi_over_rhoj_minus_1 * inv_r2_plus_eta2); 
  
  /*const float r_w_over_r2 = r * wi_dx / r2;
  const float mj_over_rhoj = mj / rho_j;
  pi->drho_dt += 2.0 * dtp.delta * hi * dtp.soundspeed * (rho_j - rho_i) * r_w_over_r2 * mj_over_rhoj;*/

}
#endif
