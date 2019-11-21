#ifndef MONAGHAN94_BOUNDARY_H
#define MONAGHAN94_BOUNDARY_H

#include "part.h"
#include "error.h"

extern struct boundary_parameters boundary_params;

struct boundary_parameters{
   float D; //Coefficient to scale the repulsive force
   float r0; //Initial seperation of the fluid particles.
   float p1; //Larger power for equation
   float p2; //Smaller power for equation
};


INLINE void init_boundary(struct boundary_parameters *bp, struct swift_params *params){
  bp->r0 = parser_get_param_float(params, "Boundary:r0");
  bp->D = parser_get_opt_param_float(params, "Boundary:D", 1.0);
  bp->p1 = parser_get_opt_param_float(params, "Boundary:p1", 4.0);
  bp->p2 = parser_get_opt_param_float(params, "Boundary:p2", 2.0);
  if(bp->p2 >= bp->p1){
    error("Illegal declaration of boundary as p2 must be < p1");
  }
}

INLINE void boundary_fluid_interaction(struct part *restrict pi, struct part * restrict pj, float r, float r2, const float *dx){
  if( r > boundary_params.r0 || (pi->is_boundary && pj->is_boundary) || !(pi->is_boundary || pj->is_boundary)) return;
  const float ir = 1.0 / r;
  const float r0_r = boundary_params.r0 * ir;
  const float term1 = pow(r0_r, boundary_params.p1);
  const float term2 = pow(r0_r, boundary_params.p2);
  const float non_directional_term = boundary_params.D*(term1 - term2);
  const float ir2 = ir*ir;
  /*if(non_directional_term != non_directional_term){
    error("Nan in non_directional_term %e %e %e %e %e %e %e %e", boundary_params.D, boundary_params.p2, boundary_params.p2, boundary_params.r0, term1, term2, r0_r, r);
  }*/
//  const float i_imass = 1.0 / pi->mass;
//  const float j_imass = 1.0 / pj->mass;

/*  if( non_directional_term * dx[0] * ir2 * i_imass > 1e7){
  struct part *p = pi;
  printf(
      "\n "
      "x=[%.9g, %.9g, %.9g], v=[%.9g, %.9g, %.9g], \n "
      " a=[%.9g, %.9g, %.9g], \n "
      "m=%.9g, a_const=[%.9g, %.9g, %.9g], P=%.9g, \n "
      "rho=%.9g, drho_dt=%.9g, h=%.9g \n "
      "time_bin=%d wakeup=%d is_boundary=%d \n",
      p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2],
      p->a_hydro[0], p->a_hydro[1], p->a_hydro[2],
      p->mass, p->a_constant[0], p->a_constant[1], p->a_constant[2], p->pressure,
      p->rho, p->drho_dt, p->h, p->time_bin, p->wakeup, p->is_boundary);

  p = pj;
  printf(
      "\n "
      "x=[%.9g, %.9g, %.9g], v=[%.9g, %.9g, %.9g], \n "
      " a=[%.9g, %.9g, %.9g], \n "
      "m=%.9g, a_const=[%.9g, %.9g, %.9g], P=%.9g, \n "
      "rho=%.9g, drho_dt=%.9g, h=%.9g \n "
      "time_bin=%d wakeup=%d is_boundary=%d \n",
      p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2],
      p->a_hydro[0], p->a_hydro[1], p->a_hydro[2],
      p->mass, p->a_constant[0], p->a_constant[1], p->a_constant[2], p->pressure,
      p->rho, p->drho_dt, p->h, p->time_bin, p->wakeup, p->is_boundary);

    error("Massive force detected, distance was r=%e\n", r);
  }*/

  pi->a_hydro[0] += non_directional_term * dx[0] * ir2;// * i_imass;
  pi->a_hydro[1] += non_directional_term * dx[1] * ir2;// * i_imass;
  pi->a_hydro[2] += non_directional_term * dx[2] * ir2;// * i_imass;

  pj->a_hydro[0] -= non_directional_term * dx[0] * ir2;// * j_imass; 
  pj->a_hydro[1] -= non_directional_term * dx[1] * ir2;// * j_imass; 
  pj->a_hydro[2] -= non_directional_term * dx[2] * ir2;// * j_imass; 
/*  if(pi->id == 2116){
    printf("x_boundary = %e, y_boundary = %e\n", non_directional_term * dx[0] * ir2,  non_directional_term * dx[1] * ir2);
  }
  if(pj->id == 2116){
    printf("x_boundary = %e, y_boundary = %e\n",-non_directional_term * dx[0] * ir2,  non_directional_term * dx[1] * ir2);
  }*/
}

INLINE void boundary_fluid_interaction_nonsym(struct part *restrict pi, const struct part * restrict pj, float r, float r2, const float *dx){
  if( r > boundary_params.r0 || (pi->is_boundary && pj->is_boundary) || !(pi->is_boundary || pj->is_boundary)) return;
  const float r0_r = boundary_params.r0 / r;
  const float term1 = pow(r0_r, boundary_params.p1);
  const float term2 = pow(r0_r, boundary_params.p2);
  const float non_directional_term = boundary_params.D*(term1 - term2);
  const float ir2 = 1.0 / r2;
  const float i_imass = 1.0 / pi->mass;

  pi->a_hydro[0] += non_directional_term * dx[0] * ir2 * i_imass;
  pi->a_hydro[1] += non_directional_term * dx[1] * ir2 * i_imass;
  pi->a_hydro[2] += non_directional_term * dx[2] * ir2 * i_imass;
}

INLINE void boundary_end_force( struct part *p ){
  if(p->is_boundary){
    p->a_hydro[0] = 0.0;
    p->a_hydro[1] = 0.0;
    p->a_hydro[2] = 0.0;
  }
} 

#endif
