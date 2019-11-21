#ifndef SWIFT_SIN_BOUNDARY_MOTION_H
#define SWIFT_SIN_BOUNDARY_MOTION_H

#include <math.h>
#include "timeline.h"
#include "common_io.h"
#include "inline.h"
#include "error.h"
#include "physical_constants.h"
#include "part.h"

extern struct boundary_motion_parameters bom;

/**
 *  * @brief The parameters of the equation of state for the gas.
 *   *
 *    * This equation of state takes a single argument, the soundspeed.
 *     */
struct boundary_motion_parameters {
  double time_base; //Time base computed from engine_init
  double period; //Time period for the sine wave
  double max_acceleration; //Maximum acceleration
  int dimension; //Dimension of the acceleration - only a single dimension of value 0 (x), 1 (y), or 2 (z) is accepted.
};

//Must be called after engine_init
INLINE void init_boundary_motion(struct boundary_motion_parameters *bm,
                            double time_base,
                            struct swift_params *params){ 
  bm->time_base = time_base;
  bm->max_acceleration = parser_get_param_float(params, "Boundary_Motion:max_acceleration");
  bm->period = parser_get_param_float(params, "Boundary_Motion:sin_period");
  bm->dimension = parser_get_param_int(params, "Boundary_Motion:dimension");

}


INLINE void apply_boundary_motion(struct part *p, integertime_t ti_start){
  if(p->is_boundary){
    const double time = ti_start * bom.time_base;
    const double period_multiple = time / bom.period;
    const double value = period_multiple * 2 * M_PI;
    const double cos_value = cos(value);
    const double acceleration = bom.max_acceleration * cos_value;
    p->a_hydro[bom.dimension] += acceleration;
    


  }
}

#endif
