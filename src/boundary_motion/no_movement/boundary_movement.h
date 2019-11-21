#ifndef SWIFT_NO_BOUNDARY_MOTION_H
#define SWIFT_NO_BOUNDARY_MOTION_H

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
 
};


INLINE void init_boundary_motion(struct boundary_motion_parameters *bdm,
                            double time_base,
                            struct swift_params *params){ 


}


INLINE void apply_boundary_motion(struct part *p, integertime_t ti_start){

}

#endif
