//
// Created by yuyttenh on 11/04/22.
//

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "space_getsid.h"
#include "timers.h"

/* Import the gradient loop functions. */
#define FUNCTION slope_estimate
#define FUNCTION_TASK_LOOP TASK_LOOP_SLOPE_ESTIMATE
#include "runner_doiact_functions_grid_hydro.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the gradient limiter loop functions. */
#define FUNCTION slope_limiter
#define FUNCTION_TASK_LOOP TASK_LOOP_SLOPE_LIMITER
#include "runner_doiact_functions_grid_hydro.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the flux loop functions. */
#define FUNCTION flux_exchange
#define FUNCTION_TASK_LOOP TASK_LOOP_FLUX_EXCHANGE
#include "runner_doiact_functions_grid_hydro.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP