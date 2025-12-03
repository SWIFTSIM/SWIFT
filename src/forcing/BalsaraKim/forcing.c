/* Config parameters */
#include <config.h>

/* Standard includes. */
#include <float.h>
#include <string.h>
#include <stdio.h> 

/* Local includes. */
#include "forcing.h"
 #include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"
#include "hydro.h"

/**  
 * @brief updates the forcing terms
 *   
 * increases the current supernova index after one has happend
 * 
 * @param terms The #forcing_temrs properties of the run
 * @param time_old The previous system time
 */
void forcing_update(struct forcing_terms *terms, const double time_old) {

  /* if the current time is later than the latest SN event */
  if (time_old >= terms->times[terms->size - 1]) {
    /* we do not want any more energy injections */
    terms->final_injection = 1;
  }

  /* else, if the next SN has happened, update the index to the next SN event */
  else if (time_old >= terms->times[terms->t_index]) {
    message("%d particles passed time condition", terms->counter);
    message("updating forcing term index at time: %f", time_old);
    terms->t_index += 1;
    terms->counter = 0;
  }

}
