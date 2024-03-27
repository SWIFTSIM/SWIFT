//
// Created by yuyttenh on 29/03/22.
//

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_branch_FUNCTION,
   runner_doself_FUNCTION, calling the pairwise interaction function
   runner_iact_FUNCTION. */

#define PASTE(x, y) x##_##y

#define _DOPAIR_BRANCH(f) PASTE(runner_dopair_branch, f)
#define DOPAIR_BRANCH _DOPAIR_BRANCH(FUNCTION)

#define _DOPAIR(f) PASTE(runner_dopair, f)
#define DOPAIR _DOPAIR(FUNCTION)

#define _DOPAIR_BOUNDARY(f) PASTE(runner_dopair_boundary, f)
#define DOPAIR_BOUNDARY _DOPAIR_BOUNDARY(FUNCTION)

#define _DOSELF(f) PASTE(runner_doself, f)
#define DOSELF _DOSELF(FUNCTION)

#define _IACT(f) PASTE(runner_iact, f)
#define IACT _IACT(FUNCTION)

#define _IACT_BOUNDARY(f) PASTE(runner_iact_boundary, f)
#define IACT_BOUNDARY _IACT_BOUNDARY(FUNCTION)

#define _IACT_BOUNDARY_PARTICLES(f) PASTE(runner_doiact_boundary_particles, f)
#define IACT_BOUNDARY_PARTICLES _IACT_BOUNDARY_PARTICLES(FUNCTION)

#define _TIMER_DOSELF(f) PASTE(timer_doself, f)
#define TIMER_DOSELF _TIMER_DOSELF(FUNCTION)

#define _TIMER_DOPAIR(f) PASTE(timer_dopair, f)
#define TIMER_DOPAIR _TIMER_DOPAIR(FUNCTION)

void DOSELF(struct runner *restrict r, struct cell *restrict c);
void DOPAIR_BRANCH(struct runner *restrict r, struct cell *restrict ci,
                   struct cell *restrict cj);
void DOPAIR_BOUNDARY(struct runner *restrict r, struct cell *restrict c);
