/* TODO: dox, license, etc */
#include "runner_doiact_rt.h"


/* TODO: TEMPORARY */
#include <stdio.h>



void DOSELF1_RT(struct runner *r, struct cell *c, int timer){
  TIMER_TIC; 
  printf("called DOSELF1_RT for cell %p\n", c);
  if (timer) TIMER_TOC(TIMER_DOSELF_RT);
}

void DOPAIR1_RT(struct runner *r, struct cell *ci, struct cell *cj, int timer){
  TIMER_TIC; 
  printf("called DOPAIR1_RT for cells %p, %p\n", ci, cj);
  if (timer) TIMER_TOC(TIMER_DOPAIR_RT);
}

void DOSELF1_BRANCH_RT(struct runner *r, struct cell *c, int timer) {
  DOSELF1_RT(r, c, timer);
}

void DOPAIR1_BRANCH_RT(struct runner *r, struct cell *ci, struct cell *cj, int timer) {
  DOPAIR1_RT(r, ci, cj, timer);
}

void DOSUB_SELF1_RT(struct runner *r, struct cell *c, int timer) {
  DOSELF1_RT(r, c, timer);
}

void DOSUB_PAIR1_RT(struct runner *r, struct cell *ci, struct cell *cj, int timer) {
  DOPAIR1_RT(r, ci, cj, timer);
}
