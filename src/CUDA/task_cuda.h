#ifndef SWIFT_CUDA_TASK_H 
#define SWIFT_CUDA_TASK_H

extern "C" {
#include "../task.h"
}

struct task_cuda {

  /* Indices of the cell (we don't want to use pointers) */
  int ci, cj;

  /* Indices of the tasks this task unlocks. */
  int *unlocks;

  /* Flags used to carry additional information (e.g. sort directions) */
  int flags;

  /* Rank of a task in the order */
  int rank;

  /* Weight of the task */
  int weight;

  /* Number of tasks unlocks by this one */
  int nr_unlock_tasks;

  /* Number of unsatisfied dependencies*/
  int wait;

  /* Type of the task */
  enum task_types type;

  /* Sub-type of the task */
  enum task_subtypes subtype;

  /* Should this task be skipped */
  char skip;

  /* Is this task implicit (may not need on GPU?) */
  char implicit;

  /* Pointer to the CPU task used for initialisation */
  struct task *task;

  /* Size of unlock array during initialisation. */
  int size_unlocks;

#ifdef REDUCED_TRANSFER
  /* Pointer to the loaded/unloaded cell for load/unload tasks*/
  struct cell *cell;
#endif

#ifdef CUDA_TASK_TIMERS
  /* Executing block*/
  int blockID;
  /* Start/end time for task */
  long long int tic , toc;
#endif
};


#endif /* SWIFT_CUDA_TASK_H */
