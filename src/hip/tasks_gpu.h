/* Config parameters. */
#include "../config.h"

struct tasks_self_gpu {
  struct task_gpu *tgpu;
};

/**
 * @brief A task to be run by the #scheduler.
 */
struct task_gpu {

  /*! Pointers to the cells this task acts upon */
  struct cell *ci, *cj;

  /*! List of tasks unlocked by this one */
  struct task_gpu **unlock_tasks;

  /*! Flags used to carry additional information (e.g. sort directions) */
  long long flags;

#ifdef WITH_MPI

  /*! Buffer for this task's communications */
  void *buff;

  /*! MPI request corresponding to this task */
  MPI_Request req;

#endif

  /*! Rank of a task in the order */
  int rank;

  /*! Weight of the task */
  float weight;

  /*! Number of tasks unlocked by this one */
  int nr_unlock_tasks;

  /*! Number of unsatisfied dependencies */
  int wait;

  /*! Type of the task */
  enum task_types type;

  /*! Sub-type of the task (for the tasks that have one */
  enum task_subtypes subtype;

  /*! Should the scheduler skip this task ? */
  char skip;

  /*! Is this task implicit (i.e. does not do anything) ? */
  char implicit;

#ifdef SWIFT_DEBUG_TASKS
  /*! ID of the queue or runner owning this task */
  short int rid;

  /*! Information about the direction of the pair task */
  short int sid;
#endif

  /*! Start and end time of this task */
  ticks tic, toc;

  /* Total time spent running this task */
  ticks total_ticks;

#ifdef SWIFT_DEBUG_CHECKS
  /* When was this task last run? */
  integertime_t ti_run;
#endif /* SWIFT_DEBUG_CHECKS */
};
