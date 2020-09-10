.. Neighbour Loop Task
   Mladen Ivkovic Sep 2020

.. _task_adding_your_own_neighbour_loop:
.. highlight:: c

Adding a Particle Interaction/Neighbour Loop Task
=================================================

There are quite a few subtle and not so subtle differences when adding tasks that include
particle neighbour loops and/or particle interactions compared to "independent" tasks, where
no particles interact with each other, but work is done on them individually.

Particle interactions are handled on a cell basis. If only particles of one cell
interact with each other, the task is referred to as a  ``self`` task. When particles of
two different cells interact with each other, we call the task type ``pair``. For a 
particle neighbour loop, one typically requires both ``self`` and ``pair`` tasks.
Furthermore, sometimes interaction tasks which have too much of a workload can get split 
into smaller units of work with tasks of the type ``sub_self`` and ``sub_pair``.

For the next paragraphs, let's assume that we want to implement the task ``new_iact``.




Adding it to the Task List
--------------------------
First you will need to add it to the task list situated in ``task.h`` and ``task.c``.
Here we aren't adding a task **type**, as the type will always be a ``pair`` or a 
``self``, but **subtypes**.

In ``task.h``, you need to provide an additional entry to the enum ``task_subtypes`` 
(e.g. ``task_subtype_new_iact``). The last entry ``task_subtype_count`` should always 
stay at the end as it is a counter of the number of elements.
For example::

    enum task_subtypes {
      task_subtype_none = 0,        task_subtype_density,           task_subtype_gradient,
      task_subtype_force,           task_subtype_limiter,           task_subtype_grav,
      task_subtype_external_grav,   task_subtype_tend_part,         task_subtype_tend_gpart,
      task_subtype_tend_spart,      task_subtype_tend_sink,         task_subtype_tend_bpart,
      task_subtype_xv,              task_subtype_rho,               task_subtype_part_swallow,
      task_subtype_bpart_merger,    task_subtype_gpart,             task_subtype_multipole,
      task_subtype_spart,           task_subtype_stars_density,     task_subtype_stars_feedback,
      task_subtype_new_iact,        task_subtype_count
    } __attribute__((packed));


In ``task.c``, you will find an array containing the name of each task and need to add your own (e.g. ``new_iact``).
Be careful with the order that should be the same than in the previous list.
For example::

    const char *subtaskID_names[task_subtype_count] = {
            "none",         "density",      "gradient",         "force", 
            "limiter",      "grav",         "external_grav",    "tend_part",
            "tend_gpart",   "tend_spart",   "tend_sink",        "tend_bpart",
            "xv",           "rho",          "part_swallow",     "bpart_merger",
            "gpart",        "multipole",    "spart",            "stars_density",
            "stars_feedback","sf_count",    "bpart_rho",        "sink",
            "new_iact"
   };



Adding it to the Cells
----------------------

Each cell contains a list to its tasks and therefore you need to provide a link for it.

In ``cell.h``, add a pointer to a task in the ``struct cell``.  In order to stay clean, 
please put the new task in the same group (e.g. ``struct hydro{...}`` inside ``struct cell``)
than the other tasks.
We won't be adding just one task though, but an entire (linked) list of them, since we're
going to need a ``self`` type task and multiple ``pair`` type tasks to have a complete
neighbour loop. So instead of pointing to a single task, we store a struct ``link`` in
the cell struct.  For example::

  struct cell {
    /* Lot of stuff before. */
    
    /*! Hydro variables */
    struct {
        /*! Pointer to the #part data. */
        struct part *parts;

        /*! Pointer to the #xpart data. */
        struct xpart *xparts;

        /* Lot of stuff */

        /*! Task for sorting the stars again after a SF event */
        struct task *stars_resort;

        /*! My new interaction task */
        struct link *new_iact;

        /* Lot of stuff after */

    } hydro;

    /* Lot of stuff after */
  }


Adding a new Timer
------------------

As SWIFT is HPC oriented, any new task need to be optimized.
It cannot be done without timing the function.

In ``timers.h``, you will find an enum that contains all the tasks.
You will need to add yours inside it.
For example::

  enum {
    timer_none = 0,
    timer_prepare,
    timer_init,
    timer_drift_part,
    timer_drift_gpart,
    timer_kick1,
    timer_kick2,
    timer_timestep,
    timer_endforce,
    timer_dosort,
    timer_doself_density,
    timer_doself_gradient,
    timer_doself_force,
    timer_dopair_density,
    timer_dopair_gradient,
    timer_dopair_force,
    timer_dosub_self_density,
    timer_dosub_self_gradient,
    timer_dosub_self_force,
    timer_dosub_pair_density,
    timer_dosub_pair_gradient,
    timer_dosub_pair_force,
    timer_doself_subset,
    timer_dopair_subset,
    timer_dopair_subset_naive,
    timer_dosub_subset,
    timer_do_ghost,
    timer_do_extra_ghost,
    timer_dorecv_part,
    timer_do_cooling,
    timer_gettask,
    timer_qget,
    timer_qsteal,
    timer_locktree,
    timer_runners,
    timer_step,
    timer_cooling,
    timer_new_iact,
    timer_count,
  };

As for ``task.h``,
you will need to give a name to your timer in ``timers.c``::

  const char* timers_names[timer_count] = {
    "none",
    "prepare",
    "init",
    "drift_part",
    "kick1",
    "kick2",
    "timestep",
    "endforce",
    "dosort",
    "doself_density",
    "doself_gradient",
    "doself_force",
    "dopair_density",
    "dopair_gradient",
    "dopair_force",
    "dosub_self_density",
    "dosub_self_gradient",
    "dosub_self_force",
    "dosub_pair_density",
    "dosub_pair_gradient",
    "dosub_pair_force",
    "doself_subset",
    "dopair_subset",
    "dopair_subset_naive",
    "dosub_subset",
    "do_ghost",
    "do_extra_ghost",
    "dorecv_part",
    "gettask",
    "qget",
    "qsteal",
    "locktree",
    "runners",
    "step",
    "cooling",
    "new_iact",
  };


You can now easily time
your functions by using::

  TIMER_TIC;
  /* Your complicated functions */
  if (timer) TIMER_TOC(timer_new_iact);


Adding your Task to the System
------------------------------

Now the tricky part happens.
SWIFT is able to deal automatically with the conflicts between tasks, but unfortunately 
cannot understand the dependencies.

To implement your new task in the task system, you will need to modify a few functions 
in ``engine_maketasks.c``.

First, you will need to add mainly three functions: ``scheduler_addtask``, ``engine_addlink``,
``scheduler_addunlocks`` in the ``engine_make_extra_hydroloop_tasks_mapper`` functions 
(depending on the type of task you implement, you will need to write it to a different 
function). The (hydro) particle interaction tasks are first created only for the density 
loop, and then replicated in ``engine_make_extra_hydroloop_tasks_mapper`` for everything
else.

In ``engine_make_extra_hydroloop_tasks_mapper``, we add the task through the following 
call::


    struct task *t_new_iact = NULL;

    /* ... lots of stuff ... */

    /* Self-interaction? */
    else if (t_type == task_type_self && t_subtype == task_subtype_density) {

      /* ... lots of stuff ... */

      t_new_iact = scheduler_addtask(sched, task_type_self, task_subtype_new_iact, 
                                     flags, 0, ci, NULL);

      /* Link the tasks to the cells */
      engine_addlink(e, &ci->new_iact, t_new_iact);

      /* Create the task dependencies */
      scheduler_addunlock(sched, ci->task_that_unlocks_this_one, t_new_iact);
      scheduler_addunlock(sched, t_new_iact, ci->task_that_will_be_unlocked_by_this_one);
    }

    /* Otherwise, pair interaction? */
    else if (t_type == task_type_pair && t_subtype == task_subtype_density) {

      /* ... lots of stuff ... */

      t_new_iact = 
        scheduler_addtask(sched, task_type_pair, task_subtype_new_iact, 
                          flags, 0, ci, cj);
      engine_addlink(e, &ci->new_iact, t_new_iact);
      engine_addlink(e, &cj->new_iact, t_new_iact);

      /* ... lots of stuff ... */

      if (ci->nodeID == nodeID) {

        /* ... lots of stuff ... */

        scheduler_addunlock(sched, ci->task_that_unlocks_this_one, t_new_iact);
        scheduler_addunlock(sched, t_new_iact, ci->task_that_will_be_unlocked_by_this_one);
      }

      if (cj->nodeID == nodeID) {

        if (ci->hydro.super != cj->hydro.super) {

          /* ... lots of stuff ... */

          scheduler_addunlock(sched, cj->task_that_unlocks_this_one, t_new_iact);
          scheduler_addunlock(sched, t_new_iact, cj->task_that_will_be_unlocked_by_this_one);

        }
      }
    }

    /* Otherwise, sub-self interaction? */
    else if (t_type == task_type_sub_self &&
             t_subtype == task_subtype_density) {

        /* You need to do the same as for task_type_self above */
    }

    /* Otherwise, sub-pair interaction? */
    else if (t_type == task_type_sub_pair &&
             t_subtype == task_subtype_density) {

        /* You need to do the same as for task_type_pair above */
    }





The next step is to activate your task
in ``engine_marktasks_mapper`` in ``engine_marktasks.c``::


  /* Single-cell task? */
  if (t_type == task_type_self || t_type == task_type_sub_self) {

    /* ... lots of stuff ...  */

    else if (t_subtype == task_subtype_new_iact) {
      scheduler_activate(s, t);
    }
  }

  /* Pair? */
  else if (t_type == task_type_pair || t_type == task_type_sub_pair) {

    /* ... lots of stuff ...  */

    else if (t_subtype == task_subtype_new_iact) {
      scheduler_activate(s, t);
    }
  }


Then you will need to update the estimate for the number of tasks in 
``engine_estimate_nr_tasks`` in ``engine.c`` by modifying ``n1`` or ``n2``.
``n1`` is the expected maximal number of tasks per top-level/super cell. ``n2``
``n2`` is the expected maximum number of tasks for all other cells, independent
of the depth of the tree. Most likely ``n2`` won't need updating, and you will
only need to update ``n1``. As to how to update ``n1``, you just need to count
the number of tasks that you will be adding, e.g. 1 self + (3^3-1)/2 = 13 pair 
tasks + 1 ghost, etc... All these numbers can be overwritten at run time by 
the user anyway in the parameter file (``Scheduler: tasks_per_cell``).

and give the task an estimate of the computational cost that it will have in 
``scheduler_reweight`` in  ``scheduler.c``::

      case task_type_self:
        if (t->subtype == task_subtype_grav) {
          cost = 1.f * (wscale * gcount_i) * gcount_i;
        /* ... lots of stuff ... */
        else if (t->subtype == task_subtype_new_iact)
          cost = 1.f * wscale * scount_i * count_i;
        else
          error("Untreated sub-type for selfs: %s",
                subtaskID_names[t->subtype]);
        break;

Similarly, you'll need to update ``case task_type_sub_self``, ``task_type_pair``, 
and ``task_type_sub_pair`` as well.



Initially, the engine will need to skip the task that updates the particles.
If this is the case for your task, you will need to add it in ``engine_skip_force_and_kick``.

Finally, you also need to initialize your new variables and pointers in 
``space_rebuild_recycle_mapper`` in ``space.c``. Additionally, you need to 
initialize the ``link`` structs in ``cell_clean_links`` in ``cell.c``.




Implementing your Task
----------------------

The last part is situated in ``runner_main.c``, where the actual functions executed
by the task are called inside the function in ``runner_main`` in the switch::

    /* Different types of tasks... */
    switch (t->type) {
      case task_type_self:
        if (t->subtype == task_subtype_density)
          runner_doself1_branch_density(r, ci);
        /* ... lots of stuff ... */
        else if (t->subtype == task_subtype_new_iact)
          runner_doself_branch_new_iact(r, ci, 1);
        else
          error("Unknown/invalid task subtype (%s).",
                subtaskID_names[t->subtype]);
        break;
        
      case task_type_pair:
        /* ... lots of stuff ... */
        else if (t->subtype == task_subtype_new_iact)
          runner_dopair_branch_new_iact(r, ci, cj, 1);
        else
          error("Unknown/invalid task subtype (%s/%s).",
                taskID_names[t->type], subtaskID_names[t->subtype]);
        break;

      case task_type_sub_self:
        /* ... lots of stuff ... */
        else if (t->subtype == task_subtype_new_iact)
          runner_dosub_self_new_iact(r, ci, 1);
        else
          error("Unknown/invalid task subtype (%s/%s).",
                taskID_names[t->type], subtaskID_names[t->subtype]);
        break;

      case task_type_sub_pair:
        /* ... lots of stuff ... */
        else if (t->subtype == task_subtype_new_iact)
          runner_dosub_pair_new_iact(r, ci, cj, 1);
        else
          error("Unknown/invalid task subtype (%s/%s).",
                taskID_names[t->type], subtaskID_names[t->subtype]);
        break;



The functions ``runner_doself1_branch_density``, ``runner_dopair_branch_new_iact``,
``runner_dosub_self_new_iact``,  and ``runner_dosub_pair_new_iact`` still need to be
implemented by you. If you only plan on doing this type of particle interaction once
per time step, you can get away with directly implementing these functions and call 
it a day. But if you intend to use the same kind of particle loop more than once, as 
it's done in e.g. the hydro density and force loops, it's better to construct the
functions using macros.
For example, you could have a file ``runner_doiact_my_stuff.h``::

    /* File runner_doiact_my_stuff.h */

    #define PASTE(x, y) x##_##y

    #define _DOSELF1_BRANCH_NEW(f) PASTE(runner_doself_branch, f)
    #define DOSELF1_BRANCH_NEW _DOSELF1_BRANCH_NEW(FUNCTION)

    #define _DOPAIR1_BRANCH_NEW(f) PASTE(runner_dopair_branch, f)
    #define DOPAIR1_BRANCH_NEW _DOPAIR1_BRANCH_NEW(FUNCTION)

    #define _DOSUB_PAIR1_NEW(f) PASTE(runner_dosub_pair, f)
    #define DOSUB_PAIR1_NEW _DOSUB_PAIR1_NEW(FUNCTION)

    #define _DOSUB_SELF1_NEW(f) PASTE(runner_dosub_self, f)
    #define DOSUB_SELF1_NEW _DOSUB_SELF1_NEW(FUNCTION)

    #define _IACT_NEW(f) PASTE(runner_iact, f)
    #define IACT_NEW _IACT_NEW(FUNCTION)

    void DOSELF1_BRANCH_NEW(struct runner *r, struct cell *c, int timer);
    void DOPAIR1_BRANCH_NEW(struct runner *r, struct cell *ci, struct cell *cj, 
                           int timer);

    void DOSUB_SELF1_NEW(struct runner *r, struct cell *ci, int timer);
    void DOSUB_PAIR1_NEW(struct runner *r, struct cell *ci, struct cell *cj,
                           int timer);



And a second file, ``runner_doiact_function_my_stuff.h``, where you define those
functions which have been declared using the macros, e.g. ::


    #include "runner_doiact_my_stuff.h"

    void DOSELF1_BRANCH_NEW(struct runner *r, struct cell *c, int timer) {
      /* do your stuff, call IACT_NEW(...) at some point...*/
    }

    void DOPAIR1_BRANCH_NEW(struct runner *r, struct cell *ci, struct cell *cj, int timer) {
      /* do your stuff, call IACT_NEW(...) at some point...*/
    }

    void DOSUB_SELF1_NEW(struct runner *r, struct cell *c, int timer) {
      /* do your stuff, call IACT_NEW(...) at some point...*/
    }

    void DOSUB_PAIR1_NEW(struct runner *r, struct cell *ci, struct cell *cj, int timer) {
      /* do your stuff, call IACT_NEW(...) at some point...*/
    }


Then we also need a ``runner_doiact_my_suff.c`` file where the functions declared in
``runner_doiact_my_suff.h`` are defined by including them with ``FUNCTION`` defined::


    #include "../config.h"
    /* other includes too... */

    /* Import the new interaction loop functions. */
    #define FUNCTION new_iact
    #include "runner_doiact_functions_my_stuff.h"
    #undef FUNCTION




Finally, we include them in ``runner_main.c`` as follows::

    /* ... lots of includes and stuff ... */

    /* Import new interaction loop functions. */
    #define FUNCTION new_iact
    #include "runner_doiact_my_suff.h"
    #undef FUNCTION

    /**
     * @brief The #runner main thread routine.
     *
     * @param data A pointer to this thread's data.
     */
    void *runner_main(void *data) {
        /* ... */
    }


The functions ``runner_doself_branch_density``, ``runner_dopair_branch_new_iact``,
``runner_dosub_self_new_iact``,  and ``runner_dosub_pair_new_iact`` will be properly
found and linked this way. All that's left for you to do is to write the function
into which ``IACT_NEW`` will expand, in the above case it would be ``runner_iact_new_iact``.





Finalizing your Task
--------------------

Now that you have done the easiest part, you can start debugging by implementing a 
test and/or an example. Before creating your merge request with your new task, do 
not forget the most funny part that consists in writing a nice and beautiful 
documentation ;)


Things to Keep in Mind
----------------------

- If you are inserting a new neighbour loop in between existing loops, or want to
  insert more than one neighbour loop, usually a new ghost task in between them is
  also needed.

- Neighbour loops may also require MPI communication tasks.
