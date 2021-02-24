.. locks
   Mladen Ivkovic Feb 2021

.. _task_locks:
.. highlight:: c


Task Locks
=============

Swift doesn't only deal with dependencies, but with data conflicts as
well.  We can't have two independent tasks working on the same data
concurrently, even if the tasks have no dependencies between them.
For this reason, each task locks the data it works on when it begins,
and unlocks the data again once it's finished. Data here refers to the
particle (gas, gravity, stars,...) content of a cell as well as of the
hierarchy of cells in the tree at levels closer to the root than the
cell itself.  By default, it is assumed that the task is doing work on
hydro particles. If your task requires other particle types, you will
need to specify that in ``src/task.c``. Suppose you have implemented a
task with type ``task_type_new`` that requires both stars and hydro
particles: ::


    /**
     * @brief Try to lock the cells associated with this task.
     *
     * @param t the #task.
     */
    int task_lock(struct task *t) {

      const enum task_types type = t->type;
      const enum task_subtypes subtype = t->subtype;
      struct cell *ci = t->ci, *cj = t->cj;

      switch (type) {
        /* lots of stuff */

        case task_type_new:

          /* is the data locked already? */
          if (ci->hydro.hold || ci->stars.hold) return 0;

          /* if it's not locked already, lock first particle type (hydro)
           * if something locked it in the meantime, exit with failure. */
          if (cell_locktree(ci) != 0) return 0;

          /* if it's not locked already, lock first particle type (stars)
           * if something locked it in the meantime, first unlock what you
           * locked, and then exit with failure. */
          if (cell_slocktree(ci) != 0) {
            cell_unlocktree(ci);
            return 0;
          }

      /* lots of other stuff */
      }

      /* lots of other stuff */

    }



Similarly, don't forget to write your unlocking routines in ``task_unlock()`` !
