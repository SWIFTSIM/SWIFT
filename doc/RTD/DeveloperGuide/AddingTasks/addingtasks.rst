.. _NewTask:

How to add a new task to SWIFT?
=================================
.. highlight:: c



.. toctree::
   :maxdepth: 0

This tutorial will step through how to add a new task to swift. First we will go through the 
idealology of adding a new task to SWIFT. This will be followed by an example of how to add a task
for an imposed external gravitational field to SWIFT and a task to include "cooling" to the gas particles.

In the simplest case adding a new tasks requires changes to five files, namely:

* task.h
* cell.h
* timers.h
* task.c
* engine.c

Further, implementation details of what the task will then do should be added to another file
(for example runner_myviptask.c) which will contain the actual task implementation.

So now lets look at what needs to change in each of the files above, starting with task.h

--------------
**task.h**
--------------
Within task.h there exists a structure of the form:: 

    /* The different task types. */
    enum task_types {
        task_type_none = 0,
        task_type_sort,
        task_type_self,
        task_type_pair,
	      .
	      .
	      .
	task_type_my_new_task,
	task_type_psort,
	task_type_split_cell,
	task_type_count
    };

Within this task structure your new task should be added. Add the task entry anywhere in the struct before the
task_type_count member. This last entry is used to count the number of tasks and must always be the last entry.

--------------
**task.c**
--------------

Within task.c the addition of the new task type must be include in a character list (which at the moment is only 
used for debugging purposes)::

     /* Task type names. */
     const char *taskID_names[task_type_count] = {
     "none",  "sort",    "self",    "pair",    "sub",
     "ghost", "kick1",   "kick2",   "send",    "recv",
     "link",  "grav_pp", "grav_mm", "grav_up", "grav_down",
     "my_new_task", "psort", "split_cell"};

The new task type should be added to this list in the same order as it was added within the task_types struct in 
task.h

--------------
**cell.h**
--------------

cell.h contains pointers to all of the tasks associated with that cell. You must include your new task type 
here e.g::

   struct task *my_new_task;

--------------
**timers.h**
--------------

Within timers.h is an enumerated list of timers associated with each task. The timers measure the time required 
to execute a given task and this information is used in improve scheduling the task in future iterations::

    /* The timers themselves. */
    enum {
      timer_none = 0,
      timer_prepare,
      timer_kick1,
           .
	   .
	   .
      timer_new_task,
      timer_step,
      timer_count,
    };

--------------
**engine.c**
--------------

Finally, in engine.c the new task is added so that the scheduler knows to include the task in the list of tasks 
to be scheduled. Knowing where to add the task in engine.c is a little bit more difficult. This will depend on 
the type of task involved and whether it is a task that acts only on a individual particle independent of other 
particles (e.g. a cooling a task) or whether the task depends on other tasks (e.g. density, force or feedback). 

If we assume that the task is a particle only task then the first place to modify is the engine_mkghosts() 
function. Within this function the new task must be added to the list of tasks 
(within the c->nodeID == e->nodeID if clause)::

       /* Generate the external gravity task*/
         c->my_new_task = scheduler_addtask(s, task_type_my_new_task, task_subtype_none, 0, 0,
                                   c, NULL, 0);


That's pretty much it - but what about dependencies and conflicts?
Remember SWIFT automatically handles conflicts (by understanding which tasks need to write to the same data) so 
you (the developer) don't need to worry about conflicts. Dependencies do however need to be managed and they will
be task specific. The following two examples, implementing cooling and an imposed external gravitational field
will illustrate how dependencies should be treated. 


Examples:

:ref:`ExternalGravityExample`

:ref:`CoolingExample`
