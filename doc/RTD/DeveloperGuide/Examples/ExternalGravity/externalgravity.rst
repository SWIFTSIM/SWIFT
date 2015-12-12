.. _ExternalGravityExample:

External Gravity Task Example
----------------------------------

An example of how to implement an external gravity task in SWIFT
=====================================================================

An external gravitational field can be imposed in SWIFT to mimic self-gravity. This is done by assigning
a gravitational force that falls as $1/ r^2$ (mathjax support to be included).

In order to do this we update the files as described in :ref:`NewTask`. For the specific case of adding an 
external graviational field the additions are as follows:


--------------
**task.h**
--------------

Code (snapshot Nov 2015)::

     /* The different task types. */
     enum task_types {
        task_type_none = 0,
     	task_type_sort,
     	task_type_self,
     	task_type_pair,
     	task_type_sub,
     	task_type_ghost,
     	task_type_kick1,
     	task_type_kick2,
     	task_type_send,
     	task_type_recv,
     	task_type_link,
     	task_type_grav_pp,
     	task_type_grav_mm,
     	task_type_grav_up,
     	task_type_grav_down,
     	**task_type_grav_external,**
     	task_type_psort,
     	task_type_split_cell,
     	task_type_count
     };

Task of type - task_type_grav_external - added to list of tasks.

--------------
**task.c**
--------------

Code (snapshot Nov 2015)::

       /* Task type names. */
       const char *taskID_names[task_type_count] = {
           "none",  "sort",    "self",    "pair",    "sub",
    	   "ghost", "kick1",   "kick2",   "send",    "recv",
    	   "link",  "grav_pp", "grav_mm", "grav_up", "grav_down", "grav_external",
    	   "psort", "split_cell"
        };

Task added to list of task names (used only for debugging purposed).


--------------
**cell.h**
--------------

Code (snapshot Nov 2015)::

     /* The ghost task to link density to interactions. */
        struct task *ghost, *kick1, *kick2, *grav_external;

Struture of type "task" declared (or pointer to a task at least). 



--------------
**timers.h**
--------------

Code (snapshot Nov 2015)::

    /* The timers themselves. */
    enum {
      timer_none = 0,
      timer_prepare,
      timer_kick1,
      timer_kick2,
      timer_dosort,
      timer_doself_density,
      timer_doself_force,
      timer_doself_grav,
      timer_dopair_density,
      timer_dopair_force,
      timer_dopair_grav,
      timer_dosub_density,
      timer_dosub_force,
      timer_dosub_grav,
      timer_dopair_subset,
      timer_doghost,
      timer_dograv_external,
      timer_gettask,
      timer_qget,
      timer_qsteal,
      timer_runners,
      timer_step,
      timer_count,
      };

The timer list is updated to include a timer task. 


--------------
**engine.c**
--------------

Code (snapshot Nov 2015)::

    void engine_mkghosts(struct engine *e, struct cell *c, struct cell *super) {

        int k;
  	struct scheduler *s = &e->sched;

  	/* Am I the super-cell? */
  	if (super == NULL && c->nr_tasks > 0) {

           /* Remember me. */
    	   super = c;

    	   /* Local tasks only... */
           if (c->nodeID == e->nodeID) {

               /* Generate the external gravity task*/
      	       c->grav_external = scheduler_addtask(s, task_type_grav_external, task_subtype_none, 0, 0,
                                   c, NULL, 0);

 	       /* Enforce gravity calculated before kick 2 */
      	       scheduler_addunlock(s, c->grav_external, c->kick2);
    	       }
	   }
     }


The first function call adds the task to the scheduler. The second function call takes care of the dependency 
involved in imposing an external gravitational field. These two functions are worth considering due to their 
obvious importance. 



The function prototype for the addtask function is (**found in scheduler.c**)::

        struct task *scheduler_addtask(struct scheduler *s, int type, int subtype,
                               int flags, int wait, struct cell *ci,
                               struct cell *cj, int tight) {

This function adds a task to the scheduler. In the call to this function in engine.c we used the actual 
parameters **s** for the scheduler, **task_type_grav_external** for the (task) type, task_subtype_none for 
the (task) subtype, zeros for the flags and wait parameters, **c** for the pointer to our cell, NULL for the
cell we interact with since there is none and 0 for the tight parameter. 

The function prototype for the addunlock function is(**found in scheduler.c**)::

        void scheduler_addunlock(struct scheduler *s, struct task *ta,
                         struct task *tb) {

This function signals when the unlock a certain task. In our case we use the external gravity task to unlock the 
kick2 task - i.e. kick2 depends on external gravity. So when calling the addunlock function the 
order is the **ta** task should be the task to unlock and **tb** should the task that does the unlocking.


--------------
**runner.c**
--------------

In runner.c the implementation of the external gravity task is taken care of. The function prototype is::

        void runner_dograv_external(struct runner *r, struct cell *c) {

The function takes a pointer to a runner struct and a pointer to the cell struct. The entire function call is::



       void runner_dograv_external(struct runner *r, struct cell *c) {

            struct part *p, *parts = c->parts;
	    float rinv;
  	    int i, ic, k, count = c->count;
  	    float dt_step = r->e->dt_step;
  	    TIMER_TIC

  	    /* Recurse? */
  	    if (c->split) {
    	       for (k = 0; k < 8; k++)
      	           if (c->progeny[k] != NULL) runner_dograv_external(r, c->progeny[k]);
    		   return;
            }		   

  	    /* Loop over the parts in this cell. */
  	    for (i = 0; i < count; i++) {

	        /* Get a direct pointer on the part. */
	 	p = &parts[i];
	 
	        /* Is this part within the time step? */
	 	if (p->dt <= dt_step) {
		   rinv = 1 / sqrtf((p->x[0])*(p->x[0]) + (p->x[1])*(p->x[1]) + (p->x[2])*(p->x[2]));
		   for(ic=0;ic<3;ic++){
		       p->grav_accel[ic] = - const_G * (p->x[ic]) * rinv * rinv * rinv;
		   }
	        }
            }
            TIMER_TOC(timer_dograv_external);
        }


The key component of this function is the calculation of **rinv** and then the imposition of the 
**grav_accel** to this particle. **rinv** is calculated assuming the centre of the gravitational 
potential lies at the origin. The acceleration of each particle then is calculated by multiplying
the graviational constant by the component of the position along one axis divided by R^3. The 
gravitational acceleration is then added to the total particle acceleration **a**.



.. toctree::
   :maxdepth: 1
