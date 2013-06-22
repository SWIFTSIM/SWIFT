/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 ******************************************************************************/

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

/* Local headers. */
#include "error.h"
#include "cycle.h"
#include "atomic.h"
#include "timers.h"
#include "const.h"
#include "vector.h"
#include "lock.h"
#include "task.h"
#include "part.h"
#include "debug.h"
#include "cell.h"
#include "space.h"
#include "queue.h"
#include "kernel.h"
#include "scheduler.h"


/**
 * @brief Mapping function to append a ghost task to each cell.
 *
 * Looks for the super cell, e.g. the highest-level cell above each
 * cell for which a pair is defined. All ghosts below this cell will
 * depend on the ghost of their parents (sounds spooky, but it isn't).
 *
 * A kick2-task is appended to each super cell.
 */

void scheduler_map_mkghosts ( struct cell *c , void *data ) {

    struct scheduler *s = (struct scheduler *)data;
    struct cell *finger;

    /* Find the super cell, i.e. the highest cell hierarchically above
       this one to still have at least one task associated with it. */
    c->super = c;
    for ( finger = c->parent ; finger != NULL ; finger = finger->parent )
        if ( finger->nr_tasks > 0 )
            c->super = finger;
            
    /* Make the ghost task */
    if ( c->super != c || c->nr_tasks > 0 )
        c->ghost = scheduler_addtask( s , task_type_ghost , task_subtype_none , 0 , 0 , c , NULL , 0 );

    /* Append a kick task if we are the active super cell. */
    if ( c->super == c && c->nr_tasks > 0 )
        c->kick2 = scheduler_addtask( s , task_type_kick2 , task_subtype_none , 0 , 0 , c , NULL , 0 );
    
    /* If we are not the super cell ourselves, make our ghost depend
       on our parent cell. */
    if ( c->super != c )
        task_addunlock( c->parent->ghost , c->ghost );
        
    }


/**
 * @brief Split tasks that may be too large.
 *
 * @param s The #scheduler we are working in.
 */
 
void scheduler_splittasks ( struct scheduler *s ) {

    int j, k, ind, sid, tid = 0, redo;
    struct cell *ci, *cj;
    double hi, hj, shift[3];
    struct task *t, *t_old;
    // float dt_step = s->dt_step;
    int pts[7][8] = { { -1 , 12 , 10 ,  9 ,  4 ,  3 ,  1 ,  0 } ,
                      { -1 , -1 , 11 , 10 ,  5 ,  4 ,  2 ,  1 } ,
                      { -1 , -1 , -1 , 12 ,  7 ,  6 ,  4 ,  3 } , 
                      { -1 , -1 , -1 , -1 ,  8 ,  7 ,  5 ,  4 } ,
                      { -1 , -1 , -1 , -1 , -1 , 12 , 10 ,  9 } ,
                      { -1 , -1 , -1 , -1 , -1 , -1 , 11 , 10 } ,
                      { -1 , -1 , -1 , -1 , -1 , -1 , -1 , 12 } };

    /* Loop through the tasks... */
    // #pragma omp parallel default(none) shared(s,tid,pts,space_subsize) private(ind,j,k,t,t_old,redo,ci,cj,hi,hj,sid,shift)
    {
    redo = 0; t_old = t = NULL;
    while ( 1 ) {
    
        /* Get a pointer on the task. */
        if ( redo ) {
            redo = 0;
            t = t_old;
            }
        else {
            if ( ( ind = atomic_inc( &tid ) ) < s->nr_tasks )
                t_old = t = &s->tasks[ s->tasks_ind[ ind ] ];
            else
                break;
            }
        
        /* Empty task? */
        if ( t->ci == NULL || ( t->type == task_type_pair && t->cj == NULL ) ) {
            t->type = task_type_none;
            t->skip = 1;
            continue;
            }
        
        /* Self-interaction? */
        if ( t->type == task_type_self ) {
        
            /* Get a handle on the cell involved. */
            ci = t->ci;
            
            /* Ingore this task? */
            /* if ( ci->dt_min > dt_step ) {
                t->skip = 1;
                continue;
                } */
            
            /* Is this cell even split? */
            if ( ci->split ) {
            
                /* Make a sub? */
                if ( scheduler_dosub && ci->count < space_subsize && ci->maxdepth - ci->depth < scheduler_maxsubdepth ) {

                    /* convert to a self-subtask. */
                    t->type = task_type_sub;

                    }

                /* Otherwise, make tasks explicitly. */
                else {

                    /* Take a step back (we're going to recycle the current task)... */
                    redo = 1;

                    /* Add the self taks. */
                    for ( k = 0 ; ci->progeny[k] == NULL ; k++ );
                    t->ci = ci->progeny[k];
                    for ( k += 1 ; k < 8 ; k++ )
                        if ( ci->progeny[k] != NULL )
                            scheduler_addtask( s , task_type_self , task_subtype_density , 0 , 0 , ci->progeny[k] , NULL , 0 );

                    /* Make a task for each pair of progeny. */
                    for ( j = 0 ; j < 8 ; j++ )
                        if ( ci->progeny[j] != NULL )
                            for ( k = j + 1 ; k < 8 ; k++ )
                                if ( ci->progeny[k] != NULL )
                                    scheduler_addtask( s , task_type_pair , task_subtype_density , pts[j][k] , 0 , ci->progeny[j] , ci->progeny[k] , 0 );
                    }

                }
        
            }
    
        /* Pair interaction? */
        else if ( t->type == task_type_pair ) {
            
            /* Get a handle on the cells involved. */
            ci = t->ci;
            cj = t->cj;
            hi = ci->dmin;
            hj = cj->dmin;

            /* Ingore this task? */
            /* if ( ci->dt_min > dt_step && cj->dt_min > dt_step ) {
                t->skip = 1;
                continue;
                } */
            
            /* Get the sort ID, use space_getsid and not t->flags
               to make sure we get ci and cj swapped if needed. */
            sid = space_getsid( s->space , &ci , &cj , shift );
                
            /* Should this task be split-up? */
            if ( ci->split && cj->split &&
                 ci->h_max*kernel_gamma*space_stretch < hi/2 &&
                 cj->h_max*kernel_gamma*space_stretch < hj/2 ) {
                 
                /* Replace by a single sub-task? */
                if ( scheduler_dosub &&
                     ci->count < space_subsize && cj->count < space_subsize &&
                     ci->maxdepth - ci->depth < scheduler_maxsubdepth && cj->maxdepth - cj->depth < scheduler_maxsubdepth &&
                     sid != 0 && sid != 2 && sid != 6 && sid != 8 ) {
                
                    /* Make this task a sub task. */
                    t->type = task_type_sub;

                    }
                    
                /* Otherwise, split it. */
                else {

                    /* Take a step back (we're going to recycle the current task)... */
                    redo = 1;

                    /* For each different sorting type... */
                    switch ( sid ) {

                        case 0: /* (  1 ,  1 ,  1 ) */
                            t->ci = ci->progeny[7]; t->cj = cj->progeny[0]; t->flags = 0;
                            break;

                        case 1: /* (  1 ,  1 ,  0 ) */
                            t->ci = ci->progeny[6]; t->cj = cj->progeny[0]; t->flags = 1; t->tight = 1;
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 1 , 0 , ci->progeny[7] , cj->progeny[1] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[6] , cj->progeny[1] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[7] , cj->progeny[0] , 1 );
                            break;

                        case 2: /* (  1 ,  1 , -1 ) */
                            t->ci = ci->progeny[6]; t->cj = cj->progeny[1]; t->flags = 2; t->tight = 1;
                            break;

                        case 3: /* (  1 ,  0 ,  1 ) */
                            t->ci = ci->progeny[5]; t->cj = cj->progeny[0]; t->flags = 3; t->tight = 1;
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 3 , 0 , ci->progeny[7] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[5] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 6 , 0 , ci->progeny[7] , cj->progeny[0] , 1 );
                            break;

                        case 4: /* (  1 ,  0 ,  0 ) */
                            t->ci = ci->progeny[4]; t->cj = cj->progeny[0]; t->flags = 4; t->tight = 1;
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 5 , 0 , ci->progeny[5] , cj->progeny[0] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[6] , cj->progeny[0] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[7] , cj->progeny[0] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 3 , 0 , ci->progeny[4] , cj->progeny[1] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 4 , 0 , ci->progeny[5] , cj->progeny[1] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 6 , 0 , ci->progeny[6] , cj->progeny[1] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[7] , cj->progeny[1] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 1 , 0 , ci->progeny[4] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[5] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 4 , 0 , ci->progeny[6] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 5 , 0 , ci->progeny[7] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[4] , cj->progeny[3] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 1 , 0 , ci->progeny[5] , cj->progeny[3] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 3 , 0 , ci->progeny[6] , cj->progeny[3] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 4 , 0 , ci->progeny[7] , cj->progeny[3] , 1 );
                            break;

                        case 5: /* (  1 ,  0 , -1 ) */
                            t->ci = ci->progeny[4]; t->cj = cj->progeny[1]; t->flags = 5; t->tight = 1;
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 5 , 0 , ci->progeny[6] , cj->progeny[3] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[4] , cj->progeny[3] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[6] , cj->progeny[1] , 1 );
                            break;

                        case 6: /* (  1 , -1 ,  1 ) */
                            t->ci = ci->progeny[5]; t->cj = cj->progeny[2]; t->flags = 6; t->tight = 1;
                            break;

                        case 7: /* (  1 , -1 ,  0 ) */
                            t->ci = ci->progeny[4]; t->cj = cj->progeny[3]; t->flags = 6; t->tight = 1;
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[5] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[4] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[5] , cj->progeny[3] , 1 );
                            break;

                        case 8: /* (  1 , -1 , -1 ) */
                            t->ci = ci->progeny[4]; t->cj = cj->progeny[3]; t->flags = 8; t->tight = 1;
                            break;

                        case 9: /* (  0 ,  1 ,  1 ) */
                            t->ci = ci->progeny[3]; t->cj = cj->progeny[0]; t->flags = 9; t->tight = 1;
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 9 , 0 , ci->progeny[7] , cj->progeny[4] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[3] , cj->progeny[4] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[7] , cj->progeny[0] , 1 );
                            break;

                        case 10: /* (  0 ,  1 ,  0 ) */
                            t->ci = ci->progeny[2]; t->cj = cj->progeny[0]; t->flags = 10; t->tight = 1;
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 11 , 0 , ci->progeny[3] , cj->progeny[0] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[6] , cj->progeny[0] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 6 , 0 , ci->progeny[7] , cj->progeny[0] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 9 , 0 , ci->progeny[2] , cj->progeny[1] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 10 , 0 , ci->progeny[3] , cj->progeny[1] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[6] , cj->progeny[1] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 7 , 0 , ci->progeny[7] , cj->progeny[1] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 1 , 0 , ci->progeny[2] , cj->progeny[4] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[3] , cj->progeny[4] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 10 , 0 , ci->progeny[6] , cj->progeny[4] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 11 , 0 , ci->progeny[7] , cj->progeny[4] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[2] , cj->progeny[5] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 1 , 0 , ci->progeny[3] , cj->progeny[5] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 9 , 0 , ci->progeny[6] , cj->progeny[5] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 10 , 0 , ci->progeny[7] , cj->progeny[5] , 1 );
                            break;

                        case 11: /* (  0 ,  1 , -1 ) */
                            t->ci = ci->progeny[2]; t->cj = cj->progeny[1]; t->flags = 11; t->tight = 1;
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 11 , 0 , ci->progeny[6] , cj->progeny[5] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[2] , cj->progeny[5] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 6 , 0 , ci->progeny[6] , cj->progeny[1] , 1 );
                            break;

                        case 12: /* (  0 ,  0 ,  1 ) */
                            t->ci = ci->progeny[1]; t->cj = cj->progeny[0]; t->flags = 12; t->tight = 1;
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 11 , 0 , ci->progeny[3] , cj->progeny[0] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 5 , 0 , ci->progeny[5] , cj->progeny[0] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 2 , 0 , ci->progeny[7] , cj->progeny[0] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 9 , 0 , ci->progeny[1] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 12 , 0 , ci->progeny[3] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 8 , 0 , ci->progeny[5] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 5 , 0 , ci->progeny[7] , cj->progeny[2] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 3 , 0 , ci->progeny[1] , cj->progeny[4] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 6 , 0 , ci->progeny[3] , cj->progeny[4] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 12 , 0 , ci->progeny[5] , cj->progeny[4] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 11 , 0 , ci->progeny[7] , cj->progeny[4] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 0 , 0 , ci->progeny[1] , cj->progeny[6] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 3 , 0 , ci->progeny[3] , cj->progeny[6] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 9 , 0 , ci->progeny[5] , cj->progeny[6] , 1 );
                            t = scheduler_addtask( s , task_type_pair , t->subtype , 12 , 0 , ci->progeny[7] , cj->progeny[6] , 1 );
                            break;

                        }
                        
                    }

                } /* split this task? */
                
            /* Otherwise, if not spilt, stitch-up the sorting. */
            else {
            
                /* Create the sort for ci. */
                // lock_lock( &ci->lock );
                if ( ci->sorts == NULL )
                    ci->sorts = scheduler_addtask( s , task_type_sort , 0 , 1 << sid , 0 , ci , NULL , 0 );
                else
                    ci->sorts->flags |= (1 << sid);
                // lock_unlock_blind( &ci->lock );
                task_addunlock( ci->sorts , t );
                
                /* Create the sort for cj. */
                // lock_lock( &cj->lock );
                if ( cj->sorts == NULL )
                    cj->sorts = scheduler_addtask( s , task_type_sort , 0 , 1 << sid , 0 , cj , NULL , 0 );
                else
                    cj->sorts->flags |= (1 << sid);
                // lock_unlock_blind( &cj->lock );
                task_addunlock( cj->sorts , t );
                
                }
                
            } /* pair interaction? */
    
        } /* loop over all tasks. */
        
        }
        
    }
    
    
/**
 * @brief Add a #task to the #scheduler.
 *
 * @param s The #scheduler we are working in.
 * @param type The type of the task.
 * @param subtype The sub-type of the task.
 * @param flags The flags of the task.
 * @param wait 
 * @param ci The first cell to interact.
 * @param cj The second cell to interact.
 * @param tight
 */
 
struct task *scheduler_addtask ( struct scheduler *s , int type , int subtype , int flags , int wait , struct cell *ci , struct cell *cj , int tight ) {

    int ind;
    struct task *t;
    
    /* Get the next free task. */
    ind = atomic_inc( &s->tasks_next );
    t = &s->tasks[ ind ];
    
    /* Copy the data. */
    t->type = type;
    t->subtype = subtype;
    t->flags = flags;
    t->wait = wait;
    t->ci = ci;
    t->cj = cj;
    t->skip = 0;
    t->tight = tight;
    t->nr_unlock_tasks = 0;
    
    /* Init the lock. */
    lock_init( &t->lock );
    
    /* Add an index for it. */
    // lock_lock( &s->lock );
    s->tasks_ind[ atomic_inc( &s->nr_tasks ) ] = ind;
    // lock_unlock_blind( &s->lock );
    
    /* Return a pointer to the new task. */
    return t;

    }



/** 
 * @brief Sort the tasks in topological order over all queues.
 *
 * @param s The #scheduler.
 */
 
void scheduler_ranktasks ( struct scheduler *s ) {

    int i, j = 0, k, temp, left = 0, rank;
    struct task *t, *tasks = s->tasks;
    int *tid = s->tasks_ind, nr_tasks = s->nr_tasks;
    
    /* Run throught the tasks and get all the waits right. */
    for ( i = 0 , k = 0 ; k < nr_tasks ; k++ ) {
        tid[k] = k;
        for ( j = 0 ; j < tasks[k].nr_unlock_tasks ; j++ )
            tasks[k].unlock_tasks[j]->wait += 1;
        }
        
    /* Main loop. */
    for ( j = 0 , rank = 0 ; left < nr_tasks ; rank++ ) {
        
        /* Load the tids of tasks with no waits. */
        for ( k = left ; k < nr_tasks ; k++ )
            if ( tasks[ tid[k] ].wait == 0 ) {
                temp = tid[j]; tid[j] = tid[k]; tid[k] = temp;
                j += 1;
                }
                
        /* Did we get anything? */
        if ( j == left )
            error( "Unsatisfiable task dependencies detected." );

        /* Unlock the next layer of tasks. */
        for ( i = left ; i < j ; i++ ) {
            t = &tasks[ tid[i] ];
            t->rank = rank;
            tid[i] = t - tasks;
            if ( tid[i] >= nr_tasks )
                error( "Task index overshoot." );
            /* printf( "scheduler_ranktasks: task %i of type %s has rank %i.\n" , i , 
                (t->type == task_type_self) ? "self" : (t->type == task_type_pair) ? "pair" : "sort" , rank ); */
            for ( k = 0 ; k < t->nr_unlock_tasks ; k++ )
                t->unlock_tasks[k]->wait -= 1;
            }
            
        /* The new left (no, not tony). */
        left = j;
            
        }
        
    /* Run throught the tasks backwards and set their weight and maxdepth. */
    for ( k = nr_tasks-1 ; k >= 0 ; k-- ) {
        t = &tasks[ tid[k] ];
        t->maxdepth = 0;
        switch ( t->type ) {
            case task_type_none:
                t->weight = 0;
                break;
            case task_type_sort:
                t->weight = t->ci->count * ( sizeof(int)*8 - __builtin_clz( t->ci->count ) );
                break;
            case task_type_self:
                t->weight = t->ci->count * t->ci->count;
                break;
            case task_type_pair:
                t->weight = t->ci->count * t->cj->count;
                break;
            case task_type_sub:
                if ( t->cj != NULL )
                    t->weight = t->ci->count * t->cj->count;
                else
                    t->weight = t->ci->count * t->ci->count;
                break;
            case task_type_ghost:
                if ( t->ci == t->ci->super )
                    t->weight = t->ci->count;
                else
                    t->weight = 0;
                break;
            case task_type_kick2:
                t->weight = t->ci->count;
                break;
            }
        for ( j = 0 ; j < t->nr_unlock_tasks ; j++ ) {
            t->weight += t->unlock_tasks[j]->weight;
            if ( t->unlock_tasks[j]->maxdepth > t->maxdepth )
                t->maxdepth = t->unlock_tasks[j]->maxdepth;
            }
        }
        
    }


/**
 * @brief (Re)allocate the task arrays.
 *
 * @param s The #scheduler.
 * @param size The maximum number of tasks in the #scheduler.
 */
 
void scheduler_reset ( struct scheduler *s , int size ) {

    int k;

    /* Do we need to re-allocate? */
    if ( size > s->size ) {

        /* Free exising task lists if necessary. */
        if ( s->tasks != NULL )
            free( s->tasks );
        if ( s->tasks_ind != NULL )
            free( s->tasks_ind );

        /* Allocate the new lists. */
        if ( ( s->tasks = (struct task *)malloc( sizeof(struct task) * size ) ) == NULL ||
             ( s->tasks_ind = (int *)malloc( sizeof(int) * size ) ) == NULL )
            error( "Failed to allocate task lists." );
            
        }
        
    /* Reset the counters. */
    s->size = size;
    s->nr_tasks = 0;
    s->tasks_next = 0;
    s->waiting = 0;
    
    /* Set the task pointers in the queues. */
    for ( k = 0 ; k < s->nr_queues ; k++ )
        s->queues[k].tasks = s->tasks;

    }


/**
 * @brief Start the scheduler, i.e. fill the queues with ready tasks.
 *
 * @param s The #scheduler.
 */
 
void scheduler_start ( struct scheduler *s ) {

    int k, j;
    struct task *t;
    
    /* Run through the tasks and get all the waits right. */
    // #pragma omp parallel for schedule(static) private(t,j)
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        t = &s->tasks[k];
        if ( !t->skip )
            for ( j = 0 ; j < t->nr_unlock_tasks ; j++ )
                atomic_inc( &t->unlock_tasks[j]->wait );
        }
        
    /* Loop over the tasks and enqueue whoever is ready. */
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        t = &s->tasks[k];
        if ( !t->skip && t->wait == 0 )
            scheduler_enqueue( s , t );
        }
        
    }


/**
 * @brief Put a task on one of the queues.
 *
 * @param s The #scheduler.
 * @param t The #task.
 */
 
void scheduler_enqueue ( struct scheduler *s , struct task *t ) {

    int k, qid = -1;
    
    /* Ignore skipped tasks. */
    if ( t->skip )
        return;
        
    /* Find the previous owner for each task type. */
    switch ( t->type ) {
        case task_type_self:
        case task_type_sort:
        case task_type_ghost:
        case task_type_kick2:
            qid = t->ci->super->owner;
            break;
        case task_type_pair:
        case task_type_sub:
            qid = t->ci->super->owner;
            if ( t->cj != NULL &&
                 ( qid < 0 || s->queues[qid].count > s->queues[t->cj->super->owner].count ) )
                qid = t->cj->super->owner;
            break;
        }
        
    /* If no previous owner, find the shortest queue. */
    if ( qid < 0 )
        for ( qid = 0 , k = 1 ; k < s->nr_queues ; k++ )
            if ( s->queues[k].count < s->queues[qid].count )
                qid = k;
                
    /* Increase the waiting counter. */
    atomic_inc( &s->waiting );
            
    /* Insert the task into that queue. */
    queue_insert( &s->queues[qid] , t );
        
    }


/**
 * @brief Take care of a tasks dependencies.
 *
 * @param s The #scheduler.
 * @param t The finished #task.
 */
 
void scheduler_done ( struct scheduler *s , struct task *t ) {

    int k;
    struct task *t2;

    /* Release whatever locks this task held. */
    switch ( t->type ) {
        case task_type_self:
        case task_type_sort:
            cell_unlocktree( t->ci );
            break;
        case task_type_pair:
        case task_type_sub:
            cell_unlocktree( t->ci );
            if ( t->cj != NULL )
                cell_unlocktree( t->cj );
            break;
        }
        
    /* Loop through the dependencies and add them to a queue if
       they are ready. */
    for ( k = 0 ; k < t->nr_unlock_tasks ; k++ ) {
        t2 = t->unlock_tasks[k];
        if ( atomic_dec( &t2->wait ) == 1 && !t2->skip )
            scheduler_enqueue( s , t2 );
        }
        
    /* Task definitely done. */
    pthread_mutex_lock( &s->sleep_mutex );
    atomic_dec( &s->waiting );
    pthread_cond_broadcast( &s->sleep_cond );
    pthread_mutex_unlock( &s->sleep_mutex );

    }


/**
 * @brief Get a task, preferably from the given queue.
 *
 * @param s The #scheduler.
 * @param qid The ID of the prefered #queue.
 *
 * @return A pointer to a #task or @c NULL if there are no available tasks.
 */
 
struct task *scheduler_gettask ( struct scheduler *s , int qid ) {

    struct task *res;
    int k, nr_queues = s->nr_queues;

    /* Loop as long as there are tasks... */
    while ( s->waiting > 0 ) {
        
        /* Try more than once before sleeping. */
        for ( int tries = 0 ; tries < scheduler_flag_maxsteal ; tries++ ) {
        
            /* Try to get a task from the suggested queue. */
            if ( ( res = queue_gettask( &s->queues[qid] , qid , 0 ) ) != NULL )
                return res;

            /* If unsucessful, try stealing from the largest queue. */
            if ( s->flags & scheduler_flag_steal ){
                int qids[ nr_queues ];
                for ( k = 0 ; k < nr_queues ; k++ )
                    qids[k] = k;
                for ( k = 0 ; k < nr_queues ; k++ ) {
                    int j = k + ( rand() % (nr_queues - k) );
                    int temp = qids[j];
                    qids[j] = qids[k];
                    if ( temp == qid )
                        continue;
                    if ( s->queues[temp].count > 0 && ( res = queue_gettask( &s->queues[temp] , qid , 0 ) ) != NULL )
                        return res;
                    }
                }
                
            }
            
        /* If we failed, take a short nap. */
        pthread_mutex_lock( &s->sleep_mutex );
        if ( s->waiting > 0 )
            pthread_cond_wait( &s->sleep_cond , &s->sleep_mutex );
        pthread_mutex_unlock( &s->sleep_mutex );
        
        }
        
    /* No milk today. */
    return NULL;

    }


/**
 * @brief Initialize the #scheduler.
 *
 * @param s The #scheduler.
 * @param nr_queues The number of queues in this scheduler.
 * @param flags The #scheduler flags.
 */
 
void scheduler_init ( struct scheduler *s , struct space *space , int nr_queues , unsigned int flags ) {
    
    int k;
    
    /* Init the lock. */
    lock_init( &s->lock );

    /* Allocate the queues. */
    if ( ( s->queues = (struct queue *)malloc( sizeof(struct queue) * nr_queues ) ) == NULL )
        error( "Failed to allocate queues." );
        
    /* Initialize each queue. */
    for ( k = 0 ; k < nr_queues ; k++ )
        queue_init( &s->queues[k] , NULL );
        
    /* Init the sleep mutex and cond. */
    if ( pthread_cond_init( &s->sleep_cond , NULL ) != 0 ||
         pthread_mutex_init( &s->sleep_mutex , NULL ) != 0 )
        error( "Failed to initialize sleep barrier." );
        
    /* Set the scheduler variables. */
    s->nr_queues = nr_queues;
    s->flags = flags;
    s->space = space;
    
    /* Init other values. */
    s->tasks = NULL;
    s->tasks_ind = NULL;
    s->waiting = 0;
    s->size = 0;
    s->nr_tasks = 0;
    s->tasks_next = 0;

    }

