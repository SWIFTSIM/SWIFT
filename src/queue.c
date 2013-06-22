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

/* Local headers. */
#include "cycle.h"
#include "lock.h"
#include "task.h"
#include "cell.h"
#include "queue.h"
#include "error.h"
#include "inline.h"

/* Define the timer macros. */
#ifdef TIMER_VERBOSE
    #ifndef TIMER
        #define TIMER
    #endif
#endif
#ifdef TIMER
    #define TIMER_TIC ticks tic = getticks();
    #define TIMER_TOC(t) timer_toc( t , tic )
    #define TIMER_TIC2 ticks tic2 = getticks();
    #define TIMER_TOC2(t) timer_toc( t , tic2 )
    INLINE static ticks timer_toc ( int t , ticks tic ) {
        ticks d = (getticks() - tic);
        __sync_add_and_fetch( &queue_timer[t] , d );
        return d;
        }
#else
    #define TIMER_TIC
    #define TIMER_TOC(t)
    #define TIMER_TIC2
    #define TIMER_TOC2(t)
#endif


/* Counter macros. */
#ifdef COUNTER
    #define COUNT(c) ( __sync_add_and_fetch( &queue_counter[ c ] , 1 ) )
#else
    #define COUNT(c)
#endif


/* The timers. */
ticks queue_timer[ queue_timer_count ];

/* The counters. */
int queue_counter[ queue_counter_count ];

        

/**
 * @brief Insert a used tasks into the given queue.
 *
 * @param q The #queue.
 * @param t The #task.
 */
 
void queue_insert ( struct queue *q , struct task *t ) {

    /* Lock the queue. */
    if ( lock_lock( &q->lock ) != 0 )
        error( "Failed to get queue lock." );
        
    /* Does the queue need to be grown? */
    if ( q->count == q->size ) {
        int *temp;
        q->size *= queue_sizegrow;
        if ( ( temp = (int *)malloc( sizeof(int) * q->size ) ) == NULL )
            error( "Failed to allocate new indices." );
        memcpy( temp , q->tid , sizeof(int) * q->count );
        free( q->tid );
        q->tid = temp;
        }
        
    /* Drop the task at the end of the queue. */
    q->tid[ q->count ] = ( t - q->tasks );
    q->count += 1;
    
    /* Shuffle up. */
    for ( int k = q->count - 1 ; k > 0 ; k /= 2 )
        if ( q->tasks[ q->tid[k] ].maxdepth > q->tasks[ q->tid[k/2] ].maxdepth ) {
            int temp = q->tid[k];
            q->tid[k] = q->tid[k/2];
            q->tid[k/2] = temp;
            }
        else
            break;
    
    /* Unlock the queue. */
    if ( lock_unlock( &q->lock ) != 0 )
        error( "Failed to unlock queue." );
    }


/** 
 * @brief Initialize the given queue.
 *
 * @param q The #queue.
 * @param size The maximum size of the queue.
 * @param tasks List of tasks to which the queue indices refer to.
 */
 
void queue_init ( struct queue *q , struct task *tasks ) {
    
    /* Allocate the task list if needed. */
    q->size = queue_sizeinit;
    if ( ( q->tid = (int *)malloc( sizeof(int) * q->size ) ) == NULL )
        error( "Failed to allocate queue tids." );
        
    /* Set the tasks pointer. */
    q->tasks = tasks;
        
    /* Init counters. */
    q->count = 0;
    
    /* Init the queue lock. */
    if ( lock_init( &q->lock ) != 0 )
        error( "Failed to init queue lock." );

    }


/**
 * @brief Get a task free of dependencies and conflicts.
 *
 * @param q The task #queue.
 * @param blocking Block until access to the queue is granted.
 * @param keep Remove the returned task from this queue.
 */
 
struct task *queue_gettask ( struct queue *q , int qid , int blocking ) {

    int k, qcount, *qtid, type;
    lock_type *qlock = &q->lock;
    struct task *qtasks, *res = NULL;
    struct cell *ci, *cj;
    TIMER_TIC
    
    /* If there are no tasks, leave immediately. */
    if ( q->count == 0 ) {
        TIMER_TOC(queue_timer_gettask);
        return NULL;
        }

    /* Main loop, while there are tasks... */
    while ( q->count > 0 ) {
    
        /* Grab the task lock. */
        if ( lock_lock( qlock ) != 0 )
            error( "Locking the qlock failed.\n" );
            
        /* Set some pointers we will use often. */
        qtid = q->tid;
        qtasks = q->tasks;
            
        /* Loop over the remaining task IDs. */
        qcount = q->count;
        for ( k = 0 ; k < qcount ; k++ ) {
        
            /* Put a finger on the task. */
            res = &qtasks[ qtid[k] ];
            ci = res->ci;
            cj = res->cj;
            type = res->type;
            
            /* Is this task blocked? */
            if ( res->wait )
                continue;
                
            /* Try to lock ci. */
            if ( type == task_type_self || 
                 type == task_type_sort || 
                 (type == task_type_sub && cj == NULL) ) {
                if ( cell_locktree( ci ) != 0 )
                    continue;
                }
            else if ( type == task_type_pair || (type == task_type_sub && cj != NULL) ) {
                if ( ci->hold || cj->hold || ci->wait || cj->wait )
                    continue;
                if ( cell_locktree( ci ) != 0 )
                    continue;
                if ( cell_locktree( cj ) != 0 ) {
                    cell_unlocktree( ci );
                    continue;
                    }
                }
            
            /* If we made it this far, we're safe. */
            break;
        
            } /* loop over the task IDs. */
            
        /* Did we get a task? */
        if ( k < qcount ) {
        
            /* Another one bites the dust. */
            q->count -= 1;
        
            /* Own the cells involved. */
            ci->super->owner = qid;
            if ( cj != NULL )
                cj->super->owner = qid;
                
            /* Swap this task with the last task and re-heap. */
            if ( k < q->count ) {
                qtid[ k ] = qtid[ q->count ];
                while ( 1 ) {
                    int i = 2*k;
                    if ( i >= q->count )
                        break;
                    if ( i+1 < q->count && qtasks[ qtid[i+1] ].maxdepth > qtasks[ qtid[i] ].maxdepth )
                        i += 1;
                    if ( qtasks[ qtid[i] ].maxdepth > qtasks[ qtid[k] ].maxdepth ) {
                        int temp = qtid[i];
                        qtid[i] = qtid[k];
                        qtid[k] = temp;
                        k = i;
                        }
                    else
                        break;
                    }
                }
                
            }
        else
            res = NULL;
    
        /* Release the task lock. */
        if ( lock_unlock( qlock ) != 0 )
            error( "Unlocking the qlock failed.\n" );
            
        /* Leave? */
        if ( res != NULL || !blocking )
            break;
    
        } /* while there are tasks. */
        
    /* Take the money and run. */
    TIMER_TOC(queue_timer_gettask);
    return res;

    }


