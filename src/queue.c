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

/* MPI headers. */
#ifdef WITH_MPI
    #include <mpi.h>
#endif

/* Local headers. */
#include "const.h"
#include "cycle.h"
#include "lock.h"
#include "task.h"
#include "timers.h"
#include "space.h"
#include "part.h"
#include "multipole.h"
#include "cell.h"
#include "queue.h"
#include "error.h"
#include "inline.h"



/* Counter macros. */
#ifdef COUNTER
    #define COUNT(c) ( __sync_add_and_fetch( &queue_counter[ c ] , 1 ) )
#else
    #define COUNT(c)
#endif


/* The counters. */
int queue_counter[ queue_counter_count ];

        

/**
 * @brief Insert a used tasks into the given queue.
 *
 * @param q The #queue.
 * @param t The #task.
 */
 
void queue_insert ( struct queue *q , struct task *t ) {

    int k, *tid;
    struct task *tasks;

    /* Lock the queue. */
    if ( lock_lock( &q->lock ) != 0 )
        error( "Failed to get queue lock." );
        
    tid = q->tid;
    tasks = q->tasks;

    /* Does the queue need to be grown? */
    if ( q->count == q->size ) {
        int *temp;
        q->size *= queue_sizegrow;
        if ( ( temp = (int *)malloc( sizeof(int) * q->size ) ) == NULL )
            error( "Failed to allocate new indices." );
        memcpy( temp , tid , sizeof(int) * q->count );
        free( tid );
        q->tid = tid = temp;
        }
        
    /* Drop the task at the end of the queue. */
    tid[ q->count ] = ( t - tasks );
    q->count += 1;
    
    /* Shuffle up. */
    for ( k = q->count - 1 ; k > 0 ; k = (k-1)/2 )
        if ( tasks[ tid[k] ].weight > tasks[ tid[(k-1)/2] ].weight ) {
            int temp = tid[k];
            tid[k] = tid[(k-1)/2];
            tid[(k-1)/2] = temp;
            }
        else
            break;
            
    /* Check the queue's consistency. */
    /* for ( k = 1 ; k < q->count ; k++ )
        if ( tasks[ tid[(k-1)/2] ].weight < tasks[ tid[k] ].weight )
            error( "Queue heap is disordered." ); */
    
    /* Unlock the queue. */
    if ( lock_unlock( &q->lock ) != 0 )
        error( "Failed to unlock queue." );
    }


/** 
 * @brief Initialize the given queue.
 *
 * @param q The #queue.
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
 * @param super The super-cell tat might conflict with the #queue
 * @param blocking Block until access to the queue is granted.
 */
 
struct task *queue_gettask ( struct queue *q , struct cell *super , int blocking ) {

    int k, qcount, *qtid, gotcha;
    lock_type *qlock = &q->lock;
    struct task *qtasks, *res = NULL;
    
    /* If there are no tasks, leave immediately. */
    if ( q->count == 0 )
        return NULL;

    /* Grab the task lock. */
    if ( blocking ) {
        if ( lock_lock( qlock ) != 0 )
            error( "Locking the qlock failed.\n" );
        }
    else {
        if ( lock_trylock( qlock ) != 0 )
            return NULL;
        }

    /* Set some pointers we will use often. */
    qtid = q->tid;
    qtasks = q->tasks;
    qcount = q->count;
    gotcha = 0;

    /* Loop over the task IDs looking for tasks with the same super-cell. */
    if ( super != NULL ) {
        for ( k = 0 ; k < qcount && k < queue_maxsuper ; k++ ) {

            /* Put a finger on the task. */
            res = &qtasks[ qtid[k] ];

            /* Try to lock the task and exit if successful. */
            if ( ( res->ci->super == super || ( res->cj != NULL && res->cj->super == super ) ) &&
                 task_lock( res ) ) {
                gotcha = 1;
                break;
                }

            } /* loop over the task IDs. */
        }

    /* Loop over the task IDs again if nothing was found, take anything. */
    if ( !gotcha ) {
        for ( k = 0 ; k < qcount ; k++ ) {

            /* Put a finger on the task. */
            res = &qtasks[ qtid[k] ];

            /* Try to lock the task and exit if successful. */
            if ( task_lock( res ) )
                break;

            } /* loop over the task IDs. */
        }

    /* Did we get a task? */
    if ( k < qcount ) {

        /* Another one bites the dust. */
        qcount = q->count -= 1;

        /* Swap this task with the last task and re-heap. */
        if ( k < qcount ) {
            qtid[ k ] = qtid[ qcount ];
            int w = qtasks[ qtid[k] ].weight;
            while ( k > 0 && w > qtasks[ qtid[(k-1)/2] ].weight ) {
                int temp = q->tid[k];
                q->tid[k] = q->tid[(k-1)/2];
                q->tid[(k-1)/2] = temp;
                k = (k-1)/2;
                }
            int i;
            while ( ( i = 2*k+1 ) < qcount ) {
                if ( i+1 < qcount && qtasks[ qtid[i+1] ].weight > qtasks[ qtid[i] ].weight )
                    i += 1;
                if ( qtasks[ qtid[i] ].weight > w ) {
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

    /* Check the queue's consistency. */
    /* for ( k = 1 ; k < q->count ; k++ )
        if ( qtasks[ qtid[(k-1)/2] ].weight < qtasks[ qtid[k] ].weight )
            error( "Queue heap is disordered." ); */
    
    /* Release the task lock. */
    if ( lock_unlock( qlock ) != 0 )
        error( "Unlocking the qlock failed.\n" );
            
    /* Take the money and run. */
    return res;

    }


