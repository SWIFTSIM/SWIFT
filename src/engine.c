/*******************************************************************************
 * This file is part of GadgetSMP.
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
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <omp.h>
#include <sched.h>

/* Local headers. */
#include "cycle.h"
#include "const.h"
#include "lock.h"
#include "task.h"
#include "part.h"
#include "cell.h"
#include "space.h"
#include "queue.h"
#include "engine.h"
#include "runner.h"
#include "runner_iact.h"

/* Error macro. */
#define error(s) { fprintf( stderr , "%s:%s:%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }

/* Convert cell location to ID. */
#define cell_getid( cdim , i , j , k ) ( (int)(k) + (cdim)[2]*( (int)(j) + (cdim)[1]*(int)(i) ) )


/**
 * @brief Prepare the #engine by re-building the cells and tasks.
 *
 * @param e The #engine to prepare.
 * @param force Flag to force re-building the cell and task structure.
 */
 
void engine_prepare ( struct engine *e , int force ) {

    int j, k, qid, changes;
    struct space *s = e->s;

    /* Rebuild the space. */
    changes = space_rebuild( e->s , force , 0 );
    // printf( "engine_prepare: space_rebuild with %i changes.\n" , changes );
    
    /* Has anything changed? */
    if ( changes ) {
    
        /* Rank the tasks in topological order. */
        engine_ranktasks( e );
    
        /* Clear the queues. */
        for ( k = 0 ; k < e->nr_queues ; k++ )
            e->queues[k].count = 0;
            
        /* Re-allocate the queue buffers? */
        for ( k = 0 ; k < e->nr_queues ; k++ )
            queue_init( &e->queues[k] , s->nr_tasks , s->tasks );
        
        /* Fill the queues (round-robin). */
        for ( k = 0 ; k < s->nr_tasks ; k++ ) {
            if ( s->tasks[ s->tasks_ind[k] ].type == task_type_none )
                continue;
            qid = k % e->nr_queues;
            e->queues[qid].tid[ e->queues[qid].count ] = s->tasks_ind[k];
            e->queues[qid].count += 1;
            }
            
        }

    /* Re-set the particle data. */
    #pragma omp parallel for
    for ( k = 0 ; k < s->nr_parts ; k++ ) {
        s->parts[k].wcount = 0.0f;
        s->parts[k].wcount_dh = 0.0f;
        s->parts[k].rho = 0.0f;
        s->parts[k].rho_dh = 0.0f;
        }
    
    /* Run throught the tasks and get all the waits right. */
    #pragma omp parallel for private(j)
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        for ( j = 0 ; j < s->tasks[k].nr_unlock_tasks ; j++ )
            __sync_add_and_fetch( &s->tasks[k].unlock_tasks[j]->wait , 1 );
        for ( j = 0 ; j < s->tasks[k].nr_unlock_cells ; j++ )
            __sync_add_and_fetch( &s->tasks[k].unlock_cells[j]->wait , 1 );
        }
    
    /* Re-set the queues.*/
    for ( k = 0 ; k < e->nr_queues ; k++ )
        e->queues[k].next = 0;
    
    }


/** 
 * @brief Sort the tasks in topological order over all queues.
 *
 * @param e The #engine.
 */
 
void engine_ranktasks ( struct engine *e ) {

    int i, j = 0, k, temp, left = 0, rank;
    struct task *t;
    struct space *s = e->s;
    int *tid = s->tasks_ind;

    /* Run throught the tasks and get all the waits right. */
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        for ( j = 0 ; j < s->tasks[k].nr_unlock_tasks ; j++ )
            s->tasks[k].unlock_tasks[j]->wait += 1;
        }
        
    /* Main loop. */
    for ( rank = 0 ; left < s->nr_tasks ; rank++ ) {
        
        /* Load the tids of tasks with no waits. */
        for ( k = left ; k < s->nr_tasks ; k++ )
            if ( s->tasks[ tid[k] ].wait == 0 ) {
                temp = tid[j]; tid[j] = tid[k]; tid[k] = temp;
                j += 1;
                }

        /* Traverse the task tree and add tasks with no weight. */
        for ( i = left ; i < j ; i++ ) {
            t = &s->tasks[ tid[i] ];
            t->rank = rank;
            s->tasks_ind[i] = t - s->tasks;
            /* printf( "engine_ranktasks: task %i of type %s has rank %i.\n" , i , 
                (t->type == task_type_self) ? "self" : (t->type == task_type_pair) ? "pair" : "sort" , rank ); */
            for ( k = 0 ; k < t->nr_unlock_tasks ; k++ )
                t->unlock_tasks[k]->wait -= 1;
            }
            
        /* The new left (no, not tony). */
        left = j;
            
        }
        
    }


/**
 * @brief Implements a barrier for the #runner threads.
 *
 * @param e The #engine.
 */
 
void engine_barrier( struct engine *e ) {

    /* First, get the barrier mutex. */
    if ( pthread_mutex_lock( &e->barrier_mutex ) != 0 )
        error( "Failed to get barrier mutex." );
        
    /* Wait for the barrier to close. */
    while ( e->barrier_count < 0 )
        if ( pthread_cond_wait( &e->barrier_cond , &e->barrier_mutex ) != 0 )
            error( "Eror waiting for barrier to close." );
        
    /* Once I'm in, increase the barrier count. */
    e->barrier_count += 1;
    
    /* If all threads are in, send a signal... */
    if ( e->barrier_count == e->nr_threads )
        if ( pthread_cond_broadcast( &e->barrier_cond ) != 0 )
            error( "Failed to broadcast barrier full condition." );
        
    /* Wait for barrier to be released. */
    while ( e->barrier_count > 0 )
        if ( pthread_cond_wait( &e->barrier_cond , &e->barrier_mutex ) != 0 )
            error( "Error waiting for barrier to be released." );
            
    /* Decrease the counter before leaving... */
    e->barrier_count += 1;
    
    /* If I'm the last one out, signal the condition again. */
    if ( e->barrier_count == 0 )
        if ( pthread_cond_broadcast( &e->barrier_cond ) != 0 )
            error( "Failed to broadcast empty barrier condition." );
            
    /* Last but not least, release the mutex. */
    if ( pthread_mutex_unlock( &e->barrier_mutex ) != 0 )
        error( "Failed to get unlock the barrier mutex." );

    }
    
    
/**
 * @brief Let the #engine loose to compute the forces.
 *
 * @param e The #engine.
 * @param sort_queues Flag to try to sort the queues topologically.
 */
 
void engine_run ( struct engine *e , int sort_queues , float dt_max ) {

    int k;

    /* Re-set the queues.*/
    if ( sort_queues ) {
        #pragma omp parallel for default(none), shared(e)
        for ( k = 0 ; k < e->nr_queues ; k++ ) {
            queue_sort( &e->queues[k] );
            e->queues[k].next = 0;
            }
        }
        
    /* Set the maximum dt. */
    e->dt_max = dt_max;
    e->s->dt_max = dt_max;
    
    /* Cry havoc and let loose the dogs of war. */
    e->barrier_count = -e->barrier_count;
    if ( pthread_cond_broadcast( &e->barrier_cond ) != 0 )
        error( "Failed to broadcast barrier open condition." );
        
    /* Sit back and wait for the runners to come home. */
    while ( e->barrier_count < e->nr_threads )
        if ( pthread_cond_wait( &e->barrier_cond , &e->barrier_mutex ) != 0 )
            error( "Error while waiting for barrier." );
    
    }
    
    
/**
 * @brief init an engine with the given number of threads, queues, and
 *      the given policy.
 *
 * @param e The #engine.
 * @param s The #space in which this #runner will run.
 * @param nr_threads The number of threads to spawn.
 * @param nr_queues The number of task queues to create.
 * @param policy The queueing policy to use.
 */
 
void engine_init ( struct engine *e , struct space *s , int nr_threads , int nr_queues , int policy ) {

    #if defined(HAVE_SETAFFINITY)
        cpu_set_t cpuset;
    #endif
    int k, qid, nrq;
    
    /* Store the values. */
    e->s = s;
    e->nr_threads = nr_threads;
    e->nr_queues = nr_queues;
    e->policy = policy;
    
    /* First of all, init the barrier and lock it. */
    if ( pthread_mutex_init( &e->barrier_mutex , NULL ) != 0 )
        error( "Failed to initialize barrier mutex." );
    if ( pthread_cond_init( &e->barrier_cond , NULL ) != 0 )
        error( "Failed to initialize barrier condition variable." );
    if ( pthread_mutex_lock( &e->barrier_mutex ) != 0 )
        error( "Failed to lock barrier mutex." );
    e->barrier_count = 0;
    
    /* Allocate the queues. */
    if ( posix_memalign( (void *)(&e->queues) , 64 , nr_queues * sizeof(struct queue) ) != 0 )
        error( "Failed to allocate queues." );
    bzero( e->queues , nr_queues * sizeof(struct queue) );
        
    /* Init the queues. */
    for ( k = 0 ; k < nr_queues ; k++ )
        queue_init( &e->queues[k] , s->nr_tasks , s->tasks );
        
    /* Rank the tasks in topological order. */
    engine_ranktasks( e );
    
    /* How many queues to fill initially? */
    for ( nrq = 0 , k = nr_queues ; k > 0 ; k = k / 2 )
        nrq += 1;
        
    /* Fill the queues (round-robin). */
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        if ( s->tasks[ s->tasks_ind[k] ].type == task_type_none )
            continue;
        // qid = 0;
        // qid = k % nrq;
        qid = k % nr_queues;
        e->queues[qid].tid[ e->queues[qid].count ] = s->tasks_ind[k];
        e->queues[qid].count += 1;
        }
        
    /* Sort the queues topologically. */
    for ( k = 0 ; k < nr_queues ; k++ )
        queue_sort( &e->queues[k] );
        
    /* Allocate and init the threads. */
    if ( ( e->runners = (struct runner *)malloc( sizeof(struct runner) * nr_threads ) ) == NULL )
        error( "Failed to allocate threads array." );
    for ( k = 0 ; k < nr_threads ; k++ ) {
        e->runners[k].id = k;
        e->runners[k].e = e;
        if ( pthread_create( &e->runners[k].thread , NULL , &runner_main , &e->runners[k] ) != 0 )
            error( "Failed to create runner thread." );
        #if defined(HAVE_SETAFFINITY)
            /* Set the cpu mask to zero | e->id. */
            CPU_ZERO( &cpuset );
            CPU_SET( e->runners[k].id , &cpuset );

            /* Apply this mask to the runner's pthread. */
            if ( pthread_setaffinity_np( e->runners[k].thread , sizeof(cpu_set_t) , &cpuset ) != 0 )
                error( "Failed to set thread affinity." );
        #endif
        }
        
    /* Wait for the runner threads to be in place. */
    while ( e->barrier_count != e->nr_threads )
        if ( pthread_cond_wait( &e->barrier_cond , &e->barrier_mutex ) != 0 )
            error( "Error while waiting for runner threads to get in place." );
    
    }
    
    
    
