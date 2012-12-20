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
#include "lock.h"
#include "task.h"
#include "cell.h"
#include "queue.h"

/* Error macro. */
#define error(s) { fprintf( stderr , "%s:%s:%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }

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
    #ifndef INLINE
    # if __GNUC__ && !__GNUC_STDC_INLINE__
    #  define INLINE extern inline
    # else
    #  define INLINE inline
    # endif
    #endif
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

    int k;

    /* Lock the queue. */
    if ( lock_lock( &q->lock ) != 0 )
        error( "Failed to get queue lock." );
        
    /* Bubble-up the tasks. */
    for ( k = q->count ; k > q->next ; k-- )
        q->tid[k] = q->tid[k-1];
    q->tid[ q->next ] = t - q->tasks;
    q->count += 1;
    q->next += 1;
    
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
 
void queue_init ( struct queue *q , int size , struct task *tasks ) {
    
    /* Allocate the task list if needed. */
    if ( q->tid == NULL || q->size < size ) {
        if ( q->tid != NULL )
            free( q->tid );
        q->size = size;
        if ( ( q->tid = (int *)malloc( sizeof(int) * size ) ) == NULL )
            error( "Failed to allocate queue tids." );
        }
        
    /* Set the tasks pointer. */
    q->tasks = tasks;
        
    /* Init counters. */
    q->count = 0;
    q->next = 0;
    
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
 
struct task *queue_gettask_old ( struct queue *q , int blocking , int keep ) {

    int k, tid = -1, qcount, *qtid = q->tid;
    lock_type *qlock = &q->lock;
    struct task *qtasks = q->tasks, *res = NULL;
    TIMER_TIC
    
    /* If there are no tasks, leave immediately. */
    if ( q->next >= q->count ) {
        TIMER_TOC(queue_timer_gettask);
        return NULL;
        }

    /* Main loop, while there are tasks... */
    while ( q->next < q->count ) {
    
        /* Grab the task lock. */
        // if ( blocking ) {
            if ( lock_lock( qlock ) != 0 )
                error( "Locking the task_lock failed.\n" );
        //     }
        // else {
        //     if ( lock_trylock( qlock ) != 0 )
        //         break;
        //     }
            
        /* Loop over the remaining task IDs. */
        qcount = q->count;
        for ( k = q->next ; k < qcount ; k++ ) {
        
            /* Put a finger on the task. */
            res = &qtasks[ qtid[k] ];
            
            /* Is this task blocked? */
            if ( res->wait )
                continue;
            
            /* Different criteria for different types. */
            if ( res->type == task_type_self || (res->type == task_type_sub && res->cj == NULL) ) {
                if ( res->ci->hold || cell_locktree( res->ci ) != 0 )
                    continue;
                }
            else if ( res->type == task_type_pair || (res->type == task_type_sub && res->cj != NULL) ) {
                if ( res->ci->hold || res->cj->hold || res->ci->wait || res->cj->wait )
                    continue;
                if ( cell_locktree( res->ci ) != 0 )
                    continue;
                if ( cell_locktree( res->cj ) != 0 ) {
                    cell_unlocktree( res->ci );
                    continue;
                    }
                }
            
            /* If we made it this far, we're safe. */
            break;
        
            } /* loop over the task IDs. */
            
        /* Did we get a task? */
        if ( k < qcount ) {
        
            /* Do we need to swap? */
            if ( k != q->next )
                COUNT(queue_counter_swap);
        
            /* get the task ID. */
            tid = qtid[k];
        
            /* Remove the task? */
            if ( keep ) {
            
                /* Bubble-up. */
                q->count = qcount - 1;
                for ( ; k < qcount - 1 ; k++ )
                    qtid[k] = qtid[k+1];
            
                }
                
            /* No, leave it in the queue. */
            else {
            
                TIMER_TIC2

                /* Bubble-down the task. */
                while ( k > q->next ) {
                    qtid[ k ] = qtid[ k-1 ];
                    k -= 1;
                    }
                qtid[ q->next ] = tid;
                
                /* up the counter. */
                q->next += 1;
                
                TIMER_TOC2(queue_timer_bubble);
            
                }
            
            }
    
        /* Release the task lock. */
        if ( lock_unlock( qlock ) != 0 )
            error( "Unlocking the task_lock failed.\n" );
            
        /* Leave? */
        if ( tid >= 0 ) {
            TIMER_TOC(queue_timer_gettask);
            return &qtasks[tid];
            }
        else if ( !blocking )
            break;
    
        } /* while there are tasks. */
        
    /* No beef. */
    TIMER_TOC(queue_timer_gettask);
    return NULL;

    }


struct task *queue_gettask ( struct queue *q , int rid , int blocking , int keep ) {

    int k, tid = -1, qcount, *qtid = q->tid, hits;
    lock_type *qlock = &q->lock;
    struct task *qtasks = q->tasks, *res = NULL;
    struct cell *ci_best = NULL, *cj_best = NULL;
    int ind_best, score_best = -1, score;
    TIMER_TIC
    
    /* If there are no tasks, leave immediately. */
    if ( q->next >= q->count ) {
        TIMER_TOC(queue_timer_gettask);
        return NULL;
        }

    /* Main loop, while there are tasks... */
    while ( q->next < q->count ) {
    
        /* Grab the task lock. */
        // if ( blocking ) {
            if ( lock_lock( qlock ) != 0 )
                error( "Locking the task_lock failed.\n" );
        //     }
        // else {
        //     if ( lock_trylock( qlock ) != 0 )
        //         break;
        //     }
            
        /* Loop over the remaining task IDs. */
        qcount = q->count; ind_best = -1; hits = 0;
        for ( k = q->next ; k < qcount && hits < queue_maxhits ; k++ ) {
        
            /* Put a finger on the task. */
            res = &qtasks[ qtid[k] ];
            
            /* Is this task blocked? */
            if ( res->wait )
                continue;
                
            /* Get the score for this task. */
            if ( res->type == task_type_self || res->type == task_type_ghost || res->type == task_type_sort || ( res->type == task_type_sub && res->cj == NULL ) )
                score = ( res->ci->owner == rid );
            else
                score = ( res->ci->owner == rid ) + ( res->cj->owner == rid );
            if ( score <= score_best )
                continue;
                
            /* Try to lock ci. */
            if ( res->type == task_type_self || (res->type == task_type_sub && res->cj == NULL) ) {
                if ( res->ci != ci_best && res->ci != cj_best && cell_locktree( res->ci ) != 0 )
                    continue;
                }
            else if ( res->type == task_type_pair || (res->type == task_type_sub && res->cj != NULL) ) {
                if ( res->ci->hold || res->cj->hold || res->ci->wait || res->cj->wait )
                    continue;
                if ( res->ci != ci_best && res->ci != cj_best && cell_locktree( res->ci ) != 0 )
                    continue;
                if ( res->cj != ci_best && res->cj != cj_best && cell_locktree( res->cj ) != 0 ) {
                    if ( res->ci != ci_best && res->ci != cj_best )
                        cell_unlocktree( res->ci );
                    continue;
                    }
                }
            
            /* If we owned a previous task, unlock it. */
            if ( ind_best >= 0 ) {
                res = &qtasks[ qtid[ ind_best ] ];
                if ( res->type == task_type_self || res->type == task_type_pair || res->type == task_type_sub )
                    if ( res->ci != ci_best && res->ci != cj_best )
                        cell_unlocktree( res->ci );
                if ( res->type == task_type_pair || (res->type == task_type_sub && res->cj != NULL) )
                    if ( res->cj != ci_best && res->cj != cj_best )
                        cell_unlocktree( res->cj );
                }
            
            /* If we made it this far, we're safe. */
            ind_best = k;
            ci_best = qtasks[ qtid[ k ] ].ci;
            cj_best = qtasks[ qtid[ k ] ].cj;
            score_best = score;
            hits += 1;
            
            /* Should we bother looking any farther? */
            if ( ( qtasks[ qtid[ ind_best ] ].cj == NULL && score_best == 1 ) ||
                score_best == 2 );
                break;
        
            } /* loop over the task IDs. */
            
        /* Did we get a task? */
        if ( ind_best >= 0 ) {
        
            /* Do we need to swap? */
            if ( ind_best != q->next )
                COUNT(queue_counter_swap);
        
            /* get the task ID. */
            tid = qtid[ ind_best ];
            
            /* Own the cells involved. */
            qtasks[ tid ].ci->owner = rid;
            if ( qtasks[ tid ].cj != NULL )
                qtasks[ tid ].cj->owner = rid;
        
            /* Remove the task? */
            if ( keep ) {
            
                /* Bubble-up. */
                q->count = qcount - 1;
                for ( k = ind_best ; k < qcount - 1 ; k++ )
                    qtid[k] = qtid[k+1];
            
                }
                
            /* No, leave it in the queue. */
            else {
            
                TIMER_TIC2

                /* Bubble-down the task. */
                for ( k = ind_best ; k > q->next ; k-- )
                    qtid[ k ] = qtid[ k-1 ];
                qtid[ q->next ] = tid;
                
                /* up the counter. */
                q->next += 1;
                
                TIMER_TOC2(queue_timer_bubble);
            
                }
            
            }
    
        /* Release the task lock. */
        if ( lock_unlock( qlock ) != 0 )
            error( "Unlocking the task_lock failed.\n" );
            
        /* Leave? */
        if ( tid >= 0 ) {
            TIMER_TOC(queue_timer_gettask);
            return &qtasks[tid];
            }
        else if ( !blocking )
            break;
    
        } /* while there are tasks. */
        
    /* No beef. */
    TIMER_TOC(queue_timer_gettask);
    return NULL;

    }


/**
 * @brief Sort the tasks IDs according to their weight and constraints.
 *
 * @param q The #queue.
 */
 
void queue_sort ( struct queue *q ) {

    struct {
        short int lo, hi;
        } qstack[20];
    int qpos, i, j, k, lo, hi, imin, temp;
    int pivot_weight, pivot_wait;
    int *weight, *wait;
    int *data = q->tid;
    struct task *t;
    
    printf( "queue_sort: sorting queue with %i tasks.\n" , q->count );
        
    /* Allocate and pre-compute each task's weight. */
    if ( ( weight = (int *)alloca( sizeof(int) * q->count ) ) == NULL ||
         ( wait = (int *)alloca( sizeof(int) * q->count ) ) == NULL )
        error( "Failed to allocate weight buffer." );
    for ( k = 0 ; k < q->count ; k++ ) {
        t = &q->tasks[ q->tid[k] ];
        switch ( t->type ) {
            case task_type_self:
                wait[k] = t->rank;
                weight[k] = 0; // t->ci->count * t->ci->count;
                break;
            case task_type_pair:
                wait[k] = t->rank;
                weight[k] = 0; // t->ci->count * t->cj->count;
                break;
            case task_type_sub:
                wait[k] = t->rank;
                weight[k] = 0; // (t->cj == NULL) ? t->ci->count * t->ci->count : t->ci->count * t->cj->count;
                break;
            case task_type_sort:
                wait[k] = t->rank;
                weight[k] = 0; // t->ci->count;
                break;
            case task_type_ghost:
                wait[k] = t->rank;
                weight[k] = 0; // t->ci->count;
                break;
            }
        }
        
    /* Sort tasks. */
    qstack[0].lo = 0; qstack[0].hi = q->count - 1; qpos = 0;
    while ( qpos >= 0 ) {
        lo = qstack[qpos].lo; hi = qstack[qpos].hi;
        qpos -= 1;
        if ( hi - lo < 15 ) {
            for ( i = lo ; i < hi ; i++ ) {
                imin = i;
                for ( j = i+1 ; j <= hi ; j++ )
                    if ( ( wait[ j ] < wait[ imin ] ) ||
                         ( wait[ j ] == wait[ imin ] && weight[ j ] > weight[ imin ] ) )
                if ( imin != i ) {
                    temp = data[imin]; data[imin] = data[i]; data[i] = temp;
                    temp = wait[imin]; wait[imin] = wait[i]; wait[i] = temp;
                    temp = weight[imin]; weight[imin] = weight[i]; weight[i] = temp;
                    }
                }
            }
        else {
            pivot_weight = weight[ ( lo + hi ) / 2 ];
            pivot_wait = wait[ ( lo + hi ) / 2 ];
            i = lo; j = hi;
            while ( i <= j ) {
                while ( ( wait[ i ] < pivot_wait ) ||
                        ( wait[ i ] == pivot_wait && weight[ i ] > pivot_weight ) )
                    i++;
                while ( ( wait[ j ] > pivot_wait ) ||
                        ( wait[ j ] == pivot_wait && weight[ j ] < pivot_weight ) )
                    j--;
                if ( i <= j ) {
                    if ( i < j ) {
                        temp = data[i]; data[i] = data[j]; data[j] = temp;
                        temp = wait[i]; wait[i] = wait[j]; wait[j] = temp;
                        temp = weight[i]; weight[i] = weight[j]; weight[j] = temp;
                        }
                    i += 1; j -= 1;
                    }
                }
            if ( j > ( lo + hi ) / 2 ) {
                if ( lo < j ) {
                    qpos += 1;
                    qstack[qpos].lo = lo;
                    qstack[qpos].hi = j;
                    }
                if ( i < hi ) {
                    qpos += 1;
                    qstack[qpos].lo = i;
                    qstack[qpos].hi = hi;
                    }
                }
            else {
                if ( i < hi ) {
                    qpos += 1;
                    qstack[qpos].lo = i;
                    qstack[qpos].hi = hi;
                    }
                if ( lo < j ) {
                    qpos += 1;
                    qstack[qpos].lo = lo;
                    qstack[qpos].hi = j;
                    }
                }
            }
        }
                
    }
    
    
