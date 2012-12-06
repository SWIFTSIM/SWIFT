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

/* The timers. */
ticks runner_timer[ runner_timer_count ];

/* The counters. */
int runner_counter[ runner_counter_count ];

        

const float runner_shift[13*3] = {
     5.773502691896258e-01 ,  5.773502691896258e-01 ,  5.773502691896258e-01 ,
     7.071067811865475e-01 ,  7.071067811865475e-01 ,  0.0                   ,
     5.773502691896258e-01 ,  5.773502691896258e-01 , -5.773502691896258e-01 ,
     7.071067811865475e-01 ,  0.0                   ,  7.071067811865475e-01 ,
     1.0                   ,  0.0                   ,  0.0                   ,
     7.071067811865475e-01 ,  0.0                   , -7.071067811865475e-01 ,
     5.773502691896258e-01 , -5.773502691896258e-01 ,  5.773502691896258e-01 ,
     7.071067811865475e-01 , -7.071067811865475e-01 ,  0.0                   ,
     5.773502691896258e-01 , -5.773502691896258e-01 , -5.773502691896258e-01 ,
     0.0                   ,  7.071067811865475e-01 ,  7.071067811865475e-01 ,
     0.0                   ,  1.0                   ,  0.0                   ,
     0.0                   ,  7.071067811865475e-01 , -7.071067811865475e-01 ,
     0.0                   ,  0.0                   ,  1.0                   ,
    };
const char runner_flip[27] = { 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 0 ,
                               0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }; 


/* Import the density functions. */
#define FUNCTION density
#include "runner_doiact.h"

#undef FUNCTION
#define FUNCTION force
#include "runner_doiact.h"


/**
 * @brief Sort the entries in ascending order using QuickSort.
 *
 * @param sort The entries
 * @param N The number of entries.
 */
 
void runner_dosort_ascending ( struct entry *sort , int N ) {

    struct {
        short int lo, hi;
        } qstack[10];
    int qpos, i, j, lo, hi, imin;
    struct entry temp;
    float pivot;
        
    /* Sort parts in cell_i in decreasing order with quicksort */
    qstack[0].lo = 0; qstack[0].hi = N - 1; qpos = 0;
    while ( qpos >= 0 ) {
        lo = qstack[qpos].lo; hi = qstack[qpos].hi;
        qpos -= 1;
        if ( hi - lo < 15 ) {
            for ( i = lo ; i < hi ; i++ ) {
                imin = i;
                for ( j = i+1 ; j <= hi ; j++ )
                    if ( sort[j].d < sort[imin].d )
                        imin = j;
                if ( imin != i ) {
                    temp = sort[imin]; sort[imin] = sort[i]; sort[i] = temp;
                    }
                }
            }
        else {
            pivot = sort[ ( lo + hi ) / 2 ].d;
            i = lo; j = hi;
            while ( i <= j ) {
                while ( sort[i].d < pivot ) i++;
                while ( sort[j].d > pivot ) j--;
                if ( i <= j ) {
                    if ( i < j ) {
                        temp = sort[i]; sort[i] = sort[j]; sort[j] = temp;
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
    
    
/**
 * @brief Sort the particles in the given cell along all cardinal directions.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
 
void runner_dosort ( struct runner *r , struct cell *c , int flags ) {

    struct entry *finger;
    struct entry *fingers[8];
    struct part *parts = c->parts;
    int j, k, count = c->count;
    int i, ind, off[8], inds[8], temp_i;
    // float shift[3];
    float buff[8], px[3];
    struct cell *temp_c;
    TIMER_TIC
    
    /* Does this cell even need to be sorted? */
    for ( temp_c = c ; temp_c != NULL && temp_c->nr_pairs == 0 ; temp_c = temp_c->parent );
    if ( temp_c == NULL )
        return;

    /* start by allocating the entry arrays. */
    if ( lock_lock( &c->lock ) != 0 )
        error( "Failed to lock cell." );
    if ( c->sort == NULL )
        if ( ( c->sort = (struct entry *)malloc( sizeof(struct entry) * (c->count + 1) * 14 ) ) == NULL )
            error( "Failed to allocate sort memory." );
    if ( lock_unlock( &c->lock ) != 0 )
        error( "Failed to unlock cell." );
        
    /* Does this cell have any progeny? */
    if ( c->split ) {
    
        /* Loop over the 13+1 different sort arrays. */
        for ( j = 0 ; j < 14 ; j++ ) {
        
            /* Has this sort array been flagged? */
            if ( !( flags & (1 << j) ) )
                continue;
                
            /* Init the particle index offsets. */
            for ( off[0] = 0 , k = 1 ; k < 8 ; k++ )
                if ( c->progeny[k-1] != NULL )
                    off[k] = off[k-1] + c->progeny[k-1]->count;
                else
                    off[k] = off[k-1];

            /* Init the entries and indices. */
            for ( k = 0 ; k < 8 ; k++ ) {
                inds[k] = k;
                if ( c->progeny[k] != NULL && c->progeny[k]->count > 0 ) {
                    fingers[k] = &c->progeny[k]->sort[ j*(c->progeny[k]->count + 1) ];
                    buff[k] = fingers[k]->d;
                    off[k] = off[k];
                    }
                else
                    buff[k] = FLT_MAX;
                }

            /* Sort the buffer. */
            for ( i = 0 ; i < 7 ; i++ )
                for ( k = i+1 ; k < 8 ; k++ )
                    if ( buff[ inds[k] ] < buff[ inds[i] ] ) {
                        temp_i = inds[i]; inds[i] = inds[k]; inds[k] = temp_i;
                        }

            /* For each entry in the new sort list. */
            finger = &c->sort[ j*(count + 1) ];
            for ( ind = 0 ; ind < count ; ind++ ) {

                /* Copy the minimum into the new sort array. */
                finger[ind].d = buff[inds[0]];
                finger[ind].i = fingers[inds[0]]->i + off[inds[0]];

                /* Update the buffer. */
                fingers[inds[0]] += 1;
                buff[inds[0]] = fingers[inds[0]]->d;

                /* Find the smallest entry. */
                for ( k = 1 ; k < 8 && buff[inds[k]] < buff[inds[k-1]] ; k++ ) {
                    temp_i = inds[k-1]; inds[k-1] = inds[k]; inds[k] = temp_i;
                    }

                } /* Merge. */
            
            /* Add a sentinel. */
            c->sort[ j*(c->count + 1) + c->count ].d = FLT_MAX;
            c->sort[ j*(c->count + 1) + c->count ].i = 0;
            
            } /* loop over sort arrays. */
    
        } /* progeny? */
        
    /* Otherwise, just sort. */
    else {
    
        /* Fill the sort array. */
        for ( k = 0 ; k < count ; k++ ) {
            px[0] = parts[k].x[0];
            px[1] = parts[k].x[1];
            px[2] = parts[k].x[2];
            for ( j = 0 ; j < 13 ; j++ )
                if ( flags & (1 << j) ) {
                    c->sort[ j*(count + 1) + k].i = k;
                    c->sort[ j*(count + 1) + k].d = px[0]*runner_shift[ 3*j + 0 ] + px[1]*runner_shift[ 3*j + 1 ] + px[2]*runner_shift[ 3*j + 2 ];
                    }
            if ( flags & (1 << 14) ) {
                c->sort[ 14*(count + 1) + k ].i = k;
                c->sort[ 14*(count + 1) + k ].d = parts[k].dt;
                }
            }

        /* Add the sentinel and sort. */
        for ( j = 0 ; j < 14 ; j++ )
            if ( flags & (1 << j) ) {
                c->sort[ j*(count + 1) + c->count ].d = FLT_MAX;
                c->sort[ j*(count + 1) + c->count ].i = 0;
                runner_dosort_ascending( &c->sort[ j*(count + 1) ] , c->count );
                }
            
        }
        
    /* Verify the sorting. */
    /* for ( j = 0 ; j < 13 ; j++ ) {
        if ( !( flags & (1 << j) ) )
            continue;
        finger = &c->sort[ j*(c->count + 1) ];
        for ( k = 1 ; k < c->count ; k++ ) {
            if ( finger[k].d < finger[k-1].d )
                error( "Sorting failed, ascending array." );
            if ( finger[k].i >= c->count )
                error( "Sorting failed, indices borked." );
            }
        } */

    #ifdef TIMER_VERBOSE
        printf( "runner_dosort[%02i]: %i parts at depth %i (flags = %i%i%i%i%i%i%i%i%i%i%i%i%i) took %.3f ms.\n" ,
            r->id , c->count , c->depth ,
            (flags & 0x1000) >> 12 , (flags & 0x800) >> 11 , (flags & 0x400) >> 10 , (flags & 0x200) >> 9 , (flags & 0x100) >> 8 , (flags & 0x80) >> 7 , (flags & 0x40) >> 6 , (flags & 0x20) >> 5 , (flags & 0x10) >> 4 , (flags & 0x8) >> 3 , (flags & 0x4) >> 2 , (flags & 0x2) >> 1 , (flags & 0x1) >> 0 , 
            ((double)TIMER_TOC(runner_timer_dosort)) / CPU_TPS * 1000 ); fflush(stdout);
    #else
        TIMER_TOC(runner_timer_dosort);
    #endif

    }
    
    
/**
 * @brief Intermediate task between density and force
 *
 * @param r The runner thread.
 * @param ci THe cell.
 */
 
void runner_doghost ( struct runner *r , struct cell *c ) {

    struct part *p;
    struct cell *finger, *finger_prev;;
    int i, k, redo, count = c->count;
    int *pid;
    float ihg, ihg2;
    TIMER_TIC
    
    /* Recurse? */
    if ( c->split ) {
        for ( k = 0 ; k < 8 ; k++ )
            if ( c->progeny[k] != NULL )
                runner_doghost( r , c->progeny[k] );
        return;
        }
    
    /* Init the IDs that have to be updated. */
    if ( ( pid = (int *)alloca( sizeof(int) * count ) ) == NULL )
        error( "Call to alloca failed." );
    for ( k = 0 ; k < count ; k++ )
        pid[k] = k;
        
    /* While there are particles that need to be updated... */
    while ( count > 0 ) {
    
        /* Reset the redo-count. */
        redo = 0;
    
        /* Loop over the parts in this cell. */
        for ( i = 0 ; i < count ; i++ ) {

            /* Get a direct pointer on the part. */
            p = &c->parts[ pid[i] ];

            /* Adjust the computed rho. */
            ihg = kernel_igamma / p->h;
            ihg2 = ihg * ihg;
            p->rho *= ihg * ihg2;
            p->rho_dh *= ihg2 * ihg2;
            
            /* Update the smoothing length. */
            p->h -= ( p->wcount + kernel_root - const_nwneigh ) / p->wcount_dh;

            /* Did we get the right number density? */
            if ( p->wcount + kernel_root > const_nwneigh + 1 ||
                 p->wcount + kernel_root < const_nwneigh - 1 ) {
                // printf( "runner_doghost: particle %lli (h=%e,depth=%i) has bad wcount=%f.\n" , p->id , p->h , c->depth , p->wcount + kernel_root ); fflush(stdout);
                // p->h += ( p->wcount + kernel_root - const_nwneigh ) / p->wcount_dh;
                pid[redo] = pid[i];
                redo += 1;
                p->wcount = 0.0;
                p->wcount_dh = 0.0;
                p->rho = 0.0;
                p->rho_dh = 0.0;
                continue;
                }

            /* Reset the acceleration. */
            for ( k = 0 ; k < 3 ; k++ )
                p->a[k] = 0.0f;

            /* Reset the time derivatives. */
            p->u_dt = 0.0f;
            p->h_dt = 0.0f;

            /* Compute this particle's time step. */
            p->dt = const_cfl * p->h / sqrtf( const_gamma * ( const_gamma - 1.0f ) * p->u );

            /* Compute the pressure. */
            // p->P = p->rho * p->u * ( const_gamma - 1.0f );

            /* Compute the P/Omega/rho2. */
            p->POrho2 = p->u * ( const_gamma - 1.0f ) / ( p->rho + p->h * p->rho_dh / 3.0f );

            }
            
        /* Re-set the counter for the next loop (potentially). */
        count = redo;
        if ( count > 0 ) {
        
            // error( "Bad smoothing length, fixing this isn't implemented yet." );
            
            /* Climb up the cell hierarchy. */
            finger_prev = c;
            for ( finger = c ; finger != NULL ; finger = finger->parent ) {
            
                /* Run through this cell's density interactions. */
                for ( k = 0 ; k < finger->nr_density ; k++ ) {
                
                    /* Self-interaction? */
                    if ( finger->density[k]->type == task_type_self )
                        runner_doself_subset_density( r , finger , c->parts , pid , count );
                        
                    /* Otherwise, pair interaction? */
                    else if ( finger->density[k]->type == task_type_pair ) {
                    
                        /* Left or right? */
                        if ( finger->density[k]->ci == finger )
                            runner_dopair_subset_density( r , finger , c->parts , pid , count , finger->density[k]->cj );
                        else
                            runner_dopair_subset_density( r , finger , c->parts , pid , count , finger->density[k]->ci );
                        
                        }
                
                    /* Otherwise, sub interaction? */
                    else if ( finger->density[k]->type == task_type_sub ) 
                        runner_dosub_subset_density( r , finger->density[k]->ci , finger->density[k]->cj , finger_prev , c->parts , pid , count , finger->density[k]->flags );
                
                    }
                    
                /* Keep a finger on the previous cell. */
                finger_prev = finger;
            
                }
        
            }
            
        }

    #ifdef TIMER_VERBOSE
        printf( "runner_doghost[%02i]: %i parts at depth %i took %.3f ms.\n" ,
            r->id , c->count , c->depth ,
            ((double)TIMER_TOC(runner_timer_doghost)) / CPU_TPS * 1000 ); fflush(stdout);
    #else
        TIMER_TOC(runner_timer_doghost);
    #endif
    
    }


/**
 * @brief The #runner main thread routine.
 *
 * @param data A pointer to this thread's data.
 */
 
void *runner_main ( void *data ) {

    struct runner *r = (struct runner *)data;
    struct engine *e = r->e;
    int threadID = r->id;
    int k, qid, naq, keep, tpq;
    struct queue *queues[ e->nr_queues ], *myq;
    struct task *t;
    struct cell *ci, *cj;
    unsigned int myseed = rand() + r->id;
    #ifdef TIMER
        ticks stalled;
    #endif
    
    /* Main loop. */
    while ( 1 ) {
    
        /* Wait at the barrier. */
        engine_barrier( e );
        
        /* Set some convenient local data. */
        keep = e->policy & engine_policy_keep;
        myq = &e->queues[ threadID % e->nr_queues ];
        tpq = ceil( ((double)e->nr_threads) / e->nr_queues );
        #ifdef TIMER
            stalled = 0;
        #endif
        
        /* Set up the local list of active queues. */
        naq = e->nr_queues;
        for ( k = 0 ; k < naq ; k++ )
            queues[k] = &e->queues[k];
    
        /* Set up the local list of active queues. */
        naq = e->nr_queues;
        for ( k = 0 ; k < naq ; k++ )
            queues[k] = &e->queues[k];
    
        /* Loop while there are tasks... */
        while ( 1 ) {
        
            /* Remove any inactive queues. */
            for ( k = 0 ; k < naq ; k++ )
                if ( queues[k]->next == queues[k]->count ) {
                    naq -= 1;
                    queues[k] = queues[naq];
                    k -= 1;
                    }
            if ( naq == 0 )
                break;
        
            /* Get a task, how and from where depends on the policy. */
            TIMER_TIC
            t = NULL;
            if ( e->nr_queues == 1 ) {
                t = queue_gettask_old( &e->queues[0] , 1 , 0 );
                }
            else if ( e->policy & engine_policy_steal ) {
                if ( ( myq->next == myq->count ) ||
                     ( t = queue_gettask( myq , r->id , 0 , 0 ) ) == NULL ) {
                    TIMER_TIC2
                    qid = rand_r( &myseed ) % naq;
                    keep = ( e->policy & engine_policy_keep ) &&
                           ( myq->count <= myq->size-tpq );
                    if ( myq->next == myq->count )
                        COUNT(runner_counter_steal_empty);
                    else
                        COUNT(runner_counter_steal_stall);
                    t = queue_gettask( queues[qid] , r->id , 0 , keep );
                    if ( t != NULL && keep )
                        queue_insert( myq , t );
                    TIMER_TOC2(runner_timer_steal);
                    }
                }
            else if ( e->policy & engine_policy_rand ) {
                qid = rand_r( &myseed ) % naq;
                t = queue_gettask( queues[qid] , r->id , e->policy & engine_policy_block , 0 );
                }
            else {
                t = queue_gettask( &e->queues[threadID] , r->id , e->policy & engine_policy_block , 0 );
                }
            TIMER_TOC(runner_timer_getpair);
            
            /* Did I get anything? */
            if ( t == NULL ) {
                COUNT(runner_counter_stall);
                #ifdef TIMER
                    if ( !stalled )
                        stalled = getticks();
                #endif
                continue;
                }
            #ifdef TIMER
            else if ( stalled ) {
                stalled = getticks() - stalled;
                __sync_add_and_fetch( &runner_timer[runner_timer_stalled] , stalled );
                #ifdef TIMER_VERBOSE
                    printf( "runner_main[%02i]: stalled %.3f ms\n" , r->id , ((double)stalled) / CPU_TPS * 1000 );
                    fflush(stdout);
                #endif
                stalled = 0;
                }
            #endif
        
            /* Get the cells. */
            ci = t->ci;
            cj = t->cj;
            
            /* Different types of tasks... */
            switch ( t->type ) {
                case task_type_self:
                    if ( t->subtype == task_subtype_density )
                        runner_doself_density( r , ci );
                    else if ( t->subtype == task_subtype_force )
                        runner_doself_force( r , ci );
                    else
                        error( "Unknown task subtype." );
                    cell_unlocktree( ci );
                    break;
                case task_type_pair:
                    if ( t->subtype == task_subtype_density )
                        runner_dopair_density( r , ci , cj );
                    else if ( t->subtype == task_subtype_force )
                        runner_dopair_force( r , ci , cj );
                    else
                        error( "Unknown task subtype." );
                    cell_unlocktree( ci );
                    cell_unlocktree( cj );
                    break;
                case task_type_sort:
                    runner_dosort( r , ci , t->flags );
                    break;
                case task_type_sub:
                    if ( t->subtype == task_subtype_density )
                        runner_dosub_density( r , ci , cj , t->flags );
                    else if ( t->subtype == task_subtype_force )
                        runner_dosub_force( r , ci , cj , t->flags );
                    else
                        error( "Unknown task subtype." );
                    cell_unlocktree( ci );
                    if ( cj != NULL )
                        cell_unlocktree( cj );
                    break;
                case task_type_ghost:
                    if ( ci->super == ci )
                        runner_doghost( r , ci );
                    break;
                default:
                    error( "Unknown task type." );
                }
                
            /* Resolve any dependencies. */
            for ( k = 0 ; k < t->nr_unlock_tasks ; k++ )
                if ( __sync_fetch_and_sub( &t->unlock_tasks[k]->wait , 1 ) == 0 )
                    abort();
            for ( k = 0 ; k < t->nr_unlock_cells ; k++ )
                __sync_fetch_and_sub( &t->unlock_cells[k]->wait , 1 );
        
            } /* main loop. */
            
    	/* Any leftover stalls? */    
        #ifdef TIMER
        if ( stalled ) {
            stalled = getticks() - stalled;
            __sync_add_and_fetch( &runner_timer[runner_timer_stalled] , stalled );
            #ifdef TIMER_VERBOSE
                printf( "runner_main[%02i]: stalled %.3f ms\n" , r->id , ((double)stalled) / CPU_TPS * 1000 );
                fflush(stdout);
            #endif
            stalled = 0;
            }
        #endif
            
        }
        
    /* Be kind, rewind. */
    return NULL;

    }
    

