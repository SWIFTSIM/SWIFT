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
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <omp.h>
#include <sched.h>

/* Local headers. */
#include "cycle.h"
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
 
void engine_prepare ( struct engine *e ) {

    int j, k, qid;
    struct space *s = e->s;
    struct queue *q;
    float dt_step = e->dt_step;
    
    TIMER_TIC

    /* Rebuild the space. */
    // tic = getticks();
    space_prepare( e->s );
    // printf( "engine_prepare: space_prepare with %i changes took %.3f ms.\n" , changes , (double)(getticks() - tic) / CPU_TPS * 1000 );
    
    // tic = getticks();
    /* Init the queues (round-robin). */
    for ( qid = 0 ; qid < e->nr_queues ; qid++ )
        queue_init( &e->queues[qid] , s->nr_tasks , s->tasks );

    /* Fill the queues (round-robin). */
    for ( qid = 0 , k = 0 ; k < s->nr_tasks ; k++ ) {
        if ( s->tasks[ s->tasks_ind[k] ].skip )
            continue;
        q = &e->queues[qid];
        qid = ( qid + 1 ) % e->nr_queues;
        q->tid[ q->count ] = s->tasks_ind[k];
        q->count += 1;
        }
    // printf( "engine_prepare: re-filling queues took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );

    /* Re-set the particle data. */
    // tic = getticks();
    #pragma omp parallel for schedule(static) 
    for ( k = 0 ; k < s->nr_parts ; k++ )
        if ( s->parts[k].dt <= dt_step ) {
            s->parts[k].wcount = 0.0f;
            s->parts[k].wcount_dh = 0.0f;
            s->parts[k].rho = 0.0f;
            s->parts[k].rho_dh = 0.0f;
	    s->parts[k].div_v = 0.0f;
	    for ( k = 0 ; k < 3 ; ++k)
	      s->parts[k].curl_v[k] = 0.0f;
            }
    // printf( "engine_prepare: re-setting particle data took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );
    
    /* Run throught the tasks and get all the waits right. */
    // tic = getticks();
    #pragma omp parallel for schedule(static) private(j)
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        for ( j = 0 ; j < s->tasks[k].nr_unlock_tasks ; j++ )
            __sync_add_and_fetch( &s->tasks[k].unlock_tasks[j]->wait , 1 );
        for ( j = 0 ; j < s->tasks[k].nr_unlock_cells ; j++ )
            __sync_add_and_fetch( &s->tasks[k].unlock_cells[j]->wait , 1 );
        }
    // printf( "engine_prepare: preparing task dependencies took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );
    
    /* Re-set the queues.*/
    for ( k = 0 ; k < e->nr_queues ; k++ )
        e->queues[k].next = 0;
        
    TIMER_TOC( timer_prepare );
    
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
 
void engine_step ( struct engine *e , int sort_queues ) {

    int k, nr_parts = e->s->nr_parts;
    struct part *restrict parts = e->s->parts, *restrict p;
    struct xpart *restrict xp;
    float dt = e->dt, hdt = 0.5*dt, dt_step, dt_max, dt_min, ldt_min, ldt_max;
    double epot = 0.0, ekin = 0.0, lepot, lekin, lmom[3], mom[3] = { 0.0 , 0.0 , 0.0 };
    int threadID, nthreads, count = 0, lcount;
    // #ifdef __SSE2__
    //     VEC_MACRO(4,float) hdtv = _mm_set1_ps( hdt );
    // #endif

    /* Get the maximum dt. */
    dt_step = 2.0f*dt;
    for ( k = 0 ; k < 32 && (e->step & (1 << k)) == 0 ; k++ )
        dt_step *= 2;
    // dt_step = FLT_MAX;
        
    /* Set the maximum dt. */
    e->dt_step = dt_step;
    e->s->dt_step = dt_step;
    printf( "engine_step: dt_step set to %.3e.\n" , dt_step ); fflush(stdout);
    
    /* First kick. */
    TIMER_TIC
    #pragma omp parallel for schedule(static) private(p,xp)
    for ( k = 0 ; k < nr_parts ; k++ ) {
        
        /* Get a handle on the part. */
        p = &parts[k];
        xp = p->xtras;
        
        /* Step and store the velocity and internal energy. */
        // #ifdef __SSE__
        //     _mm_store_ps( &v_bar[4*k] , _mm_load_ps( &p->v[0] ) + hdtv * _mm_load_ps( &p->a[0] ) );
        // #else
            xp->v_old[0] = p->v[0] + hdt * p->a[0];
            xp->v_old[1] = p->v[1] + hdt * p->a[1];
            xp->v_old[2] = p->v[2] + hdt * p->a[2];
        // #endif
        // xp->u_old = fmaxf( p->u + hdt * p->u_dt , FLT_EPSILON );
        xp->u_old = p->u + hdt * p->u_dt;
        
        /* Move the particles with the velocitie at the half-step. */
        p->x[0] += dt * xp->v_old[0];
        p->x[1] += dt * xp->v_old[1];
        p->x[2] += dt * xp->v_old[2];
        
        /* Update positions and energies at the half-step. */
        p->v[0] += dt * p->a[0];
        p->v[1] += dt * p->a[1];
        p->v[2] += dt * p->a[2];
        p->u *= expf( p->u_dt / p->u * dt );
        p->h *= expf( p->h_dt / p->h * dt );
        
        /* Integrate other values if this particle will not be updated. */
        if ( p->dt > dt_step ) {
            p->rho *= expf( -3.0f * p->h_dt / p->h * dt );
            p->POrho2 = p->u * ( const_gamma - 1.0f ) / ( p->rho + p->h * p->rho_dh / 3.0f );
            }
        
        }
    TIMER_TOC( timer_kick1 );
        
    // for(k=0; k<10; ++k)
    //   printParticle(parts, k);
    // printParticle( parts , 494849 );
 
    /* Prepare the space. */
    engine_prepare( e );
    
    /* Sort the queues?*/
    if ( sort_queues ) {
        #pragma omp parallel for default(none), shared(e)
        for ( k = 0 ; k < e->nr_queues ; k++ ) {
            queue_sort( &e->queues[k] );
            e->queues[k].next = 0;
            }
        }
        
    /* Start the clock. */
    TIMER_TIC_ND
    
    /* Cry havoc and let loose the dogs of war. */
    e->barrier_count = -e->barrier_count;
    if ( pthread_cond_broadcast( &e->barrier_cond ) != 0 )
        error( "Failed to broadcast barrier open condition." );
        
    /* Sit back and wait for the runners to come home. */
    while ( e->barrier_count < e->nr_threads )
        if ( pthread_cond_wait( &e->barrier_cond , &e->barrier_mutex ) != 0 )
            error( "Error while waiting for barrier." );
            
    /* Stop the clock. */
    TIMER_TOC(timer_step);
    
    // for(k=0; k<10; ++k)
    //   printParticle(parts, k);
    // printParticle( parts , 494849 );

    /* Second kick. */
    TIMER_TIC_ND
    dt_min = FLT_MAX; dt_max = 0.0f;
    #pragma omp parallel private(p,xp,k,ldt_min,lcount,ldt_max,lmom,lekin,lepot,threadID,nthreads)
    {
        threadID = omp_get_thread_num();
        nthreads = omp_get_num_threads();
        ldt_min = FLT_MAX; ldt_max = 0.0f;
        lmom[0] = 0.0; lmom[1] = 0.0; lmom[2] = 0.0;
        lekin = 0.0; lepot = 0.0; lcount = 0;
        for ( k = nr_parts * threadID / nthreads ; k < nr_parts * (threadID + 1) / nthreads ; k++ ) {

            /* Get a handle on the part. */
            p = &parts[k];
            xp = p->xtras;

            /* Scale the derivatives if they're freshly computed. */
            if ( p->dt <= dt_step ) {
                p->u_dt *= p->POrho2;
                p->h_dt *= p->h * 0.333333333f;
                lcount += 1;
                }
                
            /* Update the particle's time step. */
            p->dt = const_cfl * p->h / ( p->c + p->v_sig );

            /* Update positions and energies at the half-step. */
            // #ifdef __SSE__
            //     _mm_store_ps( &p->v[0] , _mm_load_ps( &v_bar[4*k] ) + hdtv * _mm_load_ps( &p->a[0] )  );
            // #else
                p->v[0] = xp->v_old[0] + hdt * p->a[0];
                p->v[1] = xp->v_old[1] + hdt * p->a[1];
                p->v[2] = xp->v_old[2] + hdt * p->a[2];
            // #endif
            // p->u = fmaxf( xp->u_old + hdt * p->u_dt , FLT_EPSILON );
            p->u = xp->u_old + hdt * p->u_dt;
            
            /* Get the smallest/largest dt. */
            ldt_min = fminf( ldt_min , p->dt );
            ldt_max = fmaxf( ldt_max , p->dt );
            
            /* Collect total energy. */
            lekin += 0.5 * p->mass * ( p->v[0]*p->v[0] + p->v[1]*p->v[1] + p->v[2]*p->v[2] );
            lepot += p->mass * p->u;

            /* Collect momentum */
            lmom[0] += p->mass * p->v[0];
            lmom[1] += p->mass * p->v[1];
            lmom[2] += p->mass * p->v[2];

            }
        #pragma omp critical
        {
            dt_min = fminf( dt_min , ldt_min );
            dt_max = fmaxf( dt_max , ldt_max );
            mom[0] += lmom[0];
            mom[1] += lmom[1];
            mom[2] += lmom[2];
            epot += lepot;
            ekin += lekin;
            count += lcount;
        }
    }
    TIMER_TOC( timer_kick2 );
    e->dt_min = dt_min;
    printf( "engine_step: dt_min/dt_max is %e/%e.\n" , dt_min , dt_max ); fflush(stdout);
    printf( "engine_step: etot is %e (ekin=%e, epot=%e).\n" , ekin+epot , ekin , epot ); fflush(stdout);
    printf( "engine_step: total momentum is [ %e , %e , %e ].\n" , mom[0] , mom[1] , mom[2] ); fflush(stdout);
    printf( "engine_step: updated %i parts (dt_step=%.3e).\n" , count , dt_step ); fflush(stdout);
        
    /* Increase the step counter. */
    e->step += 1;
    
    /* Does the time step need adjusting? */
    while ( dt_min < e->dt ) {
        e->dt *= 0.5;
        e->step *= 2;
        printf( "engine_step: dt_min dropped below time step, adjusting to dt=%e.\n" , e->dt );
        }
    while ( dt_min > 2*e->dt && (e->step & 1) == 0 ) {
        e->dt *= 2.0;
        e->step /= 2;
        printf( "engine_step: dt_min is larger than twice the time step, adjusting to dt=%e.\n" , e->dt );
        }
    
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
    e->dt_min = 0.0f;
    e->step = 0;
    
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
    // for ( k = 0 ; k < nr_queues ; k++ )
    //     queue_sort( &e->queues[k] );
        
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
    
    
    
