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
 */
 
void engine_prepare ( struct engine *e ) {

    int j, k, qid;
    struct space *s = e->s;
    struct queue *q;
    
    TIMER_TIC

    /* Rebuild the space. */
    // tic = getticks();
    space_prepare( e->s );
    // printf( "engine_prepare: space_prepare took %.3f ms.\n" , (double)(getticks() - tic) / CPU_TPS * 1000 );
    
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

    /* Run throught the tasks and get all the waits right. */
    // tic = getticks();
    #pragma omp parallel for schedule(static) private(j)
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        if ( s->tasks[k].skip )
            continue;
        for ( j = 0 ; j < s->tasks[k].nr_unlock_tasks ; j++ )
            __sync_add_and_fetch( &s->tasks[k].unlock_tasks[j]->wait , 1 );
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
 * @brief Mapping function to set dt_min and dt_max, do the first
 * kick.
 */

void engine_map_kick_first ( struct cell *c , void *data ) {

    int j, k;
    struct engine *e = (struct engine *)data;
    float pdt, dt_step = e->dt_step, dt = e->dt, hdt = 0.5f*dt;
    float dt_min, dt_max, h_max, dx, h2dx_max, dx_max, a[3], v[3], u, u_dt, h, h_dt, rho, v_old[3];
    double x[3], x_old[3];
    struct part *restrict p;
    struct xpart *restrict xp;
    struct cpart *restrict cp;

    /* No children? */
    if ( !c->split ) {
    
        /* Init the min/max counters. */
        dt_min = FLT_MAX;
        dt_max = 0.0f;
        h_max = 0.0f;
        dx_max = 0.0f;
        h2dx_max = 0.0f;
    
        /* Loop over parts. */
        for ( k = 0 ; k < c->count ; k++ ) {
            
            /* Get a handle on the kth particle. */
            p = &c->parts[k];
            xp = p->xtras;
            cp = &c->cparts[k];
            
            /* Load the data locally. */
            a[0] = p->a[0]; a[1] = p->a[1]; a[2] = p->a[2];
            v[0] = p->v[0]; v[1] = p->v[1]; v[2] = p->v[2];
            x[0] = p->x[0]; x[1] = p->x[1]; x[2] = p->x[2];
            x_old[0] = xp->x_old[0]; x_old[1] = xp->x_old[1]; x_old[2] = xp->x_old[2];
            h = p->h;
            u = p->u;
            h_dt = p->force.h_dt;
            u_dt = p->force.u_dt;
            pdt = p->dt;
            
            /* Store the min/max dt. */
            dt_min = fminf( dt_min , pdt );
            dt_max = fmaxf( dt_max , pdt );
            
            /* Step and store the velocity and internal energy. */
            xp->v_old[0] = v_old[0] = v[0] + hdt * a[0];
            xp->v_old[1] = v_old[1] = v[1] + hdt * a[1];
            xp->v_old[2] = v_old[2] = v[2] + hdt * a[2];
            xp->u_old = p->u + hdt * p->force.u_dt;

            /* Move the particles with the velocitie at the half-step. */
            p->x[0] = x[0] += dt * v_old[0];
            p->x[1] = x[1] += dt * v_old[1];
            p->x[2] = x[2] += dt * v_old[2];
            dx = sqrtf( (x[0] - x_old[0])*(x[0] - x_old[0]) +
                        (x[1] - x_old[1])*(x[1] - x_old[1]) +
                        (x[2] - x_old[2])*(x[2] - x_old[2]) );
            dx_max = fmaxf( dx_max , dx );
            h2dx_max = fmaxf( h2dx_max , dx*2 + h );

            /* Update positions and energies at the half-step. */
            p->v[0] = v[0] += dt * a[0];
            p->v[1] = v[1] += dt * a[1];
            p->v[2] = v[2] += dt * a[2];
            p->u = u += u_dt * dt;
            p->h = h += h_dt * dt;
            h_max = fmaxf( h_max , h );

        
            /* Fill the cpart. */
            cp->x[0] = x[0];
            cp->x[1] = x[1];
            cp->x[2] = x[2];
            cp->h = h;
            cp->dt = pdt;
            
            /* Integrate other values if this particle will not be updated. */
            /* Init fields for density calculation. */
            if ( pdt > dt_step ) {
                rho = p->rho *= expf( -3.0f * h_dt / h * dt );
                p->force.POrho2 = u * ( const_gamma - 1.0f ) / ( rho + h * p->rho_dh / 3.0f );
                }
            else {
                p->density.wcount = 0.0f;
                p->density.wcount_dh = 0.0f;
                p->rho = 0.0f;
                p->rho_dh = 0.0f;
	            p->density.div_v = 0.0f;
	            for ( j = 0 ; j < 3 ; ++j)
	                p->density.curl_v[j] = 0.0f;
                }
                
            }
            
        }
        
    /* Otherwise, agregate data from children. */
    else {
    
        /* Init with the first non-null child. */
        for ( k = 0 ; c->progeny[k] == NULL ; k++ );
        dt_min = c->progeny[k]->dt_min;
        dt_max = c->progeny[k]->dt_max;
        h_max = c->progeny[k]->h_max;
        dx_max = c->progeny[k]->dx_max;
        h2dx_max = c->progeny[k]->h2dx_max;
        
        /* Loop over the remaining progeny. */
        for ( k += 1 ; k < 8 ; k++ )
            if ( c->progeny[k] != NULL ) {
                dt_min = fminf( dt_min , c->progeny[k]->dt_min );
                dt_max = fmaxf( dt_max , c->progeny[k]->dt_max );
                h_max = fmaxf( h_max , c->progeny[k]->h_max );
                dx_max = fmaxf( dx_max , c->progeny[k]->dx_max );
                h2dx_max = fmaxf( h2dx_max , c->progeny[k]->h2dx_max );
                }
    
        }

    /* Store the values. */
    c->dt_min = dt_min;
    c->dt_max = dt_max;
    c->h_max = h_max;
    c->dx_max = dx_max;
    c->h2dx_max = h2dx_max;
    
    /* Clean out the task pointers. */
    // c->sorts[0] = NULL;
    // c->nr_tasks = 0;
    // c->nr_density = 0;
    
    }


/**
 * @brief Compute the force on a single particle brute-force.
 */

void engine_single ( double *dim , long long int pid , struct part *__restrict__ parts , int N , int periodic ) {

    int i, k;
    double r2, dx[3];
    float fdx[3], ihg;
    struct part p, p2;
    
    /* Find "our" part. */
    for ( k = 0 ; k < N && parts[k].id != pid ; k++ );
    if ( k == N )
        error( "Part not found." );
    p = parts[k];
    
    /* Clear accumulators. */
    ihg = kernel_igamma / p.h;
    p.a[0] = 0.0f; p.a[1] = 0.0f; p.a[2] = 0.0f;
    p.force.u_dt = 0.0f; p.force.h_dt = 0.0f; p.force.v_sig = 0.0f;
            
    /* Loop over all particle pairs (force). */
    for ( k = 0 ; k < N ; k++ ) {
        if ( parts[k].id == p.id )
            continue;
        for ( i = 0 ; i < 3 ; i++ ) {
            dx[i] = p.x[i] - parts[k].x[i];
            if ( periodic ) {
                if ( dx[i] < -dim[i]/2 )
                    dx[i] += dim[i];
                else if ( dx[i] > dim[i]/2 )
                    dx[i] -= dim[i];
                }
            fdx[i] = dx[i];
            }
        r2 = fdx[0]*fdx[0] + fdx[1]*fdx[1] + fdx[2]*fdx[2];
        if ( r2 < p.h*p.h || r2 < parts[k].h*parts[k].h ) {
            p.a[0] = 0.0f; p.a[1] = 0.0f; p.a[2] = 0.0f;
            p.force.u_dt = 0.0f; p.force.h_dt = 0.0f; p.force.v_sig = 0.0f;
            p2 = parts[k];
            runner_iact_nonsym_force( r2 , fdx , p.h , p2.h , &p , &p2 );
            printf( "pairs_simple: interacting particles %lli and %lli (rho=%.3e,r=%e): a=[%.3e,%.3e,%.3e], u_dt=%.3e, h_dt=%.3e, v_sig=%.3e.\n" ,
                pid , p2.id , p2.rho , sqrtf(r2) ,
                p.a[0] , p.a[1] , p.a[2] ,
                p.force.u_dt , p.force.h_dt , p.force.v_sig );
            }
        }
        
    /* Dump the result. */
    ihg = kernel_igamma / p.h;
    printf( "pairs_single: part %lli (h=%e) has wcount=%f, rho=%f.\n" , p.id , p.h , p.density.wcount + 32.0/3 , ihg * ihg * ihg * ( p.rho + p.mass*kernel_root ) );
    fflush(stdout);
    
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
    float pdt, h, m, v[3], u;
    double x[3];
    double epot = 0.0, ekin = 0.0, lepot, lekin, lmom[3], mom[3] = { 0.0 , 0.0 , 0.0 };
    double lang[3], ang[3] = { 0.0 , 0.0 , 0.0 };
    // double lent, ent = 0.0;
    int threadID, nthreads, count = 0, lcount;
    
    TIMER_TIC2

    /* Get the maximum dt. */
    dt_step = 2.0f*dt;
    for ( k = 0 ; k < 32 && (e->step & (1 << k)) == 0 ; k++ )
        dt_step *= 2;
    // dt_step = FLT_MAX;
        
    /* Set the maximum dt. */
    e->dt_step = dt_step;
    e->s->dt_step = dt_step;
    printf( "engine_step: dt_step set to %.3e (dt=%.3e).\n" , dt_step , e->dt ); fflush(stdout);
    
    // printParticle( parts , 432626 );
    
    /* First kick. */
    TIMER_TIC
    space_map_cells_post( e->s , 1 , &engine_map_kick_first , e );
    TIMER_TOC( timer_kick1 );
        
    // for(k=0; k<10; ++k)
    //   printParticle(parts, k);
    // printParticle( parts , 432626 );
    // printParticle( parts , 432628 );
 
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
    TIMER_TOC(timer_runners);
    
    // for(k=0; k<10; ++k)
    //   printParticle(parts, k);
    // engine_single( e->s->dim , 432626 , e->s->parts , e->s->nr_parts , e->s->periodic );
    // printParticle( parts , 432626 );
    // printParticle( parts , 432628 );

    /* Second kick. */
    TIMER_TIC_ND
    dt_min = FLT_MAX; dt_max = 0.0f;
    #pragma omp parallel private(p,xp,k,ldt_min,lcount,ldt_max,lmom,lang,lekin,lepot,threadID,nthreads,pdt,v,h,m,x,u)
    {
        threadID = omp_get_thread_num();
        nthreads = omp_get_num_threads();
        ldt_min = FLT_MAX; ldt_max = 0.0f;
        lmom[0] = 0.0; lmom[1] = 0.0; lmom[2] = 0.0;
        lang[0] = 0.0; lang[1] = 0.0; lang[2] = 0.0;
        lekin = 0.0; lepot = 0.0; /* lent=0.0; */ lcount = 0;
        for ( k = nr_parts * threadID / nthreads ; k < nr_parts * (threadID + 1) / nthreads ; k++ ) {

            /* Get a handle on the part. */
            p = &parts[k];
            xp = p->xtras;
            
            /* Get local copies of particle data. */
            pdt = p->dt;
            h = p->h;
            m = p->mass;
            x[0] = p->x[0]; x[1] = p->x[1]; x[2] = p->x[2];

            /* Scale the derivatives if they're freshly computed. */
            if ( pdt <= dt_step ) {
                p->force.h_dt *= h * 0.333333333f;
                lcount += 1;
                }
                
            /* Update the particle's time step. */
            p->dt = pdt = const_cfl * h / ( p->force.v_sig );

            /* Update positions and energies at the half-step. */
            p->v[0] = v[0] = xp->v_old[0] + hdt * p->a[0];
            p->v[1] = v[1] = xp->v_old[1] + hdt * p->a[1];
            p->v[2] = v[2] = xp->v_old[2] + hdt * p->a[2];
            p->u = u = xp->u_old + hdt * p->force.u_dt;
            
            /* Get the smallest/largest dt. */
            ldt_min = fminf( ldt_min , pdt );
            ldt_max = fmaxf( ldt_max , pdt );
            
            /* Collect total energy. */
            lekin += 0.5 * m * ( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
            lepot += m * u;

            /* Collect momentum */
            lmom[0] += m * p->v[0];
            lmom[1] += m * p->v[1];
            lmom[2] += m * p->v[2];
	    
	        /* Collect angular momentum */
	        lang[0] += m * ( x[1]*v[2] - v[2]*x[1] );
	        lang[1] += m * ( x[2]*v[0] - v[0]*x[2] );
	        lang[2] += m * ( x[0]*v[1] - v[1]*x[0] );

	        /* Collect entropic function */
	        // lent += u * pow( p->rho, 1.f-const_gamma );

            }
        #pragma omp critical
        {
            dt_min = fminf( dt_min , ldt_min );
            dt_max = fmaxf( dt_max , ldt_max );
            mom[0] += lmom[0];
            mom[1] += lmom[1];
            mom[2] += lmom[2];
            ang[0] += lang[0];
            ang[1] += lang[1];
            ang[2] += lang[2];
	        // ent += (const_gamma -1.) * lent;
            epot += lepot;
            ekin += lekin;
            count += lcount;
        }
    }
    TIMER_TOC( timer_kick2 );
    e->dt_min = dt_min;
    // printParticle( parts , 432626 );
    printf( "engine_step: dt_min/dt_max is %e/%e.\n" , dt_min , dt_max ); fflush(stdout);
    printf( "engine_step: etot is %e (ekin=%e, epot=%e).\n" , ekin+epot , ekin , epot ); fflush(stdout);
    printf( "engine_step: total momentum is [ %e , %e , %e ].\n" , mom[0] , mom[1] , mom[2] ); fflush(stdout);
    printf( "engine_step: total angular momentum is [ %e , %e , %e ].\n" , ang[0] , ang[1] , ang[2] ); fflush(stdout);
    // printf( "engine_step: total entropic function is %e .\n", ent ); fflush(stdout);
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
        
    TIMER_TOC2(timer_step);
    
    }
    
    
/**
 * @brief init an engine with the given number of threads, queues, and
 *      the given policy.
 *
 * @param e The #engine.
 * @param s The #space in which this #runner will run.
 * @param dt The initial time step to use.
 * @param nr_threads The number of threads to spawn.
 * @param nr_queues The number of task queues to create.
 * @param policy The queueing policy to use.
 */
 
void engine_init ( struct engine *e , struct space *s , float dt , int nr_threads , int nr_queues , int policy ) {

    #if defined(HAVE_SETAFFINITY)
        cpu_set_t cpuset;
    #endif
    int k;
    
    /* Store the values. */
    e->s = s;
    e->nr_threads = nr_threads;
    e->nr_queues = nr_queues;
    e->policy = policy;
    e->step = 0;
    
    /* First of all, init the barrier and lock it. */
    if ( pthread_mutex_init( &e->barrier_mutex , NULL ) != 0 )
        error( "Failed to initialize barrier mutex." );
    if ( pthread_cond_init( &e->barrier_cond , NULL ) != 0 )
        error( "Failed to initialize barrier condition variable." );
    if ( pthread_mutex_lock( &e->barrier_mutex ) != 0 )
        error( "Failed to lock barrier mutex." );
    e->barrier_count = 0;
    
    /* Run through the parts and get the minimum time step. */
    for ( k = 0 ; k < s->nr_parts ; k++ )
        if ( s->parts[k].dt > 0.0f )
            while ( s->parts[k].dt < dt )
                dt *= 0.5;
    e->dt_min = dt;
    
    /* Allocate the queues. */
    if ( posix_memalign( (void *)(&e->queues) , 64 , nr_queues * sizeof(struct queue) ) != 0 )
        error( "Failed to allocate queues." );
    bzero( e->queues , nr_queues * sizeof(struct queue) );
        
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
    
    
    
