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
#include "scheduler.h"
#include "engine.h"
#include "runner.h"
#include "runner_iact.h"
#include "error.h"

/* Convert cell location to ID. */
#define cell_getid( cdim , i , j , k ) ( (int)(k) + (cdim)[2]*( (int)(j) + (cdim)[1]*(int)(i) ) )


/**
 * @brief Fill the #space's task list.
 *
 * @param s The #space we are working in.
 */
 
void engine_maketasks ( struct engine *e ) {

    struct space *s = e->s;
    struct scheduler *sched = &e->sched;
    int i, j, k, ii, jj, kk, iii, jjj, kkk, cid, cjd, sid;
    int *cdim = s->cdim;
    struct task *t, *t2;
    struct cell *ci, *cj;

    /* Re-set the scheduler. */
    scheduler_reset( sched , s->tot_cells * engine_maxtaskspercell );
    
    /* Run through the highest level of cells and add pairs. */
    for ( i = 0 ; i < cdim[0] ; i++ )
        for ( j = 0 ; j < cdim[1] ; j++ )
            for ( k = 0 ; k < cdim[2] ; k++ ) {
                cid = cell_getid( cdim , i , j , k );
                if ( s->cells[cid].count == 0 )
                    continue;
                ci = &s->cells[cid];
                if ( ci->count == 0 )
                    continue;
                scheduler_addtask( sched , task_type_self , task_subtype_density , 0 , 0 , ci , NULL , 0 );
                for ( ii = -1 ; ii < 2 ; ii++ ) {
                    iii = i + ii;
                    if ( !s->periodic && ( iii < 0 || iii >= cdim[0] ) )
                        continue;
                    iii = ( iii + cdim[0] ) % cdim[0];
                    for ( jj = -1 ; jj < 2 ; jj++ ) {
                        jjj = j + jj;
                        if ( !s->periodic && ( jjj < 0 || jjj >= cdim[1] ) )
                            continue;
                        jjj = ( jjj + cdim[1] ) % cdim[1];
                        for ( kk = -1 ; kk < 2 ; kk++ ) {
                            kkk = k + kk;
                            if ( !s->periodic && ( kkk < 0 || kkk >= cdim[2] ) )
                                continue;
                            kkk = ( kkk + cdim[2] ) % cdim[2];
                            cjd = cell_getid( cdim , iii , jjj , kkk );
                            cj = &s->cells[cjd];
                            if ( cid >= cjd || cj->count == 0 )
                                continue;
                            sid = sortlistID[ (kk+1) + 3*( (jj+1) + 3*(ii+1) ) ];
                            t = scheduler_addtask( sched , task_type_pair , task_subtype_density , sid , 0 , ci , cj , 1 );
                            }
                        }
                    }
                }

    /* Split the tasks. */
    scheduler_splittasks( sched );
    
    /* Count the number of tasks associated with each cell and
       store the density tasks in each cell, and make each sort
       depend on the sorts of its progeny. */
    // #pragma omp parallel for private(t,j)
    for ( k = 0 ; k < sched->nr_tasks ; k++ ) {
        t = &sched->tasks[k];
        if ( t->skip )
            continue;
        if ( t->type == task_type_sort && t->ci->split )
            for ( j = 0 ; j < 8 ; j++ ) {
                if ( t->ci->progeny[j] != NULL ) {
                    if ( t->ci->progeny[j]->sorts == NULL )
                        t->ci->progeny[j]->sorts = scheduler_addtask( sched , task_type_sort , task_subtype_none , t->flags , 0 , t->ci->progeny[j] , NULL , 0 );
                    t->ci->progeny[j]->sorts->skip = 0;
                    task_addunlock( t->ci->progeny[j]->sorts , t );
                    }
                }
        if ( t->type == task_type_self ) {
            atomic_inc( &t->ci->nr_tasks );
            if ( t->subtype == task_subtype_density ) {
                t->ci->density[ atomic_inc( &t->ci->nr_density ) ] = t;
                }
            }
        else if ( t->type == task_type_pair ) {
            atomic_inc( &t->ci->nr_tasks );
            atomic_inc( &t->cj->nr_tasks );
            if ( t->subtype == task_subtype_density ) {
                t->ci->density[ atomic_inc( &t->ci->nr_density ) ] = t;
                t->cj->density[ atomic_inc( &t->cj->nr_density ) ] = t;
                }
            }
        else if ( t->type == task_type_sub ) {
            atomic_inc( &t->ci->nr_tasks );
            if ( t->cj != NULL )
                atomic_inc( &t->cj->nr_tasks );
            if ( t->subtype == task_subtype_density ) {
                t->ci->density[ atomic_inc( &t->ci->nr_density ) ] = t;
                if ( t->cj != NULL )
                    t->cj->density[ atomic_inc( &t->cj->nr_density ) ] = t;
                }
            }
        }
        
    /* Append a ghost task to each cell. */
    space_map_cells_pre( s , 1 , &scheduler_map_mkghosts , sched );
    
    /* Run through the tasks and make force tasks for each density task.
       Each force task depends on the cell ghosts and unlocks the kick2 task
       of its super-cell. */
    kk = sched->nr_tasks;
    // #pragma omp parallel for private(t,t2)
    for ( k = 0 ; k < kk ; k++ ) {
    
        /* Get a pointer to the task. */
        t = &sched->tasks[k];
        
        /* Skip? */
        if ( t->skip )
            continue;
        
        /* Self-interaction? */
        if ( t->type == task_type_self && t->subtype == task_subtype_density ) {
            task_addunlock( t , t->ci->super->ghost );
            t2 = scheduler_addtask( sched , task_type_self , task_subtype_force , 0 , 0 , t->ci , NULL , 0 );
            task_addunlock( t->ci->ghost , t2 );
            task_addunlock( t2 , t->ci->super->kick2 );
            }
            
        /* Otherwise, pair interaction? */
        else if ( t->type == task_type_pair && t->subtype == task_subtype_density ) {
            task_addunlock( t , t->ci->super->ghost );
            if ( t->ci->super != t->cj->super )
                task_addunlock( t , t->cj->super->ghost );
            t2 = scheduler_addtask( sched , task_type_pair , task_subtype_force , 0 , 0 , t->ci , t->cj , 0 );
            task_addunlock( t->ci->ghost , t2 );
            task_addunlock( t->cj->ghost , t2 );
            task_addunlock( t2 , t->ci->super->kick2 );
            if ( t->ci->super != t->cj->super )
                task_addunlock( t2 , t->cj->super->kick2 );
            }
    
        /* Otherwise, sub interaction? */
        else if ( t->type == task_type_sub && t->subtype == task_subtype_density ) {
            task_addunlock( t , t->ci->super->ghost );
            if ( t->cj != NULL && t->ci->super != t->cj->super )
                task_addunlock( t , t->cj->super->ghost );
            t2 = scheduler_addtask( sched , task_type_sub , task_subtype_force , t->flags , 0 , t->ci , t->cj , 0 );
            task_addunlock( t->ci->ghost , t2 );
            if ( t->cj != NULL )
                task_addunlock( t->cj->ghost , t2 );
            task_addunlock( t2 , t->ci->super->kick2 );
            if ( t->cj != NULL && t->ci->super != t->cj->super )
                task_addunlock( t2 , t->cj->super->kick2 );
            }
            
        }
        
    /* Rank the tasks. */
    scheduler_ranktasks( sched );
            
    /* Count the number of each task type. */
    int counts[ task_type_count+1 ];
    for ( k = 0 ; k <= task_type_count ; k++ )
        counts[k] = 0;
    for ( k = 0 ; k < sched->nr_tasks ; k++ )
        if ( !sched->tasks[k].skip )
            counts[ (int)sched->tasks[k].type ] += 1;
        else
            counts[ task_type_count ] += 1;
    printf( "engine_maketasks: task counts are [ %s=%i" , taskID_names[0] , counts[0] );
    for ( k = 1 ; k < task_type_count ; k++ )
        printf( " %s=%i" , taskID_names[k] , counts[k] );
    printf( " skipped=%i ]\n" , counts[ task_type_count ] ); fflush(stdout); 
    
    }
    
    

/**
 * @brief Mark tasks to be skipped and set the sort flags accordingly.
 * 
 * @return 1 if the space has to be rebuilt, 0 otherwise.
 */
 
int engine_marktasks ( struct engine *e ) {

    struct scheduler *s = &e->sched;
    int k, nr_tasks = s->nr_tasks, *ind = s->tasks_ind;
    struct task *t, *tasks = s->tasks;
    float dt_step = e->dt_step;
    struct cell *ci, *cj;
    
    /* Run through the tasks and mark as skip or not. */
    for ( k = 0 ; k < nr_tasks ; k++ ) {
    
        /* Get a handle on the kth task. */
        t = &tasks[ ind[k] ];
        
        /* Sort-task? Note that due to the task ranking, the sorts
           will all come before the pairs and/or subs. */
        if ( t->type == task_type_sort ) {
        
            /* Re-set the flags. */
            t->flags = 0;
            t->skip = 1;
        
            }
        
        /* Single-cell task? */
        else if ( t->type == task_type_self ||
                  t->type == task_type_ghost ||
                ( t->type == task_type_sub && t->cj == NULL ) ) {
             
            /* Set this task's skip. */
            t->skip = ( t->ci->dt_min > dt_step );
            
            }
        
        /* Pair? */
        else if ( t->type == task_type_pair || ( t->type == task_type_sub && t->cj != NULL ) ) {
            
            /* Local pointers. */
            ci = t->ci;
            cj = t->cj;
            
            /* Set this task's skip. */
            t->skip = ( ci->dt_min > dt_step && cj->dt_min > dt_step );
            
            /* Too much particle movement? */
            if ( t->tight &&
                 ( fmaxf( ci->h_max , cj->h_max ) + ci->dx_max + cj->dx_max > cj->dmin || 
                   ci->dx_max > space_maxreldx*ci->h_max || cj->dx_max > space_maxreldx*cj->h_max ) )
                return 1;
                
            /* Set the sort flags. */
            if ( !t->skip && t->type == task_type_pair ) {
                ci->sorts->flags |= (1 << t->flags);
                ci->sorts->skip = 0;
                cj->sorts->flags |= (1 << t->flags);
                cj->sorts->skip = 0;
                }
                
            }
            
        /* Kick2? */
        else if ( t->type == task_type_kick2 )
            t->skip = 0;
            
        /* None? */
        else if ( t->type == task_type_none )
            t->skip = 1;
            
        }
        
    /* All is well... */
    return 0;
    
    }


/**
 * @brief Prepare the #engine by re-building the cells and tasks.
 *
 * @param e The #engine to prepare.
 */
 
void engine_prepare ( struct engine *e ) {
    
    int rebuild;
    
    TIMER_TIC

    /* Run through the tasks and mark as skip or not. */
    // tic = getticks();
    rebuild = ( e->step == 0 || engine_marktasks( e ) );
    // printf( "space_prepare: space_marktasks took %.3f ms.\n" , (double)(getticks() - tic)/CPU_TPS*1000 );
        
    /* Did this not go through? */
    if ( rebuild ) {
    
        /* Re-build the space. */
        // tic = getticks();
        space_rebuild( e->s , 0.0 );
        // printf( "engine_prepare: space_rebuild took %.3f ms.\n" , (double)(getticks() - tic)/CPU_TPS*1000 );
    
        /* Re-build the tasks. */
        // tic = getticks();
        engine_maketasks( e );
        // printf( "engine_prepare: engine_maketasks took %.3f ms.\n" , (double)(getticks() - tic)/CPU_TPS*1000 );
    
        /* Run through the tasks and mark as skip or not. */
        // tic = getticks();
        if ( engine_marktasks( e ) )
            error( "engine_marktasks failed after space_rebuild." );
        // printf( "engine_prepare: engine_marktasks took %.3f ms.\n" , (double)(getticks() - tic)/CPU_TPS*1000 );
        
        }

    /* Start the scheduler. */
    scheduler_start( &e->sched , (1 << task_type_sort) | 
                                 (1 << task_type_self) |
                                 (1 << task_type_pair) | 
                                 (1 << task_type_sub) |
                                 (1 << task_type_ghost) | 
                                 (1 << task_type_kick2) );
    
    TIMER_TOC( timer_prepare );
    
    }


/**
 * @brief Implements a barrier for the #runner threads.
 *
 * @param e The #engine.
 */
 
void engine_barrier ( struct engine *e ) {

    /* First, get the barrier mutex. */
    if ( pthread_mutex_lock( &e->barrier_mutex ) != 0 )
        error( "Failed to get barrier mutex." );
        
    /* This thread is no longer running. */
    e->barrier_running -= 1;
        
    /* If all threads are in, send a signal... */
    if ( e->barrier_running == 0 )
        if ( pthread_cond_broadcast( &e->barrier_cond ) != 0 )
            error( "Failed to broadcast barrier full condition." );
        
    /* Wait for the barrier to open. */
    while ( e->barrier_launch == 0 )
        if ( pthread_cond_wait( &e->barrier_cond , &e->barrier_mutex ) != 0 )
            error( "Eror waiting for barrier to close." );
        
    /* This thread has been launched. */
    e->barrier_running += 1;
    e->barrier_launch -= 1;
    
    /* If I'm the last one out, signal the condition again. */
    if ( e->barrier_launch == 0 )
        if ( pthread_cond_broadcast( &e->barrier_cond ) != 0 )
            error( "Failed to broadcast empty barrier condition." );
            
    /* Last but not least, release the mutex. */
    if ( pthread_mutex_unlock( &e->barrier_mutex ) != 0 )
        error( "Failed to get unlock the barrier mutex." );

    }
    
    
/**
 * @brief Mapping function to collect the data from the second kick.
 */

void engine_collect_kick2 ( struct cell *c ) {

    int k, updated = 0;
    float dt_min = FLT_MAX, dt_max = 0.0f;
    double ekin = 0.0, epot = 0.0;
    float mom[3] = { 0.0f , 0.0f , 0.0f }, ang[3] = { 0.0f , 0.0f , 0.0f };
    struct cell *cp;
    
    /* If I am a super-cell, return immediately. */
    if ( c->kick2 != NULL || c->count == 0 )
        return;
        
    /* If this cell is not split, I'm in trouble. */
    if ( !c->split )
        error( "Cell has no super-cell." );
        
    /* Collect the values from the progeny. */
    for ( k = 0 ; k < 8 ; k++ )
        if ( ( cp = c->progeny[k] ) != NULL ) {
            engine_collect_kick2( cp );
            dt_min = fminf( dt_min , cp->dt_min );
            dt_max = fmaxf( dt_max , cp->dt_max );
            updated += cp->updated;
            ekin += cp->ekin;
            epot += cp->epot;
            mom[0] += cp->mom[0]; mom[1] += cp->mom[1]; mom[2] += cp->mom[2];
            ang[0] += cp->ang[0]; ang[1] += cp->ang[1]; ang[2] += cp->ang[2];
            }
    
    /* Store the collected values in the cell. */
    c->dt_min = dt_min;
    c->dt_max = dt_max;
    c->updated = updated;
    c->ekin = ekin;
    c->epot = epot;
    c->mom[0] = mom[0]; c->mom[1] = mom[1]; c->mom[2] = mom[2];
    c->ang[0] = ang[0]; c->ang[1] = ang[1]; c->ang[2] = ang[2];
        
    }


/**
 * @brief Compute the force on a single particle brute-force.
 */

void engine_single_density ( double *dim , long long int pid , struct part *__restrict__ parts , int N , int periodic ) {

    int i, k;
    double r2, dx[3];
    float fdx[3], ih;
    struct part p;
    
    /* Find "our" part. */
    for ( k = 0 ; k < N && parts[k].id != pid ; k++ );
    if ( k == N )
        error( "Part not found." );
    p = parts[k];
    
    /* Clear accumulators. */
    ih = 1.0f / p.h;
    p.rho = 0.0f; p.rho_dh = 0.0f;
    p.density.wcount = 0.0f; p.density.wcount_dh = 0.0f;
	p.density.div_v = 0.0;
	for ( k=0 ; k < 3 ; k++)
		p.density.curl_v[k] = 0.0;
            
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
        if ( r2 < p.h*p.h*kernel_gamma2 ) {
            runner_iact_nonsym_density( r2 , fdx , p.h , parts[k].h , &p , &parts[k] );
            }
        }
        
    /* Dump the result. */
    p.rho = ih * ih * ih * ( p.rho + p.mass*kernel_root );
    p.rho_dh = p.rho_dh * ih * ih * ih * ih;
    p.density.wcount = ( p.density.wcount + kernel_root ) * ( 4.0f / 3.0 * M_PI * kernel_gamma3 );
    printf( "pairs_single_density: part %lli (h=%e) has wcount=%e, rho=%e, rho_dh=%e.\n" , p.id , p.h , p.density.wcount , p.rho , p.rho_dh );
    fflush(stdout);
    
    }


void engine_single_force ( double *dim , long long int pid , struct part *__restrict__ parts , int N , int periodic ) {

    int i, k;
    double r2, dx[3];
    float fdx[3];
    struct part p;
    
    /* Find "our" part. */
    for ( k = 0 ; k < N && parts[k].id != pid ; k++ );
    if ( k == N )
        error( "Part not found." );
    p = parts[k];
    
    /* Clear accumulators. */
    p.a[0] = 0.0f; p.a[1] = 0.0f; p.a[2] = 0.0f;
    p.force.u_dt = 0.0f; p.force.h_dt = 0.0f; p.force.v_sig = 0.0f;
            
    /* Loop over all particle pairs (force). */
    for ( k = 0 ; k < N ; k++ ) {
    // for ( k = N-1 ; k >= 0 ; k-- ) {
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
        if ( r2 < p.h*p.h*kernel_gamma2 || r2 < parts[k].h*parts[k].h*kernel_gamma2 ) {
            p.a[0] = 0.0f; p.a[1] = 0.0f; p.a[2] = 0.0f;
            p.force.u_dt = 0.0f; p.force.h_dt = 0.0f; p.force.v_sig = 0.0f;
            runner_iact_nonsym_force( r2 , fdx , p.h , parts[k].h , &p , &parts[k] );
            double dvdr = ( (p.v[0]-parts[k].v[0])*fdx[0] + (p.v[1]-parts[k].v[1])*fdx[1] + (p.v[2]-parts[k].v[2])*fdx[2] ) / sqrt(r2);
            printf( "pairs_single_force: part %lli and %lli interact (r=%.3e,dvdr=%.3e) with a=[%.3e,%.3e,%.3e], dudt=%.3e.\n" ,
                p.id , parts[k].id , sqrt(r2) , dvdr , p.a[0] , p.a[1], p.a[2] , p.force.u_dt );
            }
        }
        
    /* Dump the result. */
    // printf( "pairs_single_force: part %lli (h=%e) has a=[%.3e,%.3e,%.3e], udt=%e.\n" , p.id , p.h , p.a[0] , p.a[1] , p.a[2] , p.force.u_dt );
    fflush(stdout);
    
    }
    
    
/**
 * @breif Launch the runners.
 *
 * @param e The #engine.
 * @param nr_runners The number of #runners to let loose.
 */
 
void engine_launch ( struct engine *e , int nr_runners ) {

    /* Cry havoc and let loose the dogs of war. */
    e->barrier_launch = nr_runners;
    if ( pthread_cond_broadcast( &e->barrier_cond ) != 0 )
        error( "Failed to broadcast barrier open condition." );
        
    /* Sit back and wait for the runners to come home. */
    while ( e->barrier_launch || e->barrier_running )
        if ( pthread_cond_wait( &e->barrier_cond , &e->barrier_mutex ) != 0 )
            error( "Error while waiting for barrier." );
            
    }


/**
 * @brief Let the #engine loose to compute the forces.
 *
 * @param e The #engine.
 * @param sort_queues Flag to try to sort the queues topologically.
 */
 
void engine_step ( struct engine *e ) {

    int k;
    float dt = e->dt, dt_step, dt_max = 0.0f, dt_min = FLT_MAX;
    double epot = 0.0, ekin = 0.0;
    float mom[3] = { 0.0 , 0.0 , 0.0 };
    float ang[3] = { 0.0 , 0.0 , 0.0 };
    int count = 0;
    struct cell *c;
    struct space *s = e->s;
    
    TIMER_TIC2

    /* Get the maximum dt. */
    if ( e->policy & engine_policy_multistep ) {
        dt_step = 2.0f*dt;
        for ( k = 0 ; k < 32 && (e->step & (1 << k)) == 0 ; k++ )
            dt_step *= 2;
        }
    else
        dt_step = FLT_MAX;
        
    /* Set the maximum dt. */
    e->dt_step = dt_step;
    e->s->dt_step = dt_step;
    // printf( "engine_step: dt_step set to %.3e (dt=%.3e).\n" , dt_step , e->dt ); fflush(stdout);
    
    // printParticle( parts , 432626 );
    
    /* First kick. */
    TIMER_TIC
    scheduler_start( &e->sched , (1 << task_type_kick1) );
    engine_launch( e , ( e->nr_threads > 16 ) ? 16 : e->nr_threads );
    TIMER_TOC( timer_kick1 );
    
    /* Check if all the kick1 threads have executed. */
    for ( k = 0 ; k < e->sched.nr_tasks ; k++ )
        if ( e->sched.tasks[k].type == task_type_kick1 &&
             e->sched.tasks[k].tic == 0 )
            error( "Not all kick1 tasks completed." );
        
    // for(k=0; k<10; ++k)
    //   printParticle(parts, k);
    // printParticle( e->s->parts , 3392063069037 , e->s->nr_parts );
 
    /* Prepare the space. */
    engine_prepare( e );
    
    // engine_single_density( e->s->dim , 3392063069037 , e->s->parts , e->s->nr_parts , e->s->periodic );

    /* Send off the runners. */
    TIMER_TIC_ND
    engine_launch( e , e->nr_threads );
    TIMER_TOC(timer_runners);
            
    // engine_single_force( e->s->dim , 8328423931905 , e->s->parts , e->s->nr_parts , e->s->periodic );
    
    // for(k=0; k<10; ++k)
    //   printParticle(parts, k);
    // printParticle( parts , 432626 );
    // printParticle( e->s->parts , 3392063069037 , e->s->nr_parts );
    // printParticle( e->s->parts , 8328423931905 , e->s->nr_parts );

    /* Collect the cell data from the second kick. */
    for ( k = 0 ; k < s->nr_cells ; k++ ) {
        c = &s->cells[k];
        engine_collect_kick2( c );
        dt_min = fminf( dt_min , c->dt_min );
        dt_max = fmaxf( dt_max , c->dt_max );
        ekin += c->ekin;
        epot += c->epot;
        count += c->updated;
        mom[0] += c->mom[0]; mom[1] += c->mom[1]; mom[2] += c->mom[2];
        ang[0] += c->ang[0]; ang[1] += c->ang[1]; ang[2] += c->ang[2];
        }
    
    e->dt_min = dt_min;
    e->dt_max = dt_max;
    e->count_step = count;
    e->ekin = ekin;
    e->epot = epot;
    // printParticle( e->s->parts , 382557 , e->s->nr_parts );
    // printf( "engine_step: dt_min/dt_max is %e/%e.\n" , dt_min , dt_max ); fflush(stdout);
    // printf( "engine_step: etot is %e (ekin=%e, epot=%e).\n" , ekin+epot , ekin , epot ); fflush(stdout);
    // printf( "engine_step: total momentum is [ %e , %e , %e ].\n" , mom[0] , mom[1] , mom[2] ); fflush(stdout);
    // printf( "engine_step: total angular momentum is [ %e , %e , %e ].\n" , ang[0] , ang[1] , ang[2] ); fflush(stdout);
    // printf( "engine_step: total entropic function is %e .\n", ent ); fflush(stdout);
    // printf( "engine_step: updated %i parts (dt_step=%.3e).\n" , count , dt_step ); fflush(stdout);
        
    /* Increase the step. */
    e->step += 1;

    /* Does the time step need adjusting? */
    if ( e->policy & engine_policy_fixdt ) {
        e->dt = e->dt_orig;
        }
    else {
        if ( e->dt == 0 ) {
            e->nullstep += 1;
            e->dt = e->dt_orig;
            while ( dt_min < e->dt )
                e->dt *= 0.5;
            while ( dt_min > 2*e->dt )
                e->dt *= 2.0;
            // printf( "engine_step: dt_min=%.3e, adjusting time step to dt=%e.\n" , dt_min , e->dt );
            }
        else {
            while ( dt_min < e->dt ) {
                e->dt *= 0.5;
                e->step *= 2;
                e->nullstep *= 2;
                // printf( "engine_step: dt_min dropped below time step, adjusting to dt=%e.\n" , e->dt );
                }
            while ( dt_min > 2*e->dt && (e->step & 1) == 0 ) {
                e->dt *= 2.0;
                e->step /= 2;
                e->nullstep /= 2;
                // printf( "engine_step: dt_min is larger than twice the time step, adjusting to dt=%e.\n" , e->dt );
                }
            }
        } 
    
    /* Set the system time. */
    e->time = e->dt * (e->step - e->nullstep);
        
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
    float dt_min = dt;
    
    /* Store the values. */
    e->s = s;
    e->nr_threads = nr_threads;
    e->policy = policy;
    e->step = 0;
    e->nullstep = 0;
    e->time = 0.0;
    
    /* First of all, init the barrier and lock it. */
    if ( pthread_mutex_init( &e->barrier_mutex , NULL ) != 0 )
        error( "Failed to initialize barrier mutex." );
    if ( pthread_cond_init( &e->barrier_cond , NULL ) != 0 )
        error( "Failed to initialize barrier condition variable." );
    if ( pthread_mutex_lock( &e->barrier_mutex ) != 0 )
        error( "Failed to lock barrier mutex." );
    e->barrier_running = 0;
    e->barrier_launch = 0;
    
    /* Run through the parts and get the minimum time step. */
    e->dt_orig = dt;
    for ( k = 0 ; k < s->nr_parts ; k++ )
        if ( s->parts[k].dt < dt_min )
            dt_min = s->parts[k].dt;
    if ( dt_min == 0.0f )
        dt = 0.0f;
    else
        while ( dt > dt_min )
            dt *= 0.5f;
    e->dt = dt;
    
    /* Init the scheduler. */
    scheduler_init( &e->sched , e->s , nr_queues , scheduler_flag_steal );
        
    /* Append a kick1 task to each cell. */
    scheduler_reset( &e->sched , s->tot_cells );
    space_map_cells_pre( e->s , 1 , &scheduler_map_mkkick1 , &e->sched );
    
    /* Allocate and init the threads. */
    if ( ( e->runners = (struct runner *)malloc( sizeof(struct runner) * nr_threads ) ) == NULL )
        error( "Failed to allocate threads array." );
    for ( k = 0 ; k < nr_threads ; k++ ) {
        e->runners[k].id = k;
        e->runners[k].e = e;
        e->barrier_running += 1;
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
    while ( e->barrier_running || e->barrier_launch )
        if ( pthread_cond_wait( &e->barrier_cond , &e->barrier_mutex ) != 0 )
            error( "Error while waiting for runner threads to get in place." );
    
    }
    
    
    
