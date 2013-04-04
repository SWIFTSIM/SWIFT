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
#include "lock.h"
#include "task.h"
#include "part.h"
#include "cell.h"
#include "space.h"
#include "queue.h"
#include "engine.h"
#include "runner.h"
#include "runner_iact.h"
#include "error.h"

/* Convert cell location to ID. */
#define cell_getid( cdim , i , j , k ) ( (int)(k) + (cdim)[2]*( (int)(j) + (cdim)[1]*(int)(i) ) )

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
 * @param flags Cell flag.
 * @param clock Flag indicating whether to record the timing or not, needed
 *      for recursive calls.
 */
 
void runner_dosort ( struct runner *r , struct cell *c , int flags , int clock ) {

    struct entry *finger;
    struct entry *fingers[8];
    struct cpart *cparts = c->cparts;
    int j, k, count = c->count;
    int i, ind, off[8], inds[8], temp_i, missing;
    // float shift[3];
    float buff[8], px[3];
    
    TIMER_TIC
    
    /* Clean-up the flags, i.e. filter out what's already been sorted. */
    flags &= ~c->sorted;
    if ( flags == 0 )
        return;
    
    /* start by allocating the entry arrays. */
    if ( lock_lock( &c->mlock ) != 0 )
        error( "Failed to lock cell." );
    if ( c->sort == NULL || c->sortsize < c->count ) {
        if ( c->sort != NULL )
            free( c->sort );
        c->sortsize = c->count * 1.1;
        if ( ( c->sort = (struct entry *)malloc( sizeof(struct entry) * (c->sortsize + 1) * 13 ) ) == NULL )
            error( "Failed to allocate sort memory." );
        }
    if ( lock_unlock( &c->mlock ) != 0 )
        error( "Failed to unlock cell." );
        
    /* Does this cell have any progeny? */
    if ( c->split ) {
    
        /* Fill in the gaps within the progeny. */
        for ( k = 0 ; k < 8 ; k++ ) {
            if ( c->progeny[k] == NULL )
                continue;
            if ( c->progeny[k]->sorts == NULL )
                missing = flags;
            else
                missing = ( c->progeny[k]->sorts->flags ^ flags ) & flags;
            if ( missing )
                runner_dosort( r , c->progeny[k] , missing , 0 );
            }
    
        /* Loop over the 13 different sort arrays. */
        for ( j = 0 ; j < 13 ; j++ ) {
        
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
            
            /* Mark as sorted. */
            c->sorted |= ( 1 << j );
            
            } /* loop over sort arrays. */
    
        } /* progeny? */
        
    /* Otherwise, just sort. */
    else {
    
        /* Fill the sort array. */
        for ( k = 0 ; k < count ; k++ ) {
            px[0] = cparts[k].x[0];
            px[1] = cparts[k].x[1];
            px[2] = cparts[k].x[2];
            for ( j = 0 ; j < 13 ; j++ )
                if ( flags & (1 << j) ) {
                    c->sort[ j*(count + 1) + k].i = k;
                    c->sort[ j*(count + 1) + k].d = px[0]*runner_shift[ 3*j + 0 ] + px[1]*runner_shift[ 3*j + 1 ] + px[2]*runner_shift[ 3*j + 2 ];
                    }
            }

        /* Add the sentinel and sort. */
        for ( j = 0 ; j < 13 ; j++ )
            if ( flags & (1 << j) ) {
                c->sort[ j*(count + 1) + c->count ].d = FLT_MAX;
                c->sort[ j*(count + 1) + c->count ].i = 0;
                runner_dosort_ascending( &c->sort[ j*(count + 1) ] , c->count );
                c->sorted |= ( 1 << j );
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
            ((double)TIMER_TOC(timer_dosort)) / CPU_TPS * 1000 ); fflush(stdout);
    #else
        if ( clock )
            TIMER_TOC(timer_dosort);
    #endif

    }
    
    
/**
 * @brief Intermediate task between density and force
 *
 * @param r The runner thread.
 * @param c The cell.
 */
 
void runner_doghost ( struct runner *r , struct cell *c ) {

    struct part *p;
    struct cpart *cp;
    struct cell *finger;
    int i, k, redo, count = c->count;
    int *pid;
    float h, hg, ihg, ihg2, ihg4, h_corr;
    float normDiv_v, normCurl_v;
    float dt_step = r->e->dt_step;
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
            cp = &c->cparts[ pid[i] ];
            
            /* Is this part within the timestep? */
            if ( cp->dt <= dt_step ) {
            
  	            /* Some smoothing length multiples. */
	            h = p->h;
                hg = kernel_gamma * h;
                ihg = 1.0f / hg;
                ihg2 = ihg * ihg;
                ihg4 = ihg2 * ihg2;

                /* Adjust the computed weighted number of neighbours. */
                p->density.wcount += kernel_wroot;

		        /* Final operation on the density. */
                p->rho = ihg * ihg2 * ( p->rho + p->mass*kernel_root );
                p->rho_dh *= ihg4;
                    
                /* Compute the smoothing length update (Newton step). */
                if ( p->density.wcount_dh == 0.0f )
                    h_corr = p->h;
                else {
                    h_corr = ( const_nwneigh - p->density.wcount ) / p->density.wcount_dh;

                    /* Truncate to the range [ -p->h/2 , p->h ]. */
                    h_corr = fminf( h_corr , h );
                    h_corr = fmaxf( h_corr , -h/2 );
                    }
                
                /* Apply the correction to p->h. */
                p->h += h_corr;
                cp->h = p->h;

                /* Did we get the right number density? */
                if ( p->density.wcount > const_nwneigh + const_delta_nwneigh ||
                     p->density.wcount < const_nwneigh - const_delta_nwneigh ) {
                    // printf( "runner_doghost: particle %lli (h=%e,depth=%i) has bad wcount=%.3f.\n" , p->id , p->h , c->depth , p->density.wcount ); fflush(stdout);
                    // p->h += ( p->density.wcount + kernel_root - const_nwneigh ) / p->density.wcount_dh;
                    pid[redo] = pid[i];
                    redo += 1;
                    p->density.wcount = 0.0;
                    p->density.wcount_dh = 0.0;
                    p->rho = 0.0;
                    p->rho_dh = 0.0;
		            p->density.div_v = 0.0;
		            for ( k=0 ; k < 3 ; k++)
		                p->density.curl_v[k] = 0.0;
                    continue;
                    }

                /* Pre-compute some stuff for the balsara switch. */
		        normDiv_v = fabs( p->density.div_v / p->rho * ihg4 );
		        normCurl_v = sqrtf( p->density.curl_v[0] * p->density.curl_v[0] + p->density.curl_v[1] * p->density.curl_v[1] + p->density.curl_v[2] * p->density.curl_v[2] ) / p->rho * ihg4;
                
                /* As of here, particle force variables will be set. Do _NOT_
                   try to read any particle density variables! */
                
                /* Compute this particle's sound speed. */
                p->force.c = sqrtf( const_gamma * ( const_gamma - 1.0f ) * p->u );

                /* Compute the P/Omega/rho2. */
                p->force.POrho2 = p->u * ( const_gamma - 1.0f ) / ( p->rho + hg * p->rho_dh / 3.0f );

		        /* Balsara switch */
		        p->force.balsara = normCurl_v / ( normDiv_v + normCurl_v + 0.0001f * p->force.c * ihg );

                /* Reset the acceleration. */
                for ( k = 0 ; k < 3 ; k++ )
                    p->a[k] = 0.0f;

                /* Reset the time derivatives. */
                p->force.u_dt = 0.0f;
                p->force.h_dt = 0.0f;
                p->force.v_sig = 0.0f;

                }

            }
            
        /* Re-set the counter for the next loop (potentially). */
        count = redo;
        if ( count > 0 ) {
        
            // error( "Bad smoothing length, fixing this isn't implemented yet." );
            
            /* Climb up the cell hierarchy. */
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
                    else if ( finger->density[k]->type == task_type_sub ) {
                    
                        /* Left or right? */
                        if ( finger->density[k]->ci == finger )
                            runner_dosub_subset_density( r , finger , c->parts , pid , count , finger->density[k]->cj , -1 );
                        else
                            runner_dosub_subset_density( r , finger , c->parts , pid , count , finger->density[k]->ci , -1 );
                        
                        }
                
                    }
                    
                }
        
            }
            
        }

    #ifdef TIMER_VERBOSE
        printf( "runner_doghost[%02i]: %i parts at depth %i took %.3f ms.\n" ,
            r->id , c->count , c->depth ,
            ((double)TIMER_TOC(timer_doghost)) / CPU_TPS * 1000 ); fflush(stdout);
    #else
        TIMER_TOC(timer_doghost);
    #endif
    
    }
    
    
/**
 * @brief Compute the second kick of the given cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 */
 
void runner_dokick2 ( struct runner *r , struct cell *c ) {

    int k, count = 0, nr_parts = c->count;
    float dt_min = FLT_MAX, dt_max = 0.0f, ekin = 0.0f, epot = 0.0f;
    float mom[3] = { 0.0f , 0.0f , 0.0f }, ang[3] = { 0.0f , 0.0f , 0.0f };
    float x[3], v[3], u, h, pdt, m;
    float dt_step = r->e->dt_step, dt = r->e->dt, hdt = 0.5f*dt;
    struct part *p, *parts = c->parts;
    struct xpart *xp;
    
    TIMER_TIC
    
    /* Loop over the particles and kick them. */
    for ( k = 0 ; k < nr_parts ; k++ ) {

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
            count += 1;
            }

        /* Update the particle's time step. */
        p->dt = pdt = const_cfl * h / ( p->force.v_sig );

        /* Update positions and energies at the half-step. */
        p->v[0] = v[0] = xp->v_old[0] + hdt * p->a[0];
        p->v[1] = v[1] = xp->v_old[1] + hdt * p->a[1];
        p->v[2] = v[2] = xp->v_old[2] + hdt * p->a[2];
        p->u = u = xp->u_old + hdt * p->force.u_dt;

        /* Get the smallest/largest dt. */
        dt_min = fminf( dt_min , pdt );
        dt_max = fmaxf( dt_max , pdt );

        /* Collect total energy. */
        ekin += 0.5 * m * ( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
        epot += m * u;

        /* Collect momentum */
        mom[0] += m * p->v[0];
        mom[1] += m * p->v[1];
        mom[2] += m * p->v[2];

	    /* Collect angular momentum */
	    ang[0] += m * ( x[1]*v[2] - v[2]*x[1] );
	    ang[1] += m * ( x[2]*v[0] - v[0]*x[2] );
	    ang[2] += m * ( x[0]*v[1] - v[1]*x[0] );

	    /* Collect entropic function */
	    // lent += u * pow( p->rho, 1.f-const_gamma );

        }

    #ifdef TIMER_VERBOSE
        printf( "runner_dokick2[%02i]: %i parts at depth %i took %.3f ms.\n" ,
            r->id , c->count , c->depth ,
            ((double)TIMER_TOC(timer_kick2)) / CPU_TPS * 1000 ); fflush(stdout);
    #else
        TIMER_TOC(timer_kick2);
    #endif
        
    /* Store the computed values in the cell. */
    c->dt_min = dt_min;
    c->dt_max = dt_max;
    c->updated = count;
    c->ekin = ekin;
    c->epot = epot;
    c->mom[0] = mom[0]; c->mom[1] = mom[1]; c->mom[2] = mom[2];
    c->ang[0] = ang[0]; c->ang[1] = ang[1]; c->ang[2] = ang[2];
        
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
        myq = &e->queues[ threadID * e->nr_queues / e->nr_threads ];
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
                    TIMER_TOC2(timer_steal);
                    }
                }
            else if ( e->policy & engine_policy_rand ) {
                qid = rand_r( &myseed ) % naq;
                t = queue_gettask( queues[qid] , r->id , e->policy & engine_policy_block , 0 );
                }
            else {
                t = queue_gettask( &e->queues[threadID] , r->id , e->policy & engine_policy_block , 0 );
                }
            TIMER_TOC(timer_getpair);
            
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
                timers_toc( timer_stalled , stalled );
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
                        runner_doself1_density( r , ci );
                    else if ( t->subtype == task_subtype_force )
                        runner_doself2_force( r , ci );
                    else
                        error( "Unknown task subtype." );
                    cell_unlocktree( ci );
                    break;
                case task_type_pair:
                    if ( t->subtype == task_subtype_density )
                        runner_dopair1_density( r , ci , cj );
                    else if ( t->subtype == task_subtype_force )
                        runner_dopair2_force( r , ci , cj );
                    else
                        error( "Unknown task subtype." );
                    cell_unlocktree( ci );
                    cell_unlocktree( cj );
                    break;
                case task_type_sort:
                    runner_dosort( r , ci , t->flags , 1 );
                    break;
                case task_type_sub:
                    if ( t->subtype == task_subtype_density )
                        runner_dosub1_density( r , ci , cj , t->flags );
                    else if ( t->subtype == task_subtype_force )
                        runner_dosub2_force( r , ci , cj , t->flags );
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
                case task_type_kick2:
                    runner_dokick2( r , ci );
                    break;
                default:
                    error( "Unknown task type." );
                }
                
            /* Resolve any dependencies. */
            for ( k = 0 ; k < t->nr_unlock_tasks ; k++ )
                if ( atomic_dec( &t->unlock_tasks[k]->wait ) == 0 )
                    error( "Task negative wait." );
        
            } /* main loop. */
            
    	/* Any leftover stalls? */    
        #ifdef TIMER
        if ( stalled ) {
            timers_toc( timer_stalled , stalled );
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
    

