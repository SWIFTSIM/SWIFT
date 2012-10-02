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
#include <float.h>
#include <limits.h>
#include <math.h>

/* Local headers. */
#include "cycle.h"
#include "lock.h"
#include "space.h"
#include "runner.h"

/* Error macro. */
#define error(s) { printf( "%s:%s:%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }

/* Convert cell location to ID. */
#define cell_getid( cdim , i , j , k ) ( (int)(k) + (cdim)[2]*( (int)(j) + (cdim)[1]*(int)(i) ) )

/* Split size. */
int space_splitsize = space_splitsize_default;

/* Task type names. */
const char *taskID_names[tid_count] = { "none" , "sort" , "self" , "pair" , "sub" };

/* Map shift vector to sortlist. */
const int sortlistID[27] = {
    /* ( -1 , -1 , -1 ) */   0 ,
    /* ( -1 , -1 ,  0 ) */   1 , 
    /* ( -1 , -1 ,  1 ) */   2 ,
    /* ( -1 ,  0 , -1 ) */   3 ,
    /* ( -1 ,  0 ,  0 ) */   4 , 
    /* ( -1 ,  0 ,  1 ) */   5 ,
    /* ( -1 ,  1 , -1 ) */   6 ,
    /* ( -1 ,  1 ,  0 ) */   7 , 
    /* ( -1 ,  1 ,  1 ) */   8 ,
    /* (  0 , -1 , -1 ) */   9 ,
    /* (  0 , -1 ,  0 ) */   10 , 
    /* (  0 , -1 ,  1 ) */   11 ,
    /* (  0 ,  0 , -1 ) */   12 ,
    /* (  0 ,  0 ,  0 ) */   0 , 
    /* (  0 ,  0 ,  1 ) */   12 ,
    /* (  0 ,  1 , -1 ) */   11 ,
    /* (  0 ,  1 ,  0 ) */   10 , 
    /* (  0 ,  1 ,  1 ) */   9 ,
    /* (  1 , -1 , -1 ) */   8 ,
    /* (  1 , -1 ,  0 ) */   7 , 
    /* (  1 , -1 ,  1 ) */   6 ,
    /* (  1 ,  0 , -1 ) */   5 ,
    /* (  1 ,  0 ,  0 ) */   4 , 
    /* (  1 ,  0 ,  1 ) */   3 ,
    /* (  1 ,  1 , -1 ) */   2 ,
    /* (  1 ,  1 ,  0 ) */   1 , 
    /* (  1 ,  1 ,  1 ) */   0 
    };
    
    
/**
 * @brief Mapping function to draw a specific cell (gnuplot).
 */

void space_map_clearsort ( struct cell *c , void *data ) {

    if ( c->sort != NULL ) {
        free( c->sort );
        c->sort = NULL;
        }

    }


/**
 * @brief Get a task free of dependencies and conflicts.
 *
 * @param s The #space.
 */
 
struct task *space_gettask ( struct space *s ) {

    int k, tid = -1;
    struct task *res = NULL;
    struct cell *c;

    /* Main loop, while there are tasks... */
    while ( s->next_task < s->nr_tasks ) {
    
        /* Grab the task lock. */
        if ( lock_lock( &s->task_lock ) != 0 )
            error( "Locking the task_lock failed.\n" );
            
        /* Loop over the remaining task IDs. */
        for ( k = s->next_task ; k < s->nr_tasks ; k++ ) {
        
            /* Put a finger on the task. */
            res = &s->tasks[ s->tasks_ind[k] ];
            
            /* Is this task blocked? */
            if ( res->wait )
                continue;
            
            /* Different criteria for different types. */
            switch ( res->type ) {
                case tid_self:
                    if ( res->ci->lock || res->ci->wait > 0 )
                        continue;
                    break;
                case tid_pair:
                    if ( res->ci->lock || res->cj->lock || res->ci->wait || res->cj->wait )
                        continue;
                    break;
                case tid_sort:
                    if ( res->ci->lock )
                        continue;
                    break;
                }
            
            /* If we made it this far, we're safe. */
            break;
        
            } /* loop over the task IDs. */
            
        /* Did we get a task? */
        if ( k < s->nr_tasks ) {
        
            // /* Swap to front. */
            // tid = s->tasks_ind[k];
            // s->tasks_ind[k] = s->tasks_ind[ s->next_task ];
            // s->tasks_ind[ s->next_task ] = tid;
            
            /* Bubble-down the task. */
            tid = s->tasks_ind[k];
            while ( k > s->next_task ) {
                s->tasks_ind[ k ] = s->tasks_ind[ k-1 ];
                k -= 1;
                }
            s->tasks_ind[ s->next_task ] = tid;
            
            /* Lock the cells, if needed. */
            if ( s->tasks[tid].type != tid_sort ) {
                for ( c = res->ci ; c != NULL ; c = c->parent )
                    __sync_fetch_and_add( &c->lock , 1 );
                for ( c = res->cj ; c != NULL ; c = c->parent )
                    __sync_fetch_and_add( &c->lock , 1 );
                }
            
            /* up the counter. */
            s->next_task += 1;
        
            }
    
        /* Release the task lock. */
        if ( lock_unlock( &s->task_lock ) != 0 )
            error( "Locking the task_lock failed.\n" );
            
        /* Leave? */
        if ( tid >= 0 )
            return &s->tasks[tid];
    
        } /* while there are tasks. */
        
    /* No beef. */
    return NULL;

    }


/**
 * @brief Map a function to all particles in a aspace.
 *
 * @param s The #space we are working in.
 */
 
void space_map_parts ( struct space *s , void (*fun)( struct part *p , struct cell *c , void *data ) , void *data ) {

    int i;

    void rec_map ( struct cell *c ) {
    
        int k;
        
        /* No progeny? */
        if ( !c->split )
            for ( k = 0 ; k < c->count ; k++ )
                fun( &c->parts[k] , c , data );
                
        /* Otherwise, recurse. */
        else
            for ( k = 0 ; k < 8 ; k++ )
                if ( c->progeny[k] != NULL )
                    rec_map( c->progeny[k] );
                
        }
        
    /* Call the recursive function on all higher-level cells. */
    for ( i = 0 ; i < s->nr_cells ; i++ )
        rec_map( &s->cells[i] );

    }


/**
 * @brief Map a function to all particles in a aspace.
 *
 * @param s The #space we are working in.
 */
 
void space_map_cells ( struct space *s , void (*fun)( struct cell *c , void *data ) , void *data ) {

    int i;

    void rec_map ( struct cell *c ) {
    
        int k;
        
        /* No progeny? */
        if ( !c->split )
            fun( c , data );
                
        /* Otherwise, recurse. */
        else
            for ( k = 0 ; k < 8 ; k++ )
                if ( c->progeny[k] != NULL )
                    rec_map( c->progeny[k] );
                
        }
        
    /* Call the recursive function on all higher-level cells. */
    for ( i = 0 ; i < s->nr_cells ; i++ )
        rec_map( &s->cells[i] );

    }


/**
 * @brief Add a #task to the #space.
 *
 * @param s The #space we are working in.
 */
 
struct task *space_addtask ( struct space *s , int type , int flags , int wait , struct cell *ci , struct cell *cj , struct task *unlock_tasks[] , int nr_unlock_tasks , struct cell *unlock_cells[] , int nr_unlock_cells ) {

    struct task *t = &s->tasks[ s->nr_tasks ];
    
    /* Copy the data. */
    t->type = type;
    t->flags = flags;
    t->wait = wait;
    t->ci = ci;
    t->cj = cj;
    if ( unlock_tasks != NULL )
        memcpy( t->unlock_tasks , unlock_tasks , sizeof(struct task *) * nr_unlock_tasks );
    t->nr_unlock_tasks = nr_unlock_tasks;
    if ( unlock_cells != NULL )
        memcpy( t->unlock_cells , unlock_cells , sizeof(struct task *) * nr_unlock_cells );
    t->nr_unlock_cells = nr_unlock_cells;
    
    /* Increase the task counter. */
    s->nr_tasks += 1;
    
    /* Return a pointer to the new task. */
    return t;

    }



/**
 * @brief Remove an unlock_task from the given task.
 *
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 */
 
void task_rmunlock( struct task *ta , struct task *tb ) {

    int k;
    
    for ( k = 0 ; k < ta->nr_unlock_tasks ; k++ )
        if ( ta->unlock_tasks[k] == tb ) {
            ta->nr_unlock_tasks -= 1;
            ta->unlock_tasks[k] = ta->unlock_tasks[ ta->nr_unlock_tasks ];
            return;
            }
    error( "Task not found." );

    }
    

/**
 * @brief Add an unlock_task to the given task.
 *
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 */
 
void task_addunlock( struct task *ta , struct task *tb ) {

    int k;
    
    /* Bogus? */
    if ( ta == NULL || tb == NULL )
        return;
    
    /* Check if ta already unlocks tb. */
    for ( k = 0 ; k < ta->nr_unlock_tasks ; k++ )
        if ( ta->unlock_tasks[k] == tb )
            return;

    if ( ta->nr_unlock_tasks == task_maxunlock )
        error( "Too many unlock_tasks in task." );
        
    ta->unlock_tasks[ ta->nr_unlock_tasks] = tb;
    ta->nr_unlock_tasks += 1;

    }
    

/**
 * @brief Split tasks that may be too large.
 *
 * @param s The #space we are working in.
 */
 
void space_splittasks ( struct space *s ) {

    int k, sid, tid;
    struct cell *ci, *cj;
    double hi, hj, shift[3];
    struct task *t;

    /* Loop through the tasks... */
    for ( tid = 0 ; tid < s->nr_tasks ; tid++ ) {
    
        /* Get a pointer on the task. */
        t = &s->tasks[tid];
    
        /* If this task isn't a pair, i'm not interested. */
        if ( t->type != tid_pair )
            continue;
            
        /* Get a handle on the cells involved. */
        ci = t->ci;
        cj = t->cj;
        hi = fmax( ci->h[0] , fmax( ci->h[1] , ci->h[2] ) );
        hj = fmax( cj->h[0] , fmax( cj->h[1] , cj->h[2] ) );
            
        /* Should this task be split-up? */
        if ( ci->split && cj->split && ci->r_max < hi/2 && cj->r_max < hj/2 ) {
        
            /* Get the relative distance between the pairs, wrapping. */
            for ( k = 0 ; k < 3 ; k++ ) {
                if ( cj->loc[k] - ci->loc[k] < -s->dim[k]/2 )
                    shift[k] = s->dim[k];
                else if ( cj->loc[k] - ci->loc[k] > s->dim[k]/2 )
                    shift[k] = -s->dim[k];
                else
                    shift[k] = 0.0;
                }

            /* Get the sorting index. */
            for ( sid = 0 , k = 0 ; k < 3 ; k++ )
                sid = 3*sid + ( (cj->loc[k] - ci->loc[k] + shift[k] < 0) ? 0 : (cj->loc[k] - ci->loc[k] + shift[k] > 0) ? 2 : 1 );

            /* Flip? */
            if ( sid < 13 ) {
                cj = t->ci;
                ci = t->cj;
                t->ci = ci; t->cj = cj;
                }
            else
                sid = 26 - sid;
                
            /* Remove the dependency of this task on the sorts of ci and cj. */
            task_rmunlock( ci->sorts[sid] , t );
            task_rmunlock( cj->sorts[sid] , t );
            ci->nr_pairs -= 1;
            cj->nr_pairs -= 1;
            t->nr_unlock_cells = 0;
                
            /* For each different sorting type... */
            switch ( sid ) {
            
                case 0: /* (  1 ,  1 ,  1 ) */
                    t->ci = ci->progeny[7]; t->cj = cj->progeny[0];
                    task_addunlock( ci->progeny[7]->sorts[0] , t ); task_addunlock( cj->progeny[0]->sorts[0] , t );
                    ci->progeny[7]->nr_pairs += 1;
                    cj->progeny[0]->nr_pairs += 1;
                    break;
                    
                case 1: /* (  1 ,  1 ,  0 ) */
                    if ( space_dosub &&
                         !ci->progeny[6]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[1]->split ) {
                        t->type = tid_sub; t->flags = 1;
                        task_addunlock( ci->progeny[6]->sorts[1] , t ); task_addunlock( cj->progeny[0]->sorts[1] , t );
                        task_addunlock( ci->progeny[7]->sorts[1] , t ); task_addunlock( cj->progeny[1]->sorts[1] , t );
                        task_addunlock( ci->progeny[6]->sorts[0] , t ); task_addunlock( cj->progeny[1]->sorts[0] , t );
                        task_addunlock( ci->progeny[7]->sorts[2] , t ); task_addunlock( cj->progeny[0]->sorts[2] , t );
                        }
                    else {
                        t->ci = ci->progeny[6]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[6]->sorts[1] , t ); task_addunlock( cj->progeny[0]->sorts[1] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[1] , t ); task_addunlock( cj->progeny[1]->sorts[1] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[0] , t ); task_addunlock( cj->progeny[1]->sorts[0] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[2] , t ); task_addunlock( cj->progeny[0]->sorts[2] , t );
                        }
                    ci->progeny[6]->nr_pairs += 2;
                    ci->progeny[7]->nr_pairs += 2;
                    cj->progeny[0]->nr_pairs += 2;
                    cj->progeny[1]->nr_pairs += 2;
                    break;
                    
                case 2: /* (  1 ,  1 , -1 ) */
                    t->ci = ci->progeny[6]; t->cj = cj->progeny[1];
                    task_addunlock( ci->progeny[6]->sorts[2] , t ); task_addunlock( cj->progeny[1]->sorts[2] , t );
                    ci->progeny[6]->nr_pairs += 1;
                    cj->progeny[1]->nr_pairs += 1;
                    break;
                    
                case 3: /* (  1 ,  0 ,  1 ) */
                    if ( space_dosub &&
                         !ci->progeny[5]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[2]->split ) {
                        t->type = tid_sub; t->flags = 3;
                        task_addunlock( ci->progeny[5]->sorts[3] , t ); task_addunlock( cj->progeny[0]->sorts[3] , t );
                        task_addunlock( ci->progeny[7]->sorts[3] , t ); task_addunlock( cj->progeny[2]->sorts[3] , t );
                        task_addunlock( ci->progeny[5]->sorts[0] , t ); task_addunlock( cj->progeny[2]->sorts[0] , t );
                        task_addunlock( ci->progeny[7]->sorts[6] , t ); task_addunlock( cj->progeny[0]->sorts[6] , t );
                        }
                    else {
                        t->ci = ci->progeny[5]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[5]->sorts[3] , t ); task_addunlock( cj->progeny[0]->sorts[3] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[3] , t ); task_addunlock( cj->progeny[2]->sorts[3] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[0] , t ); task_addunlock( cj->progeny[2]->sorts[0] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[6] , t ); task_addunlock( cj->progeny[0]->sorts[6] , t );
                        }
                    ci->progeny[5]->nr_pairs += 2;
                    ci->progeny[7]->nr_pairs += 2;
                    cj->progeny[0]->nr_pairs += 2;
                    cj->progeny[2]->nr_pairs += 2;
                    break;
                    
                case 4: /* (  1 ,  0 ,  0 ) */
                    if ( space_dosub &&
                         !ci->progeny[4]->split && !ci->progeny[5]->split && !ci->progeny[6]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[1]->split && !cj->progeny[2]->split && !cj->progeny[3]->split ) {
                        t->type = tid_sub; t->flags = 4;
                        task_addunlock( ci->progeny[4]->sorts[4] , t ); task_addunlock( cj->progeny[0]->sorts[4] , t );
                        task_addunlock( ci->progeny[5]->sorts[5] , t ); task_addunlock( cj->progeny[0]->sorts[5] , t );
                        task_addunlock( ci->progeny[6]->sorts[7] , t ); task_addunlock( cj->progeny[0]->sorts[7] , t );
                        task_addunlock( ci->progeny[7]->sorts[8] , t ); task_addunlock( cj->progeny[0]->sorts[8] , t );
                        task_addunlock( ci->progeny[4]->sorts[3] , t ); task_addunlock( cj->progeny[1]->sorts[3] , t );
                        task_addunlock( ci->progeny[5]->sorts[4] , t ); task_addunlock( cj->progeny[1]->sorts[4] , t );
                        task_addunlock( ci->progeny[6]->sorts[6] , t ); task_addunlock( cj->progeny[1]->sorts[6] , t );
                        task_addunlock( ci->progeny[7]->sorts[7] , t ); task_addunlock( cj->progeny[1]->sorts[7] , t );
                        task_addunlock( ci->progeny[4]->sorts[1] , t ); task_addunlock( cj->progeny[2]->sorts[1] , t );
                        task_addunlock( ci->progeny[5]->sorts[2] , t ); task_addunlock( cj->progeny[2]->sorts[2] , t );
                        task_addunlock( ci->progeny[6]->sorts[4] , t ); task_addunlock( cj->progeny[2]->sorts[4] , t );
                        task_addunlock( ci->progeny[7]->sorts[5] , t ); task_addunlock( cj->progeny[2]->sorts[5] , t );
                        task_addunlock( ci->progeny[4]->sorts[0] , t ); task_addunlock( cj->progeny[3]->sorts[0] , t );
                        task_addunlock( ci->progeny[5]->sorts[1] , t ); task_addunlock( cj->progeny[3]->sorts[1] , t );
                        task_addunlock( ci->progeny[6]->sorts[3] , t ); task_addunlock( cj->progeny[3]->sorts[3] , t );
                        task_addunlock( ci->progeny[7]->sorts[4] , t ); task_addunlock( cj->progeny[3]->sorts[4] , t );
                        }
                    else {
                        t->ci = ci->progeny[4]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[4]->sorts[4] , t ); task_addunlock( cj->progeny[0]->sorts[4] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[5] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[5] , t ); task_addunlock( cj->progeny[0]->sorts[5] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[7] , t ); task_addunlock( cj->progeny[0]->sorts[7] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[8] , t ); task_addunlock( cj->progeny[0]->sorts[8] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[4] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[4]->sorts[3] , t ); task_addunlock( cj->progeny[1]->sorts[3] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[5] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[4] , t ); task_addunlock( cj->progeny[1]->sorts[4] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[6] , t ); task_addunlock( cj->progeny[1]->sorts[6] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[7] , t ); task_addunlock( cj->progeny[1]->sorts[7] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[4] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[4]->sorts[1] , t ); task_addunlock( cj->progeny[2]->sorts[1] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[2] , t ); task_addunlock( cj->progeny[2]->sorts[2] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[4] , t ); task_addunlock( cj->progeny[2]->sorts[4] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[5] , t ); task_addunlock( cj->progeny[2]->sorts[5] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[4] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[4]->sorts[0] , t ); task_addunlock( cj->progeny[3]->sorts[0] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[5] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[1] , t ); task_addunlock( cj->progeny[3]->sorts[1] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[3] , t ); task_addunlock( cj->progeny[3]->sorts[3] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[4] , t ); task_addunlock( cj->progeny[3]->sorts[4] , t );
                        }
                    ci->progeny[4]->nr_pairs += 4;
                    ci->progeny[5]->nr_pairs += 4;
                    ci->progeny[6]->nr_pairs += 4;
                    ci->progeny[7]->nr_pairs += 4;
                    cj->progeny[0]->nr_pairs += 4;
                    cj->progeny[1]->nr_pairs += 4;
                    cj->progeny[2]->nr_pairs += 4;
                    cj->progeny[3]->nr_pairs += 4;
                    break;
                    
                case 5: /* (  1 ,  0 , -1 ) */
                    if ( space_dosub &&
                         !ci->progeny[4]->split && !ci->progeny[6]->split &&
                         !cj->progeny[1]->split && !cj->progeny[3]->split ) {
                        t->type = tid_sub; t->flags = 5;
                        task_addunlock( ci->progeny[4]->sorts[5] , t ); task_addunlock( cj->progeny[1]->sorts[5] , t );
                        task_addunlock( ci->progeny[6]->sorts[5] , t ); task_addunlock( cj->progeny[3]->sorts[5] , t );
                        task_addunlock( ci->progeny[4]->sorts[2] , t ); task_addunlock( cj->progeny[3]->sorts[2] , t );
                        task_addunlock( ci->progeny[6]->sorts[8] , t ); task_addunlock( cj->progeny[1]->sorts[8] , t );
                        }
                    else {
                        t->ci = ci->progeny[4]; t->cj = cj->progeny[1];
                        task_addunlock( ci->progeny[4]->sorts[5] , t ); task_addunlock( cj->progeny[1]->sorts[5] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[5] , t ); task_addunlock( cj->progeny[3]->sorts[5] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[4] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[4]->sorts[2] , t ); task_addunlock( cj->progeny[3]->sorts[2] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[8] , t ); task_addunlock( cj->progeny[1]->sorts[8] , t );
                        }
                    ci->progeny[4]->nr_pairs += 2;
                    ci->progeny[6]->nr_pairs += 2;
                    cj->progeny[1]->nr_pairs += 2;
                    cj->progeny[3]->nr_pairs += 2;
                    break;
                    
                case 6: /* (  1 , -1 ,  1 ) */
                    t->ci = ci->progeny[5]; t->cj = cj->progeny[2];
                    task_addunlock( ci->progeny[5]->sorts[6] , t ); task_addunlock( cj->progeny[2]->sorts[6] , t );
                    ci->progeny[5]->nr_pairs += 1;
                    cj->progeny[2]->nr_pairs += 1;
                    break;
                    
                case 7: /* (  1 , -1 ,  0 ) */
                    if ( space_dosub &&
                         !ci->progeny[4]->split && !ci->progeny[5]->split &&
                         !cj->progeny[2]->split && !cj->progeny[3]->split ) {
                        t->type = tid_sub; t->flags = 7;
                        task_addunlock( ci->progeny[4]->sorts[6] , t ); task_addunlock( cj->progeny[3]->sorts[6] , t );
                        task_addunlock( ci->progeny[5]->sorts[8] , t ); task_addunlock( cj->progeny[2]->sorts[8] , t );
                        task_addunlock( ci->progeny[4]->sorts[7] , t ); task_addunlock( cj->progeny[2]->sorts[7] , t );
                        task_addunlock( ci->progeny[5]->sorts[7] , t ); task_addunlock( cj->progeny[3]->sorts[7] , t );
                        }
                    else {
                        t->ci = ci->progeny[4]; t->cj = cj->progeny[3];
                        task_addunlock( ci->progeny[4]->sorts[6] , t ); task_addunlock( cj->progeny[3]->sorts[6] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[8] , t ); task_addunlock( cj->progeny[2]->sorts[8] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[4] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[4]->sorts[7] , t ); task_addunlock( cj->progeny[2]->sorts[7] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[5] , cj->progeny[3] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[7] , t ); task_addunlock( cj->progeny[3]->sorts[7] , t );
                        }
                    ci->progeny[4]->nr_pairs += 2;
                    ci->progeny[5]->nr_pairs += 2;
                    cj->progeny[2]->nr_pairs += 2;
                    cj->progeny[3]->nr_pairs += 2;
                    break;
                    
                case 8: /* (  1 , -1 , -1 ) */
                    t->ci = ci->progeny[4]; t->cj = cj->progeny[3];
                    task_addunlock( ci->progeny[4]->sorts[8] , t ); task_addunlock( cj->progeny[3]->sorts[8] , t );
                    ci->progeny[4]->nr_pairs += 1;
                    cj->progeny[3]->nr_pairs += 1;
                    break;
                    
                case 9: /* (  0 ,  1 ,  1 ) */
                    if ( space_dosub &&
                         !ci->progeny[3]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[4]->split ) {
                        t->type = tid_sub; t->flags = 9;
                        task_addunlock( ci->progeny[3]->sorts[9] , t ); task_addunlock( cj->progeny[0]->sorts[9] , t );
                        task_addunlock( ci->progeny[7]->sorts[9] , t ); task_addunlock( cj->progeny[4]->sorts[9] , t );
                        task_addunlock( ci->progeny[3]->sorts[0] , t ); task_addunlock( cj->progeny[4]->sorts[0] , t );
                        task_addunlock( ci->progeny[7]->sorts[8] , t ); task_addunlock( cj->progeny[0]->sorts[8] , t );
                        }
                    else {
                        t->ci = ci->progeny[3]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[3]->sorts[9] , t ); task_addunlock( cj->progeny[0]->sorts[9] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[9] , t ); task_addunlock( cj->progeny[4]->sorts[9] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[3] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[0] , t ); task_addunlock( cj->progeny[4]->sorts[0] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[8] , t ); task_addunlock( cj->progeny[0]->sorts[8] , t );
                        }
                    ci->progeny[3]->nr_pairs += 2;
                    ci->progeny[7]->nr_pairs += 2;
                    cj->progeny[0]->nr_pairs += 2;
                    cj->progeny[4]->nr_pairs += 2;
                    break;
                    
                case 10: /* (  0 ,  1 ,  0 ) */
                    if ( !ci->progeny[2]->split && !ci->progeny[3]->split && !ci->progeny[6]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[1]->split && !cj->progeny[4]->split && !cj->progeny[5]->split ) {
                        t->type = tid_sub; t->flags = 10;
                        task_addunlock( ci->progeny[2]->sorts[10] , t ); task_addunlock( cj->progeny[0]->sorts[10] , t );
                        task_addunlock( ci->progeny[3]->sorts[11] , t ); task_addunlock( cj->progeny[0]->sorts[11] , t );
                        task_addunlock( ci->progeny[6]->sorts[7] , t ); task_addunlock( cj->progeny[0]->sorts[7] , t );
                        task_addunlock( ci->progeny[7]->sorts[6] , t ); task_addunlock( cj->progeny[0]->sorts[6] , t );
                        task_addunlock( ci->progeny[2]->sorts[9] , t ); task_addunlock( cj->progeny[1]->sorts[9] , t );
                        task_addunlock( ci->progeny[3]->sorts[10] , t ); task_addunlock( cj->progeny[1]->sorts[10] , t );
                        task_addunlock( ci->progeny[6]->sorts[8] , t ); task_addunlock( cj->progeny[1]->sorts[8] , t );
                        task_addunlock( ci->progeny[7]->sorts[7] , t ); task_addunlock( cj->progeny[1]->sorts[7] , t );
                        task_addunlock( ci->progeny[2]->sorts[1] , t ); task_addunlock( cj->progeny[4]->sorts[1] , t );
                        task_addunlock( ci->progeny[3]->sorts[2] , t ); task_addunlock( cj->progeny[4]->sorts[2] , t );
                        task_addunlock( ci->progeny[6]->sorts[10] , t ); task_addunlock( cj->progeny[4]->sorts[10] , t );
                        task_addunlock( ci->progeny[7]->sorts[11] , t ); task_addunlock( cj->progeny[4]->sorts[11] , t );
                        task_addunlock( ci->progeny[2]->sorts[0] , t ); task_addunlock( cj->progeny[5]->sorts[0] , t );
                        task_addunlock( ci->progeny[3]->sorts[1] , t ); task_addunlock( cj->progeny[5]->sorts[1] , t );
                        task_addunlock( ci->progeny[6]->sorts[9] , t ); task_addunlock( cj->progeny[5]->sorts[9] , t );
                        task_addunlock( ci->progeny[7]->sorts[10] , t ); task_addunlock( cj->progeny[5]->sorts[10] , t );
                        }
                    else {
                        t->ci = ci->progeny[2]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[2]->sorts[10] , t ); task_addunlock( cj->progeny[0]->sorts[10] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[3] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[11] , t ); task_addunlock( cj->progeny[0]->sorts[11] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[7] , t ); task_addunlock( cj->progeny[0]->sorts[7] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[6] , t ); task_addunlock( cj->progeny[0]->sorts[6] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[2] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[2]->sorts[9] , t ); task_addunlock( cj->progeny[1]->sorts[9] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[3] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[10] , t ); task_addunlock( cj->progeny[1]->sorts[10] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[8] , t ); task_addunlock( cj->progeny[1]->sorts[8] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[7] , t ); task_addunlock( cj->progeny[1]->sorts[7] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[2] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[2]->sorts[1] , t ); task_addunlock( cj->progeny[4]->sorts[1] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[3] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[2] , t ); task_addunlock( cj->progeny[4]->sorts[2] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[10] , t ); task_addunlock( cj->progeny[4]->sorts[10] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[11] , t ); task_addunlock( cj->progeny[4]->sorts[11] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[2] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[2]->sorts[0] , t ); task_addunlock( cj->progeny[5]->sorts[0] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[3] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[1] , t ); task_addunlock( cj->progeny[5]->sorts[1] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[9] , t ); task_addunlock( cj->progeny[5]->sorts[9] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[10] , t ); task_addunlock( cj->progeny[5]->sorts[10] , t );
                        }
                    ci->progeny[2]->nr_pairs += 4;
                    ci->progeny[3]->nr_pairs += 4;
                    ci->progeny[6]->nr_pairs += 4;
                    ci->progeny[7]->nr_pairs += 4;
                    cj->progeny[0]->nr_pairs += 4;
                    cj->progeny[1]->nr_pairs += 4;
                    cj->progeny[4]->nr_pairs += 4;
                    cj->progeny[5]->nr_pairs += 4;
                    break;
                    
                case 11: /* (  0 ,  1 , -1 ) */
                    if ( space_dosub &&
                         !ci->progeny[2]->split && !ci->progeny[6]->split &&
                         !cj->progeny[1]->split && !cj->progeny[5]->split ) {
                        t->type = tid_sub; t->flags = 11;
                        task_addunlock( ci->progeny[2]->sorts[11] , t ); task_addunlock( cj->progeny[1]->sorts[11] , t );
                        task_addunlock( ci->progeny[6]->sorts[11] , t ); task_addunlock( cj->progeny[5]->sorts[11] , t );
                        task_addunlock( ci->progeny[2]->sorts[2] , t ); task_addunlock( cj->progeny[5]->sorts[2] , t );
                        task_addunlock( ci->progeny[6]->sorts[6] , t ); task_addunlock( cj->progeny[1]->sorts[6] , t );
                        }
                    else {
                        t->ci = ci->progeny[2]; t->cj = cj->progeny[1];
                        task_addunlock( ci->progeny[2]->sorts[11] , t ); task_addunlock( cj->progeny[1]->sorts[11] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[11] , t ); task_addunlock( cj->progeny[5]->sorts[11] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[2] , cj->progeny[5] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[2]->sorts[2] , t ); task_addunlock( cj->progeny[5]->sorts[2] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[6] , cj->progeny[1] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[6]->sorts[6] , t ); task_addunlock( cj->progeny[1]->sorts[6] , t );
                        }
                    ci->progeny[2]->nr_pairs += 2;
                    ci->progeny[6]->nr_pairs += 2;
                    cj->progeny[1]->nr_pairs += 2;
                    cj->progeny[5]->nr_pairs += 2;
                    break;
                    
                case 12: /* (  0 ,  0 ,  1 ) */
                    if ( space_dosub &&
                         !ci->progeny[1]->split && !ci->progeny[3]->split && !ci->progeny[5]->split && !ci->progeny[7]->split &&
                         !cj->progeny[0]->split && !cj->progeny[2]->split && !cj->progeny[4]->split && !cj->progeny[6]->split ) {
                        t->type = tid_sub; t->flags = 12;
                        task_addunlock( ci->progeny[1]->sorts[12] , t ); task_addunlock( cj->progeny[0]->sorts[12] , t );
                        task_addunlock( ci->progeny[3]->sorts[11] , t ); task_addunlock( cj->progeny[0]->sorts[11] , t );
                        task_addunlock( ci->progeny[5]->sorts[5] , t ); task_addunlock( cj->progeny[0]->sorts[5] , t );
                        task_addunlock( ci->progeny[7]->sorts[2] , t ); task_addunlock( cj->progeny[0]->sorts[2] , t );
                        task_addunlock( ci->progeny[1]->sorts[9] , t ); task_addunlock( cj->progeny[2]->sorts[9] , t );
                        task_addunlock( ci->progeny[3]->sorts[12] , t ); task_addunlock( cj->progeny[2]->sorts[12] , t );
                        task_addunlock( ci->progeny[5]->sorts[8] , t ); task_addunlock( cj->progeny[2]->sorts[8] , t );
                        task_addunlock( ci->progeny[7]->sorts[5] , t ); task_addunlock( cj->progeny[2]->sorts[5] , t );
                        task_addunlock( ci->progeny[1]->sorts[3] , t ); task_addunlock( cj->progeny[4]->sorts[3] , t );
                        task_addunlock( ci->progeny[3]->sorts[6] , t ); task_addunlock( cj->progeny[4]->sorts[6] , t );
                        task_addunlock( ci->progeny[5]->sorts[12] , t ); task_addunlock( cj->progeny[4]->sorts[12] , t );
                        task_addunlock( ci->progeny[7]->sorts[11] , t ); task_addunlock( cj->progeny[4]->sorts[11] , t );
                        task_addunlock( ci->progeny[1]->sorts[0] , t ); task_addunlock( cj->progeny[6]->sorts[0] , t );
                        task_addunlock( ci->progeny[3]->sorts[3] , t ); task_addunlock( cj->progeny[6]->sorts[3] , t );
                        task_addunlock( ci->progeny[5]->sorts[9] , t ); task_addunlock( cj->progeny[6]->sorts[9] , t );
                        task_addunlock( ci->progeny[7]->sorts[12] , t ); task_addunlock( cj->progeny[6]->sorts[12] , t );
                        }
                    else {
                        t->ci = ci->progeny[1]; t->cj = cj->progeny[0];
                        task_addunlock( ci->progeny[1]->sorts[12] , t ); task_addunlock( cj->progeny[0]->sorts[12] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[3] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[11] , t ); task_addunlock( cj->progeny[0]->sorts[11] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[5] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[5] , t ); task_addunlock( cj->progeny[0]->sorts[5] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[0] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[2] , t ); task_addunlock( cj->progeny[0]->sorts[2] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[1] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[1]->sorts[9] , t ); task_addunlock( cj->progeny[2]->sorts[9] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[3] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[12] , t ); task_addunlock( cj->progeny[2]->sorts[12] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[5] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[8] , t ); task_addunlock( cj->progeny[2]->sorts[8] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[2] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[5] , t ); task_addunlock( cj->progeny[2]->sorts[5] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[1] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[1]->sorts[3] , t ); task_addunlock( cj->progeny[4]->sorts[3] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[3] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[6] , t ); task_addunlock( cj->progeny[4]->sorts[6] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[5] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[12] , t ); task_addunlock( cj->progeny[4]->sorts[12] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[4] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[11] , t ); task_addunlock( cj->progeny[4]->sorts[11] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[1] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[1]->sorts[0] , t ); task_addunlock( cj->progeny[6]->sorts[0] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[3] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[3]->sorts[3] , t ); task_addunlock( cj->progeny[6]->sorts[3] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[5] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[5]->sorts[9] , t ); task_addunlock( cj->progeny[6]->sorts[9] , t );
                        t = space_addtask( s , tid_pair , 0 , 0 , ci->progeny[7] , cj->progeny[6] , NULL , 0 , NULL , 0 );
                        task_addunlock( ci->progeny[7]->sorts[12] , t ); task_addunlock( cj->progeny[6]->sorts[12] , t );
                        }
                    ci->progeny[1]->nr_pairs += 4;
                    ci->progeny[3]->nr_pairs += 4;
                    ci->progeny[5]->nr_pairs += 4;
                    ci->progeny[7]->nr_pairs += 4;
                    cj->progeny[0]->nr_pairs += 4;
                    cj->progeny[2]->nr_pairs += 4;
                    cj->progeny[4]->nr_pairs += 4;
                    cj->progeny[6]->nr_pairs += 4;
                    break;
            
                }
                
            /* Take a step back... */
            tid -= 1;
            
            } /* split this task? */
    
        } /* loop over all tasks. */
        
    }
    
    
/**
 * @brief Fill the #space's task list.
 *
 * @param s The #space we are working in.
 * @param do_sort Flag to add sorting tasks to the list.
 */
 
void space_maketasks ( struct space *s , int do_sort ) {

    int i, j, k, ii, jj, kk, iii, jjj, kkk, cid, cjd;
    int *cdim = s->cdim;
    int nr_tasks_old = s->nr_tasks;
    struct task *t;
    int pts[7][8] = { { -1 , 12 , 10 , 9 , 4 , 3 , 1 , 0 } ,
                      { -1 , -1 , 11 , 10 , 5 , 4 , 2 , 1 } ,
                      { -1 , -1 , -1 , 12 , 7 , 6 , 4 , 3 } , 
                      { -1 , -1 , -1 , -1 , 8 , 7 , 5 , 4 } ,
                      { -1 , -1 , -1 , -1 , -1 , 12 , 10 , 9 } ,
                      { -1 , -1 , -1 , -1 , -1 , -1 , 11 , 10 } ,
                      { -1 , -1 , -1 , -1 , -1 , -1 , -1 , 12 } };
    int counts[tid_count];

    /* Recursive function to generate tasks in the cell tree. */
    void maketasks_rec ( struct cell *c , struct task *sort_up[] , int nr_sort_up , struct cell *parent ) {

        int j, k, nr_sort = 0;
        struct task *sort[14], *t;

        /* Clear the waits on this cell. */
        c->wait = 0;
        sort[0] = NULL;
        
        /* Start by generating the sort task. */
        if ( c->count > 0 ) {
        
            if ( do_sort ) {
                if ( c->count < 1000 ) {
                    sort[0] = space_addtask( s , tid_sort , 0x1fff , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    for ( k = 0 ; k < 13 ; k++ )
                        c->sorts[k] = sort[0];
                    nr_sort = 1;
                    }
                else if ( c->count < 5000 ) {
                    sort[0] = space_addtask( s , tid_sort , 0xf , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    sort[1] = space_addtask( s , tid_sort , 0xf0 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    sort[2] = space_addtask( s , tid_sort , 0x1f00 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    for ( k = 0 ; k < 4 ; k++ )
                        c->sorts[k] = sort[0];
                    for ( k = 4 ; k < 8 ; k++ )
                        c->sorts[k] = sort[1];
                    for ( k = 8 ; k < 13 ; k++ )
                        c->sorts[k] = sort[2];
                    nr_sort = 3;
                    }
                else {
                    c->sorts[0] = sort[0] = space_addtask( s , tid_sort , 0x1 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[1] = sort[1] = space_addtask( s , tid_sort , 0x2 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[2] = sort[2] = space_addtask( s , tid_sort , 0x4 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[3] = sort[3] = space_addtask( s , tid_sort , 0x8 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[4] = sort[4] = space_addtask( s , tid_sort , 0x10 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[5] = sort[5] = space_addtask( s , tid_sort , 0x20 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[6] = sort[6] = space_addtask( s , tid_sort , 0x40 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[7] = sort[7] = space_addtask( s , tid_sort , 0x80 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[8] = sort[8] = space_addtask( s , tid_sort , 0x100 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[9] = sort[9] = space_addtask( s , tid_sort , 0x200 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[10] = sort[10] = space_addtask( s , tid_sort , 0x400 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[11] = sort[11] = space_addtask( s , tid_sort , 0x800 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    c->sorts[12] = sort[12] = space_addtask( s , tid_sort , 0x1000 , 0 , c , NULL , sort_up , nr_sort_up , NULL , 0 );
                    nr_sort = 13;
                    }
                }

            /* Generate a self-interaction if not split. */
            if ( !c->split && c->count > 1 )
                space_addtask( s , tid_self , 0 , 0 , c , NULL , NULL , 0 , NULL , 0 );
                
            }
            
        /* Otherwise, add the interactions between progeny. */
        if ( c->split ) {
        
            /* Recurse. */
            for ( k = 0 ; k < 8 ; k++ )
                if ( c->progeny[k] != NULL )
                    maketasks_rec( c->progeny[k] , sort , nr_sort , c );
                        
            /* Worth splitting into several tasks? */
            if ( !space_dosub || c->count > 2*space_splitsize ) {
            
                /* Make a task for each pair of progeny. */
                for ( j = 0 ; j < 8 ; j++ )
                    if ( c->progeny[j] != NULL && c->progeny[j]->count > 0 )
                        for ( k = j + 1 ; k < 8 ; k++ )
                            if ( c->progeny[k] != NULL && c->progeny[k]->count > 0 ) {
                                t = space_addtask( s , tid_pair , 0 , 0 , c->progeny[j] , c->progeny[k] , NULL , 0 , NULL , 0 );
                                task_addunlock( c->progeny[j]->sorts[ pts[j][k] ] , t );
                                task_addunlock( c->progeny[k]->sorts[ pts[j][k] ] , t );
                                c->progeny[k]->nr_pairs += 1;
                                c->progeny[j]->nr_pairs += 1;
                                }
                                
                }
                
            /* Otherwise, dispatch as one large task. */
            else {
            
                /* Add the task. */
                t = space_addtask( s , tid_sub , 0 , 0 , c , NULL , NULL , 0 , NULL , 0 );
                
                /* Make it depend on all the sorts of its progeny. */
                for ( k = 0 ; k < 8 ; k++ )
                    for ( j = 0 ; j < 13 ; j++ )
                        if ( c->progeny[k] != NULL )
                            task_addunlock( c->progeny[k]->sorts[j] , t );
            
                }

            }

        }
        
    /* Allocate the task-list, if needed. */
    if ( s->tasks == NULL )
        if ( posix_memalign( (void *)&s->tasks , 64 , sizeof(struct task) * s->tot_cells * 14 ) != 0 )
            error( "Failed to allocate task list." );
    s->nr_tasks = 0;
    
    /* Loop over the cells and get their sub-tasks. */
    for ( k = 0 ; k < s->nr_cells ; k++ )
        maketasks_rec( &s->cells[k] , NULL , 0 , NULL );

    /* Run through the highest level of cells and add pairs. */
    for ( i = 0 ; i < cdim[0] ; i++ )
        for ( j = 0 ; j < cdim[1] ; j++ )
            for ( k = 0 ; k < cdim[2] ; k++ ) {
                cid = cell_getid( cdim , i , j , k );
                if ( s->cells[cid].count == 0 )
                    continue;
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
                            if ( s->cells[cjd].count == 0 )
                                continue;
                            if ( cid >= cjd )
                                continue;
                            t = space_addtask( s , tid_pair , 0 , 0 , &s->cells[cid] , &s->cells[cjd] , NULL , 0 , NULL , 0 );
                            task_addunlock( s->cells[cid].sorts[ sortlistID[ (kk+1) + 3*( (jj+1) + 3*(ii+1) ) ] ] , t );
                            task_addunlock( s->cells[cjd].sorts[ sortlistID[ (kk+1) + 3*( (jj+1) + 3*(ii+1) ) ] ] , t );
                            s->cells[cid].nr_pairs += 1;
                            s->cells[cjd].nr_pairs += 1;
                            }
                        }
                    }
                }

    /* Split the tasks. */
    space_splittasks( s );
        
    /* Did we already create indices? */
    if ( s->tasks_ind == NULL )
        if ( ( s->tasks_ind = (int *)malloc( sizeof(int) * s->nr_tasks ) ) == NULL )
            error( "Failed to allocate task indices." );
    
    /* Did the number of tasks change, i.e. do we have to re-index? */
    if ( nr_tasks_old != s->nr_tasks )
        for ( k = 0 ; k < s->nr_tasks ; k++ )
            s->tasks_ind[k] = k;
            
    /* Remove sort tasks with no dependencies. */
    for ( k = 0 ; k < s->nr_tasks ; k++ ) {
        t = &s->tasks[k];
        if ( t->type == tid_sort && t->nr_unlock_tasks == 0 ) {
            t->type = tid_none;
            if ( t->ci->split )
                for ( j = 0 ; j < 8 ; j++ )
                    if ( t->ci->progeny[j] != NULL && t->flags & ( 1 << j ) )
                        task_rmunlock( t->ci->progeny[j]->sorts[j] , t );
            }
        }
            
    /* Count the number of each task type. */
    for ( k = 0 ; k < tid_count ; k++ )
        counts[k] = 0;
    for ( k = 0 ; k < s->nr_tasks ; k++ )
        counts[ s->tasks[k].type ] += 1;
    printf( "space_maketasks: task counts are [ %s=%i" , taskID_names[0] , counts[0] );
    for ( k = 1 ; k < tid_count ; k++ )
        printf( " %s=%i" , taskID_names[k] , counts[k] );
    printf( " ]\n" );
        
    /* Re-set the next task pointer. */
    s->next_task = 0;
            
    }
    
    
/**
 * @brief Sort the parts into eight bins along the given pivots.
 *
 * @param c The #cell array to be sorted.
 */
 
void cell_split ( struct cell *c  ) {

    int i, j, k, kk;
    struct part temp, *parts = c->parts;
    int left[8], right[8];
    double pivot[3];
    
    /* Init the pivot. */
    for ( k = 0 ; k < 3 ; k++ )
        pivot[k] = c->loc[k] + c->h[k]/2;
    
    /* Split along the x-axis. */
    i = 0; j = c->count - 1;
    while ( i <= j ) {
        while ( i <= c->count-1 && parts[i].x[0] <= pivot[0] )
            i += 1;
        while ( j >= 0 && parts[j].x[0] > pivot[0] )
            j -= 1;
        if ( i < j ) {
            temp = parts[i]; parts[i] = parts[j]; parts[j] = temp;
            }
        }
    for ( k = 0 ; k <= j ; k++ )
        if ( parts[k].x[0] > pivot[0] )
            error( "cell_split: sorting failed." );
    for ( k = i ; k < c->count ; k++ )
        if ( parts[k].x[0] < pivot[0] )
            error( "cell_split: sorting failed." );
    left[1] = i; right[1] = c->count - 1;
    left[0] = 0; right[0] = j;
    
    /* Split along the y axis, twice. */
    for ( k = 1 ; k >= 0 ; k-- ) {
        i = left[k]; j = right[k];
        while ( i <= j ) {
            while ( i <= right[k] && parts[i].x[1] <= pivot[1] )
                i += 1;
            while ( j >= left[k] && parts[j].x[1] > pivot[1] )
                j -= 1;
            if ( i < j ) {
                temp = parts[i]; parts[i] = parts[j]; parts[j] = temp;
                }
            }
        for ( kk = left[k] ; kk <= j ; kk++ )
            if ( parts[kk].x[1] > pivot[1] ) {
                printf( "cell_split: ival=[%i,%i], i=%i, j=%i.\n" , left[k] , right[k] , i , j );
                error( "sorting failed (left)." );
                }
        for ( kk = i ; kk <= right[k] ; kk++ )
            if ( parts[kk].x[1] < pivot[1] )
                error( "sorting failed (right)." );
        left[2*k+1] = i; right[2*k+1] = right[k];
        left[2*k] = left[k]; right[2*k] = j;
        }

    /* Split along the z axis, four times. */
    for ( k = 3 ; k >= 0 ; k-- ) {
        i = left[k]; j = right[k];
        while ( i <= j ) {
            while ( i <= right[k] && parts[i].x[2] <= pivot[2] )
                i += 1;
            while ( j >= left[k] && parts[j].x[2] > pivot[2] )
                j -= 1;
            if ( i < j ) {
                temp = parts[i]; parts[i] = parts[j]; parts[j] = temp;
                }
            }
        for ( kk = left[k] ; kk <= j ; kk++ )
            if ( parts[kk].x[2] > pivot[2] ) {
                printf( "cell_split: ival=[%i,%i], i=%i, j=%i.\n" , left[k] , right[k] , i , j );
                error( "sorting failed (left)." );
                }
        for ( kk = i ; kk <= right[k] ; kk++ )
            if ( parts[kk].x[2] < pivot[2] ) {
                printf( "cell_split: ival=[%i,%i], i=%i, j=%i.\n" , left[k] , right[k] , i , j );
                error( "sorting failed (right)." );
                }
        left[2*k+1] = i; right[2*k+1] = right[k];
        left[2*k] = left[k]; right[2*k] = j;
        }
        
    /* Store the counts and offsets. */
    for ( k = 0 ; k < 8 ; k++ ) {
        c->progeny[k]->count = right[k] - left[k] + 1;
        if ( c->progeny[k]->count < 0 )
            abort();
        c->progeny[k]->parts = &c->parts[ left[k] ];
        }
        
    /* Verify a few sub-cells. */
    /* for ( k = 0 ; k < c->progeny[0]->count ; k++ )
        if ( c->progeny[0]->parts[k].x[0] > pivot[0] ||
             c->progeny[0]->parts[k].x[1] > pivot[1] ||
             c->progeny[0]->parts[k].x[2] > pivot[2] )
            error( "Sorting failed (progeny=0)." );
    for ( k = 0 ; k < c->progeny[1]->count ; k++ )
        if ( c->progeny[1]->parts[k].x[0] > pivot[0] ||
             c->progeny[1]->parts[k].x[1] > pivot[1] ||
             c->progeny[1]->parts[k].x[2] <= pivot[2] )
            error( "Sorting failed (progeny=1)." );
    for ( k = 0 ; k < c->progeny[2]->count ; k++ )
        if ( c->progeny[2]->parts[k].x[0] > pivot[0] ||
             c->progeny[2]->parts[k].x[1] <= pivot[1] ||
             c->progeny[2]->parts[k].x[2] > pivot[2] )
            error( "Sorting failed (progeny=2)." ); */

    }


/**
 * @brief Split cells that contain too many particles.
 *
 * @param s The #space we are working in.
 * @param c The #cell under consideration.
 */
 
void space_split ( struct space *s , struct cell *c ) {

    int k, count;
    double r, r_limit, r_max = 0.0;
    struct cell *temp;
    
    /* Check the depth. */
    if ( c->depth > s->maxdepth )
        s->maxdepth = c->depth;
    
    /* Set the minimum cutoff. */
    r_limit = fmin( c->h[0] , fmin( c->h[1] , c->h[2] ) ) / 2;
    
    /* Count the particles below that. */
    for ( count = 0 , k = 0 ; k < c->count ; k++ ) {
        r = c->parts[k].r;
        if ( r <= r_limit )
            count += 1;
        if ( r > r_max )
            r_max = r;
        }
    c->r_max = r_max;
            
    /* Split or let it be? */
    if ( count > c->count*space_splitratio && c->count > space_splitsize ) {
    
        /* No longer just a leaf. */
        c->split = 1;
        
        /* Create the cell's progeny. */
        for ( k = 0 ; k < 8 ; k++ ) {
            temp = space_getcell( s );
            temp->count = 0;
            temp->loc[0] = c->loc[0];
            temp->loc[1] = c->loc[1];
            temp->loc[2] = c->loc[2];
            temp->h[0] = c->h[0]/2;
            temp->h[1] = c->h[1]/2;
            temp->h[2] = c->h[2]/2;
            if ( k & 4 )
                temp->loc[0] += temp->h[0];
            if ( k & 2 )
                temp->loc[1] += temp->h[1];
            if ( k & 1 )
                temp->loc[2] += temp->h[2];
            temp->depth = c->depth + 1;
            temp->split = 0;
            temp->r_max = 0.0;
            temp->parent = c;
            c->progeny[k] = temp;
            }
            
        /* Split the cell data. */
        cell_split( c );
            
        /* Recurse? */
        for ( k = 0 ; k < 8 ; k++ )
            space_split( s , c->progeny[k] );
            
        /* Remove any progeny with zero parts. */
        for ( k = 0 ; k < 8 ; k++ )
            if ( c->progeny[k]->count == 0 ) {
                space_recycle( s , c->progeny[k] );
                c->progeny[k] = NULL;
                }
            
        }
        
    /* Otherwise, set the progeny to null. */
    else {
        bzero( c->progeny , sizeof(struct cell *) * 8 );
        c->split = 0;
        }

    }


/**
 * @brief Return a used cell to the cell buffer.
 *
 * @param s The #space.
 * @param c The #cell.
 */
 
void space_recycle ( struct space *s , struct cell *c ) {

    /* Clear the cell. */
    if ( lock_destroy( &c->lock ) != 0 )
        error( "Failed to destroy spinlock." );
    
    /* Hook this cell into the buffer. */
    c->next = s->cells_new;
    s->cells_new = c;
    s->tot_cells -= 1;
    
    }


/**
 * @brief Get a new empty cell.
 *
 * @param s The #space.
 */
 
struct cell *space_getcell ( struct space *s ) {

    struct cell *c;
    int k;
    
    /* Is the buffer empty? */
    if ( s->cells_new == NULL ) {
        if ( posix_memalign( (void *)&s->cells_new , 64 , space_cellallocchunk * sizeof(struct cell) ) != 0 )
            error( "Failed to allocate more cells." );
        bzero( s->cells_new , space_cellallocchunk * sizeof(struct cell) );
        for ( k = 0 ; k < space_cellallocchunk-1 ; k++ )
            s->cells_new[k].next = &s->cells_new[k+1];
        s->cells_new[ space_cellallocchunk-1 ].next = NULL;
        }

    /* Pick off the next cell. */
    c = s->cells_new;
    s->cells_new = c->next;
    s->tot_cells += 1;
    
    /* Init some things in the cell. */
    bzero( c , sizeof(struct cell) );
    if ( lock_init( &c->lock ) != 0 )
        error( "Failed to initialize cell spinlock." );
        
    return c;

    }


/**
 * @brief Split the space into cells given the array of particles.
 *
 * @param The #space to initialize.
 * @param dim Spatial dimensions of the domain.
 * @param parts Pointer to an array of #part.
 * @param N The number of parts in the space.
 * @param periodic flag whether the domain is periodic or not.
 *
 * Makes a grid of edge length > r_max and fills the particles
 * into the respective cells. Cells containing more than #space_maxppc
 * parts with a cutoff below half the cell width are then split
 * recursively.
 */


void space_init ( struct space *s , double dim[3] , struct part *parts , int N , int periodic , double h_max ) {

    int i, j, k;
    int nr_cells, cdim[3];
    double r_min, r_max, h[3], ih[3];
    struct cell *c, *cells;
    struct part *parts_new, *finger;
    
    
    /* Get the minimum and maximum cutoff radii. */
    r_min = parts[0].r; r_max = r_min;
    for ( k = 1 ; k < N ; k++ )
        if ( parts[k].r < r_min )
            r_min = parts[k].r;
        else if ( parts[k].r > r_max )
            r_max = parts[k].r;
            
    /* Get the cell width. */
    if ( h_max < r_max )
        h_max = r_max;
    for ( k = 0 ; k < 3 ; k++ ) {
        cdim[k] = ceil( dim[k] / h_max );
        h[k] = dim[k] / cdim[k];
        ih[k] = 1.0 / h[k];
        }
        
    /* Allocate the highest level of cells. */
    nr_cells = cdim[0] * cdim[1] * cdim[2];
    if ( posix_memalign( (void *)&cells , 64 , nr_cells * sizeof(struct cell) ) != 0 )
        error( "Failed to allocate cells." );
    bzero( cells , nr_cells * sizeof(struct cell) );
    for ( k = 0 ; k < nr_cells ; k++ )
        if ( lock_init( &cells[k].lock ) != 0 )
            error( "Failed to init spinlock." );
        
    /* Set the cell locations. */
    for ( i = 0 ; i < cdim[0] ; i++ )
        for ( j = 0 ; j < cdim[1] ; j++ )
            for ( k = 0 ; k < cdim[2] ; k++ ) {
                c = &cells[ cell_getid( cdim , i , j , k ) ];
                c->loc[0] = i*h[0]; c->loc[1] = j*h[1]; c->loc[2] = k*h[2];
                c->h[0] = h[0]; c->h[1] = h[1]; c->h[2] = h[2];
                }
        
    /* Run through the particles and get the counts for each cell. */
    for ( k = 0 ; k < N ; k++ )
        cells[ cell_getid( cdim , parts[k].x[0]*ih[0] , parts[k].x[1]*ih[1] , parts[k].x[2]*ih[2] ) ].count += 1;
        
    /* Allocate the new part buffer and set the part pointers in each cell. */
    if ( posix_memalign( (void *)&parts_new , 64 , N * sizeof(struct part) ) != 0 )
        error( "Failed to allocate parts." );
    for ( finger = parts_new , k = 0 ; k < nr_cells ; k++ ) {
        c = &cells[ k ];
        c->parts = finger;
        finger = &finger[ c->count ];
        c->count = 0;
        }
    for ( k = 0 ; k < N ; k++ ) {
        c = &cells[ cell_getid( cdim , parts[k].x[0]*ih[0] , parts[k].x[1]*ih[1] , parts[k].x[2]*ih[2] ) ];
        c->parts[ c->count ] = parts[k];
        c->count += 1;
        }
        
    /* Store eveything in the space. */
    s->r_min = r_min; s->r_max = r_max;
    s->dim[0] = dim[0]; s->dim[1] = dim[1]; s->dim[2] = dim[2];
    s->periodic = periodic;
    s->parts = parts_new;
    s->nr_parts = N;
    s->h[0] = h[0]; s->h[1] = h[1]; s->h[2] = h[2];
    s->ih[0] = ih[0]; s->ih[1] = ih[1]; s->ih[2] = ih[2];
    s->cdim[0] = cdim[0]; s->cdim[1] = cdim[1]; s->cdim[2] = cdim[2];
    s->cells = cells;
    s->nr_cells = nr_cells;
    s->tot_cells = nr_cells;
    if ( lock_init( &s->task_lock ) != 0 )
        error( "Failed to create task spin-lock." );
    
    /* Loop over the cells and split them. */
    for ( k = 0 ; k < nr_cells ; k++ )
        space_split( s , &cells[k] );
        
    }

