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

/* MPI headers. */
#ifdef WITH_MPI
    #include <mpi.h>
#endif

/* Local headers. */
#include "const.h"
#include "cycle.h"
#include "atomic.h"
#include "lock.h"
#include "space.h"
#include "cell.h"
#include "task.h"
#include "error.h"

/* Task type names. */
const char *taskID_names[task_type_count] = {   
    "none" , "sort" , "self" , "pair" , "sub" , "ghost" , 
    "kick1" , "kick2" , "send_xv" , "recv_xv" , "send_rho" ,
    "recv_rho" };


/**
 * @brief Unlock the cell held by this task.
 * 
 * @param t The #task.
 */
 
void task_unlock ( struct task *t ) {

    /* Act based on task type. */
    switch ( t->type ) {
        case task_type_self:
        case task_type_sort:
            cell_unlocktree( t->ci );
            break;
        case task_type_pair:
        case task_type_sub:
            cell_unlocktree( t->ci );
            if ( t->cj != NULL )
                cell_unlocktree( t->cj );
            break;
        }
        
    }


/**
 * @brief Try to lock the cells associated with this task.
 *
 * @param t the #task.
 */
 
int task_lock ( struct task *t ) {

    int type = t->type;
    struct cell *ci = t->ci, *cj = t->cj;

    /* Communication task? */
    if ( type == task_type_recv_xv || type == task_type_recv_rho ||
         type == task_type_send_xv || type == task_type_send_rho ) {
    
        #ifdef WITH_MPI
            /* Check the status of the MPI request. */
            int res, err;
            MPI_Status stat;
            if ( ( err = MPI_Test( &t->req , &res , &stat ) ) != MPI_SUCCESS ) {
                char buff[ MPI_MAX_ERROR_STRING ];
                int len;
                MPI_Error_string( err , buff , &len );
                message( "MPI error: %s\n" , buff ); fflush(stdout);
                error( "Failed to test request on send/recv task." );
                }
            return res;
        #else
            error( "SWIFT was not compiled with MPI support." );
        #endif
    
        }

    /* Unary lock? */
    else if ( type == task_type_self || 
         type == task_type_sort || 
         (type == task_type_sub && cj == NULL) ) {
        if ( cell_locktree( ci ) != 0 )
            return 0;
        }
        
    /* Otherwise, binary lock. */
    else if ( type == task_type_pair || (type == task_type_sub && cj != NULL) ) {
        if ( ci->hold || cj->hold || ci->wait || cj->wait )
            return 0;
        if ( cell_locktree( ci ) != 0 )
            return 0;
        if ( cell_locktree( cj ) != 0 ) {
            cell_unlocktree( ci );
            return 0;
            }
        }
        
    /* If we made it this far, we've got a lock. */
    return 1;
            
    }


/**
 * @brief Remove all unlocks to tasks that are of the given type.
 *
 * @param t The #task.
 * @param type The task type ID to remove.
 */
 
void task_cleanunlock ( struct task *t , int type ) {

    int k;
    
    lock_lock( &t->lock );
    
    for ( k = 0 ; k < t->nr_unlock_tasks ; k++ )
        if ( t->unlock_tasks[k]->type == type ) {
            t->nr_unlock_tasks -= 1;
            t->unlock_tasks[k] = t->unlock_tasks[ t->nr_unlock_tasks ];
            }
    
    lock_unlock_blind( &t->lock );
    
    }


/**
 * @brief Remove an unlock_task from the given task.
 *
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 */
 
void task_rmunlock ( struct task *ta , struct task *tb ) {

    int k;
    
    lock_lock( &ta->lock );
    
    for ( k = 0 ; k < ta->nr_unlock_tasks ; k++ )
        if ( ta->unlock_tasks[k] == tb ) {
            ta->nr_unlock_tasks -= 1;
            ta->unlock_tasks[k] = ta->unlock_tasks[ ta->nr_unlock_tasks ];
            lock_unlock_blind( &ta->lock );
            return;
            }
    error( "Task not found." );

    }
    

/**
 * @brief Remove an unlock_task from the given task.
 *
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 *
 * Differs from #task_rmunlock in that it will not fail if
 * the task @c tb is not in the unlocks of @c ta.
 */
 
void task_rmunlock_blind ( struct task *ta , struct task *tb ) {

    int k;
    
    lock_lock( &ta->lock );
    
    for ( k = 0 ; k < ta->nr_unlock_tasks ; k++ )
        if ( ta->unlock_tasks[k] == tb ) {
            ta->nr_unlock_tasks -= 1;
            ta->unlock_tasks[k] = ta->unlock_tasks[ ta->nr_unlock_tasks ];
            break;
            }
            
    lock_unlock_blind( &ta->lock );

    }
    

/**
 * @brief Add an unlock_task to the given task.
 *
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 */
 
void task_addunlock ( struct task *ta , struct task *tb ) {

    /* Add the lock atomically. */
    ta->unlock_tasks[ atomic_inc( &ta->nr_unlock_tasks ) ] = tb;

    /* Check a posteriori if we did not overshoot. */
    if ( ta->nr_unlock_tasks > task_maxunlock )
        error( "Too many unlock_tasks in task." );
        
    }
    

void task_addunlock_old ( struct task *ta , struct task *tb ) {

    int k;
    
    lock_lock( &ta->lock );
    
    /* Check if ta already unlocks tb. */
    for ( k = 0 ; k < ta->nr_unlock_tasks ; k++ )
        if ( ta->unlock_tasks[k] == tb ) {
            error( "Duplicate unlock." );
            lock_unlock_blind( &ta->lock );
            return;
            }

    if ( ta->nr_unlock_tasks == task_maxunlock )
        error( "Too many unlock_tasks in task." );
        
    ta->unlock_tasks[ ta->nr_unlock_tasks] = tb;
    ta->nr_unlock_tasks += 1;

    lock_unlock_blind( &ta->lock );
    
    }
    

