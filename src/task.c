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
#include "lock.h"
#include "task.h"


/* Task type names. */
const char *taskID_names[task_type_count] = { "none" , "sort" , "self" , "pair" , "sub" , "ghost" };

/* Error macro. */
#define error(s) { fprintf( stderr , "%s:%s:%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }


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
 
void task_rmunlock( struct task *ta , struct task *tb ) {

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
 
void task_rmunlock_blind( struct task *ta , struct task *tb ) {

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
 
void task_addunlock( struct task *ta , struct task *tb ) {

    int k;
    
    lock_lock( &ta->lock );
    
    /* Check if ta already unlocks tb. */
    for ( k = 0 ; k < ta->nr_unlock_tasks ; k++ )
        if ( ta->unlock_tasks[k] == tb ) {
            lock_unlock_blind( &ta->lock );
            return;
            }

    if ( ta->nr_unlock_tasks == task_maxunlock )
        error( "Too many unlock_tasks in task." );
        
    ta->unlock_tasks[ ta->nr_unlock_tasks] = tb;
    ta->nr_unlock_tasks += 1;

    lock_unlock_blind( &ta->lock );
    
    }
    

