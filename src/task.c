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


/* Task type names. */
const char *taskID_names[task_type_count] = { "none" , "sort" , "self" , "pair" , "sub" , "ghost" };

/* Error macro. */
#define error(s) { printf( "%s:%s:%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }


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
    

