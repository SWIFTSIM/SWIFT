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


/* Some constants. */
#define task_maxwait                    3
#define task_maxunlock                  40


/* The different task types. */
enum task_types {
    task_type_none = 0,
    task_type_sort,
    task_type_self,
    task_type_pair,
    task_type_sub,
    task_type_ghost,
    task_type_kick2,
    task_type_count
    };
    
extern const char *taskID_names[];
    
/* The different task sub-types. */
enum task_subtypes {
    task_subtype_none = 0,
    task_subtype_density,
    task_subtype_force,
    task_subtype_count
    };
    
extern const char *taskID_names[];
    
/* Data of a task. */
struct task {

    char type, subtype, skip, tight;
    int flags, wait, rank, maxdepth;
    
    lock_type lock;
    
    struct cell *ci, *cj;
    
    int rid;
    ticks tic, toc;
    
    int nr_unlock_tasks;
    struct task *unlock_tasks[ task_maxunlock ];

    };


/* Function prototypes. */
void task_rmunlock( struct task *ta , struct task *tb );
void task_rmunlock_blind( struct task *ta , struct task *tb );
void task_cleanunlock ( struct task *t , int type );
void task_addunlock( struct task *ta , struct task *tb );
