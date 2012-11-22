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



/* The queue timers. */
enum {
    cell_timer_none = 0,
    cell_timer_tree,
    cell_timer_count,
    };
extern ticks cell_timer[ cell_timer_count ];


/* Structure to store the data of a single cell. */
struct cell {

    /* The cell location on the grid. */
    double loc[3];
    
    /* The cell dimensions. */
    double h[3];
    
    /* Max radii in this cell. */
    double r_max;
    
    /* The depth of this cell in the tree. */
    int depth, split;
    
    /* Nr of parts. */
    int count;
    
    /* Pointers to the particle data. */
    struct part *parts;
    
    /* Pointers for the sorted indices. */
    struct entry *sort;
    
    /* Number of pairs associated with this cell. */
    int nr_pairs;
    
    /* Pointers to the next level of cells. */
    struct cell *progeny[8];
    
    /* Parent cell. */
    struct cell *parent;
    
    /* Super cell, i.e. the highest-level supercell that has interactions. */
    struct cell *super;
    
    /* The tasks computing this cell's sorts. */
    struct task *sorts[14];
    
    /* The ghost task to link density to interactions. */
    struct task *ghost;
    
    /* Number of tasks that are associated with this cell. */
    int nr_tasks;
    
    /* Number of tasks this cell is waiting for and whether it is in use. */
    int wait;
    
    /* Is the data of this cell being used in a sub-cell? */
    int hold;
    
    /* Spin lock for various uses. */
    lock_type lock;
    
    /* ID of the previous owner, e.g. runner. */
    int owner;
    
    /* Linking pointer for "memory management". */
    struct cell *next;

    } __attribute__((aligned (64)));


/* Function prototypes. */
void cell_split ( struct cell *c  );
int cell_locktree( struct cell *c );
void cell_unlocktree( struct cell *c );
