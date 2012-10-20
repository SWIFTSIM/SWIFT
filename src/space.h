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




/* Some constants. */
#define space_maxdepth                  10
#define space_cellallocchunk            1000
#define space_splitratio                0.875
#define space_splitsize_default         800
#define space_dosub                     1


/* Split size. */
extern int space_splitsize;

/* Map shift vector to sortlist. */
extern const int sortlistID[27];
    
    
/* Entry in a list of sorted indices. */
struct entry {
    float d;
    int i;
    };
    
    
/* The space in which the cells reside. */
struct space {

    /* Spatial extent. */
    double dim[3];
    
    /* Cell widths. */
    double h[3], ih[3];
    
    /* The minimum and maximum cutoff radii. */
    double r_min, r_max;
    
    /* Number of cells. */
    int nr_cells, tot_cells;
    
    /* Space dimensions in number of cells. */
    int maxdepth, cdim[3];
    
    /* The (level 0) cells themselves. */
    struct cell *cells;
    
    /* Buffer of unused cells. */
    struct cell *cells_new;
    
    /* The particle data (cells have pointers to this). */
    struct part *parts;
    
    /* The sortlist data (cells hold pointers to these). */
    struct entry *sortlist;
    
    /* The total number of parts in the space. */
    int nr_parts;
    
    /* Is the space periodic? */
    int periodic;
    
    /* The list of tasks. */
    struct task *tasks;
    int nr_tasks, next_task;
    int *tasks_ind;
    lock_type task_lock;
    
    };


/* function prototypes. */
struct cell *space_getcell ( struct space *s );
struct task *space_gettask ( struct space *s );
void space_init ( struct space *s , double dim[3] , struct part *parts , int N , int periodic , double h_max );
void space_maketasks ( struct space *s , int do_sort );
void space_map_cells ( struct space *s , void (*fun)( struct cell *c , void *data ) , void *data );
void space_map_parts ( struct space *s , void (*fun)( struct part *p , struct cell *c , void *data ) , void *data );
void space_recycle ( struct space *s , struct cell *c );



