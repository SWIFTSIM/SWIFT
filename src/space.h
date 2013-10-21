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
#define space_maxdepth                  10
#define space_cellallocchunk            1000
#define space_splitratio                0.875f
#define space_splitsize_default         400
#define space_maxsize_default           8000000
#define space_subsize_default           8000000
#define space_stretch                   1.10f
#define space_maxreldx                  0.25f
#define space_qstack                    2048


/* Convert cell location to ID. */
#define cell_getid( cdim , i , j , k ) ( (int)(k) + (cdim)[2]*( (int)(j) + (cdim)[1]*(int)(i) ) )

/* Split size. */
extern int space_splitsize;
extern int space_maxsize;
extern int space_subsize;

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
    double h_max, cell_min;
    
    /* Current time step for particles. */
    float dt_step;
    
    /* Current maximum displacement for particles. */
    float dx_max;
    
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
    struct xpart *xparts;
    struct gpart *gparts;
    
    /* The total number of parts in the space. */
    int nr_parts, size_parts;
    int nr_gparts, size_gparts;
    
    /* Is the space periodic? */
    int periodic;
    
    /* General-purpose lock for this space. */
    lock_type lock;
    
    /* Number of queues in the system. */
    int nr_queues;
    
    /* The associated engine. */
    struct engine *e;
    
    /* Buffers for parts that we will receive from foreign cells. */
    struct part *parts_foreign;
    int nr_parts_foreign, size_parts_foreign;
    
    };


/* function prototypes. */
void parts_sort ( struct part *parts , struct xpart *xparts , int *ind , int N , int min , int max );
void gparts_sort ( struct gpart *gparts , int *ind , int N , int min , int max );
struct cell *space_getcell ( struct space *s );
int space_getsid ( struct space *s , struct cell **ci , struct cell **cj , double *shift );
void space_init ( struct space *s , double dim[3] , struct part *parts , int N , int periodic , double h_max );
void space_map_cells_pre ( struct space *s , int full , void (*fun)( struct cell *c , void *data ) , void *data );
void space_map_parts ( struct space *s , void (*fun)( struct part *p , struct cell *c , void *data ) , void *data );
void space_map_cells_post ( struct space *s , int full , void (*fun)( struct cell *c , void *data ) , void *data );
void space_rebuild ( struct space *s , double h_max );
void space_recycle ( struct space *s , struct cell *c );
void space_split ( struct space *s , struct cell *c );

