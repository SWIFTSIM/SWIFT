/**
 *  Generate 3D random points using a poisson disc distribution.
 *
 *  Poisson disc requires that each position is at least a given
 *  radius from all other points, so is a good way to sample a grid
 *  without introducing aliasing.
 *
 *  See: Robert Bridson. 2007. Fast Poisson disk sampling in
 *       arbitrary dimensions. In ACM SIGGRAPH 2007 sketches
 *       (SIGGRAPH '07). ACM, New York, NY, USA, , Article 22.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <values.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#define CHUNK 512
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

/*  Structs for a point and a home for the grid related data. */
struct point
{
    float x;
    float y;
    float z;
};

struct griddata
{
    struct point *cells;
    float cell_size;
    float radius;
    int depth;
    int height;
    int width;

    struct point *queue;
    int qlength;
    int qsize;

    struct point *samples;
    int nsample;
    int ssize;
};


/**
 * Get random number in range 0 to 1.
 */
static double randone()
{
    return rand() / (double) RAND_MAX;
}

/**
 *  Add a position to the active list.
 */
static void markPosition( struct griddata *grid, struct point *p )
{
    /*  Index of point on grid. */
    int index = (int) grid->depth * grid->width * floor( p->z / grid->cell_size ) +
                grid->width * floor( p->y / grid->cell_size ) +
                floor( p->x / grid->cell_size );

    /*  Check if already seen, nothing to do. */
    if ( grid->cells[index].x == -1.0f ) {

        grid->cells[index].x = p->x;
        grid->cells[index].y = p->y;
        grid->cells[index].z = p->z;

        grid->queue[grid->qlength].x = p->x;
        grid->queue[grid->qlength].y = p->y;
        grid->queue[grid->qlength].z = p->z;
        grid->qlength++;

        /*  Add more space to queue, if needed. */
        if ( grid->qlength >= grid->qsize ) {
            grid->qsize = grid->qlength + CHUNK;
            grid->queue = realloc( grid->queue, grid->qsize * sizeof(struct point) );
        }
    }
}

/**
 *  Remove a position from the active list.
 */
static void unqueuePosition( struct griddata *grid, int index )
{
    /*  Shuffle queue. */
    for ( int i = index; i < grid->qlength - 1; i++ ) {
        grid->queue[i].x = grid->queue[i+1].x;
        grid->queue[i].y = grid->queue[i+1].y;
        grid->queue[i].z = grid->queue[i+1].z;
    }
    grid->qlength--;
}

/**
 *  Test that a position is available, i.e further than the radius away from
 *  positions already marked on the grid.
 */
static int available( struct griddata *grid, struct point *p )
{
    int i = (int) p->y / grid->cell_size;
    int j = (int) p->x / grid->cell_size;
    int l = (int) p->z / grid->cell_size;
    int i0 = MAX( i - 2, 0 );
    int j0 = MAX( j - 2, 0 );
    int l0 = MAX( l - 2, 0 );
    int i1 = MIN( i + 3, grid->height );
    int j1 = MIN( j + 3, grid->width );
    int l1 = MIN( l + 3, grid->depth );

    for ( int j = j0; j < j1; j++ ) {
        for ( int i = i0; i < i1; i++ ) {
            for ( int l = l0; l < l1; l++ ) {
                int index = l * grid->width * grid->depth + i * grid->width + j;
                if ( grid->cells[index].x != -1.0 ) {
                    float dx = grid->cells[index].x - p->x;
                    float dy = grid->cells[index].y - p->y;
                    float dz = grid->cells[index].z - p->z;
                    if ( ( dx*dx + dy*dy + dz*dz ) < (grid->radius*grid->radius) ) {
                        return 0;
                    }
                }
            }
        }
    }
    return 1;
}

/**
 *  Add a selected sample to the final list.
 */
static void addSample( struct griddata *grid, struct point *p )
{
    grid->samples[grid->nsample].x = p->x;
    grid->samples[grid->nsample].y = p->y;
    grid->samples[grid->nsample].z = p->z;
    grid->nsample++;

    /*  Add more space to samples, if needed. */
    if ( grid->nsample >= grid->ssize ) {
        grid->ssize = grid->ssize + CHUNK;
        grid->samples = realloc( grid->samples, grid->ssize * sizeof(struct point) );
    }

}

void poisson_disc( struct griddata *grid, int width, int height, int depth,
                   float radius, int k )
{
    grid->radius = radius;
    grid->cell_size = radius / sqrtf( 3.0f );
    grid->width = (int) ceilf( width / grid->cell_size );
    grid->height = (int) ceilf( height / grid->cell_size );
    grid->depth = (int) ceilf( depth / grid->cell_size );

    grid->cells = (struct point *)
        malloc( sizeof(struct point) * grid->height * grid->width * grid->depth );
    printf( "# Allocated %d cells\n", grid->height * grid->width * grid->depth );
    for ( int i = 0; i < grid->height * grid->width * grid->depth; i++ ) {
        grid->cells[i].x = -1.0;
    }

    /*  Queue for active list. */
    grid->queue = (struct point *) malloc( CHUNK * sizeof(struct point) );
    grid->qsize = CHUNK;
    grid->qlength = 0;

    /*  Space for results. */
    grid->samples = (struct point *) malloc( CHUNK * sizeof(struct point) );
    grid->ssize = CHUNK;
    grid->nsample = 0;

    /* Initialise a seed point. ... at centre ... */
    struct point ip;

    ip.x = randone() * width;
    ip.y = randone() * height;
    ip.z = randone() * depth;
    markPosition( grid, &ip );

    while ( grid->qlength > 0 ) {

        /*  Grab a position from the available queue. */
        int index = (int)( randone() * grid->qlength );
        struct point *p = &grid->queue[index];
        int havenew = 0;

        /*  Look for a random position that is far enough away and not already
         *  in use. */
        for ( int j = 0; j < k; j++ ) {
            float theta = 2.0f * M_PI * randone();
            float phi = 2.0f * M_PI * randone();

            /*  Radius to this position is outside the expected radius. */
            float r = randone() * radius + radius;
            struct point np;
            np.x = p->x + r * sin( theta ) * cos( phi );
            np.y = p->y + r * sin( theta ) * sin( phi );
            np.z = p->z + r * cos( theta );

            if ( np.x >= 0.0f && np.x < width &&
                 np.y >= 0.0f && np.y < height &&
                 np.z >= 0.0f && np.z < depth &&
                 available( grid, &np ) ) {
                markPosition( grid, &np );
                havenew = 1;
                break;
            }
        }

        /*  Remove point from active list and keep as no other position is
         *  near. */
        if ( ! havenew ) {
            addSample( grid, p );
            unqueuePosition( grid, index );
        }
    }
}


int main( int argc, char *argv[] )
{
    int width = 30;
    int height = 30;
    int depth =  30;
    float radius = 0.0f;
    int k = 30;
    struct griddata grid;

    /*  Expected sample size. */
    int N = 9;

    srand( time(NULL) );

    /*  Pick radius for expected sample size. */
    radius = pow( (width * height * depth) / (N), 0.3333 );
    printf( "# Radius = %f\n", radius );

    /*  Sample is stocastic, so we may need to ask more than one to get the
     *  number of samples we require as a minimum. */
    grid.nsample = 0;
    while ( grid.nsample < N ) {
        printf( "# Sampling...\n" );
        poisson_disc( &grid, width, height, depth, radius, k );
        printf( "# Samples = %d\n", grid.nsample );
    }

    for ( int i = 0; i < grid.nsample; i++ ) {
        printf( "# %f %f %f\n",
                grid.samples[i].x,
                grid.samples[i].y,
                grid.samples[i].z );
    }


    /*  Partition the space. Slow .... */

    unsigned long int counts[N];
    for ( int i = 0; i < N; i++ ) {
        counts[i] = 0;
    }

    for ( int i = 0; i < width; i++ ) {
        for ( int j = 0; j< height; j++ ) {
            for ( int k = 0; k < depth; k++ ) {
                int select = -1;
                float rsqmax = FLT_MAX;
                for ( int l = 0; l < N; l++ ) {
                    float dx = grid.samples[l].x - (i+0.5);
                    float dy = grid.samples[l].y - (j+0.5);
                    float dz = grid.samples[l].z - (k+0.5);
                    float rsq = dx*dx + dy*dy + dz*dz;
                    if ( rsq < rsqmax ) {
                        rsqmax = rsq;
                        select = l;
                    }
                }
                counts[select]++;
                printf( "%f %f %f %d\n", i+0.5, j+0.5, k+0.5, select );
            }
        }
    }

    printf( "# Counts:\n" );
    unsigned long int total = 0;
    for ( int i = 0; i < N; i++ ) {
        printf( "#  %d %ld\n", i, counts[i] );
        total += counts[i];
    }
    printf( "# total = %ld\n", total );

    return 0;
}
