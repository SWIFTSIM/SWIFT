/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Aidan B.G. Chalk (aidan.chalk@durham.ac.uk)
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

/**
 *  @file genetic_partitioning.c
 *  @brief file with a genetic algorithm method for partitioning the space across MPI ranks.
 *
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <strings.h>
#include <values.h>
/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <fenv.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "genetic_partitioning.h"

#define CAREFUL

struct vertex {
    idx_t ci, cj;
    idx_t cost;
};

struct partitions{
    idx_t* partition;
    idx_t* costs;
};

void sort_costs(idx_t *nodes, idx_t *costs, int size)
{
/* Sort partitions in ascending order using quickshed. */
    int lo = 0, hi = size-1;
    idx_t pivot = costs[(lo+hi)/2];
    int i = 0 , j = size-1;
    int temp;
    idx_t temp2;

    if(size <= 9)
    {
        for(i = 1; i < size; i++)
        {
            j = i;
            while(j > 0 && costs[j-1] < costs[j])
            {
                temp = nodes[j];
                temp2 = costs[j];
                nodes[j] = nodes[j-1];
                costs[j] = costs[j-1];
                nodes[j-1] = temp;
                costs[j-1] = temp2;
                j--;
            }
        }
        return;
    }
    while(i < j)
    {
        while (costs[i] > pivot) i++;
        while (costs[j] < pivot) j--;
        if(i <= j){
            if( i < j )
            {
                temp = nodes[i];
                temp2 = costs[i];
                nodes[i] = nodes[j];
                costs[i] = costs[j];
                costs[j] = temp;
                costs[j] = temp2;
            }
            i += 1;
            j -= 1;
        }
    }
    /* Recurse on the sublists.*/
        sort_costs(&nodes[0], &costs[0], j);
        sort_costs(&nodes[j+1], &costs[j+1], size-j-1);

}

void sort_partitions(struct partitions *parts, idx_t *costs, int size)
{
    /* Sort partitions in ascending order using quickshed. */
    int lo = 0, hi = size-1;
    double pivot = costs[(lo+hi)/2];
    int i = 0 , j = size-1;
    idx_t* temp;
    idx_t* temp3;
    idx_t temp2;

    /*if(size == 1)
        return;*/

    /* If the partition is small enough then just do insertion sort. */
    if(size <= 9)
    {
        for(i = 1; i < size; i++)
        {
            j = i;
            while(j > 0 && costs[j-1] > costs[j])
            {
                temp = parts[j].partition;
                temp3 = parts[j].costs;
                temp2 = costs[j];
                parts[j].partition = parts[j-1].partition;
                parts[j].costs = parts[j-1].costs;
                costs[j] = costs[j-1];
                parts[j-1].partition = temp;
                parts[j-1].costs = temp3;
                costs[j-1] = temp2;
                j--;
            }
        }
        return;
    }
    while(i < j)
    {
        while (costs[i] < pivot) i++;
        while (costs[j] > pivot) j--;
        if(i <= j){
            if( i < j )
            {
                temp = parts[i].partition;
                temp3 = parts[i].costs;
                temp2 = costs[i];
                parts[i].partition = parts[j].partition;
                parts[i].costs = parts[j].costs;
                costs[i] = costs[j];
                parts[j].partition = temp;
                parts[j].costs = temp3;
                costs[j] = temp2;
            }
            i += 1;
            j -= 1;
        }
    }
    /* Recurse on the sublists.*/
        sort_partitions(&parts[0], &costs[0], j);
        sort_partitions(&parts[j+1], &costs[j+1], size-j-1);
}

void check_symmetry(struct vertex *graph, int N)
{
    for(int i = 0; i < N*27; i++)
    {
        struct vertex *v = &graph[i];
        if(v->ci == v->cj)
            continue;
        idx_t cost = v->cost;
        idx_t ci = v->ci;
        idx_t cj = v->cj;
        int app = 0;
        for(int j = cj*27; j < cj*27 + 27; j++)
        {
            struct vertex *v2 = &graph[j];
            if(v2->cj == ci)
            {
                if(app)
                    message("cj, ci appeared twice for ci = %lli cj = %lli", ci, cj);
                app = 1;
                if(v2->cost != cost)
                    message("v2->cost = %lli v->cost = %lli", v2->cost, cost);
            }
        }
    }
}

/**N == nr_cells in the partition. **/
struct vertex* partition_create_graph( struct space *s, struct task *tasks, int nr_tasks)
{
#if defined(WITH_MPI) && defined(HAVE_METIS)
    int N = s->cdim[0] * s->cdim[1] * s->cdim[2];
    struct vertex *graph;
    graph = calloc(N*27, sizeof(struct vertex));
    if(graph == NULL)
        error("Failed to allocate graph");
    int i;
    //float wscale = 1e-3;


  /* Loop over all cells in the space. */
  int cid = 0;
  for (int l = 0; l < s->cdim[0]; l++) {
    for (int m = 0; m < s->cdim[1]; m++) {
      for (int n = 0; n < s->cdim[2]; n++) {
       graph[cid*27].ci = cell_getid(s->cdim, l, m, n);
        /* Visit all neighbours of this cell, wrapping space at edges. */
        int p = 0;
        for (int i = -1; i <= 1; i++) {
          int ii = l + i;
          if (ii < 0)
            ii += s->cdim[0];
          else if (ii >= s->cdim[0])
            ii -= s->cdim[0];
          for (int j = -1; j <= 1; j++) {
            int jj = m + j;
            if (jj < 0)
              jj += s->cdim[1];
            else if (jj >= s->cdim[1])
              jj -= s->cdim[1];
            for (int k = -1; k <= 1; k++) {
              int kk = n + k;
              if (kk < 0)
                kk += s->cdim[2];
              else if (kk >= s->cdim[2])
                kk -= s->cdim[2];

              /* If not self, record id of neighbour. */
              if (i || j || k) {
                graph[cid*27+p].ci = graph[27*cid].ci;
                graph[cid*27+p].cj = cell_getid(s->cdim, ii, jj, kk);
                p++;
              }else{
                /* If self record cj == ci */
                graph[cid*27+p].ci = graph[27*cid].ci;
                graph[cid*27+p].cj = graph[27*cid].ci;
                p++;
              }
            }
          }
        }

        /* Next cell. */
        cid++;
      }
    }
  }


//Is there a good way to store cj for a given ci? Also we store each pair twice! Oh well.
    for(int tid = 0; tid < nr_tasks; tid++){
        struct task *t = &tasks[tid];
        int ci, cj;
        struct cell *cell_i, *cell_j;
           /* Skip un-interesting tasks. */
        if (t->type != task_type_self && t->type != task_type_pair &&
            t->type != task_type_sub_self && t->type != task_type_sub_pair && t->type != task_type_ghost &&
            t->type != task_type_drift && t->type != task_type_kick1 && t->type != task_type_kick2 &&
            t->type != task_type_init)
            continue;

        /* Find the upper level cells corresponding to this task.*/
        for (cell_i = t->ci; cell_i->parent != NULL; cell_i = cell_i->parent)
          ;
        if (t->cj != NULL)
          for (cell_j = t->cj; cell_j->parent != NULL; cell_j = cell_j->parent)
            ;
        else
          cell_j = NULL;

        /* If we have only 1 cell, we only need to update a single node.*/
        {
            ci = cell_i - s->cells_top;
            if(cell_j != NULL)
                cj = cell_j - s->cells_top;
            else
                cj = cell_i - s->cells_top;            

            //Find where cj is in ci subarray.
            for(i = 0; i < 27; i++)
            {
                if(graph[ci*27+i].cj == cj)
                {
                    graph[ci*27+i].cost += (t->cost)/1000;
                    break;
                }
            }
            if(cj != ci)
            for(i = 0; i < 27; i++)
            {
                if(graph[cj*27+i].cj == ci)
                {
                    graph[cj*27+i].cost += (t->cost)/1000;
                    break;
                }
            }
           
        }
    }

    check_symmetry(graph, N);
    message("Symmetry checked");
    /* Graph should be now setup!*/
    return graph;
#else
printf("Not compiled with METIS or MPI support\n");
return 0;
#endif
}

void evaluate_partition(idx_t *partition, int N, int nr_nodes, idx_t *cost, struct vertex *graph, idx_t* cost_per_rank)
{

    idx_t max_cost = -1;
    for(int i = 0; i < nr_nodes; i++)
    {
        cost_per_rank[i] = 0.;
    }

    for(int i = 0; i < N*27; i++)
    {
        int pi = partition[graph[i].ci];
        int pj = partition[graph[i].cj];
        if(graph[i].cj < graph[i].ci)
            continue;
        if(pi == pj)
        {
            cost_per_rank[pi] += graph[i].cost;
            if(cost_per_rank[pi] > max_cost)
                max_cost = cost_per_rank[pi];
        }else{
            cost_per_rank[pi] += graph[i].cost;
            if(cost_per_rank[pi] > max_cost)
                max_cost = cost_per_rank[pi];
            cost_per_rank[pj] += graph[i].cost;
            if(cost_per_rank[pj] > max_cost)
                max_cost = cost_per_rank[pj];
        }
    }


    for(int i = 0; i < nr_nodes; i++)
    {
        if(cost_per_rank[i] == 0)
        {
            //message("cost per rank was 0 somehow");
            max_cost = LLONG_MAX;
            break;
        }
    }
    fflush(stdout);
    *cost = max_cost;
}

//#define M_LOG2E 1.44269504088896340736 //log2(e)

inline double log2(const double x){
    return  log(x) * M_LOG2E;
}

void mutate_partition(struct partitions *partition, /* TODO restrict*/struct vertex *graph, int N, int nr_nodes, struct partitions *new, idx_t *temp, idx_t *temp_costs)
{
    /* For now just find the max cost node*/
    idx_t *parray = partition->partition;
    idx_t* costs = partition->costs;
    for(int i = 0; i < nr_nodes; i++)
    {
        temp_costs[i] = costs[i];
        temp[i] = i;
    } 
    /* Sort nodes by cost.*/
    sort_costs(temp, temp_costs, nr_nodes);
    idx_t maxindex = 0;
    idx_t worstNode = -1;        
    /* Pick node psuedo randomly. First with chance 50%, second with chance 25%, so on....*/
    double r = (double)rand() / (double) RAND_MAX;
    maxindex = (int)(-log2(r));
    if(maxindex >= nr_nodes) maxindex = 0;
    worstNode = temp[maxindex];
   
    idx_t guess = (2*N)/nr_nodes;
    if(guess < 27)
        guess = 27;
    idx_t cells[guess*2];
    int nr_on_node = 0;

    /* Find all of the cells allocated to the worst node.*/    
    for(int i = 0; i < N && nr_on_node < guess*2; i++)
    {
        if(parray[i] == worstNode)
        {
            cells[nr_on_node++] = i;
        }
    }

    /* Pick one randomly.*/
    idx_t choice = -1;;
    if(nr_on_node >1)
{
    choice = rand() % nr_on_node;

    choice = cells[choice];

    /* Find the neighbouring cells that aren't on the same rank.*/
    nr_on_node = 0;
    for(int i = choice*27; i < (choice+1)*27 && nr_on_node < 27; i++)
    {
        if(graph[i].ci == choice && graph[i].cj != choice)
        {
            if(parray[graph[i].cj] != parray[graph[i].ci])
            {
                cells[nr_on_node++] = graph[i].cj;
            }
        }
    }
}else{
    nr_on_node = 0;
}
    if(nr_on_node == 0)
    {
        //TODO Make a decision to do something if no neighbours are on another node! Currently just creates the original.
//Craete a new partition that is the same as the original.
        memcpy( (void*)new->partition,(const void*)partition->partition, N*sizeof(idx_t));
        memcpy( (void*)new->costs, (const void*) partition->costs, nr_nodes*sizeof(idx_t));
       /* TODO Never seems to get here so not gonna worry.*/

    }else{
        idx_t choice2 = rand() % nr_on_node;
        //Craete a new partition that is the same as the original.
        memcpy( (void*)new->partition,(const void*) partition->partition, N*sizeof(idx_t));
        memcpy( (void*)new->costs, (const void*) partition->costs, nr_nodes*sizeof(idx_t));
        for(int i = choice*27; i < (choice+1)*27; i++)
        {
            idx_t ci = graph[i].ci;
            if(ci != choice)
                error("We dumb");
            idx_t cj = graph[i].cj;
            idx_t ck = cells[choice2];
            idx_t pi = new->partition[ci];
            idx_t pj = new->partition[cj];
            idx_t pk = new->partition[ck];
            if(ci == cj)
            {
                new->costs[pi] -= graph[i].cost;
                new->costs[pk] += graph[i].cost;
                continue;
            }
            if(pi != pj )
            {
                new->costs[pi] -= graph[i].cost;
            }
            if( pj != pk)
            {
                new->costs[pk] += graph[i].cost;
            }
        }
        new->partition[choice] = new->partition[cells[choice2]];
    }
}

idx_t quick_eval(idx_t* cost_per_rank, int nr_nodes)
{
    double max = -1.;
    for(int i = 0; i < nr_nodes; i++)
    {
        if(cost_per_rank[i] > max)
        {
            max = cost_per_rank[i];
        }
    }
    return max;
}

int hash_partition(idx_t* partition, int nr_nodes, int nr_cells, struct vertex *graph)
{
    idx_t hash=0;
    for(int i = 0; i < nr_cells; i++)
    {
        hash += partition[i] * graph[i].cost;
    }
    return hash;
}

void random_partition(idx_t* partition, int nr_nodes, int nr_cells )
{
    for(int i = 0; i < nr_cells; i++)
    {
        partition[i] = -1;
    }
    for(int i = 0; i < nr_nodes; i++)
    {
        idx_t r = rand() % nr_cells;
        while(partition[r] >= 0)
        {
            r = rand() % nr_cells;
        }
        partition[r] = i;
    }
    for(int i = 0; i < nr_cells; i++)
    {
        if(partition[i] < 0)
        {
            partition[i] = rand() % nr_nodes;
        }
    }
}

void genetic_algorithm(struct space *s, struct task *tasks, int nr_tasks, idx_t *initial_partition, int nr_nodes)
{
#if defined(WITH_MPI) && defined(HAVE_METIS)
ticks tic = getticks();

struct vertex *graph = partition_create_graph(s, tasks, nr_tasks);

ticks toc = getticks();

message("Graph creation took %.3f %s.", clocks_from_ticks(toc - tic),clocks_getunit());
/*for(int i = 0; i < 27*s->nr_cells; i++)
{
    printf("Node[%i] ci = %i cj %i cost = %f\n", i, graph[i].ci, graph[i].cj, graph[i].cost);
}*/
//#define OUTPUT
#ifdef OUTPUT
FILE *file = fopen("/cosma/home/Virgo/d74ksy/swiftsim/src/partout.out", "w");
#endif

tic = getticks();
struct partitions initial;
initial.partition = initial_partition;
idx_t initial_cost;


//int test[]={ 2, 2, 1, 0, 4, 2, 3, 2, 4, 3, 3, 4, 8, 2, 1, 0, 6, 5, 5, 6, 9, 7, 8, 9, 9, 7, 8, 8, 6, 5, 5, 6, 7, 7, 5, 6, 9, 8, 8, 9, 9, 8, 8, 9, 7, 7, 5, 6, 0, 1, 1, 0, 4, 3, 3, 4, 4, 2, 3, 4, 0, 1, 1, 0 };

//12th element


//int best_cost;
//struct partitions *best_partition;
int since_improve = 0;

srand(time(NULL));

idx_t *temp = malloc(sizeof(idx_t) * nr_nodes);
if(temp == NULL)
    error("failed to allocate temp");

idx_t *temp_costs = malloc(sizeof(idx_t) * nr_nodes);
if(temp_costs == NULL)
    error("failed to allocate temp_costs");

initial.costs = malloc(nr_nodes*sizeof(idx_t));
if(initial.costs == NULL)
    error("Failed to allocate initial.costs");

/*evaluate_partition((int*)test, s->nr_cells, nr_nodes, &initial_cost, graph, (initial.costs));
message("test cost is %f", initial_cost);

printf("costs per node = [ ");
    for(int i = 0; i < nr_nodes; i++)
    {
        printf("%f ", initial.costs[i]);
    }
    printf("];\n");*/
//costs per node = [ 17397.833836 17907.148839 21065.156989 15360.329743 18486.952945 15420.812744 14943.531688 15977.754752 21729.249023 18747.356897 ];
//temp costs per node = [ 17397.833836 17907.148839 21658.795009 15360.329743 18486.952945 15420.812744 14943.531688 15977.754752 21923.904029 18747.356897 ];
evaluate_partition(initial.partition, s->nr_cells, nr_nodes, &initial_cost, graph, (initial.costs));
idx_t sum =0;
message("METIS cost is %lli", initial_cost);
for(int i = 0; i < nr_nodes; i++)
{
    sum += initial.costs[i];
}
message("METIS sum cost is %lli", sum);

/*for(int i = 0; i < s->nr_cells; i++)
{
        initial.partition[i] = i / (s->nr_cells/nr_nodes);
        if(initial.partition[i] >= nr_nodes)
            initial.partition[i] = nr_nodes-1;
     //   printf("%i ", initial[i]);
}
evaluate_partition(initial.partition, s->nr_cells, nr_nodes, &initial_cost, graph, (initial.costs));
message("Initial cost is %f", initial_cost);
*/

/*printf("costs per node = [ ");
    for(int i = 0; i < nr_nodes; i++)
    {
        sum += initial.costs[i];
        printf("%f ", initial.costs[i]);
    }
    printf("];\n");
    message("Sum of costs = %f", sum);*/

//Setup the partition array.
int pop_size = 100;
struct partitions *population;
population = malloc(sizeof(struct partitions) * pop_size);
if(population == NULL)
    error("Failed to allocate population");
idx_t *population_costs = malloc(sizeof(idx_t) * pop_size);
if(population_costs == NULL)
    error("Failed to allocate population costs array");

//Create initial population
//double *tempcheck = malloc(sizeof(double)*nr_nodes);
//double tempcost;
for(int i = 0; i < pop_size; i++)
{
    population[i].partition = malloc(sizeof(idx_t) * s->nr_cells);
    population[i].costs = malloc(nr_nodes*sizeof(idx_t));
    if(population[i].partition == NULL || population[i].costs == NULL)
        error("Failed to allocate population arrays");
    mutate_partition(&initial, graph, s->nr_cells, nr_nodes, &population[i], temp, temp_costs);
    population_costs[i] = quick_eval(population[i].costs, nr_nodes);
/*    random_partition(population[i].partition, nr_nodes, s->nr_cells);*/
    evaluate_partition(population[i].partition, s->nr_cells, nr_nodes, &population_costs[i], graph, population[i].costs);
    if(i < 10)
{
/*    printf("Hash = %i, cost = %lli, {", hash_partition(population[i].partition,nr_nodes, s->nr_cells, graph), population_costs[i]);
    for(int k = 0; k < s->nr_cells; k++)
    {
        printf("%lli ", population[i].partition[k]);
    }
    printf("}\n");*/
}
/*    evaluate_partition(population[i].partition, s->nr_cells, nr_nodes, &tempcost, graph, tempcheck);
    if(tempcost != population_costs[i])
        error("Failed at %i %f != %f", i,tempcost, population_costs[i]);*/
}
random_partition(population[pop_size-1].partition, nr_nodes, s->nr_cells);
evaluate_partition(population[pop_size-1].partition, s->nr_cells, nr_nodes, &population_costs[pop_size-1], graph, population[pop_size-1].costs);
/*printf("Hash = %i, cost = %lli, {", hash_partition(population[pop_size-1].partition,nr_nodes, s->nr_cells, graph), population_costs[pop_size-1]);
    for(int k = 0; k < s->nr_cells; k++)
    {
        printf("%lli ", population[pop_size-1].partition[k]);
    }
    printf("}\n");*/
random_partition(population[pop_size-1].partition, nr_nodes, s->nr_cells);
evaluate_partition(population[pop_size-1].partition, s->nr_cells, nr_nodes, &population_costs[pop_size-1], graph, population[pop_size-1].costs);
/*printf("Hash = %i, cost = %lli, {", hash_partition(population[pop_size-1].partition,nr_nodes, s->nr_cells, graph), population_costs[pop_size-1]);
    for(int k = 0; k < s->nr_cells; k++)
    {
        printf("%lli ", population[pop_size-1].partition[k]);
    }
    printf("}\n");*/
random_partition(population[pop_size-1].partition, nr_nodes, s->nr_cells);
evaluate_partition(population[pop_size-1].partition, s->nr_cells, nr_nodes, &population_costs[pop_size-1], graph, population[pop_size-1].costs);
/*printf("Hash = %i, cost = %lli, {", hash_partition(population[pop_size-1].partition,nr_nodes, s->nr_cells, graph), population_costs[pop_size-1]);
    for(int k = 0; k < s->nr_cells; k++)
    {
        printf("%lli ", population[pop_size-1].partition[k]);
    }
    printf("}\n");*/
//message("Trying to sort population");
sort_partitions(population, population_costs, pop_size);

#ifdef OUTPUT
for(int i = 0; i < pop_size-1; i++)
{
    fprintf(file, "%.3f, ", population_costs[i]);
}
fprintf(file, "%.3f\n", population_costs[pop_size-1]);
#endif

//message("sorted population");
//message("population cost 0 i.e. best is %f", population_costs[0]);
//if(population_costs[0] < initial_cost)
//{
    initial_cost = population_costs[0];
    memcpy(initial.costs, population[0].costs, sizeof(idx_t) * nr_nodes);
    memcpy(initial.partition, population[0].partition, sizeof(idx_t) * s->nr_cells);
//}

//message("Created initial population");
double c = 50.0;
int j = 0;
while((since_improve < 250 && j < 5000) || j < 1000) 
{
    double max_diff = (double) (population_costs[pop_size-1] - population_costs[10]);
  /* Go through the last 90 and choose 40 to keep, use probabilistic based shuffle */
    for(int i = 10; i < pop_size; i++)
    {
        int swap = rand() % (pop_size-10);
        swap += 10;
        while(swap == i)
        {
            swap = rand() % (pop_size-10);
            swap += 10;
        }
        if(population_costs[i] == LLONG_MAX || population_costs[swap] == LLONG_MAX)
            continue;
        double diff = (double)(population_costs[i] - population_costs[swap]);
        /*if(diff < 0)
        {*/
            /* Swap them. */
/*            int *temp4 = population[i].partition;
            double *temp2 = population[i].costs;
            double temp3 = population_costs[i];
            population[i].partition = population[swap].partition;
            population[i].costs = population[swap].costs;
            population_costs[i] = population_costs[swap];
            population[swap].partition = temp4;
            population[swap].costs = temp2;
            population_costs[swap] = temp3;
            diff *= -1;        
        }*/
        /* Work out difference relative to max_diff. Max diff is 5 (c - X)% chance, no diff is 50% chance.*/
        double m = (-45.0) / max_diff;
        
        double p = m * diff + c;
        int r = rand() % 10001;
        double r2 = (double)r;
        r2 = r2 / 100.0;
        if(r2 > p)
        {
            idx_t *temp4 = population[i].partition;
            idx_t *temp2 = population[i].costs;
            idx_t temp3 = population_costs[i];
            population[i].partition = population[swap].partition;
            population[i].costs = population[swap].costs;
            population_costs[i] = population_costs[swap];
            population[swap].partition = temp4;
            population[swap].costs = temp2;
            population_costs[swap] = temp3;
        }

    }

    /* Avoid all the best partitions being the same?*/
    for(int i = 1; i < 10; i++)
    {
        //TODO Compare using hash.
        if(population_costs[i] == population_costs[0])
        {
            int rande = rand() % 40 +10;
            mutate_partition(&population[rande], graph, s->nr_cells, nr_nodes, &population[i], temp, temp_costs);
            population_costs[i] = quick_eval(population[i].costs, nr_nodes);
        }
    }
    /* Generate new partitions. */
    for(int i = 50; i < pop_size; i++)
    {
        int rande = rand() % 50;
        mutate_partition(&population[rande], graph, s->nr_cells, nr_nodes, &population[i], temp, temp_costs);
        population_costs[i] = quick_eval(population[i].costs, nr_nodes);
/*    evaluate_partition(population[i].partition, s->nr_cells, nr_nodes, &tempcost, graph, tempcheck);
    if(tempcost != population_costs[i])
        error("Failed at %i %f != %f", i,tempcost, population_costs[i]);*/
    }

    if(j%100 == 0)
    {
        for(int i = 0; i < pop_size; i++)
        {
            evaluate_partition(population[i].partition, s->nr_cells, nr_nodes, &population_costs[i], graph, population[i].costs);
        }
    }
    /* Resort the partitions. */ 
sort_partitions(population, population_costs, pop_size);
#ifdef OUTPUT
for(int i = 0; i < pop_size-1; i++)
{
    fprintf(file, "%.3f, ", population_costs[i]);
}
fprintf(file, "%.3f\n", population_costs[pop_size-1]);
#endif
    if(population_costs[0] < initial_cost)
    {
        if(j % 100)
            evaluate_partition(population[0].partition, s->nr_cells, nr_nodes, &population_costs[0], graph, population[0].costs);
        if(population_costs[0] < initial_cost)
        {
        initial_cost = population_costs[0];
        memcpy(initial.costs, population[0].costs, sizeof(idx_t) *  nr_nodes);
        memcpy(initial.partition, population[0].partition, sizeof(idx_t) * s->nr_cells);
        since_improve = -1;
        }
    }
since_improve++;
j++;
}
    message("Step count = %i, since imrpove = %i", j, since_improve);
    message("Best result = %lli", initial_cost);
/*    printf("final = [ ");
    for(int i = 0; i < s->nr_cells; i++)
    {
        printf("%lli ", initial.partition[i]);
    }
    printf("];\n");*/
    sum =0.;
printf("costs per node = [ ");
    for(int i = 0; i < nr_nodes; i++)
    {
        sum += initial.costs[i];
        printf("%lli ", initial.costs[i]);
    }
    printf("];\n");
    message("Sum of costs = %lli", sum);

toc = getticks();

message("Genetic algorithm took %.3f %s.", clocks_from_ticks(toc - tic), clocks_getunit());

    for(int i = 0; i < s->nr_cells; i++)
    {
        initial_partition[i] = initial.partition[i];
    }

#ifdef OUTPUT
fclose(file);
#endif
#else
  printf("Not compiled with MPI or METIS support\n");
#endif
}

/*

Graph partitioning for SWIFT
Init data
First we want to create a graph of the costs at the top-level cells.
struct {
  int ci, cj;
  Double cost;
} graph[N*27];
We fill this list by iterating over all the tasks:
for (int tid = 0; tid < nr_tasks; tid++) {
  Struct task *t = …;
  Int ci = highest-level parent of t->ci;
  Int cj = highest-level parent of t->cj;
  Int index = offset in graph[] for ci/cj.
  Graph[index].cost += t->cost;
}
To evaluate the cost for each node in a partition, do the following: Can store max cost while doing this too!
Int part[N]
double cost_per_rank[nr_nodes];
For (int i = 0; i < N * 27; i++) {
  Int pi = part[graph[i].ci];
  Int pj = part[graph[i].cj];
  If (pi == pj)
    Cost_per_rank[pi] += graph[i].cost;
  Else {
    Cost_per_rank[pi] += graph[i].cost;
    Cost_per_rank[pj] += graph[i].cost;
  }
}


Fast lookup of neighbours of ci using graph. Large population - don’t be too strict on removing things. Pick solutions to keep probabilistically - don’t just prune “the worst”.

Mutation - probabilistically pick a node, 50% worst, 25% 2nd worst etc. move a cell from that to be adjacent to one of its non-local neighbours.
All nodes must have cells - can be enforced. They have INT_MAX solution value.

Run bigcosmovolume for 10 steps, output tasks as csv with higher level cells only. Partition it and check results - check max cost, load imbalance, reliability - modify costs see how badly it goes, lots of runs and check for divergence, feed in result and see if it diverges. How many iterations do we need to wait with no change before we stop etc.

TODO Try to convert the costs to integers? Might save slightly on recomputing various results.
*/
#ifdef SKIPQWR
int main()
{
    double *costs = malloc(100*sizeof(double));
    struct partitions *p = malloc(100*sizeof(struct partitions));
    if(p == NULL)
        error("Failed to allocate partitions");
    srand(6);
    for(int i = 0; i < 100; i++)
    {
        costs[i] = DBL_MAX;
        p[i].partition = (int*)malloc(sizeof(int));
        p[i].costs = (double*)malloc(sizeof(double));
    }

    sort_partitions(p, costs, 100);
    for(int i = 0; i < 100; i++)
    {
       printf("%f ", costs[i]);
        free(p[i].partition);
        free(p[i].costs);
    }   
    printf("\n");
    free( costs);
    free(p);
}
#endif

#ifdef SKIPQWR
int main( int argc, char *argv[])
{
struct part *parts = NULL;
struct gpart *gparts = NULL;
size_t Ngas = 0, Ngpart = 0;
int periodic = 1;
double dim[3] = {1.0, 1.0, 1.0};
struct space s;
struct engine e;
//double h_max = 0.0625;
double h_max = 0.3;
float dt_max = 0.0002f, dt_min = 0.0001f;
double time_end = 0.5;

#ifdef WITH_MPI
int res = 0, prov = 0;
  if ((res = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov)) !=
      MPI_SUCCESS)
    error("Call to MPI_Init failed with error %i.", res);
  if (prov != MPI_THREAD_MULTIPLE)
    error(
        "MPI does not provide the level of threading required "
        "(MPI_THREAD_MULTIPLE).");
#endif
#ifdef WITH_MPI
#if defined(HAVE_PARALLEL_HDF5)
  read_ic_parallel("/cosma/home/Virgo/d74ksy/swiftsim/examples/UniformBox/uniformBox.hdf5", dim, &parts, &gparts, &Ngas, &Ngpart, &periodic,
		   0, 1, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
  read_ic_serial("/cosma/home/Virgo/d74ksy/swiftsim/examples/UniformBox/uniformBox.hdf5", dim, &parts, &gparts, &Ngas, &Ngpart, &periodic,
                 0, 1, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
#else
read_ic_single("/cosma/home/Virgo/d74ksy/swiftsim/examples/UniformBox/uniformBox.hdf5", dim, &parts, &gparts, &Ngas, &Ngpart, &periodic);
#endif


message("Read %lld gas particles and %lld DM particles from the ICs", (long long int)Ngas, (long long int)(Ngpart-Ngas));

  /* Init the space. */
  bzero(&s, sizeof(struct space));
space_init(&s, dim, parts, gparts, Ngas, Ngpart, periodic, h_max, 1);

  /* Say a few nice things about the space we just created. */
    message("space dimensions are [ %.3f %.3f %.3f ].", s.dim[0], s.dim[1],
            s.dim[2]);
    message("space %s periodic.", s.periodic ? "is" : "isn't");
    message("highest-level cell dimensions are [ %i %i %i ].", s.cdim[0],
            s.cdim[1], s.cdim[2]);
    message("%zi parts in %i cells.", s.nr_parts, s.tot_cells);
    message("%zi gparts in %i cells.", s.nr_gparts, s.tot_cells);
    message("maximum depth is %d.", s.maxdepth);
    // message( "cutoffs in [ %g %g ]." , s.h_min , s.h_max ); fflush(stdout);
s.nr_gparts = 0;
  int engine_policies = 0 | engine_policy_steal | engine_policy_hydro;

    engine_init(&e, &s, dt_max, 1, 1, 1, 0,
              engine_policies, 0, time_end, dt_min, dt_max, 0);

#define nodes 20
#ifdef WITH_MPI
  struct partition initial_partition;

//  initial_partition.type = INITPART_GRID;
    initial_partition.type = INITPART_METIS_WEIGHT;
  initial_partition.grid[0] = nodes;
  initial_partition.grid[1] = 1;
  initial_partition.grid[2] = 1;
#endif

/* Initialise the particles */
  engine_init_particles(&e);
    message("Done engine_init");
    srand(7);
    /* Add random costs to all of the tasks! */
    for(int i = 0; i < e.sched.nr_tasks; i++)
    {
        e.sched.tasks[i].tic = (ticks)12345;
        e.sched.tasks[i].toc = (ticks)12345 + (rand() % 100000);
    }

#ifdef WITH_MPI
ticks tic = getticks();
partition_initial_partition(&initial_partition,
                                 0, nodes, &s);
int initial[s.nr_cells];

printf("initial = { ");
for(int i = 0; i < s.nr_cells; i++)
{
            initial[i] = s.cells[i].nodeID;
            printf("%lli, ", s.cells[i].nodeID);    
}
printf("}\n");

ticks toc = getticks();
message("Initial partition took %.3f %s.", clocks_from_ticks(toc - tic), clocks_getunit());
#else

    int initial[s.nr_cells];
    printf("initial = [ ");
    for(int i = 0; i < s.nr_cells; i++)
    {
        initial[i] = i / (s.nr_cells/10);
        if(initial[i] >= 10)
            initial[i] = 9;
        printf("%i ", initial[i]);
    }
    printf("];\n");
#endif    
    genetic_algorithm(&s, &e.sched, initial, nodes);

}
#endif
