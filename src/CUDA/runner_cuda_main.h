#ifndef SWIFT_CUDA_RUNNER_MAIN_H
#define SWIFT_CUDA_RUNNER_MAIN_H

//We want 100% occupancy I think.
#define num_blocks 128
#define num_cuda_threads 128

void create_tasks(struct engine *e);
void update_tasks(struct engine *e);
//void create_cells_and_data_tasks(struct engine *e);
void run_cuda();

#ifdef __cplusplus
extern "C" {
#endif
void test_27_cells(struct cell **cells, struct cell *main_cell, struct part *parts );
void allocate_cells(void **parts, int particles, int cells);
void allocate_cell( void **cell );
void free_parts(void *parts);
void free_cell( void *cell );
void test_125_cells( struct cell **cells, struct cell *main_cell, struct part *parts, struct engine *e);
#ifdef __cplusplus
}
#endif


#endif
