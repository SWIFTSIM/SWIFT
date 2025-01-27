#ifndef CUDA_HEADERS_H
#define CUDA_HEADERS_H
#define n_streams 1024

#ifdef __cplusplus
extern "C" {
#endif
#include "part_gpu.h"
void launch_density_kernel(struct part_soa parts_soa, int *d_task_first_part,
                           int *d_task_last_part, int *d_bundle_first_part,
                           int *d_bundle_last_part, float d_a, float d_H,
                           const char *loop_type, hipStream_t stream, int bid,
                           int block_size, int count_tasks, int tasksperbundle,
                           int numBlocks_x, int numBlocks_y, int tid,
                           int offset, int bundle_first_task, int max_parts,
                           int max_active_bin);

#ifdef __cplusplus
}
#endif

#endif  // CUDA_HEADER_H
