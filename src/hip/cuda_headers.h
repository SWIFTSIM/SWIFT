#ifndef CUDA_HEADERS_H
#define CUDA_HEADERS_H
#define n_streams 1024

#ifdef WITH_CUDA
extern "C" {
#endif

void GPU_runner_doself1_branch_gradient(struct cell_gpu *ci_gpu,
                                        struct part_gpu *parts_gpu);
void cuda_tester(struct cell **ci_list_mgd, int numBlocksTest,
                 int block_size_test, int count_tasks);
void launch_cuda_kernel(struct cell_gpu *ci_gpu, struct part_gpu *parts,
                        int numBlocks, float d_a, float d_H,
                        const char *loop_type);
void launch_cuda_kernel_streams(struct part_gpu *d_parts, int numBlocks,
                                float d_a, float d_H, const char *loop_type,
                                cudaStream_t stream, int tid, int count,
                                int max_count, float cellx, float celly,
                                float cellz, int first_part, int last_part);
void launch_cuda_kernel_bundles(struct cell_gpu *d_all_cells,
                                struct part_gpu **d_all_parts, int numBlocks,
                                float d_a, float d_H, const char *loop_type,
                                cudaStream_t stream, int bid, int block_size,
                                int count_tasks, int tasksperbundle,
                                int numBlocks_x, int numBlocks_y, int tid,
                                int offset);
void launch_cuda_kernel_bundles_revised(
    struct part_gpu *d_all_parts, int *d_task_first_part, int *d_task_last_part,
    int *d_bundle_first_part, int *d_bundle_last_part, int numBlocks, float d_a,
    float d_H, const char *loop_type, cudaStream_t stream, int bid,
    int block_size, int count_tasks, int tasksperbundle, int numBlocks_x,
    int numBlocks_y, int tid, int offset);
void launch_cuda_kernel_bundles_revised_soa(
    struct part_soa parts_gpu_soa, int *d_task_first_part,
    int *d_task_last_part, int *d_bundle_first_part, int *d_bundle_last_part,
    int numBlocks, float d_a, float d_H, const char *loop_type,
    cudaStream_t stream, int bid, int block_size, int count_tasks,
    int tasksperbundle, int numBlocks_x, int numBlocks_y, int tid, int offset,
    int bundle_first_task, int max_parts);
void launch_cuda_print_streams(int numBlocks, cudaStream_t stream, int tid);
void launch_cuda_kernel_tester(struct cell_gpu *d_ci_gpu,
                               struct part_gpu **d_parts, int numBlocks,
                               float d_a, float d_H, const char *loop_type,
                               cudaStream_t stream, int bid, int block_size,
                               int count_tasks, int tasksperbundle,
                               int numBlocks_x, int numBlocks_y, int tid);
void launch_cuda_kernel_bundles_test(struct cell_gpu *d_all_cells,
                                     struct part_gpu **d_all_parts,
                                     int numBlocks, float d_a, float d_H,
                                     int count_tasks);
void mgd_mem_cuda_kernel_bundles(struct part_gpu **parts_gpu_list,
                                 int numBlocks, float d_a, float d_H,
                                 const char *loop_type, cudaStream_t stream,
                                 int bid, int block_size, int count_tasks,
                                 int tasksperbundle, int numBlocks_x,
                                 int numBlocks_y, int tid, int offset);

#ifdef WITH_CUDA
}
#endif

#endif // CUDA_HEADER_H
