#ifndef CUDA_HEADERS_H
#define CUDA_HEADERS_H
#define n_streams 1024

#ifdef __cplusplus
extern "C" {
#endif
#include "part_gpu.h"
void launch_density_kernel(struct part_soa parts_soa, int *d_task_first_part,
                           int *d_task_last_part, float d_a, float d_H,
                           const char *loop_type, cudaStream_t stream,
                           int block_size, int count_tasks, int tasksperbundle,
                           int numBlocks_x, int numBlocks_y,
                           int bundle_first_task, int max_parts);
void launch_density_aos(struct part_aos *parts_aos, int *d_task_first_part,
                           int *d_task_last_part, float d_a, float d_H,
                           const char *loop_type, cudaStream_t stream,
                           int block_size, int count_tasks, int tasksperbundle,
                           int numBlocks_x, int numBlocks_y,
                           int bundle_first_task, int max_parts,
						   double * d_cell_x,
						   double * d_cell_y, double * d_cell_z);
void launch_density_aos_f4(struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv, float d_a, float d_H,
        cudaStream_t stream, int numBlocks_x, int numBlocks_y,
        int bundle_first_task, int2 *d_task_first_part_f4);
void launch_gradient_aos(struct part_aos_g *parts_aos, int *d_task_first_part,
                           int *d_task_last_part, float d_a, float d_H,
                           const char *loop_type, cudaStream_t stream,
                           int block_size, int count_tasks, int tasksperbundle,
                           int numBlocks_x, int numBlocks_y,
                           int bundle_first_task, int max_parts,
                           double * d_cell_x,
							double * d_cell_y, double * d_cell_z);
void launch_gradient_aos_f4(struct part_aos_f4_g_send *parts_send, struct part_aos_f4_g_recv *parts_recv, float d_a, float d_H,
        cudaStream_t stream, int numBlocks_x, int numBlocks_y,
        int bundle_first_task, int2 * d_task_first_part_f4);
void launch_force_aos(struct part_aos_f *parts_aos, int *d_task_first_part,
                           int *d_task_last_part, float d_a, float d_H,
                           const char *loop_type, cudaStream_t stream,
                           int block_size, int count_tasks, int tasksperbundle,
                           int numBlocks_x, int numBlocks_y,
                           int bundle_first_task, int max_parts,
                           double * d_cell_x,
							double * d_cell_y, double * d_cell_z);
void launch_force_aos_f4(struct part_aos_f4_f_send *parts_send, struct part_aos_f4_f_recv *parts_recv, float d_a, float d_H,
        cudaStream_t stream, int numBlocks_x, int numBlocks_y,
        int bundle_first_task, int2 * d_task_first_part_f4);
void launch_density_pair_two_kernels(struct part_soa parts_soa_ci, struct part_soa parts_soa_cj, int *d_task_first_part_ci,
	          int *d_task_first_part_cj, int *d_task_last_part_ci, int *d_task_last_part_cj, float d_a, float d_H,
              const char *loop_type, cudaStream_t stream, int bid, int block_size, int count_tasks, int tasksperbundle,
              int max_parts_i, int max_parts_j, int numBlocks_y, int tid, int offset, int bundle_first_task, int max_active_bin);
void runner_dopair1_branch_density_gpu(struct part_soa parts_soa, int *d_task_first_parts_pair,
        int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
		  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
		  int numBlocks_y, int tid, int offset, int bundle_first_task, int max_active_bin, double * d_shift_x,
		  double * d_shift_y, double * d_shift_z);
void runner_dopairci_branch_density_gpu(struct part_soa parts_soa, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z);
void runner_dopaircj_branch_density_gpu(struct part_soa parts_soa, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z);
void runner_dopairci_branch_density_gpu_aos(struct part_aos *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z);
void runner_dopairci_branch_density_gpu_aos_f4(struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv,
		      float d_a, float d_H, cudaStream_t stream,
			  int numBlocks_x, int numBlocks_y, int bundle_first_task, int4 *fparti_fpartj_lparti_lpartj_dens);
void runner_dopaircj_branch_density_gpu_aos_f4(struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv,
		      float d_a, float d_H, cudaStream_t stream,
			  int numBlocks_x, int numBlocks_y, int bundle_first_task, int4 *fparti_fpartj_lparti_lpartj_dens);
void runner_dopair_branch_density_gpu_aos_f4(struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv,
		      float d_a, float d_H, cudaStream_t stream,
			  int numBlocks_x, int numBlocks_y, int bundle_first_part, int bundle_n_parts);
void runner_dopaircj_branch_density_gpu_aos(struct part_aos *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z);
void runner_dopairci_branch_density_gpu_aos_g(struct part_aos_g *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z);
void runner_dopaircj_branch_density_gpu_aos_g(struct part_aos_g *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z);
void runner_dopair_branch_gradient_gpu_aos_f4(struct part_aos_f4_g_send *parts_send, struct part_aos_f4_g_recv *parts_recv,
		      float d_a, float d_H, cudaStream_t stream,
			  int numBlocks_x, int numBlocks_y, int bundle_first_part, int bundle_n_parts);
void runner_dopairci_branch_density_gpu_aos_f(struct part_aos_f *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z);
void runner_dopaircj_branch_density_gpu_aos_f(struct part_aos_f *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z);
void runner_dopair_branch_force_gpu_aos_f4(struct part_aos_f4_f_send *parts_send, struct part_aos_f4_f_recv *parts_recv,
		      float d_a, float d_H, cudaStream_t stream,
			  int numBlocks_x, int numBlocks_y, int bundle_first_part, int bundle_n_parts);
#ifdef __cplusplus
}
#endif

#endif // CUDA_HEADER_H
