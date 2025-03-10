#include "scheduler.h"
#include "runner_doiact_hydro.h"
#include "active.h"
#include <atomic.h>
struct pack_vars_self {
  /*List of tasks and respective cells to be packed*/
  struct task **task_list;
  struct task **top_task_list;
  struct cell **cell_list;
  /*List of cell positions*/
  double *cellx;
  double *celly;
  double *cellz;
  /*List of cell positions*/
  double *d_cellx;
  double *d_celly;
  double *d_cellz;
  int bundle_size;
  /*How many particles in a bundle*/
  int count_parts;
  /**/
  int tasks_packed;
  int top_tasks_packed;
  int *task_first_part;
  int *task_last_part;
  int *d_task_first_part;
  int *d_task_last_part;
  int *bundle_first_part;
  int *bundle_last_part;
  int *bundle_first_task_list;
  int count_max_parts;
  int launch;
  int launch_leftovers;
  int target_n_tasks;
  int nBundles;
  int tasksperbundle;

} pack_vars_self;
struct leaf_cell_list{
  struct cell **ci;
  struct cell **cj;
  int n_leaves;
  int n_packed;
};
struct pack_vars_pair {
  /*List of tasks and respective cells to be packed*/
  struct task **task_list;
  struct task **top_task_list;
  struct leaf_cell_list * leaf_list;
  struct cell **ci_list;
  struct cell **cj_list;
  /*List of cell shifts*/
  double *shiftx;
  double *shifty;
  double *shiftz;
  /*List of cell shifts*/
  double *d_shiftx;
  double *d_shifty;
  double *d_shiftz;
  int bundle_size;
  /*How many particles in a bundle*/
  int count_parts;
  /**/
  int tasks_packed;
  int top_tasks_packed;
  int *task_first_part;
  int *task_last_part;
  int *d_task_first_part;
  int *d_task_last_part;
  int *bundle_first_part;
  int *bundle_last_part;
  int *bundle_first_task_list;
  int count_max_parts;
  int launch;
  int launch_leftovers;
  int target_n_tasks;
  int nBundles;
  int tasksperbundle;
  int task_locked;

} pack_vars_pair;

struct pack_vars_pair_f4 {
  /*List of tasks and respective cells to be packed*/
  struct task **task_list;
  struct cell **ci_list;
  struct cell **cj_list;
  /*List of cell shifts*/
  float3 *shift;
  /*List of cell shifts*/
  float3 *d_shift;
  int bundle_size;
  /*How many particles in a bundle*/
  int count_parts;
  /**/
  int tasks_packed;
  int4 *fparti_fpartj_lparti_lpartj;
  int4 *d_fparti_fpartj_lparti_lpartj;
  int *bundle_first_part;
  int *bundle_last_part;
  int *bundle_first_task_list;
  int count_max_parts;
  int launch;
  int launch_leftovers;
  int target_n_tasks;
  int nBundles;
  int tasksperbundle;

} pack_vars_pair_f4;

#include "cuda/BLOCK_SIZE.h"
#include "cuda/GPU_runner_functions.h"
#include "runner_gpu_pack_functions.h"
#include "task.h"
#define CUDA_DEBUG

double runner_doself1_pack_f4(struct runner *r, struct scheduler *s,
                              struct pack_vars_self *pack_vars, struct cell *ci,
                              struct task *t,
                              struct part_aos_f4_send *parts_send,
                              int2 *task_first_part_f4) {
  /* Timers for how long this all takes.
   * t0 and t1 are from start to finish including GPU calcs
   * tp0 and tp1 only time packing and unpacking*/
  struct timespec t0, t1;  //
  clock_gettime(CLOCK_REALTIME, &t0);
  /* Find my queue for use later*/
  int qid = r->qid;
  /*Place pointers to the task and cells packed in an array for use later
   * when unpacking after the GPU offload*/
  int tasks_packed = pack_vars->tasks_packed;
  pack_vars->task_list[tasks_packed] = t;
  pack_vars->cell_list[tasks_packed] = ci;
  /* Identify row in particle arrays where this task starts*/
  task_first_part_f4[tasks_packed].x = pack_vars->count_parts;
  int *count_parts_self = &pack_vars->count_parts;
  /* This re-arranges the particle data from cell->hydro->parts into a
  long array of part structs*/
  runner_doself1_gpu_pack_neat_aos_f4(
      r, ci, parts_send, 0 /*timer. 0 no timing, 1 for timing*/,
      count_parts_self, tasks_packed, pack_vars->count_max_parts);
  /* Identify the row in the array where this task ends (row id of its
     last particle)*/
  task_first_part_f4[tasks_packed].y = pack_vars->count_parts;
  /* Identify first particle for each bundle of tasks */
  const int bundle_size = pack_vars->bundle_size;
  if (tasks_packed % bundle_size == 0) {
    int bid = tasks_packed / bundle_size;
    pack_vars->bundle_first_part[bid] = task_first_part_f4[tasks_packed].x;
    pack_vars->bundle_first_task_list[bid] = tasks_packed;
  }
  /* Tell the cell it has been packed */
  ci->pack_done++;
  /* Record that we have now done a packing (self) */
  t->done = 1;
  pack_vars->tasks_packed++;
  pack_vars->launch = 0;
  pack_vars->launch_leftovers = 0;

  /*Get a lock to the queue so we can safely decrement counter and check for launch leftover condition*/
  lock_lock(&s->queues[qid].lock);
  s->queues[qid].n_packs_self_left_d--;
  if (s->queues[qid].n_packs_self_left_d < 1) pack_vars->launch_leftovers = 1;
  lock_unlock(&s->queues[qid].lock);
  /*Have we packed enough tasks to offload to GPU?*/
  if (pack_vars->tasks_packed == pack_vars->target_n_tasks)
    pack_vars->launch = 1;

  /*Record the end of packing time*/
  clock_gettime(CLOCK_REALTIME, &t1);
  /* Release the lock on the cell */
  cell_unlocktree(ci);
  t->gpu_done = 1;
  /*Calculate time spent packing and return to runner_main*/
  return (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}

double runner_doself1_pack_f4_g(struct runner *r, struct scheduler *s,
                                struct pack_vars_self *pack_vars,
                                struct cell *ci, struct task *t,
                                struct part_aos_f4_g_send *parts_send,
                                int2 *task_first_part_f4) {

  /* Timers for how long this all takes.
   * t0 and t1 are from start to finish including GPU calcs
   * tp0 and tp1 only time packing and unpacking*/
  struct timespec t0, t1;  //
  clock_gettime(CLOCK_REALTIME, &t0);
  /* Find my queue for use later*/
  int qid = r->qid;
  /*Place pointers to the task and cells packed in an array for use later
   * when unpacking after the GPU offload*/
  int tasks_packed = pack_vars->tasks_packed;
  pack_vars->task_list[tasks_packed] = t;
  pack_vars->cell_list[tasks_packed] = ci;
  /* Identify row in particle arrays where this task starts*/
  task_first_part_f4[tasks_packed].x = pack_vars->count_parts;
  int *count_parts_self = &pack_vars->count_parts;
  /* This re-arranges the particle data from cell->hydro->parts into a
  long array of part structs*/
  runner_doself1_gpu_pack_neat_aos_f4_g(
      r, ci, parts_send, 0 /*timer. 0 no timing, 1 for timing*/,
      count_parts_self, tasks_packed, pack_vars->count_max_parts);
  /* identify the row in the array where this task ends (row id of its
     last particle)*/
  task_first_part_f4[tasks_packed].y = pack_vars->count_parts;
  /* Identify first particle for each bundle of tasks */
  const int bundle_size = pack_vars->bundle_size;
  if (tasks_packed % bundle_size == 0) {
    int bid = tasks_packed / bundle_size;
    pack_vars->bundle_first_part[bid] = task_first_part_f4[tasks_packed].x;
    pack_vars->bundle_first_task_list[bid] = tasks_packed;
  }
  /* Tell the cell it has been packed */
  ci->pack_done_g++;
  /* Record that we have now done a packing (self) */
  t->done = 1;
  pack_vars->tasks_packed++;
  pack_vars->launch = 0;
  pack_vars->launch_leftovers = 0;
  /*Get a lock to the queue so we can safely decrement counter and check for launch leftover condition*/
  lock_lock(&s->queues[qid].lock);
  s->queues[qid].n_packs_self_left_g--;
  if (s->queues[qid].n_packs_self_left_g < 1) pack_vars->launch_leftovers = 1;
  lock_unlock(&s->queues[qid].lock);

  if (pack_vars->tasks_packed == pack_vars->target_n_tasks)
    pack_vars->launch = 1;
  /*Add time to packing_time. Timer for end of GPU work after the if(launch ||
   * launch_leftovers statement)*/
  clock_gettime(CLOCK_REALTIME, &t1);
  /* Release the lock on the cell */
  cell_unlocktree(ci);
  /*Calculate time spent packing and return to runner_main*/
  return (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}

double runner_doself1_pack_f4_f(struct runner *r, struct scheduler *s,
                                struct pack_vars_self *pack_vars,
                                struct cell *ci, struct task *t,
                                struct part_aos_f4_f_send *parts_send,
                                int2 *task_first_part_f4) {

  /* Timers for how long this all takes.
   * t0 and t1 are from start to finish including GPU calcs
   * tp0 and tp1 only time packing and unpacking*/
  struct timespec t0, t1;  //
  clock_gettime(CLOCK_REALTIME, &t0);
  /* Find my queue for use later*/
  int qid = r->qid;
  /*Place pointers to the task and cells packed in an array for use later
   * when unpacking after the GPU offload*/
  int tasks_packed = pack_vars->tasks_packed;
  pack_vars->task_list[tasks_packed] = t;
  pack_vars->cell_list[tasks_packed] = ci;
  /* Identify row in particle arrays where this task starts*/
  task_first_part_f4[tasks_packed].x = pack_vars->count_parts;
  int *count_parts_self = &pack_vars->count_parts;
  /* This re-arranges the particle data from cell->hydro->parts into a
  long array of part structs*/
  runner_doself1_gpu_pack_neat_aos_f4_f(
      r, ci, parts_send, 0 /*timer. 0 no timing, 1 for timing*/,
      count_parts_self, tasks_packed, pack_vars->count_max_parts);
  /* Identify the row in the array where this task ends (row id of its
     last particle) */
  task_first_part_f4[tasks_packed].y = pack_vars->count_parts;
  /* Identify first particle for each bundle of tasks */
  const int bundle_size = pack_vars->bundle_size;
  if (tasks_packed % bundle_size == 0) {
    int bid = tasks_packed / bundle_size;
    pack_vars->bundle_first_part[bid] = task_first_part_f4[tasks_packed].x;
    pack_vars->bundle_first_task_list[bid] = tasks_packed;
  }
  /* Tell the cell it has been packed */
  ci->pack_done_f++;
  /* Record that we have now done a packing (self) */
  t->done = 1;
  pack_vars->tasks_packed++;
  pack_vars->launch = 0;
  pack_vars->launch_leftovers = 0;
  /*Get a lock to the queue so we can safely decrement counter and check for launch leftover condition*/
  lock_lock(&s->queues[qid].lock);
  s->queues[qid].n_packs_self_left_f--;
  if (s->queues[qid].n_packs_self_left_f < 1) pack_vars->launch_leftovers = 1;
  lock_unlock(&s->queues[qid].lock);
  /*Have we packed enough tasks to offload to GPU?*/
  if (pack_vars->tasks_packed == pack_vars->target_n_tasks)
    pack_vars->launch = 1;

  /*Record the end of packing time*/
  clock_gettime(CLOCK_REALTIME, &t1);
  /* Release the lock on the cell */
  cell_unlocktree(ci);
  /*Calculate time spent packing and return to runner_main*/
  return (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}

void runner_recurse_gpu(struct runner *r, struct scheduler *s,
                              struct pack_vars_pair *restrict pack_vars,
                              struct cell *ci, struct cell *cj, struct task *t,
                              struct part_aos_f4_send *parts_send,
                              struct engine *e,
                              int4 *fparti_fpartj_lparti_lpartj, int *n_leafs_found,
							  struct cell ** cells_left, struct cell ** cells_right, int depth, int n_expected_tasks) {

	/* Should we even bother? A. Nasar: For GPU code we need to be clever about this */
  if (!CELL_IS_ACTIVE(ci, e) && !CELL_IS_ACTIVE(cj, e)) return;
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_hydro_task(ci) &&
	  cell_can_recurse_in_pair_hydro_task(cj)) {
	struct cell_split_pair *csp = &cell_split_pairs[sid];
	for (int k = 0; k < csp->count; k++) {
	  const int pid = csp->pairs[k].pid;
	  const int pjd = csp->pairs[k].pjd;
	  /*Do we want to do anything before we recurse?*/

	  /*We probably want to record */
	  if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL){
		runner_recurse_gpu(r, s, pack_vars, ci->progeny[pid], cj->progeny[pjd], t, parts_send, e, fparti_fpartj_lparti_lpartj,
				n_leafs_found, cells_left, cells_right, depth + 1, n_expected_tasks);
//	        message("recursing to depth %i", depth + 1);
	  }
	}
  }
  else if (CELL_IS_ACTIVE(ci, e) || CELL_IS_ACTIVE(cj, e)) {
	/* if any cell empty: skip */
	if(ci->hydro.count == 0 || cj->hydro.count == 0) return;
	/*for all leafs to be sent add to cell list */
	cells_left[*n_leafs_found] = ci;
	cells_right[*n_leafs_found] = cj;
	/*Add leaf cells to list for each top_level task*/
	pack_vars->leaf_list[pack_vars->top_tasks_packed].ci[*n_leafs_found] = ci;
	pack_vars->leaf_list[pack_vars->top_tasks_packed].cj[*n_leafs_found] = cj;
	pack_vars->leaf_list[pack_vars->top_tasks_packed].n_leaves++;
	*n_leafs_found = *n_leafs_found + 1;
	if(*n_leafs_found >= n_expected_tasks)
		error("Created %i more than expected leaf cells. depth %i", *n_leafs_found, depth);
  }

};

double runner_dopair1_pack_f4(struct runner *r, struct scheduler *s,
                              struct pack_vars_pair *restrict pack_vars,
                              struct cell *ci, struct cell *cj, struct task *t,
                              struct part_aos_f4_send *parts_send,
                              struct engine *e,
                              int4 *fparti_fpartj_lparti_lpartj, int leaves_packed) {
  /* Timers for how long this all takes.
   * t0 and t1 are from start to finish including GPU calcs
   * tp0 and tp1 only time packing and unpacking*/
  struct timespec t0, t1;  //
  clock_gettime(CLOCK_REALTIME, &t0);
  int tasks_packed = pack_vars->tasks_packed;
  int qid = r->qid;

  double x_tmp = 0.0, y_tmp = 0.0, z_tmp = 0.0;
  struct cell *citmp, *cjtmp;
  citmp=ci;
  cjtmp=cj;
  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &citmp, &cjtmp, shift);
  if(citmp != ci) error("I'm flipped");
  /*Get the shifts in case of periodics*/
  space_getsid_GPU(e->s, &ci, &cj, &x_tmp, &y_tmp, &z_tmp);

  /*Get pointers to the list of tasks and cells packed*/
//  pack_vars->task_list[tasks_packed] = t;
  pack_vars->ci_list[tasks_packed] = ci;
  pack_vars->cj_list[tasks_packed] = cj;

  float3 shift_tmp = {x_tmp, y_tmp, z_tmp};

  const int count_ci = ci->hydro.count;
  const int count_cj = cj->hydro.count;

  /*Assign an id for this task*/
  const int tid = tasks_packed;

  /* Find first parts in task for ci and cj. Packed_tmp is index for cell i.
   * packed_tmp+1 is index for cell j */
  fparti_fpartj_lparti_lpartj[tasks_packed].x = pack_vars->count_parts;
  fparti_fpartj_lparti_lpartj[tasks_packed].y =
      pack_vars->count_parts + count_ci;

  int *count_parts = &pack_vars->count_parts;
  /* This re-arranges the particle data from cell->hydro->parts into a
  long array of part structs*/
  runner_do_ci_cj_gpu_pack_neat_aos_f4(
      r, ci, cj, parts_send, 0 /*timer. 0 no timing, 1 for timing*/,
      count_parts, tid, pack_vars->count_max_parts, count_ci, count_cj,
      shift_tmp);
  /* Find last parts in task for ci and cj*/
  fparti_fpartj_lparti_lpartj[tasks_packed].z =
      pack_vars->count_parts - count_cj;
  fparti_fpartj_lparti_lpartj[tasks_packed].w = pack_vars->count_parts;

  /* Tell the cells they have been packed */
  ci->pack_done++;
  cj->pack_done++;

  /* Identify first particle for each bundle of tasks */
  const int bundle_size = pack_vars->bundle_size;
  if (tasks_packed % bundle_size == 0) {
    int bid = tasks_packed / bundle_size;
    pack_vars->bundle_first_part[bid] =
        fparti_fpartj_lparti_lpartj[tasks_packed].x;
    pack_vars->bundle_first_task_list[bid] = tasks_packed;
  }
  /* Record that we have now done a packing (self) */
  t->done = 1;
  pack_vars->tasks_packed++;
  pack_vars->launch = 0;
  pack_vars->launch_leftovers = 0;
  pack_vars->leaf_list[pack_vars->top_tasks_packed - 1].n_packed++;

  //A. Nasar: Need to come back to this at some point!
  lock_lock(&s->queues[qid].lock);
  s->queues[qid].n_packs_pair_left_d--;
  if (s->queues[qid].n_packs_pair_left_d < 1) pack_vars->launch_leftovers = 1;
  lock_unlock(&s->queues[qid].lock);
  if (pack_vars->tasks_packed == pack_vars->target_n_tasks){
    pack_vars->launch = 1;
  }
  /*Add time to packing_time. Timer for end of GPU work after the if(launch ||
   * launch_leftovers statement)*/
  clock_gettime(CLOCK_REALTIME, &t1);
  return (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
};

double runner_dopair1_pack_f4_g(struct runner *r, struct scheduler *s,
                                struct pack_vars_pair *restrict pack_vars,
                                struct cell *ci, struct cell *cj,
                                struct task *t,
                                struct part_aos_f4_g_send *parts_send,
                                struct engine *e,
                                int4 *fparti_fpartj_lparti_lpartj) {

  /* Timers for how long this all takes.
   * t0 and t1 are from start to finish including GPU calcs
   * tp0 and tp1 only time packing and unpacking*/
  struct timespec t0, t1;  //
  clock_gettime(CLOCK_REALTIME, &t0);
  int tasks_packed = pack_vars->tasks_packed;

  int qid = r->qid;
  //  pthread_mutex_lock(&s->sleep_mutex);
  //  atomic_dec(&(s->p_g_left[qid]));
  //  pthread_cond_broadcast(&s->sleep_cond);
  //  pthread_mutex_unlock(&s->sleep_mutex);

  double x_tmp = 0.0, y_tmp = 0.0, z_tmp = 0.0;
  /*Get the shifts in case of periodics*/
  space_getsid_GPU(e->s, &ci, &cj, &x_tmp, &y_tmp, &z_tmp);

  /*Get pointers to the list of tasks and cells packed*/
  pack_vars->task_list[tasks_packed] = t;
  pack_vars->ci_list[tasks_packed] = ci;
  pack_vars->cj_list[tasks_packed] = cj;

  float3 shift_tmp = {x_tmp, y_tmp, z_tmp};

  const int count_ci = ci->hydro.count;
  const int count_cj = cj->hydro.count;

  /*Assign an id for this task*/
  const int tid = tasks_packed;

  /* Find first parts in task for ci and cj. Packed_tmp is index for cell i.
   * packed_tmp+1 is index for cell j */
  //    pack_vars->task_first_part[packed_tmp] = pack_vars->count_parts;
  //    pack_vars->task_first_part[packed_tmp + 1] = pack_vars->count_parts +
  //    count_ci;

  fparti_fpartj_lparti_lpartj[tasks_packed].x = pack_vars->count_parts;
  fparti_fpartj_lparti_lpartj[tasks_packed].y =
      pack_vars->count_parts + count_ci;

  int *count_parts = &pack_vars->count_parts;
  //    if(r->cpuid == 0)fprintf(stderr, "cpu %i before count %i\n", r->cpuid,
  //    pack_vars->count_parts);
  /* This re-arranges the particle data from cell->hydro->parts into a
  long array of part structs*/
  runner_do_ci_cj_gpu_pack_neat_aos_f4_g(
      r, ci, cj, parts_send, 0 /*timer. 0 no timing, 1 for timing*/,
      count_parts, tid, pack_vars->count_max_parts, count_ci, count_cj,
      shift_tmp);
  //	runner_doself1_gpu_pack_neat_aos(r, ci, parts_aos, 0/*timer. 0 no
  // timing, 1 for timing*/, 		  count_parts, tasks_packed,
  // pack_vars->count_max_parts); //This may cause an issue. Be sure to test
  // that
  // pack_vars->count_parts is actually increment here
  /* Find last parts in task for ci and cj. Packed_tmp is index for cell i.
   * packed_tmp+1 is index for cell j */

  //    if(r->cpuid == 0)fprintf(stderr, "cpu %i after count %i pack_vars_count
  //    %i\n", r->cpuid, *count_parts, 		pack_vars->count_parts);
  fparti_fpartj_lparti_lpartj[tasks_packed].z =
      pack_vars->count_parts - count_cj;
  fparti_fpartj_lparti_lpartj[tasks_packed].w = pack_vars->count_parts;
  //    pack_vars->task_last_part[packed_tmp] = pack_vars->count_parts -
  //    count_cj; pack_vars->task_last_part[packed_tmp + 1] =
  //    pack_vars->count_parts;

  /* Tell the cells they have been packed */
  ci->pack_done_g++;
  cj->pack_done_g++;

  /* Identify first particle for each bundle of tasks */
  const int bundle_size = pack_vars->bundle_size;
  if (tasks_packed % bundle_size == 0) {
    int bid = tasks_packed / bundle_size;
    pack_vars->bundle_first_part[bid] =
        fparti_fpartj_lparti_lpartj[tasks_packed].x;
    pack_vars->bundle_first_task_list[bid] = tasks_packed;
  }

  /* Record that we have now done a packing (self) */
  t->done = 1;
  /* Copies done. Release the lock ! */
  cell_unlocktree(ci);
  cell_unlocktree(cj);
  pack_vars->tasks_packed++;
  pack_vars->launch = 0;
  pack_vars->launch_leftovers = 0;
  /* Record that we have now done a packing (self) */
  //  int qid = r->qid;
  //  atomic_dec(&(s->queues[qid].n_packs_pair_left_g));

  lock_lock(&s->queues[qid].lock);

  s->queues[qid].n_packs_pair_left_g--;

  if (s->queues[qid].n_packs_pair_left_g < 1) pack_vars->launch_leftovers = 1;

  lock_unlock(&s->queues[qid].lock);

  //  if ((s->p_g_left[qid] < 1))
  //    pack_vars->launch_leftovers = 1;
  if (pack_vars->tasks_packed == pack_vars->target_n_tasks)
    pack_vars->launch = 1;
  /*Add time to packing_time. Timer for end of GPU work after the if(launch ||
   * launch_leftovers statement)*/
  clock_gettime(CLOCK_REALTIME, &t1);
  return (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}

double runner_dopair1_pack_f4_f(struct runner *r, struct scheduler *s,
                                struct pack_vars_pair *restrict pack_vars,
                                struct cell *ci, struct cell *cj,
                                struct task *t,
                                struct part_aos_f4_f_send *parts_send,
                                struct engine *e,
                                int4 *fparti_fpartj_lparti_lpartj) {

  /* Timers for how long this all takes.
   * t0 and t1 are from start to finish including GPU calcs
   * tp0 and tp1 only time packing and unpacking*/
  struct timespec t0, t1;  //
  clock_gettime(CLOCK_REALTIME, &t0);
  int tasks_packed = pack_vars->tasks_packed;

  /* Record that we have now done a packing (self) */
  int qid = r->qid;
  //  atomic_dec(&(s->queues[qid].n_packs_pair_left_f));
  //  pthread_mutex_lock(&s->sleep_mutex);
  atomic_dec(&(s->p_f_left[qid]));
  //  pthread_cond_broadcast(&s->sleep_cond);
  //  pthread_mutex_unlock(&s->sleep_mutex);

  double x_tmp = 0.0, y_tmp = 0.0, z_tmp = 0.0;
  /*Get the shifts in case of periodics*/
  space_getsid_GPU(e->s, &ci, &cj, &x_tmp, &y_tmp, &z_tmp);

  /*Get pointers to the list of tasks and cells packed*/
  pack_vars->task_list[tasks_packed] = t;
  pack_vars->ci_list[tasks_packed] = ci;
  pack_vars->cj_list[tasks_packed] = cj;

  float3 shift_tmp = {x_tmp, y_tmp, z_tmp};

  const int count_ci = ci->hydro.count;
  const int count_cj = cj->hydro.count;

  /*Assign an id for this task*/
  const int tid = tasks_packed;

  /* Find first parts in task for ci and cj. Packed_tmp is index for cell i.
   * packed_tmp+1 is index for cell j */
  //    pack_vars->task_first_part[packed_tmp] = pack_vars->count_parts;
  //    pack_vars->task_first_part[packed_tmp + 1] = pack_vars->count_parts +
  //    count_ci;

  fparti_fpartj_lparti_lpartj[tasks_packed].x = pack_vars->count_parts;
  fparti_fpartj_lparti_lpartj[tasks_packed].y =
      pack_vars->count_parts + count_ci;

  int *count_parts = &pack_vars->count_parts;
  //    if(r->cpuid == 0)fprintf(stderr, "cpu %i before count %i\n", r->cpuid,
  //    pack_vars->count_parts);
  /* This re-arranges the particle data from cell->hydro->parts into a
  long array of part structs*/
  runner_do_ci_cj_gpu_pack_neat_aos_f4_f(
      r, ci, cj, parts_send, 0 /*timer. 0 no timing, 1 for timing*/,
      count_parts, tid, pack_vars->count_max_parts, count_ci, count_cj,
      shift_tmp);
  //	runner_doself1_gpu_pack_neat_aos(r, ci, parts_aos, 0/*timer. 0 no
  // timing, 1 for timing*/, 		  count_parts, tasks_packed,
  // pack_vars->count_max_parts); //This may cause an issue. Be sure to test
  // that
  // pack_vars->count_parts is actually increment here
  /* Find last parts in task for ci and cj. Packed_tmp is index for cell i.
   * packed_tmp+1 is index for cell j */

  //    if(r->cpuid == 0)fprintf(stderr, "cpu %i after count %i pack_vars_count
  //    %i\n", r->cpuid, *count_parts, 		pack_vars->count_parts);
  fparti_fpartj_lparti_lpartj[tasks_packed].z =
      pack_vars->count_parts - count_cj;
  fparti_fpartj_lparti_lpartj[tasks_packed].w = pack_vars->count_parts;
  //    pack_vars->task_last_part[packed_tmp] = pack_vars->count_parts -
  //    count_cj; pack_vars->task_last_part[packed_tmp + 1] =
  //    pack_vars->count_parts;

  /* Tell the cells they have been packed */
  ci->pack_done_f++;
  cj->pack_done_f++;

  /* Identify first particle for each bundle of tasks */
  const int bundle_size = pack_vars->bundle_size;
  if (tasks_packed % bundle_size == 0) {
    int bid = tasks_packed / bundle_size;
    pack_vars->bundle_first_part[bid] =
        fparti_fpartj_lparti_lpartj[tasks_packed].x;
    pack_vars->bundle_first_task_list[bid] = tasks_packed;
  }

  /* Record that we have now done a packing (self) */
  t->done = 1;
  /* Copies done. Release the lock ! */
  cell_unlocktree(ci);
  cell_unlocktree(cj);
  pack_vars->tasks_packed++;
  pack_vars->launch = 0;
  pack_vars->launch_leftovers = 0;

  lock_lock(&s->queues[qid].lock);

  s->queues[qid].n_packs_pair_left_f--;

  if (s->queues[qid].n_packs_pair_left_f < 1) pack_vars->launch_leftovers = 1;

  lock_unlock(&s->queues[qid].lock);

  //  if ((s->p_f_left[qid] < 1))
  //    pack_vars->launch_leftovers = 1;
  if (pack_vars->tasks_packed == pack_vars->target_n_tasks)
    pack_vars->launch = 1;
  /*Add time to packing_time. Timer for end of GPU work after the if(launch ||
   * launch_leftovers statement)*/
  clock_gettime(CLOCK_REALTIME, &t1);
  return (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}

void runner_doself1_launch_f4(
    struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars,
    struct cell *ci, struct task *t, struct part_aos_f4_send *parts_send,
    struct part_aos_f4_recv *parts_recv, struct part_aos_f4_send *d_parts_send,
    struct part_aos_f4_recv *d_parts_recv, cudaStream_t *stream, float d_a,
    float d_H, struct engine *e, double *packing_time, double *gpu_time,
    double *unpack_time, int devId,
    int2 *task_first_part_f4, int2 *d_task_first_part_f4,
    cudaEvent_t *self_end) {

  struct timespec t0, t1, tp0, tp1;  //
  clock_gettime(CLOCK_REALTIME, &t0);

  /* Identify the number of GPU bundles to run in ideal case*/
  int nBundles_temp = pack_vars->nBundles;

  /*How many tasks have we packed?*/
  const int tasks_packed = pack_vars->tasks_packed;

  /*How many tasks should be in a bundle?*/
  const int bundle_size = pack_vars->bundle_size;

  /* Special case for incomplete bundles (when having leftover tasks not enough
   * to fill a bundle) */
  if (pack_vars->launch_leftovers) {
    nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
    if (tasks_packed == 0)
      error("zero tasks packed but somehow got into GPU loop");
    //	  pack_vars->bundle_first_part[nBundles_temp] =
    // pack_vars->task_first_part[tasks_packed - 1];
    pack_vars->bundle_first_part[nBundles_temp] =
        task_first_part_f4[tasks_packed - 1].x;
  }
  /* Identify the last particle for each bundle of tasks */
  for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
  }
  /* special treatment for the last bundle */
  if (nBundles_temp > 1)
    pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
  else
    pack_vars->bundle_last_part[0] = pack_vars->count_parts;
  //    clock_gettime(CLOCK_REALTIME, &t0hmemcpy);
  /*Copy arrays containing first and last part for each task to GPU*/
  //    cudaMemcpy(pack_vars->d_task_first_part, pack_vars->task_first_part,
  //               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
  //    cudaMemcpy(pack_vars->d_task_last_part, pack_vars->task_last_part,
  //               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
  //    cudaMemPrefetchAsync(d_task_first_part_self_dens_f4, tasks_packed *
  //    sizeof(int2), devId, NULL);
  /*Copy cell shifts to device*/
  //    cudaMemcpy(pack_vars->d_cellx, pack_vars->cellx,
  //               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
  //    cudaMemcpy(pack_vars->d_celly, pack_vars->celly,
  //               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
  //    cudaMemcpy(pack_vars->d_cellz, pack_vars->cellz,
  //               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
  //    clock_gettime(CLOCK_REALTIME, &t1hmemcpy);
  //    *hmemcpy_time += (t1hmemcpy.tv_sec - t0hmemcpy.tv_sec) +
  //			(t1hmemcpy.tv_nsec - t0hmemcpy.tv_nsec) / 1000000000.0;
  /* Launch the copies for each bundle and run the GPU kernel */
  /*We don't go into this loop if tasks_left_self == 1 as
   nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
  int max_parts;
  for (int bid = 0; bid < nBundles_temp; bid++) {

    max_parts = 0;
    int parts_in_bundle = 0;
    const int first_task = bid * bundle_size;
    int last_task = (bid + 1) * bundle_size;
    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {
      if (tid < tasks_packed) {
        /*Get an estimate for the max number of parts per cell in the bundle.
         *  Used for determining the number of GPU CUDA blocks*/
        int count = task_first_part_f4[tid].y - task_first_part_f4[tid].x;
        parts_in_bundle += count;
        max_parts = max(max_parts, count);
        last_task = tid;
      }
    }
    //	  const int n_tasks = last_task - first_task;

    const int first_part_tmp = pack_vars->bundle_first_part[bid];
    const int bundle_n_parts =
        pack_vars->bundle_last_part[bid] - first_part_tmp;
    //	  clock_gettime(CLOCK_REALTIME, &t0hmemcpy);
    //      cudaMemPrefetchAsync(&d_task_first_part_self_dens_f4[first_task],
    //      (last_task - first_task) * sizeof(int2),
    //    		  devId, stream[bid]);
    cudaMemcpyAsync(&d_task_first_part_f4[first_task],
                    &task_first_part_f4[first_task],
                    (last_task + 1 - first_task) * sizeof(int2),
                    cudaMemcpyHostToDevice, stream[bid]);
    //	  cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();
    //// 	  if (cu_error != cudaSuccess) { 		fprintf(
    /// stderr, 			"CUDA error in density
    // self host 2 device memcpy: %s cpuid id is: %i\n ",
    //			cudaGetErrorString(cu_error), r->cpuid);
    //		exit(0);
    //	  }
    //       clock_gettime(CLOCK_REALTIME, &t1hmemcpy);
    //       *hmemcpy_time += (t1hmemcpy.tv_sec - t0hmemcpy.tv_sec) +
    //   			(t1hmemcpy.tv_nsec - t0hmemcpy.tv_nsec) /
    //   1000000000.0;
    cudaMemcpyAsync(&d_parts_send[first_part_tmp], &parts_send[first_part_tmp],
                    bundle_n_parts * sizeof(struct part_aos_f4_send),
                    cudaMemcpyHostToDevice, stream[bid]);

    // #ifdef CUDA_DEBUG
    //	  cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();
    ////
    //										//
    // Get error code 	  if (cu_error != cudaSuccess) { 		fprintf(
    // stderr, 			"CUDA error in density self host 2 device
    // memcpy: %s cpuid id is: %i\n ",
    // cudaGetErrorString(cu_error), r->cpuid);
    //		exit(0);
    //	  }
    // #endif
    const int tasksperbundle = pack_vars->tasksperbundle;
    int tasks_left = tasksperbundle;
    if (bid == nBundles_temp - 1) {
      tasks_left = tasks_packed - (nBundles_temp - 1) * tasksperbundle;
    }
    // Will launch a 2d grid of GPU thread blocks (number of tasks is
    // the y dimension and max_parts is the x dimension
    int numBlocks_y = tasks_left;
    int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int bundle_first_task = pack_vars->bundle_first_task_list[bid];
    //	  const char *loop_type = "density";
    //	  struct first_part first_parts;
    //	  for(int i = 0; i < numBlocks_y; i++) first_parts.list[i] =
    // pack_vars->task_first_part[i]; 	  fprintf(stderr, "Launching kernel with
    // %i tasks leftovers %i\n", 			  tasks_packed,
    // pack_vars->launch_leftovers);
    // Launch the kernel
    launch_density_aos_f4(d_parts_send, d_parts_recv, d_a, d_H, stream[bid],
                          numBlocks_x, numBlocks_y, bundle_first_task,
                          d_task_first_part_f4);
    // #ifdef CUDA_DEBUG
    //	  cu_error = cudaPeekAtLastError(); // Get error code
    //	  if (cu_error != cudaSuccess) {
    //		fprintf(stderr,
    //				"CUDA error with self density kernel launch: %s
    // cpuid id is: %i\n ",
    // cudaGetErrorString(cu_error), r->cpuid); 		exit(0);
    //	  }
    // #endif
    cudaMemcpyAsync(&parts_recv[first_part_tmp], &d_parts_recv[first_part_tmp],
                    bundle_n_parts * sizeof(struct part_aos_f4_recv),
                    cudaMemcpyDeviceToHost, stream[bid]);
    cudaEventRecord(self_end[bid], stream[bid]);
    // #ifdef CUDA_DEBUG
    //	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
    //										//
    // Get error code 	  if (cu_error != cudaSuccess) {
    // fprintf(stderr, 				"CUDA error with self density
    // D2H memcpy: %s cpuid id is: %i\n ",
    // cudaGetErrorString(cu_error),
    // r->cpuid); 		error("Something's up with your cuda code");
    //	  }
    // #endif
  } /*End of looping over bundles to launch in streams*/
  /* Make sure all the kernels and copies back are finished */
  //	cudaDeviceSynchronize();

  /*Time end of GPU work*/
  clock_gettime(CLOCK_REALTIME, &t1);
  *gpu_time +=
      (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

  /* Now copy the data back from the CPU thread-local buffers to the cells */
  /* Pack length counter for use in unpacking */
  int pack_length_unpack = 0;
  ticks total_cpu_unpack_ticks = 0.;
  for (int bid = 0; bid < nBundles_temp; bid++) {

    clock_gettime(CLOCK_REALTIME, &t0);

    //		cudaStreamSynchronize(stream[bid]);
    cudaEventSynchronize(self_end[bid]);

    clock_gettime(CLOCK_REALTIME, &t1);
    *gpu_time +=
        (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

    /*Time unpacking*/
    //		clock_gettime(CLOCK_REALTIME, &tp0);

    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

      if (tid < tasks_packed) {
        struct cell *cii = pack_vars->cell_list[tid];
        struct task *tii = pack_vars->task_list[tid];

        //              struct cell *cii = ci_list_self_dens[tid];
        //              struct task *tii = task_list_self_dens[tid];

        clock_gettime(CLOCK_REALTIME, &tp0);

        //			  clock_gettime(CLOCK_REALTIME, &t0hmemcpy);
        while (cell_locktree(cii)) {
          ; /* spin until we acquire the lock */
        }
        //			  clock_gettime(CLOCK_REALTIME, &t1hmemcpy);
        //				*hmemcpy_time += (t1hmemcpy.tv_sec -
        // t0hmemcpy.tv_sec) + 				(t1hmemcpy.tv_nsec -
        // t0hmemcpy.tv_nsec) / 1000000000.0;
        const ticks tic = getticks();
        /* Do the copy */
        runner_doself1_gpu_unpack_neat_aos_f4(r, cii, parts_recv, 0,
                                              &pack_length_unpack, tid,
                                              pack_vars->count_max_parts, e);
        const ticks toc = getticks();

        total_cpu_unpack_ticks += toc - tic;
        /* Record things for debugging */
        cii->gpu_done++;
        /*Time end of unpacking*/
        clock_gettime(CLOCK_REALTIME, &tp1);
        *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
                        (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
        pthread_mutex_lock(&s->sleep_mutex);
        atomic_dec(&s->waiting);
        pthread_cond_broadcast(&s->sleep_cond);
        pthread_mutex_unlock(&s->sleep_mutex);
        /* Release the lock */
        cell_unlocktree(cii);

        /*schedule my dependencies (Only unpacks really)*/
        enqueue_dependencies(s, tii);
        /*Signal sleeping runners*/
        // MATTHIEU signal_sleeping_runners(s, tii);

        tii->gpu_done = 1;
      }
    }
    /*Time end of unpacking*/
    //		clock_gettime(CLOCK_REALTIME, &tp1);
    //		*hmemcpy_time += (tp1.tv_sec - tp0.tv_sec) +
    //		(tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
    //		*packing_time += (tp1.tv_sec - tp0.tv_sec) +
    //		(tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
  }
  /* Zero counters for the next pack operations */
  pack_vars->count_parts = 0;
  pack_vars->tasks_packed = 0;

  t->total_cpu_unpack_ticks += total_cpu_unpack_ticks;

} /*End of GPU work Self*/

void runner_doself1_launch_f4_g(
    struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars,
    struct cell *ci, struct task *t, struct part_aos_f4_g_send *parts_send,
    struct part_aos_f4_g_recv *parts_recv,
    struct part_aos_f4_g_send *d_parts_send,
    struct part_aos_f4_g_recv *d_parts_recv, cudaStream_t *stream, float d_a,
    float d_H, struct engine *e, double *packing_time, double *gpu_time,
    int2 *task_first_part_f4, int2 *d_task_first_part_f4, cudaEvent_t *self_end,
    double *unpack_time) {

  struct timespec t0, t1, tp0, tp1;
  clock_gettime(CLOCK_REALTIME, &t0);

  /* Identify the number of GPU bundles to run in ideal case*/
  int nBundles_temp = pack_vars->nBundles;

  /*How many tasks have we packed?*/
  const int tasks_packed = pack_vars->tasks_packed;

  /*How many tasks should be in a bundle?*/
  const int bundle_size = pack_vars->bundle_size;

  /* Special case for incomplete bundles (when having leftover tasks not enough
   * to fill a bundle) */
  if (pack_vars->launch_leftovers) {
    nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
    //	  if(tasks_packed == 0) error("zero tasks packed but somehow got into
    // GPU loop");
    pack_vars->bundle_first_part[nBundles_temp] =
        task_first_part_f4[tasks_packed - 1].x;
  }
  /* Identify the last particle for each bundle of tasks */
  for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
  }
  /* special treatment for the last bundle */
  if (nBundles_temp > 1)
    pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
  else
    pack_vars->bundle_last_part[0] = pack_vars->count_parts;

  /* Launch the copies for each bundle and run the GPU kernel */
  /*We don't go into this loop if tasks_left_self == 1 as
   nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
  int max_parts;
  for (int bid = 0; bid < nBundles_temp; bid++) {

    max_parts = 0;
    int parts_in_bundle = 0;
    const int first_task = bid * bundle_size;
    int last_task = (bid + 1) * bundle_size;
    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {
      if (tid < tasks_packed) {
        /*Get an estimate for the max number of parts per cell in the bundle.
         *  Used for determining the number of GPU CUDA blocks*/
        int count = task_first_part_f4[tid].y - task_first_part_f4[tid].x;
        parts_in_bundle += count;
        max_parts = max(max_parts, count);
        last_task = tid;
      }
    }

    const int first_part_tmp = pack_vars->bundle_first_part[bid];
    const int bundle_n_parts =
        pack_vars->bundle_last_part[bid] - first_part_tmp;

    cudaMemcpyAsync(&d_task_first_part_f4[first_task],
                    &task_first_part_f4[first_task],
                    (last_task + 1 - first_task) * sizeof(int2),
                    cudaMemcpyHostToDevice, stream[bid]);

    cudaMemcpyAsync(&d_parts_send[first_part_tmp], &parts_send[first_part_tmp],
                    bundle_n_parts * sizeof(struct part_aos_f4_g_send),
                    cudaMemcpyHostToDevice, stream[bid]);
    //	  fprintf(stderr, "bid %i first_part %i nparts %i\n", bid,
    // first_part_tmp, bundle_n_parts);

#ifdef CUDA_DEBUG
    cudaError_t cu_error =
        cudaPeekAtLastError();  // cudaGetLastError();        //
                                // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(stderr,
              "CUDA error in gradient self host 2 device memcpy: %s cpuid id "
              "is: %i\n ",
              cudaGetErrorString(cu_error), r->cpuid);
      exit(0);
    }
#endif
    const int tasksperbundle = pack_vars->tasksperbundle;
    int tasks_left = tasksperbundle;
    if (bid == nBundles_temp - 1) {
      tasks_left = tasks_packed - (nBundles_temp - 1) * tasksperbundle;
    }
    // Will launch a 2d grid of GPU thread blocks (number of tasks is
    // the y dimension and max_parts is the x dimension
    int numBlocks_y = tasks_left;
    int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int bundle_first_task = pack_vars->bundle_first_task_list[bid];
    //	  const char *loop_type = "density";
    // Launch the kernel
    launch_gradient_aos_f4(d_parts_send, d_parts_recv, d_a, d_H, stream[bid],
                           numBlocks_x, numBlocks_y, bundle_first_task,
                           d_task_first_part_f4);
#ifdef CUDA_DEBUG
    cu_error = cudaPeekAtLastError();  // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(
          stderr,
          "CUDA error with self gradient kernel launch: %s cpuid id is: %i\n ",
          cudaGetErrorString(cu_error), r->cpuid);
      exit(0);
    }
#endif
    cudaMemcpyAsync(&parts_recv[first_part_tmp], &d_parts_recv[first_part_tmp],
                    bundle_n_parts * sizeof(struct part_aos_f4_g_recv),
                    cudaMemcpyDeviceToHost, stream[bid]);
    cudaEventRecord(self_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
    cu_error = cudaPeekAtLastError();  // cudaGetLastError();        //
                                       // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(stderr,
              "CUDA error with self gradient D2H memcpy: %s cpuid id is: %i\n ",
              cudaGetErrorString(cu_error), r->cpuid);
      error("Something's up with your cuda code");
    }
#endif
  } /*End of looping over bundles to launch in streams*/
  //	exit(0);
  /* Make sure all the kernels and copies back are finished */
  //	cudaDeviceSynchronize();

  /*Time end of GPU work*/
  clock_gettime(CLOCK_REALTIME, &t1);
  *gpu_time +=
      (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

  /* Now copy the data back from the CPU thread-local buffers to the cells */
  /* Pack length counter for use in unpacking */
  int pack_length_unpack = 0;
  ticks total_cpu_unpack_ticks = 0.;
  for (int bid = 0; bid < nBundles_temp; bid++) {

    clock_gettime(CLOCK_REALTIME, &t0);

    //		cudaStreamSynchronize(stream[bid]);
    cudaEventSynchronize(self_end[bid]);

    clock_gettime(CLOCK_REALTIME, &t1);
    *gpu_time +=
        (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

    /*Time unpacking*/
    //		clock_gettime(CLOCK_REALTIME, &tp0);

    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

      if (tid < tasks_packed) {

        struct cell *cii = pack_vars->cell_list[tid];
        struct task *tii = pack_vars->task_list[tid];

        //              struct cell *cii = ci_list_self_dens[tid];
        //              struct task *tii = task_list_self_dens[tid];

        while (cell_locktree(cii)) {
          ; /* spin until we acquire the lock */
        }
        /*Time unpacking*/
        clock_gettime(CLOCK_REALTIME, &tp0);
        const ticks tic = getticks();

        /* Do the copy */
        runner_doself1_gpu_unpack_neat_aos_f4_g(r, cii, parts_recv, 0,
                                                &pack_length_unpack, tid,
                                                pack_vars->count_max_parts, e);
        const ticks toc = getticks();

        total_cpu_unpack_ticks += toc - tic;
        /*Time end of unpacking*/
        clock_gettime(CLOCK_REALTIME, &tp1);
        *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
                        (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;

        /* Record things for debugging */
        cii->gpu_done_g++;
        pthread_mutex_lock(&s->sleep_mutex);
        atomic_dec(&s->waiting);
        pthread_cond_broadcast(&s->sleep_cond);
        pthread_mutex_unlock(&s->sleep_mutex);
        /* Release the lock */
        cell_unlocktree(cii);

        /*schedule my dependencies (Only unpacks really)*/
        enqueue_dependencies(s, tii);
        /*Signal sleeping runners*/
        // MATTHIEU signal_sleeping_runners(s, tii);

        tii->gpu_done = 1;
      }
    }
    /*Time end of unpacking*/
    //		clock_gettime(CLOCK_REALTIME, &tp1);
    //		*unpack_time += (tp1.tv_sec - tp0.tv_sec) +
    //		(tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
    //		*packing_time += (tp1.tv_sec - tp0.tv_sec) +
    //		(tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
  }
  /* Zero counters for the next pack operations */
  pack_vars->count_parts = 0;
  pack_vars->tasks_packed = 0;

  t->total_cpu_unpack_ticks += total_cpu_unpack_ticks;

} /*End of GPU work Self Gradient*/

void runner_doself1_launch_f4_f(
    struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars,
    struct cell *ci, struct task *t, struct part_aos_f4_f_send *parts_send,
    struct part_aos_f4_f_recv *parts_recv,
    struct part_aos_f4_f_send *d_parts_send,
    struct part_aos_f4_f_recv *d_parts_recv, cudaStream_t *stream, float d_a,
    float d_H, struct engine *e, double *packing_time, double *gpu_time,
    int2 *task_first_part_f4_f, int2 *d_task_first_part_f4_f,
    cudaEvent_t *self_end, double *unpack_time) {

  struct timespec t0, t1, tp0, tp1;  //
  clock_gettime(CLOCK_REALTIME, &t0);

  /* Identify the number of GPU bundles to run in ideal case*/
  int nBundles_temp = pack_vars->nBundles;

  /*How many tasks have we packed?*/
  const int tasks_packed = pack_vars->tasks_packed;

  /*How many tasks should be in a bundle?*/
  const int bundle_size = pack_vars->bundle_size;

  /* Special case for incomplete bundles (when having leftover tasks not enough
   * to fill a bundle) */
  if (pack_vars->launch_leftovers) {
    nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
    if (tasks_packed == 0)
      error("zero tasks packed but somehow got into GPU loop");
    pack_vars->bundle_first_part[nBundles_temp] =
        task_first_part_f4_f[tasks_packed - 1].x;
  }
  /* Identify the last particle for each bundle of tasks */
  for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
  }
  /* special treatment for the last bundle */
  if (nBundles_temp > 1)
    pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
  else
    pack_vars->bundle_last_part[0] = pack_vars->count_parts;
  /*Copy arrays containing first and last part for each task to GPU*/
  //    cudaMemcpy(pack_vars->d_task_first_part, pack_vars->task_first_part,
  //               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
  //    cudaMemcpy(pack_vars->d_task_last_part, pack_vars->task_last_part,
  //               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);

  /*Copy cell shifts to device*/
  //    cudaMemcpy(pack_vars->d_cellx, pack_vars->cellx,
  //               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
  //    cudaMemcpy(pack_vars->d_celly, pack_vars->celly,
  //               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
  //    cudaMemcpy(pack_vars->d_cellz, pack_vars->cellz,
  //               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);

  /* Launch the copies for each bundle and run the GPU kernel */
  /*We don't go into this loop if tasks_left_self == 1 as
   nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
  int max_parts = 0;
  for (int bid = 0; bid < nBundles_temp; bid++) {

    max_parts = 0;
    int parts_in_bundle = 0;
    const int first_task = bid * bundle_size;
    int last_task = (bid + 1) * bundle_size;
    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {
      if (tid < tasks_packed) {
        /*Get an estimate for the max number of parts per cell in the bundle.
         *  Used for determining the number of GPU CUDA blocks*/
        int count = task_first_part_f4_f[tid].y - task_first_part_f4_f[tid].x;
        parts_in_bundle += count;
        max_parts = max(max_parts, count);
        last_task = tid;
      }
    }

    const int first_part_tmp = pack_vars->bundle_first_part[bid];
    const int bundle_n_parts =
        pack_vars->bundle_last_part[bid] - first_part_tmp;
    cudaMemcpyAsync(&d_task_first_part_f4_f[first_task],
                    &task_first_part_f4_f[first_task],
                    (last_task + 1 - first_task) * sizeof(int2),
                    cudaMemcpyHostToDevice, stream[bid]);

    cudaMemcpyAsync(&d_parts_send[first_part_tmp], &parts_send[first_part_tmp],
                    bundle_n_parts * sizeof(struct part_aos_f4_f_send),
                    cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
    cudaError_t cu_error =
        cudaPeekAtLastError();  // cudaGetLastError();        //
                                // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(stderr,
              "CUDA error in density self host 2 device memcpy: %s cpuid id "
              "is: %i\n ",
              cudaGetErrorString(cu_error), r->cpuid);
      exit(0);
    }
#endif
    const int tasksperbundle = pack_vars->tasksperbundle;
    int tasks_left = tasksperbundle;
    if (bid == nBundles_temp - 1) {
      tasks_left = tasks_packed - (nBundles_temp - 1) * tasksperbundle;
    }
    // Will launch a 2d grid of GPU thread blocks (number of tasks is
    // the y dimension and max_parts is the x dimension
    int numBlocks_y = tasks_left;
    int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int bundle_first_task = pack_vars->bundle_first_task_list[bid];
    // Launch the kernel
    launch_force_aos_f4(d_parts_send, d_parts_recv, d_a, d_H, stream[bid],
                        numBlocks_x, numBlocks_y, bundle_first_task,
                        d_task_first_part_f4_f);
#ifdef CUDA_DEBUG
    cu_error = cudaPeekAtLastError();  // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(stderr,
              "CUDA error with self force kernel launch: %s cpuid id is: %i\n ",
              cudaGetErrorString(cu_error), r->cpuid);
      exit(0);
    }
#endif
    cudaMemcpyAsync(&parts_recv[first_part_tmp], &d_parts_recv[first_part_tmp],
                    bundle_n_parts * sizeof(struct part_aos_f4_f_recv),
                    cudaMemcpyDeviceToHost, stream[bid]);
    cudaEventRecord(self_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
    cu_error = cudaPeekAtLastError();  // cudaGetLastError();        //
                                       // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(stderr,
              "CUDA error with self firce D2H memcpy: %s cpuid id is: %i\n ",
              cudaGetErrorString(cu_error), r->cpuid);
      error("Something's up with your cuda code");
    }
#endif
  } /*End of looping over bundles to launch in streams*/

  /* Make sure all the kernels and copies back are finished */
  //	cudaDeviceSynchronize();

  /*Time end of GPU work*/
  clock_gettime(CLOCK_REALTIME, &t1);
  *gpu_time +=
      (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

  /* Now copy the data back from the CPU thread-local buffers to the cells */
  /* Pack length counter for use in unpacking */
  int pack_length_unpack = 0;
  ticks total_cpu_unpack_ticks = 0.;
  for (int bid = 0; bid < nBundles_temp; bid++) {

    clock_gettime(CLOCK_REALTIME, &t0);

    //		cudaStreamSynchronize(stream[bid]);
    cudaEventSynchronize(self_end[bid]);

    clock_gettime(CLOCK_REALTIME, &t1);
    *gpu_time +=
        (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

    /*Time unpacking*/
    //		clock_gettime(CLOCK_REALTIME, &tp0);

    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

      if (tid < tasks_packed) {
        struct cell *cii = pack_vars->cell_list[tid];
        struct task *tii = pack_vars->task_list[tid];

        //              struct cell *cii = ci_list_self_dens[tid];
        //              struct task *tii = task_list_self_dens[tid];

        while (cell_locktree(cii)) {
          ; /* spin until we acquire the lock */
        }
        clock_gettime(CLOCK_REALTIME, &tp0);
        const ticks tic = getticks();

        /* Do the copy */
        runner_doself1_gpu_unpack_neat_aos_f4_f(r, cii, parts_recv, 0,
                                                &pack_length_unpack, tid,
                                                pack_vars->count_max_parts, e);
        const ticks toc = getticks();

        total_cpu_unpack_ticks += toc - tic;
        /* Record things for debugging */
        cii->gpu_done_f++;
        clock_gettime(CLOCK_REALTIME, &tp1);
        *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
                        (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
        pthread_mutex_lock(&s->sleep_mutex);
        atomic_dec(&s->waiting);
        pthread_cond_broadcast(&s->sleep_cond);
        pthread_mutex_unlock(&s->sleep_mutex);
        /* Release the lock */
        cell_unlocktree(cii);

        /*schedule my dependencies (Only unpacks really)*/
        enqueue_dependencies(s, tii);
        /*Signal sleeping runners*/
        // MATTHIEU signal_sleeping_runners(s, tii);

        tii->gpu_done = 1;
      }
    }
    /*Time end of unpacking*/
    //		clock_gettime(CLOCK_REALTIME, &tp1);
    //		*unpack_time += (tp1.tv_sec - tp0.tv_sec) +
    //		(tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
  }

  /* Zero counters for the next pack operations */
  pack_vars->count_parts = 0;
  pack_vars->tasks_packed = 0;

  t->total_cpu_unpack_ticks += total_cpu_unpack_ticks;
} /*End of GPU work Self Gradient*/

void runner_dopair1_launch_f4_one_memcpy(
    struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars,
    struct task *t, struct part_aos_f4_send *parts_send,
    struct part_aos_f4_recv *parts_recv, struct part_aos_f4_send *d_parts_send,
    struct part_aos_f4_recv *d_parts_recv, cudaStream_t *stream, float d_a,
    float d_H, struct engine *e, double *packing_time, double *gpu_time,
    double *unpack_time, int4 *fparti_fpartj_lparti_lpartj_dens,
    cudaEvent_t *pair_end) {

  struct timespec t0, t1, tp0, tp1;  //
  clock_gettime(CLOCK_REALTIME, &t0);

  /* Identify the number of GPU bundles to run in ideal case*/
  int nBundles_temp = pack_vars->nBundles;
  /*How many tasks have we packed?*/
  const int tasks_packed = pack_vars->tasks_packed;

  /*How many tasks should be in a bundle?*/
  const int bundle_size = pack_vars->bundle_size;

  /* Special case for incomplete bundles (when having leftover tasks not enough
   * to fill a bundle) */
  if (pack_vars->launch_leftovers) {
    nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
    if (tasks_packed == 0)
      error("zero pair tasks packed but somehow got into GPU loop");
    //	  pack_vars->bundle_first_part[nBundles_temp] =
    // pack_vars->task_first_part[packed_tmp - 2];
    pack_vars->bundle_first_part[nBundles_temp] =
        fparti_fpartj_lparti_lpartj_dens[tasks_packed - 1].x;
  }
  /* Identify the last particle for each bundle of tasks */
  for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
  }
  /* special treatment for the last bundle */
  if (nBundles_temp > 1)
    pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
  else
    pack_vars->bundle_last_part[0] = pack_vars->count_parts;

  /* Launch the copies for each bundle and run the GPU kernel */
  /*We don't go into this loop if tasks_left_self == 1 as
   nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
  for (int bid = 0; bid < nBundles_temp; bid++) {

    int max_parts_i = 0;
    int max_parts_j = 0;
    int parts_in_bundle_ci = 0;
    int parts_in_bundle_cj = 0;
    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {
      if (tid < tasks_packed) {
        /*Get an estimate for the max number of parts per cell in each bundle.
         *  Used for determining the number of GPU CUDA blocks*/
        int count_i = fparti_fpartj_lparti_lpartj_dens[tid].z -
                      fparti_fpartj_lparti_lpartj_dens[tid].x;
        parts_in_bundle_ci += count_i;
        max_parts_i = max(max_parts_i, count_i);
        int count_j = fparti_fpartj_lparti_lpartj_dens[tid].w -
                      fparti_fpartj_lparti_lpartj_dens[tid].y;
        parts_in_bundle_cj += count_j;
        max_parts_j = max(max_parts_j, count_j);
        //        if(count_i > 100 || count_j > 100)
        //        	error("Sending data for excessive n parts %i %i",
        //        count_i, count_j);
      }
    }
    const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
    const int bundle_n_parts =
        pack_vars->bundle_last_part[bid] - first_part_tmp_i;

    cudaMemcpyAsync(&d_parts_send[first_part_tmp_i],
                    &parts_send[first_part_tmp_i],
                    bundle_n_parts * sizeof(struct part_aos_f4_send),
                    cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
    cudaError_t cu_error =
        cudaPeekAtLastError();  // cudaGetLastError();        //
                                // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(stderr,
              "CUDA error with pair density H2D async  memcpy ci: %s cpuid id "
              "is: %i\n ",
              cudaGetErrorString(cu_error), r->cpuid);
      error("Something's up with your cuda code first_part %i bundle size %i",
            first_part_tmp_i, bundle_n_parts);
    }
#endif
    /* LAUNCH THE GPU KERNELS for ci & cj */
    // Setup 2d grid of GPU thread blocks for ci (number of tasks is
    // the y dimension and max_parts is the x dimension
    int numBlocks_y = 0;  // tasks_left;
    int numBlocks_x = (bundle_n_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int bundle_part_0 = pack_vars->bundle_first_part[bid];
    /* Launch the kernel for ci using data for ci and cj */
    runner_dopair_branch_density_gpu_aos_f4(
        d_parts_send, d_parts_recv, d_a, d_H, stream[bid], numBlocks_x,
        numBlocks_y, bundle_part_0, bundle_n_parts);

#ifdef CUDA_DEBUG
    cu_error = cudaPeekAtLastError();  // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(
          stderr,
          "CUDA error with pair density kernel launch: %s cpuid id is: %i\n "
          "nbx %i nby %i max_parts_i %i max_parts_j %i\n",
          cudaGetErrorString(cu_error), r->cpuid, numBlocks_x, numBlocks_y,
          max_parts_i, max_parts_j);
      error("Something's up with kernel launch.");
    }
#endif

    // Copy results back to CPU BUFFERS
    cudaMemcpyAsync(&parts_recv[first_part_tmp_i],
                    &d_parts_recv[first_part_tmp_i],
                    bundle_n_parts * sizeof(struct part_aos_f4_recv),
                    cudaMemcpyDeviceToHost, stream[bid]);
    cudaEventRecord(pair_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
    cu_error = cudaPeekAtLastError();  // cudaGetLastError();        //
                                       // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(stderr,
              "CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
              cudaGetErrorString(cu_error), r->cpuid);
      error("Something's up with your cuda code");
    }
#endif
  } /*End of looping over bundles to launch in streams*/

  /* Make sure all the kernels and copies back are finished */
  //	cudaDeviceSynchronize();

  /*Time end of GPU work*/
  clock_gettime(CLOCK_REALTIME, &t1);
  *gpu_time +=
      (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

  /* Now copy the data back from the CPU thread-local buffers to the cells */
  /* Pack length counter for use in unpacking */

  int pack_length_unpack = 0;
  ticks total_cpu_unpack_ticks = 0;

  for (int bid = 0; bid < nBundles_temp; bid++) {
    /*Time unpacking*/
    clock_gettime(CLOCK_REALTIME, &t0);

    //		cudaStreamSynchronize(stream[bid]);
    cudaEventSynchronize(pair_end[bid]);

    clock_gettime(CLOCK_REALTIME, &t1);
    *gpu_time +=
        (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

    ////////////

    /*Time unpacking*/
    //		clock_gettime(CLOCK_REALTIME, &tp0);

//    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {
//
//      if (tid < tasks_packed) {
//        clock_gettime(CLOCK_REALTIME, &tp0);
//        /*grab cell and task pointers*/
//        struct cell *cii = pack_vars->ci_list[tid];
//        struct cell *cjj = pack_vars->cj_list[tid];
//        struct task *tii = pack_vars->task_list[tid];
//
////        if(!pack_vars->task_locked){
////          /*Let's lock ci*/
////          while (cell_locktree(cii)) {
////            ; /* spin until we acquire the lock */
////          }
////          /*Let's lock cj*/
////          while (cell_locktree(cjj)) {
////            ; /* spin until we acquire the lock */
////          }
////          pack_vars->task_locked = 1;
////        }
//
//        const ticks tic = getticks();
//
//        /* Do the copy */
//        runner_do_ci_cj_gpu_unpack_neat_aos_f4(
//            r, cii, cjj, parts_recv, 0, &pack_length_unpack, tid,
//            2 * pack_vars->count_max_parts, e);
//
//        const ticks toc = getticks();
//
//        total_cpu_unpack_ticks += toc - tic;
//
//        /* Record things for debugging */
//        cii->gpu_done_pair++;
//        cjj->gpu_done_pair++;
//
////        if(pack_vars->task_locked){
////          /* Release the locks */
////          cell_unlocktree(cii);
////          /* Release the locks */
////          cell_unlocktree(cjj);
//          pack_vars->task_locked = 0;
////        }
//
//        /*Time end of unpacking*/
//        clock_gettime(CLOCK_REALTIME, &tp1);
//        *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
//                        (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
//        /*Signal sleeping runners*/
//        // MATTHIEU signal_sleeping_runners(s, tii);
//
//        tii->gpu_done = 1;
//      }
//    }
  }

  /* Zero counters for the next pack operations */
  pack_vars->count_parts = 0;
  pack_vars->tasks_packed = 0;

  //	/*Time end of unpacking*/
  //	clock_gettime(CLOCK_REALTIME, &t1);
  //	*packing_time += (t1.tv_sec - t0.tv_sec) +
  //	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

  /* Write the timers back to the task */
  t->total_cpu_unpack_ticks += total_cpu_unpack_ticks;

} /*End of GPU work*/

void runner_dopair1_unpack_f4(
    struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars,
    struct task *t, struct part_aos_f4_send *parts_send,
    struct part_aos_f4_recv *parts_recv, struct part_aos_f4_send *d_parts_send,
    struct part_aos_f4_recv *d_parts_recv, cudaStream_t *stream, float d_a,
    float d_H, struct engine *e, double *packing_time, double *gpu_time,
    double *unpack_time, int4 *fparti_fpartj_lparti_lpartj_dens,
    cudaEvent_t *pair_end, int cstart, int n_leaves_found){

  int topid;
  int pack_length_unpack = 0;
  ticks total_cpu_unpack_ticks = 0;
  for (topid = 0; topid < pack_vars->top_tasks_packed; topid++) {
	//lock top level cell here
	struct cell * cii = pack_vars->top_task_list[topid]->ci;
	struct cell * cjj = pack_vars->top_task_list[topid]->cj;
	const ticks tic = getticks();
	/* Do the copy */

	int n_leaves_in_task = pack_vars->leaf_list[topid].n_packed;
	for(int tid = 0; tid < n_leaves_in_task; tid++){
	  //Get pointers to the leaf cells. SEEMS I'm NOT GETTING A CORRECT POINTER
	  struct cell * cii_l = pack_vars->leaf_list[topid].ci[tid];
	  struct cell * cjj_l = pack_vars->leaf_list[topid].cj[tid];
	  runner_do_ci_cj_gpu_unpack_neat_aos_f4(
			r, cii_l, cjj_l, parts_recv, 0, &pack_length_unpack, tid,
			2 * pack_vars->count_max_parts, e);
	}

	const ticks toc = getticks();
	total_cpu_unpack_ticks += toc - tic;
	/*For some reason the code fails if we get a leaf pair task
	 *this if statement stops the code from trying to unlock same cells twice*/
	if(topid == pack_vars->top_tasks_packed -1 && cstart != n_leaves_found)
		continue;
    enqueue_dependencies(s, pack_vars->top_task_list[topid]);
    pthread_mutex_lock(&s->sleep_mutex);
    atomic_dec(&s->waiting);
    pthread_cond_broadcast(&s->sleep_cond);
    pthread_mutex_unlock(&s->sleep_mutex);
    pack_vars->task_locked = 0;
  }
}
void runner_dopair1_launch_f4_g_one_memcpy(
    struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars,
    struct task *t, struct part_aos_f4_g_send *parts_send,
    struct part_aos_f4_g_recv *parts_recv,
    struct part_aos_f4_g_send *d_parts_send,
    struct part_aos_f4_g_recv *d_parts_recv, cudaStream_t *stream, float d_a,
    float d_H, struct engine *e, double *packing_time, double *gpu_time,
    double *unpack_time, int4 *fparti_fpartj_lparti_lpartj,
    cudaEvent_t *pair_end) {

  struct timespec t0, t1, tp0, tp1;  //
  clock_gettime(CLOCK_REALTIME, &t0);

  /* Identify the number of GPU bundles to run in ideal case*/
  int nBundles_temp = pack_vars->nBundles;
  /*How many tasks have we packed?*/
  const int tasks_packed = pack_vars->tasks_packed;

  /*How many tasks should be in a bundle?*/
  const int bundle_size = pack_vars->bundle_size;

  /*tasks-packed needs decrementing before calculating packed_tmp as it was
   * incremented in runner_dopair1_pack*/
  //	const int packed_tmp = 2 * (tasks_packed - 1);

  /* Special case for incomplete bundles (when having leftover tasks not enough
   * to fill a bundle) */
  if (pack_vars->launch_leftovers) {
    nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
    if (tasks_packed == 0)
      error("zero pair tasks packed but somehow got into GPU loop");
    //	  pack_vars->bundle_first_part[nBundles_temp] =
    // pack_vars->task_first_part[packed_tmp - 2];
    pack_vars->bundle_first_part[nBundles_temp] =
        fparti_fpartj_lparti_lpartj[tasks_packed - 1].x;
  }
  /* Identify the last particle for each bundle of tasks */
  for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
  }
  /* special treatment for the last bundle */
  if (nBundles_temp > 1)
    pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
  else
    pack_vars->bundle_last_part[0] = pack_vars->count_parts;

  /* Launch the copies for each bundle and run the GPU kernel */
  /*We don't go into this loop if tasks_left_self == 1 as
   nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
  //	int max_parts = 0;
  for (int bid = 0; bid < nBundles_temp; bid++) {

    int max_parts_i = 0;
    int max_parts_j = 0;
    int parts_in_bundle_ci = 0;
    int parts_in_bundle_cj = 0;
    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {
      if (tid < tasks_packed) {
        /*Get an estimate for the max number of parts per cell in each bundle.
         *  Used for determining the number of GPU CUDA blocks*/
        int count_i = fparti_fpartj_lparti_lpartj[tid].z -
                      fparti_fpartj_lparti_lpartj[tid].x;
        parts_in_bundle_ci += count_i;
        max_parts_i = max(max_parts_i, count_i);
        int count_j = fparti_fpartj_lparti_lpartj[tid].w -
                      fparti_fpartj_lparti_lpartj[tid].y;
        parts_in_bundle_cj += count_j;
        max_parts_j = max(max_parts_j, count_j);
      }
    }
    const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
    const int bundle_n_parts =
        pack_vars->bundle_last_part[bid] - first_part_tmp_i;

    cudaMemcpyAsync(&d_parts_send[first_part_tmp_i],
                    &parts_send[first_part_tmp_i],
                    bundle_n_parts * sizeof(struct part_aos_f4_g_send),
                    cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
    cudaError_t cu_error =
        cudaPeekAtLastError();  // cudaGetLastError();        //
                                // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(stderr,
              "CUDA error with pair density H2D async  memcpy ci: %s cpuid id "
              "is: %i\n ",
              cudaGetErrorString(cu_error), r->cpuid);
      error("Something's up with your cuda code");
    }
#endif

    //	  const int tasksperbundle = pack_vars->tasksperbundle;
    /* LAUNCH THE GPU KERNELS for ci & cj */
    // Setup 2d grid of GPU thread blocks for ci (number of tasks is
    // the y dimension and max_parts is the x dimension
    int numBlocks_y = 0;  // tasks_left;
    int numBlocks_x = (bundle_n_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int bundle_part_0 = pack_vars->bundle_first_part[bid];
    //              fprintf(stderr, "bundle_part_0 %i bundle_first_task %i\n",
    //              bundle_part_0, bundle_first_task);

    /* Launch the kernel for ci using data for ci and cj */
    runner_dopair_branch_gradient_gpu_aos_f4(
        d_parts_send, d_parts_recv, d_a, d_H, stream[bid], numBlocks_x,
        numBlocks_y, bundle_part_0, bundle_n_parts);

#ifdef CUDA_DEBUG
    cu_error = cudaPeekAtLastError();  // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(
          stderr,
          "CUDA error with pair density kernel launch: %s cpuid id is: %i\n "
          "nbx %i nby %i max_parts_i %i max_parts_j %i\n",
          cudaGetErrorString(cu_error), r->cpuid, numBlocks_x, numBlocks_y,
          max_parts_i, max_parts_j);
      exit(0);
    }
#endif

    // Copy results back to CPU BUFFERS
    cudaMemcpyAsync(&parts_recv[first_part_tmp_i],
                    &d_parts_recv[first_part_tmp_i],
                    bundle_n_parts * sizeof(struct part_aos_f4_g_recv),
                    cudaMemcpyDeviceToHost, stream[bid]);
    cudaEventRecord(pair_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
    cu_error = cudaPeekAtLastError();  // cudaGetLastError();        //
                                       // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(stderr,
              "CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
              cudaGetErrorString(cu_error), r->cpuid);
      error("Something's up with your cuda code");
    }
#endif
  } /*End of looping over bundles to launch in streams*/

  /* Make sure all the kernels and copies back are finished */
  //	cudaDeviceSynchronize();

  /*Time end of GPU work*/
  clock_gettime(CLOCK_REALTIME, &t1);
  *gpu_time +=
      (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  /* Now copy the data back from the CPU thread-local buffers to the cells */
  /* Pack length counter for use in unpacking */
  int pack_length_unpack = 0;

  ticks total_cpu_unpack_ticks = 0.;

  for (int bid = 0; bid < nBundles_temp; bid++) {
    /*Time unpacking*/
    clock_gettime(CLOCK_REALTIME, &t0);

    //		cudaStreamSynchronize(stream[bid]);
    cudaEventSynchronize(pair_end[bid]);

    clock_gettime(CLOCK_REALTIME, &t1);
    *gpu_time +=
        (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

    /*Time unpacking*/
    //		clock_gettime(CLOCK_REALTIME, &tp0);
    //		int bundle_first_task = pack_vars->bundle_first_task_list[bid];

    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

      if (tid < tasks_packed) {
        clock_gettime(CLOCK_REALTIME, &tp0);
        /*grab cell and task pointers*/
        struct cell *cii = pack_vars->ci_list[tid];
        struct cell *cjj = pack_vars->cj_list[tid];
        struct task *tii = pack_vars->task_list[tid];
        /*Let's lock ci*/
        while (cell_locktree(cii)) {
          ; /* spin until we acquire the lock */
        }
        /*Let's lock cj*/
        while (cell_locktree(cjj)) {
          ; /* spin until we acquire the lock */
        }

        const ticks tic = getticks();

        /* Do the copy */
        runner_do_ci_cj_gpu_unpack_neat_aos_f4_g(
            r, cii, cjj, parts_recv, 0, &pack_length_unpack, tid,
            2 * pack_vars->count_max_parts, e);

        const ticks toc = getticks();

        total_cpu_unpack_ticks += toc - tic;

        /* Record things for debugging */
        cii->gpu_done_pair_g++;
        cjj->gpu_done_pair_g++;
        pthread_mutex_lock(&s->sleep_mutex);
        atomic_dec(&s->waiting);
        pthread_cond_broadcast(&s->sleep_cond);
        pthread_mutex_unlock(&s->sleep_mutex);
        /* Release the locks */
        cell_unlocktree(cii);
        /* Release the locks */
        cell_unlocktree(cjj);

        /*Time end of unpacking*/
        clock_gettime(CLOCK_REALTIME, &tp1);
        *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
                        (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;

        /*schedule my dependencies (Only unpacks really)*/
        enqueue_dependencies(s, tii);
        /*Signal sleeping runners*/
        // MATTHIEU signal_sleeping_runners(s, tii);

        tii->gpu_done = 1;
      }
    }
  }
  /* Zero counters for the next pack operations */
  pack_vars->count_parts = 0;
  pack_vars->tasks_packed = 0;

  /* Write the timers back to the task */
  t->total_cpu_unpack_ticks += total_cpu_unpack_ticks;
  //	/*Time end of unpacking*/
  //	clock_gettime(CLOCK_REALTIME, &t1);
  //	*packing_time += (t1.tv_sec - t0.tv_sec) +
  //	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
} /*End of GPU work*/

void runner_dopair1_launch_f4_f_one_memcpy(
    struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars,
    struct task *t, struct part_aos_f4_f_send *parts_send,
    struct part_aos_f4_f_recv *parts_recv,
    struct part_aos_f4_f_send *d_parts_send,
    struct part_aos_f4_f_recv *d_parts_recv, cudaStream_t *stream, float d_a,
    float d_H, struct engine *e, double *packing_time, double *gpu_time,
    double *unpack_time, int4 *fparti_fpartj_lparti_lpartj,
    cudaEvent_t *pair_end) {

  struct timespec t0, t1, tp0, tp1;  //
  clock_gettime(CLOCK_REALTIME, &t0);

  /* Identify the number of GPU bundles to run in ideal case*/
  int nBundles_temp = pack_vars->nBundles;
  /*How many tasks have we packed?*/
  const int tasks_packed = pack_vars->tasks_packed;

  /*How many tasks should be in a bundle?*/
  const int bundle_size = pack_vars->bundle_size;

  /*tasks-packed needs decrementing before calculating packed_tmp as it was
   * incremented in runner_dopair1_pack*/
  //	const int packed_tmp = 2 * (tasks_packed - 1);

  /* Special case for incomplete bundles (when having leftover tasks not enough
   * to fill a bundle) */
  if (pack_vars->launch_leftovers) {
    nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
    if (tasks_packed == 0)
      error("zero pair tasks packed but somehow got into GPU loop");
    //	  pack_vars->bundle_first_part[nBundles_temp] =
    // pack_vars->task_first_part[packed_tmp - 2];
    pack_vars->bundle_first_part[nBundles_temp] =
        fparti_fpartj_lparti_lpartj[tasks_packed - 1].x;
  }
  /* Identify the last particle for each bundle of tasks */
  for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
  }
  /* special treatment for the last bundle */
  if (nBundles_temp > 1)
    pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
  else
    pack_vars->bundle_last_part[0] = pack_vars->count_parts;

  /* Launch the copies for each bundle and run the GPU kernel */
  /*We don't go into this loop if tasks_left_self == 1 as
   nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
  //	int max_parts = 0;
  for (int bid = 0; bid < nBundles_temp; bid++) {

    int max_parts_i = 0;
    int max_parts_j = 0;
    int parts_in_bundle_ci = 0;
    int parts_in_bundle_cj = 0;
    //      const int first_task = bid * pack_vars->bundle_size;
    //	  int last_task = (bid + 1) * bundle_size;
    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {
      if (tid < tasks_packed) {
        /*Get an estimate for the max number of parts per cell in each bundle.
         *  Used for determining the number of GPU CUDA blocks*/
        int count_i = fparti_fpartj_lparti_lpartj[tid].z -
                      fparti_fpartj_lparti_lpartj[tid].x;
        parts_in_bundle_ci += count_i;
        max_parts_i = max(max_parts_i, count_i);
        int count_j = fparti_fpartj_lparti_lpartj[tid].w -
                      fparti_fpartj_lparti_lpartj[tid].y;
        parts_in_bundle_cj += count_j;
        max_parts_j = max(max_parts_j, count_j);

        //		  last_task = tid;
      }
    }
    const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
    const int bundle_n_parts =
        pack_vars->bundle_last_part[bid] - first_part_tmp_i;

    cudaMemcpyAsync(&d_parts_send[first_part_tmp_i],
                    &parts_send[first_part_tmp_i],
                    bundle_n_parts * sizeof(struct part_aos_f4_f_send),
                    cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
    cudaError_t cu_error =
        cudaPeekAtLastError();  // cudaGetLastError();        //
                                // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(stderr,
              "CUDA error with pair density H2D async  memcpy ci: %s cpuid id "
              "is: %i\n ",
              cudaGetErrorString(cu_error), r->cpuid);
      error("Something's up with your cuda code");
    }
#endif

    //	  const int tasksperbundle = pack_vars->tasksperbundle;
    /* LAUNCH THE GPU KERNELS for ci & cj */
    //      int tid = 0;
    //      int offset = bid * tasksperbundle;
    //      int tasks_left = tasksperbundle;
    //      if (bid == nBundles_temp - 1) {
    //        tasks_left =
    //        		tasks_packed - (nBundles_temp - 1) * tasksperbundle;
    //      }

    // Setup 2d grid of GPU thread blocks for ci (number of tasks is
    // the y dimension and max_parts is the x dimension
    int numBlocks_y = 0;  // tasks_left;
    int numBlocks_x = (bundle_n_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int bundle_part_0 = pack_vars->bundle_first_part[bid];
    //      int bundle_first_task = pack_vars->bundle_first_task_list[bid];
    //              fprintf(stderr, "bundle_part_0 %i bundle_first_task %i\n",
    //              bundle_part_0, bundle_first_task);

    /* Launch the kernel for ci using data for ci and cj */
    runner_dopair_branch_force_gpu_aos_f4(d_parts_send, d_parts_recv, d_a, d_H,
                                          stream[bid], numBlocks_x, numBlocks_y,
                                          bundle_part_0, bundle_n_parts);

#ifdef CUDA_DEBUG
    cu_error = cudaPeekAtLastError();  // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(
          stderr,
          "CUDA error with pair density kernel launch: %s cpuid id is: %i\n "
          "nbx %i nby %i max_parts_i %i max_parts_j %i\n",
          cudaGetErrorString(cu_error), r->cpuid, numBlocks_x, numBlocks_y,
          max_parts_i, max_parts_j);
      exit(0);
    }
#endif

    // Copy results back to CPU BUFFERS
    cudaMemcpyAsync(&parts_recv[first_part_tmp_i],
                    &d_parts_recv[first_part_tmp_i],
                    bundle_n_parts * sizeof(struct part_aos_f4_f_recv),
                    cudaMemcpyDeviceToHost, stream[bid]);
    cudaEventRecord(pair_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
    cu_error = cudaPeekAtLastError();  // cudaGetLastError();        //
                                       // Get error code
    if (cu_error != cudaSuccess) {
      fprintf(stderr,
              "CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
              cudaGetErrorString(cu_error), r->cpuid);
      error("Something's up with your cuda code");
    }
#endif
  } /*End of looping over bundles to launch in streams*/

  /* Make sure all the kernels and copies back are finished */
  //	cudaDeviceSynchronize();

  /*Time end of GPU work*/
  clock_gettime(CLOCK_REALTIME, &t1);
  *gpu_time +=
      (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  /* Now copy the data back from the CPU thread-local buffers to the cells */
  /* Pack length counter for use in unpacking */
  int pack_length_unpack = 0;
  ticks total_cpu_unpack_ticks = 0.;
  for (int bid = 0; bid < nBundles_temp; bid++) {
    /*Time unpacking*/
    clock_gettime(CLOCK_REALTIME, &t0);

    //		cudaStreamSynchronize(stream[bid]);
    cudaEventSynchronize(pair_end[bid]);

    clock_gettime(CLOCK_REALTIME, &t1);
    *gpu_time +=
        (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

    /*Time unpacking*/
    //		clock_gettime(CLOCK_REALTIME, &tp0);
    //		int bundle_first_task = pack_vars->bundle_first_task_list[bid];

    for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

      if (tid < tasks_packed) {
        clock_gettime(CLOCK_REALTIME, &tp0);
        /*grab cell and task pointers*/
        struct cell *cii = pack_vars->ci_list[tid];
        struct cell *cjj = pack_vars->cj_list[tid];
        struct task *tii = pack_vars->task_list[tid];
        /*Let's lock ci*/
        while (cell_locktree(cii)) {
          ; /* spin until we acquire the lock */
        }
        /*Let's lock cj*/
        while (cell_locktree(cjj)) {
          ; /* spin until we acquire the lock */
        }

        const ticks tic = getticks();

        /* Do the copy */
        runner_do_ci_cj_gpu_unpack_neat_aos_f4_f(
            r, cii, cjj, parts_recv, 0, &pack_length_unpack, tid,
            2 * pack_vars->count_max_parts, e);

        const ticks toc = getticks();

        total_cpu_unpack_ticks += toc - tic;

        /* Record things for debugging */
        cii->gpu_done_pair_f++;
        cjj->gpu_done_pair_f++;
        pthread_mutex_lock(&s->sleep_mutex);
        atomic_dec(&s->waiting);
        pthread_cond_broadcast(&s->sleep_cond);
        pthread_mutex_unlock(&s->sleep_mutex);
        //		  /* Release the locks */
        cell_unlocktree(cii);
        //		  /* Release the locks */
        cell_unlocktree(cjj);

        /*Time end of unpacking*/
        clock_gettime(CLOCK_REALTIME, &tp1);
        *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
                        (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;

        /*schedule my dependencies (Only unpacks really)*/
        enqueue_dependencies(s, tii);
        /*Signal sleeping runners*/
        // MATTHIEU signal_sleeping_runners(s, tii);

        tii->gpu_done = 1;
      }
    }
  }
  /* Zero counters for the next pack operations */
  pack_vars->count_parts = 0;
  pack_vars->tasks_packed = 0;

  /* Write the timers back to the task */
  t->total_cpu_unpack_ticks += total_cpu_unpack_ticks;
  //	/*Time end of unpacking*/
  //	clock_gettime(CLOCK_REALTIME, &t1);
  //	*packing_time += (t1.tv_sec - t0.tv_sec) +
  //	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
} /*End of GPU work*/
