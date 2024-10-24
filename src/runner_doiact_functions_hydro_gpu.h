#include "scheduler.h"
struct pack_vars_self{
    /*List of tasks and respective cells to be packed*/
	struct task **task_list;
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
	int *task_first_part;
	int *task_last_part;
	int *d_task_first_part;
	int *d_task_last_part;
	int * bundle_first_part;
	int * bundle_last_part;
	int * bundle_first_task_list;
	int count_max_parts;
	int launch;
	int launch_leftovers;
	int target_n_tasks;
	int nBundles;
	int tasksperbundle;

}pack_vars_self;

struct pack_vars_pair{
    /*List of tasks and respective cells to be packed*/
	struct task **task_list;
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
	int *task_first_part;
	int *task_last_part;
	int *d_task_first_part;
	int *d_task_last_part;
	int * bundle_first_part;
	int * bundle_last_part;
	int * bundle_first_task_list;
	int count_max_parts;
	int launch;
	int launch_leftovers;
	int target_n_tasks;
	int nBundles;
	int tasksperbundle;

}pack_vars_pair;

struct pack_vars_pair_f4{
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
	int * bundle_first_part;
	int * bundle_last_part;
	int * bundle_first_task_list;
	int count_max_parts;
	int launch;
	int launch_leftovers;
	int target_n_tasks;
	int nBundles;
	int tasksperbundle;

}pack_vars_pair_f4;

#include "task.h"
#include "runner_gpu_pack_functions.h"
#include "cuda/BLOCK_SIZE.h"
#include "cuda/GPU_runner_functions.h"
#define CUDA_DEBUG
void runner_doself1_pack(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos *parts_aos, int * packing_time){
    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	int tasks_packed = pack_vars->tasks_packed;
	pack_vars->cellx[tasks_packed] = ci->loc[0];
	pack_vars->celly[tasks_packed] = ci->loc[1];
	pack_vars->cellz[tasks_packed] = ci->loc[2];
	/*Get pointers to the list of tasks and cells packed*/
	pack_vars->task_list[tasks_packed] = t;
	pack_vars->cell_list[tasks_packed] = ci;
	//    /* Identify row in particle arrays where this task starts*/
	pack_vars->task_first_part[tasks_packed] = pack_vars->count_parts;
	int *count_parts_self = &pack_vars->count_parts;
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
	runner_doself1_gpu_pack_neat_aos(r, ci, parts_aos, 0/*timer. 0 no timing, 1 for timing*/,
		  count_parts_self, tasks_packed, pack_vars->count_max_parts);
	//    // identify the row in the array where this task ends (row id of its last particle)
	pack_vars->task_last_part[tasks_packed] = pack_vars->count_parts;
	/* Identify first particle for each bundle of tasks */
	const int bundle_size = pack_vars->bundle_size;
	if (tasks_packed % bundle_size == 0) {
	  int bid = tasks_packed / bundle_size;
	  pack_vars->bundle_first_part[bid] = pack_vars->task_first_part[tasks_packed];
	  pack_vars->bundle_first_task_list[bid] = tasks_packed;
	}
	/* Tell the cell it has been packed */
	ci->pack_done++;
	/* Record that we have now done a packing (self) */
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_self_left));
	t->done = 1;
	/* Release the lock on the cell */
	task_unlock(t);
	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_self_left == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
    clock_gettime(CLOCK_REALTIME, &t1);
    *packing_time += (t1.tv_sec - t0.tv_sec) +
                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

}

double runner_doself1_pack_f4(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_f4_send *parts_send, int2 *task_first_part_f4){
    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	int tasks_packed = pack_vars->tasks_packed;
//	pack_vars->cellx[tasks_packed] = ci->loc[0];
//	pack_vars->celly[tasks_packed] = ci->loc[1];
//	pack_vars->cellz[tasks_packed] = ci->loc[2];
	/*Get pointers to the list of tasks and cells packed*/
	pack_vars->task_list[tasks_packed] = t;
	pack_vars->cell_list[tasks_packed] = ci;
	//    /* Identify row in particle arrays where this task starts*/
//	pack_vars->task_first_part[tasks_packed] = pack_vars->count_parts;
//	d_task_first_part_self_dens_f4[tasks_packed].x = pack_vars->count_parts;
	task_first_part_f4[tasks_packed].x = pack_vars->count_parts;
	int *count_parts_self = &pack_vars->count_parts;
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
	runner_doself1_gpu_pack_neat_aos_f4(r, ci, parts_send, 0/*timer. 0 no timing, 1 for timing*/,
		  count_parts_self, tasks_packed, pack_vars->count_max_parts);
	//    // identify the row in the array where this task ends (row id of its last particle)
//	pack_vars->task_last_part[tasks_packed] = pack_vars->count_parts;
//	d_task_first_part_self_dens_f4[tasks_packed].y = pack_vars->count_parts;
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
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_self_left));
	t->done = 1;
	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_self_left == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
	  clock_gettime(CLOCK_REALTIME, &t1);
		/* Release the lock on the cell */
//		task_unlock(t);
		cell_unlocktree(ci);
	  return (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

}

void runner_doself1_pack_g(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_g *parts_aos, double *packing_time){

    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	int tasks_packed = pack_vars->tasks_packed;
	pack_vars->cellx[tasks_packed] = ci->loc[0];
	pack_vars->celly[tasks_packed] = ci->loc[1];
	pack_vars->cellz[tasks_packed] = ci->loc[2];
	/*Get pointers to the list of tasks and cells packed*/
	pack_vars->task_list[tasks_packed] = t;
	pack_vars->cell_list[tasks_packed] = ci;
	//    /* Identify row in particle arrays where this task starts*/
	pack_vars->task_first_part[tasks_packed] = pack_vars->count_parts;
	int *count_parts_self = &pack_vars->count_parts;
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
	runner_doself1_gpu_pack_neat_aos_g(r, ci, parts_aos, 0/*timer. 0 no timing, 1 for timing*/,
		  count_parts_self, tasks_packed, pack_vars->count_max_parts);
	//    // identify the row in the array where this task ends (row id of its last particle)
	pack_vars->task_last_part[tasks_packed] = pack_vars->count_parts;
	/* Identify first particle for each bundle of tasks */
	const int bundle_size = pack_vars->bundle_size;
	if (tasks_packed % bundle_size == 0) {
	  int bid = tasks_packed / bundle_size;
	  pack_vars->bundle_first_part[bid] = pack_vars->task_first_part[tasks_packed];
	  pack_vars->bundle_first_task_list[bid] = tasks_packed;
	}
	/* Tell the cell it has been packed */
	ci->pack_done_g++;
	/* Record that we have now done a packing (self) */
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_self_left_g));
	t->done = 1;
	/* Release the lock on the cell */
	task_unlock(t);
	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_self_left_g == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
    clock_gettime(CLOCK_REALTIME, &t1);
    *packing_time += (t1.tv_sec - t0.tv_sec) +
                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

}

double runner_doself1_pack_f4_g(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_f4_g_send *parts_send, int2 * task_first_part_f4){

    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	int tasks_packed = pack_vars->tasks_packed;
//	pack_vars->cellx[tasks_packed] = ci->loc[0];
//	pack_vars->celly[tasks_packed] = ci->loc[1];
//	pack_vars->cellz[tasks_packed] = ci->loc[2];
	/*Get pointers to the list of tasks and cells packed*/
	pack_vars->task_list[tasks_packed] = t;
	pack_vars->cell_list[tasks_packed] = ci;
	//    /* Identify row in particle arrays where this task starts*/
//	pack_vars->task_first_part[tasks_packed] = pack_vars->count_parts;
	task_first_part_f4[tasks_packed].x = pack_vars->count_parts;
	int *count_parts_self = &pack_vars->count_parts;
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
	runner_doself1_gpu_pack_neat_aos_f4_g(r, ci, parts_send, 0/*timer. 0 no timing, 1 for timing*/,
		  count_parts_self, tasks_packed, pack_vars->count_max_parts);
	//    // identify the row in the array where this task ends (row id of its last particle)
//	pack_vars->task_last_part[tasks_packed] = pack_vars->count_parts;
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
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_self_left_g));
	t->done = 1;
	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_self_left_g == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
    clock_gettime(CLOCK_REALTIME, &t1);
	/* Release the lock on the cell */
//	task_unlock(t);
	cell_unlocktree(ci);
    return (t1.tv_sec - t0.tv_sec) +
                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}

void runner_doself1_pack_f(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_f *parts_aos, double *packing_time){

    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	int tasks_packed = pack_vars->tasks_packed;
	pack_vars->cellx[tasks_packed] = ci->loc[0];
	pack_vars->celly[tasks_packed] = ci->loc[1];
	pack_vars->cellz[tasks_packed] = ci->loc[2];
	/*Get pointers to the list of tasks and cells packed*/
	pack_vars->task_list[tasks_packed] = t;
	pack_vars->cell_list[tasks_packed] = ci;
	//    /* Identify row in particle arrays where this task starts*/
	pack_vars->task_first_part[tasks_packed] = pack_vars->count_parts;
	int *count_parts_self = &pack_vars->count_parts;
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
	runner_doself1_gpu_pack_neat_aos_f(r, ci, parts_aos, 0/*timer. 0 no timing, 1 for timing*/,
		  count_parts_self, tasks_packed, pack_vars->count_max_parts);
	//    // identify the row in the array where this task ends (row id of its last particle)
	pack_vars->task_last_part[tasks_packed] = pack_vars->count_parts;
	/* Identify first particle for each bundle of tasks */
	const int bundle_size = pack_vars->bundle_size;
	if (tasks_packed % bundle_size == 0) {
	  int bid = tasks_packed / bundle_size;
	  pack_vars->bundle_first_part[bid] = pack_vars->task_first_part[tasks_packed];
	  pack_vars->bundle_first_task_list[bid] = tasks_packed;
	}
	/* Tell the cell it has been packed */
	ci->pack_done_f++;
	/* Record that we have now done a packing (self) */
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_self_left_f));
	t->done = 1;
	/* Release the lock on the cell */
	task_unlock(t);
	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_self_left_f == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
    clock_gettime(CLOCK_REALTIME, &t1);
    *packing_time += (t1.tv_sec - t0.tv_sec) +
                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}

double runner_doself1_pack_f4_f(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_f4_f_send *parts_send, int2 * task_first_part_f4){

    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	int tasks_packed = pack_vars->tasks_packed;
//	pack_vars->cellx[tasks_packed] = ci->loc[0];
//	pack_vars->celly[tasks_packed] = ci->loc[1];
//	pack_vars->cellz[tasks_packed] = ci->loc[2];
	/*Get pointers to the list of tasks and cells packed*/
	pack_vars->task_list[tasks_packed] = t;
	pack_vars->cell_list[tasks_packed] = ci;
	//    /* Identify row in particle arrays where this task starts*/
//	pack_vars->task_first_part[tasks_packed] = pack_vars->count_parts;
	task_first_part_f4[tasks_packed].x = pack_vars->count_parts;
	int *count_parts_self = &pack_vars->count_parts;
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
	runner_doself1_gpu_pack_neat_aos_f4_f(r, ci, parts_send, 0/*timer. 0 no timing, 1 for timing*/,
		  count_parts_self, tasks_packed, pack_vars->count_max_parts);
	//    // identify the row in the array where this task ends (row id of its last particle)
//	pack_vars->task_last_part[tasks_packed] = pack_vars->count_parts;
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
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_self_left_f));
	t->done = 1;
	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_self_left_f == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
    clock_gettime(CLOCK_REALTIME, &t1);
	/* Release the lock on the cell */
//	task_unlock(t);
	cell_unlocktree(ci);
    return (t1.tv_sec - t0.tv_sec) +
                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}

void runner_dopair1_pack(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars, struct cell *ci,
		 struct cell *cj, struct task *t, struct part_aos *parts_aos, struct engine *e, double *packing_time){

    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);
	int tasks_packed = pack_vars->tasks_packed;

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
    /*Indexing increment per task is 2 fot these arrays*/
    const int packed_tmp = tasks_packed * 2;

    /* Find first parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */
    pack_vars->task_first_part[packed_tmp] = pack_vars->count_parts;
    pack_vars->task_first_part[packed_tmp + 1] = pack_vars->count_parts + count_ci;

    int *count_parts = &pack_vars->count_parts;
//    if(r->cpuid == 0)fprintf(stderr, "cpu %i before count %i\n", r->cpuid, pack_vars->count_parts);
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
    runner_do_ci_cj_gpu_pack_neat_aos(r, ci, cj, parts_aos, 0/*timer. 0 no timing, 1 for timing*/,
    		count_parts, tid, pack_vars->count_max_parts, count_ci, count_cj, shift_tmp);
//	runner_doself1_gpu_pack_neat_aos(r, ci, parts_aos, 0/*timer. 0 no timing, 1 for timing*/,
//		  count_parts, tasks_packed, pack_vars->count_max_parts); //This may cause an issue.
	                                                              //Be sure to test that pack_vars->count_parts
	                                                              //is actually increment here
	/* Find last parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */

//    if(r->cpuid == 0)fprintf(stderr, "cpu %i after count %i pack_vars_count %i\n", r->cpuid, *count_parts,
//    		pack_vars->count_parts);
    pack_vars->task_last_part[packed_tmp] = pack_vars->count_parts - count_cj;
    pack_vars->task_last_part[packed_tmp + 1] = pack_vars->count_parts;

	/* Tell the cells they have been packed */
	ci->pack_done++;
	cj->pack_done++;

	/* Identify first particle for each bundle of tasks */
	const int bundle_size = pack_vars->bundle_size;
	if (tasks_packed % bundle_size == 0) {
	  int bid = tasks_packed / bundle_size;
	  pack_vars->bundle_first_part[bid] = pack_vars->task_first_part[packed_tmp];
	  pack_vars->bundle_first_task_list[bid] = tasks_packed;
	}

	/* Record that we have now done a packing (self) */
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_pair_left));
	t->done = 1;

	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_pair_left == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
    clock_gettime(CLOCK_REALTIME, &t1);
    *packing_time += (t1.tv_sec - t0.tv_sec) +
                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
	/* Copies done. Release the lock ! */
//	task_unlock(t);
	cell_unlocktree(ci);
	cell_unlocktree(cj);
}

double runner_dopair1_pack_f4(struct runner *r, struct scheduler *s, struct pack_vars_pair * restrict pack_vars, struct cell * restrict ci,
		 struct cell * restrict cj, struct task *t, struct part_aos_f4_send *parts_send, struct engine *e, int4 *fparti_fpartj_lparti_lpartj){

    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);
	int tasks_packed = pack_vars->tasks_packed;

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

    /* Find first parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */
//    pack_vars->task_first_part[packed_tmp] = pack_vars->count_parts;
//    pack_vars->task_first_part[packed_tmp + 1] = pack_vars->count_parts + count_ci;

    fparti_fpartj_lparti_lpartj[tasks_packed].x = pack_vars->count_parts;
    fparti_fpartj_lparti_lpartj[tasks_packed].y = pack_vars->count_parts + count_ci;

    int *count_parts = &pack_vars->count_parts;
//    if(r->cpuid == 0)fprintf(stderr, "cpu %i before count %i\n", r->cpuid, pack_vars->count_parts);
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
    runner_do_ci_cj_gpu_pack_neat_aos_f4(r, ci, cj, parts_send, 0/*timer. 0 no timing, 1 for timing*/,
    		count_parts, tid, pack_vars->count_max_parts, count_ci, count_cj, shift_tmp);
//	runner_doself1_gpu_pack_neat_aos(r, ci, parts_aos, 0/*timer. 0 no timing, 1 for timing*/,
//		  count_parts, tasks_packed, pack_vars->count_max_parts); //This may cause an issue.
	                                                              //Be sure to test that pack_vars->count_parts
	                                                              //is actually increment here
	/* Find last parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */

//    if(r->cpuid == 0)fprintf(stderr, "cpu %i after count %i pack_vars_count %i\n", r->cpuid, *count_parts,
//    		pack_vars->count_parts);
    fparti_fpartj_lparti_lpartj[tasks_packed].z = pack_vars->count_parts - count_cj;
    fparti_fpartj_lparti_lpartj[tasks_packed].w = pack_vars->count_parts;
//    pack_vars->task_last_part[packed_tmp] = pack_vars->count_parts - count_cj;
//    pack_vars->task_last_part[packed_tmp + 1] = pack_vars->count_parts;

	/* Tell the cells they have been packed */
	ci->pack_done++;
	cj->pack_done++;

	/* Identify first particle for each bundle of tasks */
	const int bundle_size = pack_vars->bundle_size;
	if (tasks_packed % bundle_size == 0) {
	  int bid = tasks_packed / bundle_size;
	  pack_vars->bundle_first_part[bid] = fparti_fpartj_lparti_lpartj[tasks_packed].x;
	  pack_vars->bundle_first_task_list[bid] = tasks_packed;
	}

	/* Record that we have now done a packing (self) */
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_pair_left));
	t->done = 1;
	/* Copies done. Release the lock ! */
	task_unlock(t);
	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_pair_left == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
    clock_gettime(CLOCK_REALTIME, &t1);
    return (t1.tv_sec - t0.tv_sec) +
                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}

void runner_dopair1_pack_g(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars, struct cell *ci,
		 struct cell *cj, struct task *t, struct part_aos_g *parts_aos, struct engine *e, double * packing_time){

    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);
	int tasks_packed = pack_vars->tasks_packed;
	const int tid_tmp = 2 * tasks_packed;
	/*shifts for ci*/
	pack_vars->shiftx[tid_tmp] = 0.0;
	pack_vars->shifty[tid_tmp] = 0.0;
	pack_vars->shiftz[tid_tmp] = 0.0;
	/*shifts for cj. Stored using strided indexing (stride of two per task)*/
	pack_vars->shiftx[tid_tmp + 1] = 0.0;
	pack_vars->shifty[tid_tmp + 1] = 0.0;
	pack_vars->shiftz[tid_tmp + 1] = 0.0;

    double x_tmp = 0.0, y_tmp = 0.0, z_tmp = 0.0;
    /*Get the shifts in case of periodics*/
    space_getsid_GPU(e->s, &ci, &cj, &x_tmp, &y_tmp, &z_tmp);

	/*Get pointers to the list of tasks and cells packed*/
	pack_vars->task_list[tasks_packed] = t;
	pack_vars->ci_list[tasks_packed] = ci;
	pack_vars->cj_list[tasks_packed] = cj;

    const double cjx = cj->loc[0];
    const double cjy = cj->loc[1];
    const double cjz = cj->loc[2];

    /*Correct the shifts for cell i*/
	pack_vars->shiftx[tid_tmp] = x_tmp + cjx;
	pack_vars->shifty[tid_tmp] = y_tmp + cjy;
	pack_vars->shiftz[tid_tmp] = z_tmp + cjz;
	/*Shift for cell j is it's position. Stored using strided indexing (stride of two per task)*/
	pack_vars->shiftx[tid_tmp + 1] = cjx;
	pack_vars->shifty[tid_tmp + 1] = cjy;
	pack_vars->shiftz[tid_tmp + 1] = cjz;

    const int count_ci = ci->hydro.count;
    const int count_cj = cj->hydro.count;

    /*Assign an id for this task*/
    const int tid = tasks_packed;
    /*Indexing increment per task is 2 fot these arrays*/
    const int packed_tmp = tasks_packed * 2;

    /* Find first parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */
    pack_vars->task_first_part[packed_tmp] = pack_vars->count_parts;
    pack_vars->task_first_part[packed_tmp + 1] = pack_vars->count_parts + count_ci;

    int *count_parts = &pack_vars->count_parts;
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
    runner_do_ci_cj_gpu_pack_neat_aos_g(r, ci, cj, parts_aos, 0/*timer. 0 no timing, 1 for timing*/,
    		count_parts, tid, pack_vars->count_max_parts, count_ci, count_cj);
	/* Find last parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */

    pack_vars->task_last_part[packed_tmp] = pack_vars->count_parts - count_cj;
    pack_vars->task_last_part[packed_tmp + 1] = pack_vars->count_parts;

	/* Tell the cells they have been packed */
	ci->pack_done_g++;
	cj->pack_done_g++;

	/* Identify first particle for each bundle of tasks */
	const int bundle_size = pack_vars->bundle_size;
	if (tasks_packed % bundle_size == 0) {
	  int bid = tasks_packed / bundle_size;
	  pack_vars->bundle_first_part[bid] = pack_vars->task_first_part[packed_tmp];
	  pack_vars->bundle_first_task_list[bid] = tasks_packed;
	}

	/* Record that we have now done a packing (self) */
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_pair_left_g));
	t->done = 1;
	/* Copies done. Release the lock ! */
	task_unlock(t);
	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_pair_left_g == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
    clock_gettime(CLOCK_REALTIME, &t1);
    *packing_time += (t1.tv_sec - t0.tv_sec) +
                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}


double runner_dopair1_pack_f4_g(struct runner *r, struct scheduler *s, struct pack_vars_pair * restrict pack_vars, struct cell * restrict ci,
		 struct cell * restrict cj, struct task *t, struct part_aos_f4_g_send *parts_send, struct engine *e, int4 *fparti_fpartj_lparti_lpartj){

    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);
	int tasks_packed = pack_vars->tasks_packed;

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

    /* Find first parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */
//    pack_vars->task_first_part[packed_tmp] = pack_vars->count_parts;
//    pack_vars->task_first_part[packed_tmp + 1] = pack_vars->count_parts + count_ci;

    fparti_fpartj_lparti_lpartj[tasks_packed].x = pack_vars->count_parts;
    fparti_fpartj_lparti_lpartj[tasks_packed].y = pack_vars->count_parts + count_ci;

    int *count_parts = &pack_vars->count_parts;
//    if(r->cpuid == 0)fprintf(stderr, "cpu %i before count %i\n", r->cpuid, pack_vars->count_parts);
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
    runner_do_ci_cj_gpu_pack_neat_aos_f4_g(r, ci, cj, parts_send, 0/*timer. 0 no timing, 1 for timing*/,
    		count_parts, tid, pack_vars->count_max_parts, count_ci, count_cj, shift_tmp);
//	runner_doself1_gpu_pack_neat_aos(r, ci, parts_aos, 0/*timer. 0 no timing, 1 for timing*/,
//		  count_parts, tasks_packed, pack_vars->count_max_parts); //This may cause an issue.
	                                                              //Be sure to test that pack_vars->count_parts
	                                                              //is actually increment here
	/* Find last parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */

//    if(r->cpuid == 0)fprintf(stderr, "cpu %i after count %i pack_vars_count %i\n", r->cpuid, *count_parts,
//    		pack_vars->count_parts);
    fparti_fpartj_lparti_lpartj[tasks_packed].z = pack_vars->count_parts - count_cj;
    fparti_fpartj_lparti_lpartj[tasks_packed].w = pack_vars->count_parts;
//    pack_vars->task_last_part[packed_tmp] = pack_vars->count_parts - count_cj;
//    pack_vars->task_last_part[packed_tmp + 1] = pack_vars->count_parts;

	/* Tell the cells they have been packed */
	ci->pack_done_g++;
	cj->pack_done_g++;

	/* Identify first particle for each bundle of tasks */
	const int bundle_size = pack_vars->bundle_size;
	if (tasks_packed % bundle_size == 0) {
	  int bid = tasks_packed / bundle_size;
	  pack_vars->bundle_first_part[bid] = fparti_fpartj_lparti_lpartj[tasks_packed].x;
	  pack_vars->bundle_first_task_list[bid] = tasks_packed;
	}

	/* Record that we have now done a packing (self) */
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_pair_left_g));
	t->done = 1;
	/* Copies done. Release the lock ! */
	task_unlock(t);
	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_pair_left_g == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
    clock_gettime(CLOCK_REALTIME, &t1);
    return (t1.tv_sec - t0.tv_sec) +
                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}


void runner_dopair1_pack_f(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars, struct cell *ci,
		 struct cell *cj, struct task *t, struct part_aos_f *parts_aos, struct engine *e, double *packing_time){

    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);
	int tasks_packed = pack_vars->tasks_packed; //Copy pasted this code again. Issue isn't here
	const int tid_tmp = 2 * tasks_packed;
	/*shifts for ci*/
	pack_vars->shiftx[tid_tmp] = 0.0;
	pack_vars->shifty[tid_tmp] = 0.0;
	pack_vars->shiftz[tid_tmp] = 0.0;
	/*shifts for cj. Stored using strided indexing (stride of two per task)*/
	pack_vars->shiftx[tid_tmp + 1] = 0.0;
	pack_vars->shifty[tid_tmp + 1] = 0.0;
	pack_vars->shiftz[tid_tmp + 1] = 0.0;

    double x_tmp = 0.0, y_tmp = 0.0, z_tmp = 0.0;
    /*Get the shifts in case of periodics*/
    space_getsid_GPU(e->s, &ci, &cj, &x_tmp, &y_tmp, &z_tmp);

	/*Get pointers to the list of tasks and cells packed*/
	pack_vars->task_list[tasks_packed] = t;
	pack_vars->ci_list[tasks_packed] = ci;
	pack_vars->cj_list[tasks_packed] = cj;

    const double cjx = cj->loc[0];
    const double cjy = cj->loc[1];
    const double cjz = cj->loc[2];

    /*Correct the shifts for cell i*/
	pack_vars->shiftx[tid_tmp] = x_tmp + cjx;
	pack_vars->shifty[tid_tmp] = y_tmp + cjy;
	pack_vars->shiftz[tid_tmp] = z_tmp + cjz;
	/*Shift for cell j is it's position. Stored using strided indexing (stride of two per task)*/
	pack_vars->shiftx[tid_tmp + 1] = cjx;
	pack_vars->shifty[tid_tmp + 1] = cjy;
	pack_vars->shiftz[tid_tmp + 1] = cjz;

    const int count_ci = ci->hydro.count;
    const int count_cj = cj->hydro.count;

    /*Assign an id for this task*/
    const int tid = tasks_packed;
    /*Indexing increment per task is 2 fot these arrays*/
    const int packed_tmp = tasks_packed * 2;

    /* Find first parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */
    pack_vars->task_first_part[packed_tmp] = pack_vars->count_parts;
    pack_vars->task_first_part[packed_tmp + 1] = pack_vars->count_parts + count_ci;

    int *count_parts = &pack_vars->count_parts;
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
    runner_do_ci_cj_gpu_pack_neat_aos_f(r, ci, cj, parts_aos, 0/*timer. 0 no timing, 1 for timing*/,
    		count_parts, tid, pack_vars->count_max_parts, count_ci, count_cj);
	/* Find last parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */

    pack_vars->task_last_part[packed_tmp] = pack_vars->count_parts - count_cj;
    pack_vars->task_last_part[packed_tmp + 1] = pack_vars->count_parts;

	/* Tell the cells they have been packed */
	ci->pack_done_f++;
	cj->pack_done_f++;

	/* Identify first particle for each bundle of tasks */
	const int bundle_size = pack_vars->bundle_size;
	if (tasks_packed % bundle_size == 0) {
	  int bid = tasks_packed / bundle_size;
	  pack_vars->bundle_first_part[bid] = pack_vars->task_first_part[packed_tmp];
	  pack_vars->bundle_first_task_list[bid] = tasks_packed;
	}

	/* Record that we have now done a packing (self) */
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_pair_left_f));
	t->done = 1;
	/* Copies done. Release the lock ! */
	task_unlock(t);
	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_pair_left_f == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
    clock_gettime(CLOCK_REALTIME, &t1);
    *packing_time += (t1.tv_sec - t0.tv_sec) +
                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}

double runner_dopair1_pack_f4_f(struct runner *r, struct scheduler *s, struct pack_vars_pair * restrict pack_vars, struct cell * restrict ci,
		 struct cell * restrict cj, struct task *t, struct part_aos_f4_f_send *parts_send, struct engine *e, int4 *fparti_fpartj_lparti_lpartj){

    /* Timers for how long this all takes.
    * t0 and t1 are from start to finish including GPU calcs
    * tp0 and tp1 only time packing and unpacking*/
    struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);
	int tasks_packed = pack_vars->tasks_packed;

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

    /* Find first parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */
//    pack_vars->task_first_part[packed_tmp] = pack_vars->count_parts;
//    pack_vars->task_first_part[packed_tmp + 1] = pack_vars->count_parts + count_ci;

    fparti_fpartj_lparti_lpartj[tasks_packed].x = pack_vars->count_parts;
    fparti_fpartj_lparti_lpartj[tasks_packed].y = pack_vars->count_parts + count_ci;

    int *count_parts = &pack_vars->count_parts;
//    if(r->cpuid == 0)fprintf(stderr, "cpu %i before count %i\n", r->cpuid, pack_vars->count_parts);
	/* This re-arranges the particle data from cell->hydro->parts into a
	long array of part structs*/
    runner_do_ci_cj_gpu_pack_neat_aos_f4_f(r, ci, cj, parts_send, 0/*timer. 0 no timing, 1 for timing*/,
    		count_parts, tid, pack_vars->count_max_parts, count_ci, count_cj, shift_tmp);
//	runner_doself1_gpu_pack_neat_aos(r, ci, parts_aos, 0/*timer. 0 no timing, 1 for timing*/,
//		  count_parts, tasks_packed, pack_vars->count_max_parts); //This may cause an issue.
	                                                              //Be sure to test that pack_vars->count_parts
	                                                              //is actually increment here
	/* Find last parts in task for ci and cj. Packed_tmp is index for cell i. packed_tmp+1 is index for cell j */

//    if(r->cpuid == 0)fprintf(stderr, "cpu %i after count %i pack_vars_count %i\n", r->cpuid, *count_parts,
//    		pack_vars->count_parts);
    fparti_fpartj_lparti_lpartj[tasks_packed].z = pack_vars->count_parts - count_cj;
    fparti_fpartj_lparti_lpartj[tasks_packed].w = pack_vars->count_parts;
//    pack_vars->task_last_part[packed_tmp] = pack_vars->count_parts - count_cj;
//    pack_vars->task_last_part[packed_tmp + 1] = pack_vars->count_parts;

	/* Tell the cells they have been packed */
	ci->pack_done_f++;
	cj->pack_done_f++;

	/* Identify first particle for each bundle of tasks */
	const int bundle_size = pack_vars->bundle_size;
	if (tasks_packed % bundle_size == 0) {
	  int bid = tasks_packed / bundle_size;
	  pack_vars->bundle_first_part[bid] = fparti_fpartj_lparti_lpartj[tasks_packed].x;
	  pack_vars->bundle_first_task_list[bid] = tasks_packed;
	}

	/* Record that we have now done a packing (self) */
	int qid = r->qid;
	atomic_dec(&(s->queues[qid].n_packs_pair_left_f));
	t->done = 1;
	/* Copies done. Release the lock ! */
	task_unlock(t);
	pack_vars->tasks_packed++;
	pack_vars->launch = 0;
	pack_vars->launch_leftovers = 0;
	if ((s->queues[qid].n_packs_pair_left_f == 0))
		pack_vars->launch_leftovers = 1;
	if(pack_vars->tasks_packed == pack_vars->target_n_tasks)
		pack_vars->launch = 1;
    /*Add time to packing_time. Timer for end of GPU work after the if(launch || launch_leftovers statement)*/
    clock_gettime(CLOCK_REALTIME, &t1);
    return (t1.tv_sec - t0.tv_sec) +
                    (t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
}

void runner_doself1_launch(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos *parts_aos, struct part_aos *d_parts_aos, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time, double *hmemcpy_time){

	struct timespec t0, t1, t0hmemcpy, t1hmemcpy, tp0, tp1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;

	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero tasks packed but somehow got into GPU loop");
	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[tasks_packed - 1];
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;
    clock_gettime(CLOCK_REALTIME, &t0hmemcpy);
    /*Copy arrays containing first and last part for each task to GPU*/
    cudaMemcpy(pack_vars->d_task_first_part, pack_vars->task_first_part,
               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_task_last_part, pack_vars->task_last_part,
               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);

    /*Copy cell shifts to device*/
    cudaMemcpy(pack_vars->d_cellx, pack_vars->cellx,
               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_celly, pack_vars->celly,
               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_cellz, pack_vars->cellz,
               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
    clock_gettime(CLOCK_REALTIME, &t1hmemcpy);
    *hmemcpy_time += (t1hmemcpy.tv_sec - t0hmemcpy.tv_sec) +
			(t1hmemcpy.tv_nsec - t0hmemcpy.tv_nsec) / 1000000000.0;
	/* Launch the copies for each bundle and run the GPU kernel */
	/*We don't go into this loop if tasks_left_self == 1 as
	 nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
	int max_parts = 0;
	for (int bid = 0; bid < nBundles_temp; bid++) {

	  max_parts = 0;
	  int parts_in_bundle = 0;
	  for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
		   tid++) {
		if (tid < tasks_packed) {
		  /*Get an estimate for the max number of parts per cell in the bundle.
		   *  Used for determining the number of GPU CUDA blocks*/
		  int count = pack_vars->task_last_part[tid] - pack_vars->task_first_part[tid];
		  parts_in_bundle += count;
		  max_parts = max(max_parts, count);
		}
	  }

	  const int first_part_tmp = pack_vars->bundle_first_part[bid];
	  const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp;

	  cudaMemcpyAsync(&d_parts_aos[first_part_tmp], &parts_aos[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos), cudaMemcpyHostToDevice, stream[bid]);

//#ifdef CUDA_DEBUG
//	  cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
//										// Get error code
//	  if (cu_error != cudaSuccess) {
//		fprintf(
//			stderr,
//			"CUDA error in density self host 2 device memcpy: %s cpuid id is: %i\n ",
//			cudaGetErrorString(cu_error), r->cpuid);
//		exit(0);
//	  }
//#endif
	  const int tasksperbundle = pack_vars->tasksperbundle;
	  int tasks_left = tasksperbundle;
	  if (bid == nBundles_temp - 1) {
		tasks_left =
        tasks_packed - (nBundles_temp - 1) * tasksperbundle;
	  }
	  // Will launch a 2d grid of GPU thread blocks (number of tasks is
	  // the y dimension and max_parts is the x dimension
	  int numBlocks_y = tasks_left;
	  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
	  int bundle_first_task = pack_vars->bundle_first_task_list[bid];
	  const char *loop_type = "density";
//	  fprintf(stderr, "Launching kernel with %i tasks leftovers %i\n",
//			  tasks_packed, pack_vars->launch_leftovers);
	  // Launch the kernel
	  launch_density_aos(
		  d_parts_aos, pack_vars->d_task_first_part, pack_vars->d_task_last_part, d_a, d_H, loop_type,
		  stream[bid], BLOCK_SIZE, tasks_packed, tasksperbundle,
		  numBlocks_x, numBlocks_y, bundle_first_task,
		  max_parts, pack_vars->d_cellx, pack_vars->d_celly, pack_vars->d_cellz);
//#ifdef CUDA_DEBUG
//	  cu_error = cudaPeekAtLastError(); // Get error code
//	  if (cu_error != cudaSuccess) {
//		fprintf(stderr,
//				"CUDA error with self density kernel launch: %s cpuid id is: %i\n ",
//				cudaGetErrorString(cu_error), r->cpuid);
//		exit(0);
//	  }
//#endif
	  cudaMemcpyAsync(&parts_aos[first_part_tmp], &d_parts_aos[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos), cudaMemcpyDeviceToHost, stream[bid]);

//#ifdef CUDA_DEBUG
//	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
//										// Get error code
//	  if (cu_error != cudaSuccess) {
//		fprintf(stderr,
//				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
//				cudaGetErrorString(cu_error), r->cpuid);
//		error("Something's up with your cuda code");
//	  }
//#endif
	}/*End of looping over bundles to launch in streams*/
	/* Make sure all the kernels and copies back are finished */
	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
    /*Time unpacking*/
	clock_gettime(CLOCK_REALTIME, &tp0);
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int tid = 0; tid < tasks_packed; tid++) {

	  struct cell *cii = pack_vars->cell_list[tid];
	  struct task *tii = pack_vars->task_list[tid];

//              struct cell *cii = ci_list_self_dens[tid];
//              struct task *tii = task_list_self_dens[tid];

	  while(cell_locktree(cii)) {
		;  /* spin until we acquire the lock */
	  }
	  /* Do the copy */
	 runner_doself1_gpu_unpack_neat_aos(r, cii, parts_aos, 0, &pack_length_unpack, tid, pack_vars->count_max_parts, e);

	  /* Record things for debugging */
	  cii->gpu_done++;

	  /* Release the lock */
	  cell_unlocktree(cii);

	  /*schedule my dependencies (Only unpacks really)*/
	  enqueue_dependencies(s, tii);
	  /*Signal sleeping runners*/
	  signal_sleeping_runners(s, tii);

	  tii->gpu_done = 1;
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
	/*Time end of unpacking*/
	clock_gettime(CLOCK_REALTIME, &tp1);
	*packing_time += (tp1.tv_sec - tp0.tv_sec) +
	(tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
} /*End of GPU work Self*/

void runner_doself1_launch_f4(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv, struct part_aos_f4_send *d_parts_send,
		struct part_aos_f4_recv *d_parts_recv, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time, double *unpack_time, int2 * d_task_first_part_self_dens_f4,
		int devId, int2 * task_first_part_f4, int2 * d_task_first_part_f4, cudaEvent_t * self_end){

	struct timespec t0, t1, tp0, tp1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;

	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero tasks packed but somehow got into GPU loop");
//	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[tasks_packed - 1];
	  pack_vars->bundle_first_part[nBundles_temp] = task_first_part_f4[tasks_packed - 1].x;
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;
//    clock_gettime(CLOCK_REALTIME, &t0hmemcpy);
    /*Copy arrays containing first and last part for each task to GPU*/
//    cudaMemcpy(pack_vars->d_task_first_part, pack_vars->task_first_part,
//               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
//    cudaMemcpy(pack_vars->d_task_last_part, pack_vars->task_last_part,
//               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
//    cudaMemPrefetchAsync(d_task_first_part_self_dens_f4, tasks_packed * sizeof(int2), devId, NULL);
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
	  for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
		   tid++) {
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
	  const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp;
//	  clock_gettime(CLOCK_REALTIME, &t0hmemcpy);
//      cudaMemPrefetchAsync(&d_task_first_part_self_dens_f4[first_task], (last_task - first_task) * sizeof(int2),
//    		  devId, stream[bid]);
      cudaMemcpyAsync(&d_task_first_part_f4[first_task], &task_first_part_f4[first_task],
    		  (last_task + 1  - first_task) * sizeof(int2), cudaMemcpyHostToDevice, stream[bid]);
//	  cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
//	  if (cu_error != cudaSuccess) {
//		fprintf(
//			stderr,
//			"CUDA error in density self host 2 device memcpy: %s cpuid id is: %i\n ",
//			cudaGetErrorString(cu_error), r->cpuid);
//		exit(0);
//	  }
//      clock_gettime(CLOCK_REALTIME, &t1hmemcpy);
//      *hmemcpy_time += (t1hmemcpy.tv_sec - t0hmemcpy.tv_sec) +
//  			(t1hmemcpy.tv_nsec - t0hmemcpy.tv_nsec) / 1000000000.0;
	  cudaMemcpyAsync(&d_parts_send[first_part_tmp], &parts_send[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos_f4_send), cudaMemcpyHostToDevice, stream[bid]);

//#ifdef CUDA_DEBUG
//	  cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
//										// Get error code
//	  if (cu_error != cudaSuccess) {
//		fprintf(
//			stderr,
//			"CUDA error in density self host 2 device memcpy: %s cpuid id is: %i\n ",
//			cudaGetErrorString(cu_error), r->cpuid);
//		exit(0);
//	  }
//#endif
	  const int tasksperbundle = pack_vars->tasksperbundle;
	  int tasks_left = tasksperbundle;
	  if (bid == nBundles_temp - 1) {
		tasks_left =
            tasks_packed - (nBundles_temp - 1) * tasksperbundle;
	  }
	  // Will launch a 2d grid of GPU thread blocks (number of tasks is
	  // the y dimension and max_parts is the x dimension
	  int numBlocks_y = tasks_left;
	  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
	  int bundle_first_task = pack_vars->bundle_first_task_list[bid];
//	  const char *loop_type = "density";
//	  struct first_part first_parts;
//	  for(int i = 0; i < numBlocks_y; i++) first_parts.list[i] = pack_vars->task_first_part[i];
//	  fprintf(stderr, "Launching kernel with %i tasks leftovers %i\n",
//			  tasks_packed, pack_vars->launch_leftovers);
	  // Launch the kernel
	  launch_density_aos_f4(
			  d_parts_send, d_parts_recv, d_a, d_H, stream[bid], numBlocks_x, numBlocks_y, bundle_first_task,
		  d_task_first_part_f4);
//#ifdef CUDA_DEBUG
//	  cu_error = cudaPeekAtLastError(); // Get error code
//	  if (cu_error != cudaSuccess) {
//		fprintf(stderr,
//				"CUDA error with self density kernel launch: %s cpuid id is: %i\n ",
//				cudaGetErrorString(cu_error), r->cpuid);
//		exit(0);
//	  }
//#endif
	  cudaMemcpyAsync(&parts_recv[first_part_tmp], &d_parts_recv[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos_f4_recv), cudaMemcpyDeviceToHost, stream[bid]);
      cudaEventRecord(self_end[bid], stream[bid]);
//#ifdef CUDA_DEBUG
//	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
//										// Get error code
//	  if (cu_error != cudaSuccess) {
//		fprintf(stderr,
//				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
//				cudaGetErrorString(cu_error), r->cpuid);
//		error("Something's up with your cuda code");
//	  }
//#endif
	}/*End of looping over bundles to launch in streams*/
	/* Make sure all the kernels and copies back are finished */
//	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int bid = 0; bid < nBundles_temp; bid++){

		clock_gettime(CLOCK_REALTIME, &t0);

//		cudaStreamSynchronize(stream[bid]);
		cudaEventSynchronize(self_end[bid]);

		clock_gettime(CLOCK_REALTIME, &t1);
		*gpu_time += (t1.tv_sec - t0.tv_sec) +
				(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	    /*Time unpacking*/
//		clock_gettime(CLOCK_REALTIME, &tp0);

		for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

		  if(tid < tasks_packed){
			  struct cell *cii = pack_vars->cell_list[tid];
			  struct task *tii = pack_vars->task_list[tid];

		//              struct cell *cii = ci_list_self_dens[tid];
		//              struct task *tii = task_list_self_dens[tid];

				clock_gettime(CLOCK_REALTIME, &tp0);

//			  clock_gettime(CLOCK_REALTIME, &t0hmemcpy);
			  while(cell_locktree(cii)) {
				;  /* spin until we acquire the lock */
			  }
//			  clock_gettime(CLOCK_REALTIME, &t1hmemcpy);
//				*hmemcpy_time += (t1hmemcpy.tv_sec - t0hmemcpy.tv_sec) +
//				(t1hmemcpy.tv_nsec - t0hmemcpy.tv_nsec) / 1000000000.0;
			  /* Do the copy */
			 runner_doself1_gpu_unpack_neat_aos_f4(r, cii, parts_recv, 0, &pack_length_unpack, tid, pack_vars->count_max_parts, e);

			  /* Record things for debugging */
			  cii->gpu_done++;
			  /*Time end of unpacking*/
			  clock_gettime(CLOCK_REALTIME, &tp1);
			  *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
			    (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;

			  /* Release the lock */
			  cell_unlocktree(cii);

			  /*schedule my dependencies (Only unpacks really)*/
			  enqueue_dependencies(s, tii);
			  /*Signal sleeping runners*/
			  signal_sleeping_runners(s, tii);

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

} /*End of GPU work Self*/

void runner_doself1_launch_g(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_g *parts_aos, struct part_aos_g *d_parts_aos, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time){


	struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;

	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero tasks packed but somehow got into GPU loop");
	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[tasks_packed - 1];
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;
    /*Copy arrays containing first and last part for each task to GPU*/
    cudaMemcpy(pack_vars->d_task_first_part, pack_vars->task_first_part,
               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_task_last_part, pack_vars->task_last_part,
               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);

    /*Copy cell shifts to device*/
    cudaMemcpy(pack_vars->d_cellx, pack_vars->cellx,
               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_celly, pack_vars->celly,
               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_cellz, pack_vars->cellz,
               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);

	/* Launch the copies for each bundle and run the GPU kernel */
	/*We don't go into this loop if tasks_left_self == 1 as
	 nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
	int max_parts = 0;
	for (int bid = 0; bid < nBundles_temp; bid++) {

	  max_parts = 0;
	  int parts_in_bundle = 0;
	  for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
		   tid++) {
		if (tid < tasks_packed) {
		  /*Get an estimate for the max number of parts per cell in the bundle.
		   *  Used for determining the number of GPU CUDA blocks*/
		  int count = pack_vars->task_last_part[tid] - pack_vars->task_first_part[tid];
		  parts_in_bundle += count;
		  max_parts = max(max_parts, count);
		}
	  }

	  const int first_part_tmp = pack_vars->bundle_first_part[bid];
	  const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp;

	  cudaMemcpyAsync(&d_parts_aos[first_part_tmp], &parts_aos[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos_g), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
	  cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(
			stderr,
			"CUDA error in density self host 2 device memcpy: %s cpuid id is: %i\n ",
			cudaGetErrorString(cu_error), r->cpuid);
		exit(0);
	  }
#endif
	  const int tasksperbundle = pack_vars->tasksperbundle;
	  int tasks_left = tasksperbundle;
	  if (bid == nBundles_temp - 1) {
		tasks_left =
  tasks_packed - (nBundles_temp - 1) * tasksperbundle;
	  }
	  // Will launch a 2d grid of GPU thread blocks (number of tasks is
	  // the y dimension and max_parts is the x dimension
	  int numBlocks_y = tasks_left;
	  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
	  int bundle_first_task = pack_vars->bundle_first_task_list[bid];
	  const char *loop_type = "density";
	  // Launch the kernel
	  launch_gradient_aos(
		  d_parts_aos, pack_vars->d_task_first_part, pack_vars->d_task_last_part, d_a, d_H, loop_type,
		  stream[bid], BLOCK_SIZE, tasks_packed, tasksperbundle,
		  numBlocks_x, numBlocks_y, bundle_first_task,
		  max_parts, pack_vars->d_cellx, pack_vars->d_celly, pack_vars->d_cellz);
#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density kernel launch: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		exit(0);
	  }
#endif
	  cudaMemcpyAsync(&parts_aos[first_part_tmp], &d_parts_aos[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos_g), cudaMemcpyDeviceToHost, stream[bid]);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
    /*Time unpacking*/
	clock_gettime(CLOCK_REALTIME, &t0);
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int tid = 0; tid < tasks_packed; tid++) {

	  struct cell *cii = pack_vars->cell_list[tid];
	  struct task *tii = pack_vars->task_list[tid];

//              struct cell *cii = ci_list_self_dens[tid];
//              struct task *tii = task_list_self_dens[tid];

	  while(cell_locktree(cii)) {
		;  /* spin until we acquire the lock */
	  }
	  /* Do the copy */
	 runner_doself1_gpu_unpack_neat_aos_g(r, cii, parts_aos, 0, &pack_length_unpack, tid, pack_vars->count_max_parts, e);

	  /* Record things for debugging */
	  cii->gpu_done_g++;

	  /* Release the lock */
	  cell_unlocktree(cii);

	  /*schedule my dependencies (Only unpacks really)*/
	  enqueue_dependencies(s, tii);
	  /*Signal sleeping runners*/
	  signal_sleeping_runners(s, tii);

	  tii->gpu_done = 1;
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
	/*Time end of unpacking*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*packing_time += (t1.tv_sec - t0.tv_sec) +
	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  } /*End of GPU work Self Gradient*/

void runner_doself1_launch_f4_g(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_f4_g_send *parts_send, struct part_aos_f4_g_recv *parts_recv, struct part_aos_f4_g_send *d_parts_send,
		struct part_aos_f4_g_recv *d_parts_recv, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time, int2 *task_first_part_f4, int2 *d_task_first_part_f4,
		cudaEvent_t *self_end, double *unpack_time){


	struct timespec t0, t1, tp0, tp1;
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;

	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
//	  if(tasks_packed == 0) error("zero tasks packed but somehow got into GPU loop");
	  pack_vars->bundle_first_part[nBundles_temp] = task_first_part_f4[tasks_packed - 1].x;
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;

	/* Launch the copies for each bundle and run the GPU kernel */
	/*We don't go into this loop if tasks_left_self == 1 as
	 nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
	int max_parts;
	for (int bid = 0; bid < nBundles_temp; bid++) {

	  max_parts = 0;
	  int parts_in_bundle = 0;
	  const int first_task = bid * bundle_size;
	  int last_task = (bid + 1) * bundle_size;
	  for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
		   tid++) {
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
	  const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp;

	  cudaMemcpyAsync(&d_task_first_part_f4[first_task], &task_first_part_f4[first_task],
			  (last_task + 1  - first_task) * sizeof(int2), cudaMemcpyHostToDevice, stream[bid]);

	  cudaMemcpyAsync(&d_parts_send[first_part_tmp], &parts_send[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos_f4_g_send), cudaMemcpyHostToDevice, stream[bid]);
//	  fprintf(stderr, "bid %i first_part %i nparts %i\n", bid, first_part_tmp, bundle_n_parts);

#ifdef CUDA_DEBUG
	  cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(
			stderr,
			"CUDA error in gradient self host 2 device memcpy: %s cpuid id is: %i\n ",
			cudaGetErrorString(cu_error), r->cpuid);
		exit(0);
	  }
#endif
	  const int tasksperbundle = pack_vars->tasksperbundle;
	  int tasks_left = tasksperbundle;
	  if (bid == nBundles_temp - 1) {
		tasks_left =
            tasks_packed - (nBundles_temp - 1) * tasksperbundle;
	  }
	  // Will launch a 2d grid of GPU thread blocks (number of tasks is
	  // the y dimension and max_parts is the x dimension
	  int numBlocks_y = tasks_left;
	  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
	  int bundle_first_task = pack_vars->bundle_first_task_list[bid];
//	  const char *loop_type = "density";
	  // Launch the kernel
	  launch_gradient_aos_f4(
		  d_parts_send, d_parts_recv, d_a, d_H, stream[bid], numBlocks_x, numBlocks_y, bundle_first_task,
		  d_task_first_part_f4);
#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self gradient kernel launch: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		exit(0);
	  }
#endif
	  cudaMemcpyAsync(&parts_recv[first_part_tmp], &d_parts_recv[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos_f4_g_recv), cudaMemcpyDeviceToHost, stream[bid]);
      cudaEventRecord(self_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self gradient D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/
//	exit(0);
	/* Make sure all the kernels and copies back are finished */
//	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int bid = 0; bid < nBundles_temp; bid++){

		clock_gettime(CLOCK_REALTIME, &t0);

//		cudaStreamSynchronize(stream[bid]);
		cudaEventSynchronize(self_end[bid]);

		clock_gettime(CLOCK_REALTIME, &t1);
		*gpu_time += (t1.tv_sec - t0.tv_sec) +
				(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	    /*Time unpacking*/
//		clock_gettime(CLOCK_REALTIME, &tp0);

		for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

			  if(tid < tasks_packed){

			  struct cell *cii = pack_vars->cell_list[tid];
			  struct task *tii = pack_vars->task_list[tid];

		//              struct cell *cii = ci_list_self_dens[tid];
		//              struct task *tii = task_list_self_dens[tid];

			  while(cell_locktree(cii)) {
				;  /* spin until we acquire the lock */
			  }
			    /*Time unpacking*/
				clock_gettime(CLOCK_REALTIME, &tp0);
			  /* Do the copy */
			 runner_doself1_gpu_unpack_neat_aos_f4_g(r, cii, parts_recv, 0, &pack_length_unpack, tid, pack_vars->count_max_parts, e);
				/*Time end of unpacking*/
				clock_gettime(CLOCK_REALTIME, &tp1);
				*unpack_time += (tp1.tv_sec - tp0.tv_sec) +
				(tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;

			  /* Record things for debugging */
			  cii->gpu_done_g++;

			  /* Release the lock */
			  cell_unlocktree(cii);

			  /*schedule my dependencies (Only unpacks really)*/
			  enqueue_dependencies(s, tii);
			  /*Signal sleeping runners*/
			  signal_sleeping_runners(s, tii);

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

} /*End of GPU work Self Gradient*/


void runner_doself1_launch_f(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_f *parts_aos, struct part_aos_f *d_parts_aos, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time){

	struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;

	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero tasks packed but somehow got into GPU loop");
	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[tasks_packed - 1];
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;
    /*Copy arrays containing first and last part for each task to GPU*/
    cudaMemcpy(pack_vars->d_task_first_part, pack_vars->task_first_part,
               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_task_last_part, pack_vars->task_last_part,
               tasks_packed * sizeof(int), cudaMemcpyHostToDevice);

    /*Copy cell shifts to device*/
    cudaMemcpy(pack_vars->d_cellx, pack_vars->cellx,
               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_celly, pack_vars->celly,
               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_cellz, pack_vars->cellz,
               tasks_packed * sizeof(double), cudaMemcpyHostToDevice);

	/* Launch the copies for each bundle and run the GPU kernel */
	/*We don't go into this loop if tasks_left_self == 1 as
	 nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
	int max_parts = 0;
	for (int bid = 0; bid < nBundles_temp; bid++) {

	  max_parts = 0;
	  int parts_in_bundle = 0;
	  for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
		   tid++) {
		if (tid < tasks_packed) {
		  /*Get an estimate for the max number of parts per cell in the bundle.
		   *  Used for determining the number of GPU CUDA blocks*/
		  int count = pack_vars->task_last_part[tid] - pack_vars->task_first_part[tid];
		  parts_in_bundle += count;
		  max_parts = max(max_parts, count);
		}
	  }

	  const int first_part_tmp = pack_vars->bundle_first_part[bid];
	  const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp;

	  cudaMemcpyAsync(&d_parts_aos[first_part_tmp], &parts_aos[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos_f), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
	  cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(
			stderr,
			"CUDA error in density self host 2 device memcpy: %s cpuid id is: %i\n ",
			cudaGetErrorString(cu_error), r->cpuid);
		exit(0);
	  }
#endif
	  const int tasksperbundle = pack_vars->tasksperbundle;
	  int tasks_left = tasksperbundle;
	  if (bid == nBundles_temp - 1) {
		tasks_left =
  tasks_packed - (nBundles_temp - 1) * tasksperbundle;
	  }
	  // Will launch a 2d grid of GPU thread blocks (number of tasks is
	  // the y dimension and max_parts is the x dimension
	  int numBlocks_y = tasks_left;
	  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
	  int bundle_first_task = pack_vars->bundle_first_task_list[bid];
	  const char *loop_type = "density";
	  // Launch the kernel
	  launch_force_aos(
		  d_parts_aos, pack_vars->d_task_first_part, pack_vars->d_task_last_part, d_a, d_H, loop_type,
		  stream[bid], BLOCK_SIZE, tasks_packed, tasksperbundle,
		  numBlocks_x, numBlocks_y, bundle_first_task,
		  max_parts, pack_vars->d_cellx, pack_vars->d_celly, pack_vars->d_cellz);
#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self force kernel launch: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		exit(0);
	  }
#endif
	  cudaMemcpyAsync(&parts_aos[first_part_tmp], &d_parts_aos[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos_f), cudaMemcpyDeviceToHost, stream[bid]);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self firce D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
    /*Time unpacking*/
	clock_gettime(CLOCK_REALTIME, &t0);
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int tid = 0; tid < tasks_packed; tid++) {

	  struct cell *cii = pack_vars->cell_list[tid];
	  struct task *tii = pack_vars->task_list[tid];

//              struct cell *cii = ci_list_self_dens[tid];
//              struct task *tii = task_list_self_dens[tid];

	  while(cell_locktree(cii)) {
		;  /* spin until we acquire the lock */
	  }
	  /* Do the copy */
	 runner_doself1_gpu_unpack_neat_aos_f(r, cii, parts_aos, 0, &pack_length_unpack, tid, pack_vars->count_max_parts, e);

	  /* Record things for debugging */
	  cii->gpu_done_f++;

	  /* Release the lock */
	  cell_unlocktree(cii);

	  /*schedule my dependencies (Only unpacks really)*/
	  enqueue_dependencies(s, tii);
	  /*Signal sleeping runners*/
	  signal_sleeping_runners(s, tii);

	  tii->gpu_done = 1;
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
	/*Time end of unpacking*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*packing_time += (t1.tv_sec - t0.tv_sec) +
	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  } /*End of GPU work Self Gradient*/

void runner_doself1_launch_f4_f(struct runner *r, struct scheduler *s, struct pack_vars_self *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_f4_f_send *parts_send, struct part_aos_f4_f_recv *parts_recv,
		struct part_aos_f4_f_send *d_parts_send, struct part_aos_f4_f_recv *d_parts_recv, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time, int2 *task_first_part_f4_f, int2 *d_task_first_part_f4_f, cudaEvent_t * self_end,
		double *unpack_time){

	struct timespec t0, t1, tp0, tp1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;

	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero tasks packed but somehow got into GPU loop");
	  pack_vars->bundle_first_part[nBundles_temp] = task_first_part_f4_f[tasks_packed - 1].x;
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;
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
	  for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
		   tid++) {
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
	  const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp;
	  cudaMemcpyAsync(&d_task_first_part_f4_f[first_task], &task_first_part_f4_f[first_task],
			  (last_task + 1  - first_task) * sizeof(int2), cudaMemcpyHostToDevice, stream[bid]);

	  cudaMemcpyAsync(&d_parts_send[first_part_tmp], &parts_send[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos_f4_f_send), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
	  cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(
			stderr,
			"CUDA error in density self host 2 device memcpy: %s cpuid id is: %i\n ",
			cudaGetErrorString(cu_error), r->cpuid);
		exit(0);
	  }
#endif
	  const int tasksperbundle = pack_vars->tasksperbundle;
	  int tasks_left = tasksperbundle;
	  if (bid == nBundles_temp - 1) {
		tasks_left =
             tasks_packed - (nBundles_temp - 1) * tasksperbundle;
	  }
	  // Will launch a 2d grid of GPU thread blocks (number of tasks is
	  // the y dimension and max_parts is the x dimension
	  int numBlocks_y = tasks_left;
	  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
	  int bundle_first_task = pack_vars->bundle_first_task_list[bid];
	  // Launch the kernel
	  launch_force_aos_f4(
		  d_parts_send, d_parts_recv, d_a, d_H, stream[bid], numBlocks_x, numBlocks_y, bundle_first_task,
		  d_task_first_part_f4_f);
#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self force kernel launch: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		exit(0);
	  }
#endif
	  cudaMemcpyAsync(&parts_recv[first_part_tmp], &d_parts_recv[first_part_tmp],
			bundle_n_parts * sizeof(struct part_aos_f4_f_recv), cudaMemcpyDeviceToHost, stream[bid]);
      cudaEventRecord(self_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self firce D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
//	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int bid = 0; bid < nBundles_temp; bid++){

		clock_gettime(CLOCK_REALTIME, &t0);

//		cudaStreamSynchronize(stream[bid]);
		cudaEventSynchronize(self_end[bid]);

		clock_gettime(CLOCK_REALTIME, &t1);
		*gpu_time += (t1.tv_sec - t0.tv_sec) +
				(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	    /*Time unpacking*/
//		clock_gettime(CLOCK_REALTIME, &tp0);

		for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

		  if(tid < tasks_packed){
			  struct cell *cii = pack_vars->cell_list[tid];
			  struct task *tii = pack_vars->task_list[tid];

		//              struct cell *cii = ci_list_self_dens[tid];
		//              struct task *tii = task_list_self_dens[tid];

			  while(cell_locktree(cii)) {
				;  /* spin until we acquire the lock */
			  }
			 clock_gettime(CLOCK_REALTIME, &tp0);
			  /* Do the copy */
			 runner_doself1_gpu_unpack_neat_aos_f4_f(r, cii, parts_recv, 0, &pack_length_unpack, tid, pack_vars->count_max_parts, e);

			  /* Record things for debugging */
			  cii->gpu_done_f++;
			  clock_gettime(CLOCK_REALTIME, &tp1);
			  *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
			  (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
			  /* Release the lock */
			  cell_unlocktree(cii);

			  /*schedule my dependencies (Only unpacks really)*/
			  enqueue_dependencies(s, tii);
			  /*Signal sleeping runners*/
			  signal_sleeping_runners(s, tii);

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
  } /*End of GPU work Self Gradient*/

void runner_dopair1_launch(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars, struct cell *ci,
		struct task *t, struct part_aos *parts_aos, struct part_aos *d_parts_aos, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time){

	struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;
	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/*tasks-packed needs decrementing before calculating packed_tmp as it was incremented in runner_dopair1_pack*/
	const int packed_tmp = 2 * (tasks_packed - 1);

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero pair tasks packed but somehow got into GPU loop");
	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[packed_tmp - 2];
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;
    /*Copy arrays containing first and last part for each task to GPU*/
    cudaMemcpy(pack_vars->d_task_first_part, pack_vars->task_first_part,
               2 * tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_task_last_part, pack_vars->task_last_part,
               2 * tasks_packed * sizeof(int), cudaMemcpyHostToDevice);

    /*Copy cell shifts to device*/
//    cudaMemcpy(pack_vars->d_shiftx, pack_vars->shiftx,
//    		2 * tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
//    cudaMemcpy(pack_vars->d_shifty, pack_vars->shifty,
//    		2 * tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
//    cudaMemcpy(pack_vars->d_shiftz, pack_vars->shiftz,
//    		2 * tasks_packed * sizeof(double), cudaMemcpyHostToDevice);

	/* Launch the copies for each bundle and run the GPU kernel */
	/*We don't go into this loop if tasks_left_self == 1 as
	 nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
	for (int bid = 0; bid < nBundles_temp; bid++) {

      int max_parts_i = 0;
      int max_parts_j = 0;
      int parts_in_bundle_ci = 0;
      int parts_in_bundle_cj = 0;
      for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
           tid++) {
        if (tid < tasks_packed) {
          /*Get an estimate for the max number of parts per cell in each bundle.
           *  Used for determining the number of GPU CUDA blocks*/
      	const int tid_tmp = 2 * tid;
          int count_i = pack_vars->task_last_part[tid_tmp]
															  - pack_vars->task_first_part[tid_tmp];
          parts_in_bundle_ci += count_i;
          max_parts_i = max(max_parts_i, count_i);
          int count_j = pack_vars->task_last_part[tid_tmp + 1]
															  - pack_vars->task_first_part[tid_tmp + 1];
          parts_in_bundle_cj += count_j;
          max_parts_j = max(max_parts_j, count_j);
        }
      }
      const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
      const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;
      cudaMemcpyAsync(&d_parts_aos[first_part_tmp_i], &parts_aos[first_part_tmp_i],
        bundle_n_parts * sizeof(struct part_aos), cudaMemcpyHostToDevice, stream[bid]);

//#ifdef CUDA_DEBUG
//      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
//                                                // Get error code
//      if (cu_error != cudaSuccess) {
//        fprintf(stderr,
//        "CUDA error with pair density H2D async  memcpy ci: %s cpuid id is: %i\n ",
//        cudaGetErrorString(cu_error), r->cpuid);
//        error("Something's up with your cuda code");
//      }
//#endif

	  const int tasksperbundle = pack_vars->tasksperbundle;
	  /* LAUNCH THE GPU KERNELS for ci & cj */
      int tid = 0;
      int offset = bid * tasksperbundle;
      int tasks_left = tasksperbundle;
      if (bid == nBundles_temp - 1) {
        tasks_left =
        		tasks_packed - (nBundles_temp - 1) * tasksperbundle;
      }

      // Setup 2d grid of GPU thread blocks for ci (number of tasks is
      // the y dimension and max_parts is the x dimension
      int numBlocks_y = tasks_left;
      int bundle_first_task = pack_vars->bundle_first_task_list[bid];
      const char *loop_type = "density";

        /* Launch the kernel for ci using data for ci and cj */
        runner_dopairci_branch_density_gpu_aos(d_parts_aos,
        pack_vars->d_task_first_part, pack_vars->d_task_last_part,
	    d_a, d_H, loop_type, stream[bid], bid, BLOCK_SIZE,
	    tasks_packed, tasksperbundle, max_parts_i, max_parts_j, numBlocks_y, tid,
	    offset, bundle_first_task, pack_vars->d_shiftx, pack_vars->d_shifty, pack_vars->d_shiftz);

        /* Launch the kernel for ci using data for ci and cj */
        runner_dopaircj_branch_density_gpu_aos(d_parts_aos,
        pack_vars->d_task_first_part, pack_vars->d_task_last_part,
	    d_a, d_H, loop_type, stream[bid], bid, BLOCK_SIZE,
	    tasks_packed, tasksperbundle, max_parts_i, max_parts_j, numBlocks_y, tid,
	    offset, bundle_first_task, pack_vars->d_shiftx, pack_vars->d_shifty, pack_vars->d_shiftz);

//#ifdef CUDA_DEBUG
//	  cu_error = cudaPeekAtLastError(); // Get error code
//	  if (cu_error != cudaSuccess) {
//		fprintf(stderr,
//				"CUDA error with self density kernel launch: %s cpuid id is: %i\n ",
//				cudaGetErrorString(cu_error), r->cpuid);
//		exit(0);
//	  }
//#endif

        // Copy results back to CPU BUFFERS
      cudaMemcpyAsync(&parts_aos[first_part_tmp_i], &d_parts_aos[first_part_tmp_i],
    	  	  bundle_n_parts * sizeof(struct part_aos), cudaMemcpyDeviceToHost, stream[bid]);

//#ifdef CUDA_DEBUG
//	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
//										// Get error code
//	  if (cu_error != cudaSuccess) {
//		fprintf(stderr,
//				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
//				cudaGetErrorString(cu_error), r->cpuid);
//		error("Something's up with your cuda code");
//	  }
//#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
    /*Time unpacking*/
	clock_gettime(CLOCK_REALTIME, &t0);
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int tid = 0; tid < tasks_packed; tid++) {

	  /*grab cell and task pointers*/
	  struct cell *cii = pack_vars->ci_list[tid];
      struct cell *cjj = pack_vars->cj_list[tid];
	  struct task *tii = pack_vars->task_list[tid];

	  /*Let's lock ci*/
	  while(cell_locktree(cii)) {
		;  /* spin until we acquire the lock */
	  }
	  /*Let's lock cj*/
	  while(cell_locktree(cjj)) {
	    ;  /* spin until we acquire the lock */
	  }
	  /* Do the copy */
	  runner_do_ci_cj_gpu_unpack_neat_aos(r, cii, cjj, parts_aos, 0,
		  	  &pack_length_unpack, tid, 2 * pack_vars->count_max_parts, e);

	  /* Record things for debugging */
	  cii->gpu_done_pair++;
	  cjj->gpu_done_pair++;

	  tii->gpu_done = 1;

      /*schedule my dependencies (Only unpacks really)*/
      enqueue_dependencies(s, tii);
      /*Signal sleeping runners*/
      signal_sleeping_runners(s, tii);

	  /* Release the locks */
	  cell_unlocktree(cii);
	  /* Release the locks */
	  cell_unlocktree(cjj);
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
	/*Time end of unpacking*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*packing_time += (t1.tv_sec - t0.tv_sec) +
	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  } /*End of GPU work*/

void runner_dopair1_launch_f4(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars,
		struct task *t, struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv, struct part_aos_f4_send *d_parts_send,
		struct part_aos_f4_recv *d_parts_recv, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time, double *unpack_time, int4 *fparti_fpartj_lparti_lpartj_dens,
		int4 *d_fparti_fpartj_lparti_lpartj_dens, cudaEvent_t * pair_end){

	struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;
	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero pair tasks packed but somehow got into GPU loop");
//	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[packed_tmp - 2];
	  pack_vars->bundle_first_part[nBundles_temp] = fparti_fpartj_lparti_lpartj_dens[tasks_packed - 1].x;
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;

	/* Launch the copies for each bundle and run the GPU kernel */
	/*We don't go into this loop if tasks_left_self == 1 as
	 nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
	for (int bid = 0; bid < nBundles_temp; bid++) {

      int max_parts_i = 0;
      int max_parts_j = 0;
      int parts_in_bundle_ci = 0;
      int parts_in_bundle_cj = 0;
      const int first_task = bid * pack_vars->bundle_size;
	  int last_task = (bid + 1) * bundle_size;
      for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
           tid++) {
        if (tid < tasks_packed) {
          /*Get an estimate for the max number of parts per cell in each bundle.
           *  Used for determining the number of GPU CUDA blocks*/
          int count_i = fparti_fpartj_lparti_lpartj_dens[tid].z
															  - fparti_fpartj_lparti_lpartj_dens[tid].x;
          parts_in_bundle_ci += count_i;
          max_parts_i = max(max_parts_i, count_i);
          int count_j = fparti_fpartj_lparti_lpartj_dens[tid].w
															  - fparti_fpartj_lparti_lpartj_dens[tid].y;
          parts_in_bundle_cj += count_j;
          max_parts_j = max(max_parts_j, count_j);

		  last_task = tid;
        }
      }
      const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
      const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;

      cudaMemcpyAsync(&d_fparti_fpartj_lparti_lpartj_dens[first_task], &fparti_fpartj_lparti_lpartj_dens[first_task],
    		  (last_task + 1  - first_task) * sizeof(int4), cudaMemcpyHostToDevice, stream[bid]);

      cudaMemcpyAsync(&d_parts_send[first_part_tmp_i], &parts_send[first_part_tmp_i],
        bundle_n_parts * sizeof(struct part_aos_f4_send), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
                                                // Get error code
      if (cu_error != cudaSuccess) {
        fprintf(stderr,
        "CUDA error with pair density H2D async  memcpy ci: %s cpuid id is: %i\n ",
        cudaGetErrorString(cu_error), r->cpuid);
        error("Something's up with your cuda code");
      }
#endif

	  const int tasksperbundle = pack_vars->tasksperbundle;
	  /* LAUNCH THE GPU KERNELS for ci & cj */
      int tasks_left = tasksperbundle;
      if (bid == nBundles_temp - 1) {
        tasks_left =
        		tasks_packed - (nBundles_temp - 1) * tasksperbundle;
      }

      // Setup 2d grid of GPU thread blocks for ci (number of tasks is
      // the y dimension and max_parts is the x dimension
      int numBlocks_y = tasks_left;
      int numBlocks_x = (max_parts_i + BLOCK_SIZE - 1) / BLOCK_SIZE;
      int bundle_first_task = pack_vars->bundle_first_task_list[bid];

      /* Launch the kernel for ci using data for ci and cj */
      runner_dopairci_branch_density_gpu_aos_f4(d_parts_send, d_parts_recv,
		      d_a, d_H, stream[bid], numBlocks_x, numBlocks_y, bundle_first_task, d_fparti_fpartj_lparti_lpartj_dens);

      numBlocks_x = (max_parts_j + BLOCK_SIZE - 1) / BLOCK_SIZE;

      /* Launch the kernel for ci using data for ci and cj */
      runner_dopaircj_branch_density_gpu_aos_f4(d_parts_send, d_parts_recv,
		      d_a, d_H, stream[bid], numBlocks_x, numBlocks_y, bundle_first_task, d_fparti_fpartj_lparti_lpartj_dens);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with pair density kernel launch: %s cpuid id is: %i\n "
				"nbx %i nby %i max_parts_i %i max_parts_j %i\n",
				cudaGetErrorString(cu_error), r->cpuid, numBlocks_x, numBlocks_y, max_parts_i, max_parts_j);
		exit(0);
	  }
#endif

        // Copy results back to CPU BUFFERS
      cudaMemcpyAsync(&parts_recv[first_part_tmp_i], &d_parts_recv[first_part_tmp_i],
    	  	  bundle_n_parts * sizeof(struct part_aos_f4_recv), cudaMemcpyDeviceToHost, stream[bid]);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
    /*Time unpacking*/
	clock_gettime(CLOCK_REALTIME, &t0);
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int tid = 0; tid < tasks_packed; tid++) {

	  /*grab cell and task pointers*/
	  struct cell *cii = pack_vars->ci_list[tid];
      struct cell *cjj = pack_vars->cj_list[tid];
	  struct task *tii = pack_vars->task_list[tid];

	  /*Let's lock ci*/
	  while(cell_locktree(cii)) {
		;  /* spin until we acquire the lock */
	  }
	  /*Let's lock cj*/
	  while(cell_locktree(cjj)) {
	    ;  /* spin until we acquire the lock */
	  }
	  /* Do the copy */
	  runner_do_ci_cj_gpu_unpack_neat_aos_f4(r, cii, cjj, parts_recv, 0,
		  	  &pack_length_unpack, tid, 2 * pack_vars->count_max_parts, e);

	  /* Record things for debugging */
	  cii->gpu_done_pair++;
	  cjj->gpu_done_pair++;

	  tii->gpu_done = 1;

      /*schedule my dependencies (Only unpacks really)*/
      enqueue_dependencies(s, tii);
      /*Signal sleeping runners*/
      signal_sleeping_runners(s, tii);

	  /* Release the locks */
	  cell_unlocktree(cii);
	  /* Release the locks */
	  cell_unlocktree(cjj);
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
	/*Time end of unpacking*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*packing_time += (t1.tv_sec - t0.tv_sec) +
	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  } /*End of GPU work*/
void runner_dopair1_launch_f4_one_memcpy(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars,
		struct task *t, struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv, struct part_aos_f4_send *d_parts_send,
		struct part_aos_f4_recv *d_parts_recv, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time, double *unpack_time, int4 *fparti_fpartj_lparti_lpartj_dens,
		cudaEvent_t * pair_end){

	struct timespec t0, t1, tp0, tp1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;
	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero pair tasks packed but somehow got into GPU loop");
//	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[packed_tmp - 2];
	  pack_vars->bundle_first_part[nBundles_temp] = fparti_fpartj_lparti_lpartj_dens[tasks_packed - 1].x;
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;

	/* Launch the copies for each bundle and run the GPU kernel */
	/*We don't go into this loop if tasks_left_self == 1 as
	 nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
	for (int bid = 0; bid < nBundles_temp; bid++) {

      int max_parts_i = 0;
      int max_parts_j = 0;
      int parts_in_bundle_ci = 0;
      int parts_in_bundle_cj = 0;
      for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
           tid++) {
        if (tid < tasks_packed) {
          /*Get an estimate for the max number of parts per cell in each bundle.
           *  Used for determining the number of GPU CUDA blocks*/
          int count_i = fparti_fpartj_lparti_lpartj_dens[tid].z
															  - fparti_fpartj_lparti_lpartj_dens[tid].x;
          parts_in_bundle_ci += count_i;
          max_parts_i = max(max_parts_i, count_i);
          int count_j = fparti_fpartj_lparti_lpartj_dens[tid].w
															  - fparti_fpartj_lparti_lpartj_dens[tid].y;
          parts_in_bundle_cj += count_j;
          max_parts_j = max(max_parts_j, count_j);
        }
      }
      const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
      const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;

      cudaMemcpyAsync(&d_parts_send[first_part_tmp_i], &parts_send[first_part_tmp_i],
        bundle_n_parts * sizeof(struct part_aos_f4_send), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
                                                // Get error code
      if (cu_error != cudaSuccess) {
        fprintf(stderr,
        "CUDA error with pair density H2D async  memcpy ci: %s cpuid id is: %i\n ",
        cudaGetErrorString(cu_error), r->cpuid);
        error("Something's up with your cuda code");
      }
#endif
	  /* LAUNCH THE GPU KERNELS for ci & cj */
      // Setup 2d grid of GPU thread blocks for ci (number of tasks is
      // the y dimension and max_parts is the x dimension
      int numBlocks_y = 0;//tasks_left;
      int numBlocks_x = (bundle_n_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
      int bundle_part_0 = pack_vars->bundle_first_part[bid];
      /* Launch the kernel for ci using data for ci and cj */
      runner_dopair_branch_density_gpu_aos_f4(d_parts_send, d_parts_recv,
		      d_a, d_H, stream[bid], numBlocks_x, numBlocks_y, bundle_part_0, bundle_n_parts);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with pair density kernel launch: %s cpuid id is: %i\n "
				"nbx %i nby %i max_parts_i %i max_parts_j %i\n",
				cudaGetErrorString(cu_error), r->cpuid, numBlocks_x, numBlocks_y, max_parts_i, max_parts_j);
		exit(0);
	  }
#endif

        // Copy results back to CPU BUFFERS
      cudaMemcpyAsync(&parts_recv[first_part_tmp_i], &d_parts_recv[first_part_tmp_i],
    	  	  bundle_n_parts * sizeof(struct part_aos_f4_recv), cudaMemcpyDeviceToHost, stream[bid]);
      cudaEventRecord(pair_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
//	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int bid = 0; bid < nBundles_temp; bid++){
		/*Time unpacking*/
		clock_gettime(CLOCK_REALTIME, &t0);

//		cudaStreamSynchronize(stream[bid]);
        cudaEventSynchronize(pair_end[bid]);

		clock_gettime(CLOCK_REALTIME, &t1);
		*gpu_time += (t1.tv_sec - t0.tv_sec) +
				(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	    /*Time unpacking*/
//		clock_gettime(CLOCK_REALTIME, &tp0);

		for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

		  if(tid < tasks_packed){
		  clock_gettime(CLOCK_REALTIME, &tp0);
		  /*grab cell and task pointers*/
		  struct cell *cii = pack_vars->ci_list[tid];
		  struct cell *cjj = pack_vars->cj_list[tid];
		  struct task *tii = pack_vars->task_list[tid];

		  /*Let's lock ci*/
		  if(tii->corner_pair == 1)fprintf(stderr, "Corner task\n");
		  while(cell_locktree(cii)) {
			;  /* spin until we acquire the lock */
		  }
		  /*Let's lock cj*/
		  while(cell_locktree(cjj)) {
			;  /* spin until we acquire the lock */
		  }
		  /* Do the copy */
		  runner_do_ci_cj_gpu_unpack_neat_aos_f4(r, cii, cjj, parts_recv, 0,
				  &pack_length_unpack, tid, 2 * pack_vars->count_max_parts, e);

		  /* Record things for debugging */
		  cii->gpu_done_pair++;
		  cjj->gpu_done_pair++;

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
		  signal_sleeping_runners(s, tii);

		  tii->gpu_done = 1;


		  }
		}
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
//	/*Time end of unpacking*/
//	clock_gettime(CLOCK_REALTIME, &t1);
//	*packing_time += (t1.tv_sec - t0.tv_sec) +
//	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  } /*End of GPU work*/

void runner_dopair1_launch_f4_mcpy_Ker_mcpy(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars,
		struct task *t, struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv, struct part_aos_f4_send *d_parts_send,
		struct part_aos_f4_recv *d_parts_recv, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time, double *unpack_time, int4 *fparti_fpartj_lparti_lpartj_dens,
		cudaEvent_t * pair_end){

	struct timespec t0, t1, tp0, tp1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;
	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/*tasks-packed needs decrementing before calculating packed_tmp as it was incremented in runner_dopair1_pack*/
//	const int packed_tmp = 2 * (tasks_packed - 1);

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero pair tasks packed but somehow got into GPU loop");
//	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[packed_tmp - 2];
	  pack_vars->bundle_first_part[nBundles_temp] = fparti_fpartj_lparti_lpartj_dens[tasks_packed - 1].x;
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;

	/* Launch the copies for each bundle and run the GPU kernel */
	/*We don't go into this loop if tasks_left_self == 1 as
	 nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
//	int max_parts = 0;
	for (int bid = 0; bid < nBundles_temp; bid++) {

      int max_parts_i = 0;
      int max_parts_j = 0;
      int parts_in_bundle_ci = 0;
      int parts_in_bundle_cj = 0;
      for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
           tid++) {
        if (tid < tasks_packed) {
          /*Get an estimate for the max number of parts per cell in each bundle.
           *  Used for determining the number of GPU CUDA blocks*/
          int count_i = fparti_fpartj_lparti_lpartj_dens[tid].z
															  - fparti_fpartj_lparti_lpartj_dens[tid].x;
          parts_in_bundle_ci += count_i;
          max_parts_i = max(max_parts_i, count_i);
          int count_j = fparti_fpartj_lparti_lpartj_dens[tid].w
															  - fparti_fpartj_lparti_lpartj_dens[tid].y;
          parts_in_bundle_cj += count_j;
          max_parts_j = max(max_parts_j, count_j);
        }
      }
      const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
      const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;

      cudaMemcpyAsync(&d_parts_send[first_part_tmp_i], &parts_send[first_part_tmp_i],
        bundle_n_parts * sizeof(struct part_aos_f4_send), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
                                                // Get error code
      if (cu_error != cudaSuccess) {
        fprintf(stderr,
        "CUDA error with pair density H2D async  memcpy ci: %s cpuid id is: %i\n ",
        cudaGetErrorString(cu_error), r->cpuid);
        error("Something's up with your cuda code");
      }
#endif
	}
      	for (int bid = 0; bid < nBundles_temp; bid++) {

            int max_parts_i = 0;
            int max_parts_j = 0;
            int parts_in_bundle_ci = 0;
            int parts_in_bundle_cj = 0;
            for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
                 tid++) {
              if (tid < tasks_packed) {
                /*Get an estimate for the max number of parts per cell in each bundle.
                 *  Used for determining the number of GPU CUDA blocks*/
                int count_i = fparti_fpartj_lparti_lpartj_dens[tid].z
      															  - fparti_fpartj_lparti_lpartj_dens[tid].x;
                parts_in_bundle_ci += count_i;
                max_parts_i = max(max_parts_i, count_i);
                int count_j = fparti_fpartj_lparti_lpartj_dens[tid].w
      															  - fparti_fpartj_lparti_lpartj_dens[tid].y;
                parts_in_bundle_cj += count_j;
                max_parts_j = max(max_parts_j, count_j);
              }
            }
            const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
            const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;
//////////////////////////////////
//	  const int tasksperbundle = pack_vars->tasksperbundle;
	  /* LAUNCH THE GPU KERNELS for ci & cj */
      // Setup 2d grid of GPU thread blocks for ci (number of tasks is
      // the y dimension and max_parts is the x dimension
      int numBlocks_y = 0;//tasks_left;
      int numBlocks_x = (bundle_n_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
      int bundle_part_0 = pack_vars->bundle_first_part[bid];
//      int bundle_first_task = pack_vars->bundle_first_task_list[bid];
//              fprintf(stderr, "bundle_part_0 %i bundle_first_task %i\n", bundle_part_0, bundle_first_task);

      /* Launch the kernel for ci using data for ci and cj */
      runner_dopair_branch_density_gpu_aos_f4(d_parts_send, d_parts_recv,
		      d_a, d_H, stream[bid], numBlocks_x, numBlocks_y, bundle_part_0, bundle_n_parts);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with pair density kernel launch: %s cpuid id is: %i\n "
				"nbx %i nby %i max_parts_i %i max_parts_j %i\n",
				cudaGetErrorString(cu_error), r->cpuid, numBlocks_x, numBlocks_y, max_parts_i, max_parts_j);
		exit(0);
	  }
#endif
      	}

    	for (int bid = 0; bid < nBundles_temp; bid++) {

          int max_parts_i = 0;
          int max_parts_j = 0;
          int parts_in_bundle_ci = 0;
          int parts_in_bundle_cj = 0;
          for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
               tid++) {
            if (tid < tasks_packed) {
              /*Get an estimate for the max number of parts per cell in each bundle.
               *  Used for determining the number of GPU CUDA blocks*/
              int count_i = fparti_fpartj_lparti_lpartj_dens[tid].z
    															  - fparti_fpartj_lparti_lpartj_dens[tid].x;
              parts_in_bundle_ci += count_i;
              max_parts_i = max(max_parts_i, count_i);
              int count_j = fparti_fpartj_lparti_lpartj_dens[tid].w
    															  - fparti_fpartj_lparti_lpartj_dens[tid].y;
              parts_in_bundle_cj += count_j;
              max_parts_j = max(max_parts_j, count_j);
            }
          }
          const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
          const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;
///////////////////////////////////////////////////////////////////
        // Copy results back to CPU BUFFERS
      cudaMemcpyAsync(&parts_recv[first_part_tmp_i], &d_parts_recv[first_part_tmp_i],
    	  	  bundle_n_parts * sizeof(struct part_aos_f4_recv), cudaMemcpyDeviceToHost, stream[bid]);
      cudaEventRecord(pair_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
//	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
    /*Time unpacking*/
//	clock_gettime(CLOCK_REALTIME, &t0);
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int bid = 0; bid < nBundles_temp; bid++){
	  clock_gettime(CLOCK_REALTIME, &t0);

//		cudaStreamSynchronize(stream[bid]);
	  cudaEventSynchronize(pair_end[bid]);

	  clock_gettime(CLOCK_REALTIME, &t1);
	  *gpu_time += (t1.tv_sec - t0.tv_sec) +
				(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	    /*Time unpacking*/
//	clock_gettime(CLOCK_REALTIME, &tp0);
//	  int bundle_first_task = pack_vars->bundle_first_task_list[bid];

	  for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

		  if(tid < tasks_packed){

		    /*grab cell and task pointers*/
		    struct cell *cii = pack_vars->ci_list[tid];
		    struct cell *cjj = pack_vars->cj_list[tid];
		    struct task *tii = pack_vars->task_list[tid];

		  /*Let's lock ci*/
		    while(cell_locktree(cii)) {
			  ;  /* spin until we acquire the lock */
		    }
		  /*Let's lock cj*/
		    while(cell_locktree(cjj)) {
			  ;  /* spin until we acquire the lock */
		    }
		  /* Do the copy */
		    /*Time unpacking*/
		    clock_gettime(CLOCK_REALTIME, &tp0);
		    runner_do_ci_cj_gpu_unpack_neat_aos_f4(r, cii, cjj, parts_recv, 0,
				  &pack_length_unpack, tid, 2 * pack_vars->count_max_parts, e);

		    /* Record things for debugging */
		    cii->gpu_done_pair++;
		    cjj->gpu_done_pair++;

		    tii->gpu_done = 1;
		    /*Time end of unpacking*/
		    clock_gettime(CLOCK_REALTIME, &tp1);
		    *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
		    (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
		    /*schedule my dependencies (Only unpacks really)*/
		    enqueue_dependencies(s, tii);
		    /*Signal sleeping runners*/
		    signal_sleeping_runners(s, tii);

		    /* Release the locks */
		    cell_unlocktree(cii);
		    /* Release the locks */
		    cell_unlocktree(cjj);

		  }
	  }
	  /*Time end of unpacking*/
//	  clock_gettime(CLOCK_REALTIME, &tp1);
//	  *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
//	  (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
//	  *packing_time += (tp1.tv_sec - tp0.tv_sec) +
//	  (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
//		*packing_time += (tp1.tv_sec - tp0.tv_sec) +
//		(tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
	/*Time end of unpacking*/
//	clock_gettime(CLOCK_REALTIME, &t1);
//	*packing_time += (t1.tv_sec - t0.tv_sec) +
//	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  } /*End of GPU work*/

void runner_dopair1_launch_g(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_g *parts_aos, struct part_aos_g *d_parts_aos, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time){

	struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;
	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/*tasks-packed needs decrementing before calculating packed_tmp as it was incremented in runner_dopair1_pack*/
	const int packed_tmp = 2 * (tasks_packed - 1);

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero pair tasks packed but somehow got into GPU loop");
	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[packed_tmp - 2];
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;
    /*Copy arrays containing first and last part for each task to GPU*/
    cudaMemcpy(pack_vars->d_task_first_part, pack_vars->task_first_part,
               2 * tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_task_last_part, pack_vars->task_last_part,
               2 * tasks_packed * sizeof(int), cudaMemcpyHostToDevice);

    /*Copy cell shifts to device*/
    cudaMemcpy(pack_vars->d_shiftx, pack_vars->shiftx,
    		2 * tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_shifty, pack_vars->shifty,
    		2 * tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_shiftz, pack_vars->shiftz,
    		2 * tasks_packed * sizeof(double), cudaMemcpyHostToDevice);

	/* Launch the copies for each bundle and run the GPU kernel */
	/*We don't go into this loop if tasks_left_self == 1 as
	 nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
	for (int bid = 0; bid < nBundles_temp; bid++) {

      int max_parts_i = 0;
      int max_parts_j = 0;
      int parts_in_bundle_ci = 0;
      int parts_in_bundle_cj = 0;
      for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
           tid++) {
        if (tid < tasks_packed) {
          /*Get an estimate for the max number of parts per cell in each bundle.
           *  Used for determining the number of GPU CUDA blocks*/
      	const int tid_tmp = 2 * tid;
          int count_i = pack_vars->task_last_part[tid_tmp]
															  - pack_vars->task_first_part[tid_tmp];
          parts_in_bundle_ci += count_i;
          max_parts_i = max(max_parts_i, count_i);
          int count_j = pack_vars->task_last_part[tid_tmp + 1]
															  - pack_vars->task_first_part[tid_tmp + 1];
          parts_in_bundle_cj += count_j;
          max_parts_j = max(max_parts_j, count_j);
        }
      }
      const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
      const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;
      cudaMemcpyAsync(&d_parts_aos[first_part_tmp_i], &parts_aos[first_part_tmp_i],
        bundle_n_parts * sizeof(struct part_aos_g), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
                                                // Get error code
      if (cu_error != cudaSuccess) {
        fprintf(stderr,
        "CUDA error with pair density H2D async  memcpy ci: %s cpuid id is: %i\n ",
        cudaGetErrorString(cu_error), r->cpuid);
        error("Something's up with your cuda code");
      }
#endif

	  const int tasksperbundle = pack_vars->tasksperbundle;
	  /* LAUNCH THE GPU KERNELS for ci & cj */
      int tid = 0;
      int offset = bid * tasksperbundle;
      int tasks_left = tasksperbundle;
      if (bid == nBundles_temp - 1) {
        tasks_left =
        		tasks_packed - (nBundles_temp - 1) * tasksperbundle;
      }

      // Setup 2d grid of GPU thread blocks for ci (number of tasks is
      // the y dimension and max_parts is the x dimension
      int numBlocks_y = tasks_left;
      int bundle_first_task = pack_vars->bundle_first_task_list[bid];
      const char *loop_type = "density";

        /* Launch the kernel for ci using data for ci and cj */
        runner_dopairci_branch_density_gpu_aos_g(d_parts_aos,
        pack_vars->d_task_first_part, pack_vars->d_task_last_part,
	    d_a, d_H, loop_type, stream[bid], bid, BLOCK_SIZE,
	    tasks_packed, tasksperbundle, max_parts_i, max_parts_j, numBlocks_y, tid,
	    offset, bundle_first_task, pack_vars->d_shiftx, pack_vars->d_shifty, pack_vars->d_shiftz);

        /* Launch the kernel for ci using data for ci and cj */
        runner_dopaircj_branch_density_gpu_aos_g(d_parts_aos,
        pack_vars->d_task_first_part, pack_vars->d_task_last_part,
	    d_a, d_H, loop_type, stream[bid], bid, BLOCK_SIZE,
	    tasks_packed, tasksperbundle, max_parts_i, max_parts_j, numBlocks_y, tid,
	    offset, bundle_first_task, pack_vars->d_shiftx, pack_vars->d_shifty, pack_vars->d_shiftz);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density kernel launch: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		exit(0);
	  }
#endif

        // Copy results back to CPU BUFFERS
      cudaMemcpyAsync(&parts_aos[first_part_tmp_i], &d_parts_aos[first_part_tmp_i],
    	  	  bundle_n_parts * sizeof(struct part_aos_g), cudaMemcpyDeviceToHost, stream[bid]);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
    /*Time unpacking*/
	clock_gettime(CLOCK_REALTIME, &t0);
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int tid = 0; tid < tasks_packed; tid++) {

	  /*grab cell and task pointers*/
	  struct cell *cii = pack_vars->ci_list[tid];
      struct cell *cjj = pack_vars->cj_list[tid];
	  struct task *tii = pack_vars->task_list[tid];

	  /*Let's lock ci*/
	  while(cell_locktree(cii)) {
		;  /* spin until we acquire the lock */
	  }
	  /*Let's lock cj*/
	  while(cell_locktree(cjj)) {
	    ;  /* spin until we acquire the lock */
	  }
	  /* Do the copy */
	  runner_do_ci_cj_gpu_unpack_neat_aos_g(r, cii, cjj, parts_aos, 0,
		  	  &pack_length_unpack, tid, 2 * pack_vars->count_max_parts, e);

	  /* Record things for debugging */
	  cii->gpu_done_pair_g++;
	  cjj->gpu_done_pair_g++;

	  tii->gpu_done = 1;

      /*schedule my dependencies (Only unpacks really)*/
      enqueue_dependencies(s, tii);
      /*Signal sleeping runners*/
      signal_sleeping_runners(s, tii);

	  /* Release the locks */
	  cell_unlocktree(cii);
	  /* Release the locks */
	  cell_unlocktree(cjj);
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
	/*Time end of unpacking*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*packing_time += (t1.tv_sec - t0.tv_sec) +
	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
} /*End of GPU work*/

void runner_dopair1_launch_f4_g_one_memcpy(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars,
		struct task *t, struct part_aos_f4_g_send *parts_send, struct part_aos_f4_g_recv *parts_recv, struct part_aos_f4_g_send *d_parts_send,
		struct part_aos_f4_g_recv *d_parts_recv, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time, double *unpack_time, int4 *fparti_fpartj_lparti_lpartj,
		cudaEvent_t * pair_end){

	struct timespec t0, t1, tp0, tp1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;
	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/*tasks-packed needs decrementing before calculating packed_tmp as it was incremented in runner_dopair1_pack*/
//	const int packed_tmp = 2 * (tasks_packed - 1);

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero pair tasks packed but somehow got into GPU loop");
//	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[packed_tmp - 2];
	  pack_vars->bundle_first_part[nBundles_temp] = fparti_fpartj_lparti_lpartj[tasks_packed - 1].x;
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;

	/* Launch the copies for each bundle and run the GPU kernel */
	/*We don't go into this loop if tasks_left_self == 1 as
	 nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
//	int max_parts = 0;
	for (int bid = 0; bid < nBundles_temp; bid++) {

      int max_parts_i = 0;
      int max_parts_j = 0;
      int parts_in_bundle_ci = 0;
      int parts_in_bundle_cj = 0;
      for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
           tid++) {
        if (tid < tasks_packed) {
          /*Get an estimate for the max number of parts per cell in each bundle.
           *  Used for determining the number of GPU CUDA blocks*/
          int count_i = fparti_fpartj_lparti_lpartj[tid].z
															  - fparti_fpartj_lparti_lpartj[tid].x;
          parts_in_bundle_ci += count_i;
          max_parts_i = max(max_parts_i, count_i);
          int count_j = fparti_fpartj_lparti_lpartj[tid].w
															  - fparti_fpartj_lparti_lpartj[tid].y;
          parts_in_bundle_cj += count_j;
          max_parts_j = max(max_parts_j, count_j);
        }
      }
      const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
      const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;

      cudaMemcpyAsync(&d_parts_send[first_part_tmp_i], &parts_send[first_part_tmp_i],
        bundle_n_parts * sizeof(struct part_aos_f4_g_send), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
                                                // Get error code
      if (cu_error != cudaSuccess) {
        fprintf(stderr,
        "CUDA error with pair density H2D async  memcpy ci: %s cpuid id is: %i\n ",
        cudaGetErrorString(cu_error), r->cpuid);
        error("Something's up with your cuda code");
      }
#endif

//	  const int tasksperbundle = pack_vars->tasksperbundle;
	  /* LAUNCH THE GPU KERNELS for ci & cj */
      // Setup 2d grid of GPU thread blocks for ci (number of tasks is
      // the y dimension and max_parts is the x dimension
      int numBlocks_y = 0;//tasks_left;
      int numBlocks_x = (bundle_n_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
      int bundle_part_0 = pack_vars->bundle_first_part[bid];
//              fprintf(stderr, "bundle_part_0 %i bundle_first_task %i\n", bundle_part_0, bundle_first_task);

      /* Launch the kernel for ci using data for ci and cj */
      runner_dopair_branch_gradient_gpu_aos_f4(d_parts_send, d_parts_recv,
		      d_a, d_H, stream[bid], numBlocks_x, numBlocks_y, bundle_part_0, bundle_n_parts);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with pair density kernel launch: %s cpuid id is: %i\n "
				"nbx %i nby %i max_parts_i %i max_parts_j %i\n",
				cudaGetErrorString(cu_error), r->cpuid, numBlocks_x, numBlocks_y, max_parts_i, max_parts_j);
		exit(0);
	  }
#endif

        // Copy results back to CPU BUFFERS
      cudaMemcpyAsync(&parts_recv[first_part_tmp_i], &d_parts_recv[first_part_tmp_i],
    	  	  bundle_n_parts * sizeof(struct part_aos_f4_g_recv), cudaMemcpyDeviceToHost, stream[bid]);
      cudaEventRecord(pair_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
//	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int bid = 0; bid < nBundles_temp; bid++){
		/*Time unpacking*/
		clock_gettime(CLOCK_REALTIME, &t0);

//		cudaStreamSynchronize(stream[bid]);
        cudaEventSynchronize(pair_end[bid]);

		clock_gettime(CLOCK_REALTIME, &t1);
		*gpu_time += (t1.tv_sec - t0.tv_sec) +
				(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	    /*Time unpacking*/
//		clock_gettime(CLOCK_REALTIME, &tp0);
//		int bundle_first_task = pack_vars->bundle_first_task_list[bid];

		for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

		  if(tid < tasks_packed){
		  clock_gettime(CLOCK_REALTIME, &tp0);
		  /*grab cell and task pointers*/
		  struct cell *cii = pack_vars->ci_list[tid];
		  struct cell *cjj = pack_vars->cj_list[tid];
		  struct task *tii = pack_vars->task_list[tid];
		  /*Let's lock ci*/
		  while(cell_locktree(cii)) {
			;  /* spin until we acquire the lock */
		  }
		  /*Let's lock cj*/
		  while(cell_locktree(cjj)) {
			;  /* spin until we acquire the lock */
		  }
		  /* Do the copy */
		  runner_do_ci_cj_gpu_unpack_neat_aos_f4_g(r, cii, cjj, parts_recv, 0,
				  &pack_length_unpack, tid, 2 * pack_vars->count_max_parts, e);

		  /* Record things for debugging */
		  cii->gpu_done_pair_g++;
		  cjj->gpu_done_pair_g++;

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
		  signal_sleeping_runners(s, tii);

		  tii->gpu_done = 1;


		  }
		}
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
//	/*Time end of unpacking*/
//	clock_gettime(CLOCK_REALTIME, &t1);
//	*packing_time += (t1.tv_sec - t0.tv_sec) +
//	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  } /*End of GPU work*/

void runner_dopair1_launch_f4_g_mcpy_Ker_mcpy(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars,
		struct task *t, struct part_aos_f4_g_send *parts_send, struct part_aos_f4_g_recv *parts_recv, struct part_aos_f4_g_send *d_parts_send,
		struct part_aos_f4_g_recv *d_parts_recv, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time, double *unpack_time, int4 *fparti_fpartj_lparti_lpartj,
		cudaEvent_t * pair_end){

	struct timespec t0, t1, tp0, tp1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;
	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/*tasks-packed needs decrementing before calculating packed_tmp as it was incremented in runner_dopair1_pack*/
//	const int packed_tmp = 2 * (tasks_packed - 1);

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero pair tasks packed but somehow got into GPU loop");
//	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[packed_tmp - 2];
	  pack_vars->bundle_first_part[nBundles_temp] = fparti_fpartj_lparti_lpartj[tasks_packed - 1].x;
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;

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
      for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
           tid++) {
        if (tid < tasks_packed) {
          /*Get an estimate for the max number of parts per cell in each bundle.
           *  Used for determining the number of GPU CUDA blocks*/
          int count_i = fparti_fpartj_lparti_lpartj[tid].z
															  - fparti_fpartj_lparti_lpartj[tid].x;
          parts_in_bundle_ci += count_i;
          max_parts_i = max(max_parts_i, count_i);
          int count_j = fparti_fpartj_lparti_lpartj[tid].w
															  - fparti_fpartj_lparti_lpartj[tid].y;
          parts_in_bundle_cj += count_j;
          max_parts_j = max(max_parts_j, count_j);

//		  last_task = tid;
        }
      }
      const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
      const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;

      cudaMemcpyAsync(&d_parts_send[first_part_tmp_i], &parts_send[first_part_tmp_i],
        bundle_n_parts * sizeof(struct part_aos_f4_g_send), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
                                                // Get error code
      if (cu_error != cudaSuccess) {
        fprintf(stderr,
        "CUDA error with pair density H2D async  memcpy ci: %s cpuid id is: %i\n ",
        cudaGetErrorString(cu_error), r->cpuid);
        error("Something's up with your cuda code");
      }
#endif
	}
      	for (int bid = 0; bid < nBundles_temp; bid++) {

          int max_parts_i = 0;
          int max_parts_j = 0;
          int parts_in_bundle_ci = 0;
          int parts_in_bundle_cj = 0;
//          const int first_task = bid * pack_vars->bundle_size;
//      	  int last_task = (bid + 1) * bundle_size;
          for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
                 tid++) {
              if (tid < tasks_packed) {
                /*Get an estimate for the max number of parts per cell in each bundle.
                 *  Used for determining the number of GPU CUDA blocks*/
                int count_i = fparti_fpartj_lparti_lpartj[tid].z
      															  - fparti_fpartj_lparti_lpartj[tid].x;
                parts_in_bundle_ci += count_i;
                max_parts_i = max(max_parts_i, count_i);
                int count_j = fparti_fpartj_lparti_lpartj[tid].w
      															  - fparti_fpartj_lparti_lpartj[tid].y;
                parts_in_bundle_cj += count_j;
                max_parts_j = max(max_parts_j, count_j);

//      		  last_task = tid;
              }
          }
          const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
          const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;
//////////////////////////////////
//	      const int tasksperbundle = pack_vars->tasksperbundle;
	      /* LAUNCH THE GPU KERNELS for ci & cj */
          // Setup 2d grid of GPU thread blocks for ci (number of tasks is
          // the y dimension and max_parts is the x dimension
          int numBlocks_y = 0;//tasks_left;
          int numBlocks_x = (bundle_n_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
          int bundle_part_0 = pack_vars->bundle_first_part[bid];
//          int bundle_first_task = pack_vars->bundle_first_task_list[bid];
//              fprintf(stderr, "bundle_part_0 %i bundle_first_task %i\n", bundle_part_0, bundle_first_task);

          /* Launch the kernel for ci using data for ci and cj */
          runner_dopair_branch_gradient_gpu_aos_f4(d_parts_send, d_parts_recv,
		        d_a, d_H, stream[bid], numBlocks_x, numBlocks_y, bundle_part_0, bundle_n_parts);

#ifdef CUDA_DEBUG
          cudaError_t cu_error = cudaPeekAtLastError(); // Get error code
	      if (cu_error != cudaSuccess) {
		    fprintf(stderr,
				"CUDA error with pair density kernel launch: %s cpuid id is: %i\n "
				"nbx %i nby %i max_parts_i %i max_parts_j %i\n",
				cudaGetErrorString(cu_error), r->cpuid, numBlocks_x, numBlocks_y, max_parts_i, max_parts_j);
		    exit(0);
	      }
#endif
      	}

    	for (int bid = 0; bid < nBundles_temp; bid++) {

          int max_parts_i = 0;
          int max_parts_j = 0;
          int parts_in_bundle_ci = 0;
          int parts_in_bundle_cj = 0;
//          const int first_task = bid * pack_vars->bundle_size;
//    	  int last_task = (bid + 1) * bundle_size;
          for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
               tid++) {
            if (tid < tasks_packed) {
              /*Get an estimate for the max number of parts per cell in each bundle.
               *  Used for determining the number of GPU CUDA blocks*/
              int count_i = fparti_fpartj_lparti_lpartj[tid].z
    															  - fparti_fpartj_lparti_lpartj[tid].x;
              parts_in_bundle_ci += count_i;
              max_parts_i = max(max_parts_i, count_i);
              int count_j = fparti_fpartj_lparti_lpartj[tid].w
    															  - fparti_fpartj_lparti_lpartj[tid].y;
              parts_in_bundle_cj += count_j;
              max_parts_j = max(max_parts_j, count_j);

//    		  last_task = tid;
            }
          }
          const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
          const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;
///////////////////////////////////////////////////////////////////
        // Copy results back to CPU BUFFERS
      cudaMemcpyAsync(&parts_recv[first_part_tmp_i], &d_parts_recv[first_part_tmp_i],
    	  	  bundle_n_parts * sizeof(struct part_aos_f4_g_recv), cudaMemcpyDeviceToHost, stream[bid]);
      cudaEventRecord(pair_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
//	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
    /*Time unpacking*/
//	clock_gettime(CLOCK_REALTIME, &t0);
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int bid = 0; bid < nBundles_temp; bid++){
	  clock_gettime(CLOCK_REALTIME, &t0);

//		cudaStreamSynchronize(stream[bid]);
	  cudaEventSynchronize(pair_end[bid]);

	  clock_gettime(CLOCK_REALTIME, &t1);
	  *gpu_time += (t1.tv_sec - t0.tv_sec) +
				(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	    /*Time unpacking*/
//	clock_gettime(CLOCK_REALTIME, &tp0);
//	  int bundle_first_task = pack_vars->bundle_first_task_list[bid];

	  for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

		  if(tid < tasks_packed){

		    /*grab cell and task pointers*/
		    struct cell *cii = pack_vars->ci_list[tid];
		    struct cell *cjj = pack_vars->cj_list[tid];
		    struct task *tii = pack_vars->task_list[tid];

		  /*Let's lock ci*/
		    while(cell_locktree(cii)) {
			  ;  /* spin until we acquire the lock */
		    }
		  /*Let's lock cj*/
		    while(cell_locktree(cjj)) {
			  ;  /* spin until we acquire the lock */
		    }
		  /* Do the copy */
		    /*Time unpacking*/
		    clock_gettime(CLOCK_REALTIME, &tp0);
		    runner_do_ci_cj_gpu_unpack_neat_aos_f4_g(r, cii, cjj, parts_recv, 0,
				  &pack_length_unpack, tid, 2 * pack_vars->count_max_parts, e);

		    /* Record things for debugging */
		    cii->gpu_done_pair_g++;
		    cjj->gpu_done_pair_g++;

		    tii->gpu_done = 1;
		    /*Time end of unpacking*/
		    clock_gettime(CLOCK_REALTIME, &tp1);
		    *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
		    (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
		    /*schedule my dependencies (Only unpacks really)*/
		    enqueue_dependencies(s, tii);
		    /*Signal sleeping runners*/
		    signal_sleeping_runners(s, tii);

		    /* Release the locks */
		    cell_unlocktree(cii);
		    /* Release the locks */
		    cell_unlocktree(cjj);

		  }
	  }
	  /*Time end of unpacking*/
//	  clock_gettime(CLOCK_REALTIME, &tp1);
//	  *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
//	  (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
//	  *packing_time += (tp1.tv_sec - tp0.tv_sec) +
//	  (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
//		*packing_time += (tp1.tv_sec - tp0.tv_sec) +
//		(tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
	/*Time end of unpacking*/
//	clock_gettime(CLOCK_REALTIME, &t1);
//	*packing_time += (t1.tv_sec - t0.tv_sec) +
//	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  } /*End of GPU work*/

void runner_dopair1_launch_f(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars, struct cell *ci,
		struct task *t, struct part_aos_f *parts_aos, struct part_aos_f *d_parts_aos, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time){

	struct timespec t0, t1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;
	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/*tasks-packed needs decrementing before calculating packed_tmp as it was incremented in runner_dopair1_pack*/
	const int packed_tmp = 2 * (tasks_packed - 1);

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero pair tasks packed but somehow got into GPU loop");
	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[packed_tmp - 2];
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;
    /*Copy arrays containing first and last part for each task to GPU*/
    cudaMemcpy(pack_vars->d_task_first_part, pack_vars->task_first_part,
               2 * tasks_packed * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_task_last_part, pack_vars->task_last_part,
               2 * tasks_packed * sizeof(int), cudaMemcpyHostToDevice);

    /*Copy cell shifts to device*/
    cudaMemcpy(pack_vars->d_shiftx, pack_vars->shiftx,
    		2 * tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_shifty, pack_vars->shifty,
    		2 * tasks_packed * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(pack_vars->d_shiftz, pack_vars->shiftz,
    		2 * tasks_packed * sizeof(double), cudaMemcpyHostToDevice);

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
      for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
           tid++) {
        if (tid < tasks_packed) {
          /*Get an estimate for the max number of parts per cell in each bundle.
           *  Used for determining the number of GPU CUDA blocks*/
      	const int tid_tmp = 2 * tid;
          int count_i = pack_vars->task_last_part[tid_tmp]
															  - pack_vars->task_first_part[tid_tmp];
          parts_in_bundle_ci += count_i;
          max_parts_i = max(max_parts_i, count_i);
          int count_j = pack_vars->task_last_part[tid_tmp + 1]
															  - pack_vars->task_first_part[tid_tmp + 1];
          parts_in_bundle_cj += count_j;
          max_parts_j = max(max_parts_j, count_j);
        }
      }
      const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
      const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;
      cudaMemcpyAsync(&d_parts_aos[first_part_tmp_i], &parts_aos[first_part_tmp_i],
        bundle_n_parts * sizeof(struct part_aos_f), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
                                                // Get error code
      if (cu_error != cudaSuccess) {
        fprintf(stderr,
        "CUDA error with pair density H2D async  memcpy ci: %s cpuid id is: %i\n ",
        cudaGetErrorString(cu_error), r->cpuid);
        error("Something's up with your cuda code");
      }
#endif

	  const int tasksperbundle = pack_vars->tasksperbundle;
	  /* LAUNCH THE GPU KERNELS for ci & cj */
      int tid = 0;
      int offset = bid * tasksperbundle;
      int tasks_left = tasksperbundle;
      if (bid == nBundles_temp - 1) {
        tasks_left =
        		tasks_packed - (nBundles_temp - 1) * tasksperbundle;
      }

      // Setup 2d grid of GPU thread blocks for ci (number of tasks is
      // the y dimension and max_parts is the x dimension
      int numBlocks_y = tasks_left;
//      int numBlocks_x = (max_parts_i + BLOCK_SIZE - 1) / BLOCK_SIZE;
      int bundle_first_task = pack_vars->bundle_first_task_list[bid];
      const char *loop_type = "density";

        /* Launch the kernel for ci using data for ci and cj */
        runner_dopairci_branch_density_gpu_aos_f(d_parts_aos,
        pack_vars->d_task_first_part, pack_vars->d_task_last_part,
	    d_a, d_H, loop_type, stream[bid], bid, BLOCK_SIZE,
	    tasks_packed, tasksperbundle, max_parts_i, max_parts_j, numBlocks_y, tid,
	    offset, bundle_first_task, pack_vars->d_shiftx, pack_vars->d_shifty, pack_vars->d_shiftz);

        /* Launch the kernel for ci using data for ci and cj */
        runner_dopaircj_branch_density_gpu_aos_f(d_parts_aos,
        pack_vars->d_task_first_part, pack_vars->d_task_last_part,
	    d_a, d_H, loop_type, stream[bid], bid, BLOCK_SIZE,
	    tasks_packed, tasksperbundle, max_parts_i, max_parts_j, numBlocks_y, tid,
	    offset, bundle_first_task, pack_vars->d_shiftx, pack_vars->d_shifty, pack_vars->d_shiftz);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density kernel launch: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		exit(0);
	  }
#endif

        // Copy results back to CPU BUFFERS
      cudaMemcpyAsync(&parts_aos[first_part_tmp_i], &d_parts_aos[first_part_tmp_i],
    	  	  bundle_n_parts * sizeof(struct part_aos_f), cudaMemcpyDeviceToHost, stream[bid]);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
    /*Time unpacking*/
	clock_gettime(CLOCK_REALTIME, &t0);
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int tid = 0; tid < tasks_packed; tid++) {

	  /*grab cell and task pointers*/
	  struct cell *cii = pack_vars->ci_list[tid];
      struct cell *cjj = pack_vars->cj_list[tid];
	  struct task *tii = pack_vars->task_list[tid];

	  /*Let's lock ci*/
	  while(cell_locktree(cii)) {
		;  /* spin until we acquire the lock */
	  }
	  /*Let's lock cj*/
	  while(cell_locktree(cjj)) {
	    ;  /* spin until we acquire the lock */
	  }
	  /* Do the copy */
	  runner_do_ci_cj_gpu_unpack_neat_aos_f(r, cii, cjj, parts_aos, 0,
		  	  &pack_length_unpack, tid, 2 * pack_vars->count_max_parts, e);

	  /* Record things for debugging */
	  cii->gpu_done_pair_f++;
	  cjj->gpu_done_pair_f++;

	  tii->gpu_done = 1;

      /*schedule my dependencies (Only unpacks really)*/
      enqueue_dependencies(s, tii);
      /*Signal sleeping runners*/
      signal_sleeping_runners(s, tii);

	  /* Release the locks */
	  cell_unlocktree(cii);
	  /* Release the locks */
	  cell_unlocktree(cjj);
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
	/*Time end of unpacking*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*packing_time += (t1.tv_sec - t0.tv_sec) +
	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
} /*End of GPU work*/

void runner_dopair1_launch_f4_f_one_memcpy(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars,
		struct task *t, struct part_aos_f4_f_send *parts_send, struct part_aos_f4_f_recv *parts_recv, struct part_aos_f4_f_send *d_parts_send,
		struct part_aos_f4_f_recv *d_parts_recv, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time, double *unpack_time, int4 *fparti_fpartj_lparti_lpartj,
		cudaEvent_t * pair_end){

	struct timespec t0, t1, tp0, tp1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;
	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/*tasks-packed needs decrementing before calculating packed_tmp as it was incremented in runner_dopair1_pack*/
//	const int packed_tmp = 2 * (tasks_packed - 1);

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero pair tasks packed but somehow got into GPU loop");
//	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[packed_tmp - 2];
	  pack_vars->bundle_first_part[nBundles_temp] = fparti_fpartj_lparti_lpartj[tasks_packed - 1].x;
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;

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
      for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
           tid++) {
        if (tid < tasks_packed) {
          /*Get an estimate for the max number of parts per cell in each bundle.
           *  Used for determining the number of GPU CUDA blocks*/
          int count_i = fparti_fpartj_lparti_lpartj[tid].z
															  - fparti_fpartj_lparti_lpartj[tid].x;
          parts_in_bundle_ci += count_i;
          max_parts_i = max(max_parts_i, count_i);
          int count_j = fparti_fpartj_lparti_lpartj[tid].w
															  - fparti_fpartj_lparti_lpartj[tid].y;
          parts_in_bundle_cj += count_j;
          max_parts_j = max(max_parts_j, count_j);

//		  last_task = tid;
        }
      }
      const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
      const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;

      cudaMemcpyAsync(&d_parts_send[first_part_tmp_i], &parts_send[first_part_tmp_i],
        bundle_n_parts * sizeof(struct part_aos_f4_f_send), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
                                                // Get error code
      if (cu_error != cudaSuccess) {
        fprintf(stderr,
        "CUDA error with pair density H2D async  memcpy ci: %s cpuid id is: %i\n ",
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
      int numBlocks_y = 0;//tasks_left;
      int numBlocks_x = (bundle_n_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
      int bundle_part_0 = pack_vars->bundle_first_part[bid];
//      int bundle_first_task = pack_vars->bundle_first_task_list[bid];
//              fprintf(stderr, "bundle_part_0 %i bundle_first_task %i\n", bundle_part_0, bundle_first_task);

      /* Launch the kernel for ci using data for ci and cj */
      runner_dopair_branch_force_gpu_aos_f4(d_parts_send, d_parts_recv,
		      d_a, d_H, stream[bid], numBlocks_x, numBlocks_y, bundle_part_0, bundle_n_parts);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with pair density kernel launch: %s cpuid id is: %i\n "
				"nbx %i nby %i max_parts_i %i max_parts_j %i\n",
				cudaGetErrorString(cu_error), r->cpuid, numBlocks_x, numBlocks_y, max_parts_i, max_parts_j);
		exit(0);
	  }
#endif

        // Copy results back to CPU BUFFERS
      cudaMemcpyAsync(&parts_recv[first_part_tmp_i], &d_parts_recv[first_part_tmp_i],
    	  	  bundle_n_parts * sizeof(struct part_aos_f4_f_recv), cudaMemcpyDeviceToHost, stream[bid]);
      cudaEventRecord(pair_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
	  cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
//	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int bid = 0; bid < nBundles_temp; bid++){
		/*Time unpacking*/
		clock_gettime(CLOCK_REALTIME, &t0);

//		cudaStreamSynchronize(stream[bid]);
        cudaEventSynchronize(pair_end[bid]);

		clock_gettime(CLOCK_REALTIME, &t1);
		*gpu_time += (t1.tv_sec - t0.tv_sec) +
				(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	    /*Time unpacking*/
//		clock_gettime(CLOCK_REALTIME, &tp0);
//		int bundle_first_task = pack_vars->bundle_first_task_list[bid];

		for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

		  if(tid < tasks_packed){
		  clock_gettime(CLOCK_REALTIME, &tp0);
		  /*grab cell and task pointers*/
		  struct cell *cii = pack_vars->ci_list[tid];
		  struct cell *cjj = pack_vars->cj_list[tid];
		  struct task *tii = pack_vars->task_list[tid];
		  /*Let's lock ci*/
		  while(cell_locktree(cii)) {
			;  /* spin until we acquire the lock */
		  }
		  /*Let's lock cj*/
		  while(cell_locktree(cjj)) {
			;  /* spin until we acquire the lock */
		  }
		  /* Do the copy */
		  runner_do_ci_cj_gpu_unpack_neat_aos_f4_f(r, cii, cjj, parts_recv, 0,
				  &pack_length_unpack, tid, 2 * pack_vars->count_max_parts, e);

		  /* Record things for debugging */
		  cii->gpu_done_pair_f++;
		  cjj->gpu_done_pair_f++;

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
		  signal_sleeping_runners(s, tii);

		  tii->gpu_done = 1;


		  }
		}
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
//	/*Time end of unpacking*/
//	clock_gettime(CLOCK_REALTIME, &t1);
//	*packing_time += (t1.tv_sec - t0.tv_sec) +
//	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  } /*End of GPU work*/

void runner_dopair1_launch_f4_f_mcpy_Ker_mcpy(struct runner *r, struct scheduler *s, struct pack_vars_pair *pack_vars,
		struct task *t, struct part_aos_f4_f_send *parts_send, struct part_aos_f4_f_recv *parts_recv, struct part_aos_f4_f_send *d_parts_send,
		struct part_aos_f4_f_recv *d_parts_recv, cudaStream_t * stream, float d_a, float d_H,
		struct engine *e, double *packing_time, double *gpu_time, double *unpack_time, int4 *fparti_fpartj_lparti_lpartj,
		cudaEvent_t * pair_end){

	struct timespec t0, t1, tp0, tp1; //
    clock_gettime(CLOCK_REALTIME, &t0);

	/* Identify the number of GPU bundles to run in ideal case*/
	int nBundles_temp = pack_vars->nBundles;
	/*How many tasks have we packed?*/
	const int tasks_packed = pack_vars->tasks_packed;

	/*How many tasks should be in a bundle?*/
	const int bundle_size = pack_vars->bundle_size;

	/* Special case for incomplete bundles (when having leftover tasks not enough to fill a bundle) */
	if (pack_vars->launch_leftovers) {
	  nBundles_temp = (tasks_packed + bundle_size - 1) / bundle_size;
	  if(tasks_packed == 0) error("zero pair tasks packed but somehow got into GPU loop");
//	  pack_vars->bundle_first_part[nBundles_temp] = pack_vars->task_first_part[packed_tmp - 2];
	  pack_vars->bundle_first_part[nBundles_temp] = fparti_fpartj_lparti_lpartj[tasks_packed - 1].x;
	}
    /* Identify the last particle for each bundle of tasks */
    for (int bid = 0; bid < nBundles_temp - 1; bid++) {
    	pack_vars->bundle_last_part[bid] = pack_vars->bundle_first_part[bid + 1];
	}
    /* special treatment for the last bundle */
    if(nBundles_temp > 1)pack_vars->bundle_last_part[nBundles_temp - 1] = pack_vars->count_parts;
    else pack_vars->bundle_last_part[0] = pack_vars->count_parts;

	/* Launch the copies for each bundle and run the GPU kernel */
	/*We don't go into this loop if tasks_left_self == 1 as
	 nBundles_temp will be zero DUHDUHDUHDUHHHHHH!!!!!*/
	for (int bid = 0; bid < nBundles_temp; bid++) {

      int max_parts_i = 0;
      int max_parts_j = 0;
      int parts_in_bundle_ci = 0;
      int parts_in_bundle_cj = 0;
//      const int first_task = bid * pack_vars->bundle_size;
//	  int last_task = (bid + 1) * bundle_size;
      for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
           tid++) {
        if (tid < tasks_packed) {
          /*Get an estimate for the max number of parts per cell in each bundle.
           *  Used for determining the number of GPU CUDA blocks*/
          int count_i = fparti_fpartj_lparti_lpartj[tid].z
															  - fparti_fpartj_lparti_lpartj[tid].x;
          parts_in_bundle_ci += count_i;
          max_parts_i = max(max_parts_i, count_i);
          int count_j = fparti_fpartj_lparti_lpartj[tid].w
															  - fparti_fpartj_lparti_lpartj[tid].y;
          parts_in_bundle_cj += count_j;
          max_parts_j = max(max_parts_j, count_j);

//		  last_task = tid;
        }
      }
      const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
      const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;

      cudaMemcpyAsync(&d_parts_send[first_part_tmp_i], &parts_send[first_part_tmp_i],
        bundle_n_parts * sizeof(struct part_aos_f4_f_send), cudaMemcpyHostToDevice, stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
                                                // Get error code
      if (cu_error != cudaSuccess) {
        fprintf(stderr,
        "CUDA error with pair density H2D async  memcpy ci: %s cpuid id is: %i\n ",
        cudaGetErrorString(cu_error), r->cpuid);
        error("Something's up with your cuda code");
      }
#endif
	}
      	for (int bid = 0; bid < nBundles_temp; bid++) {

          int max_parts_i = 0;
          int max_parts_j = 0;
          int parts_in_bundle_ci = 0;
          int parts_in_bundle_cj = 0;
//          const int first_task = bid * pack_vars->bundle_size;
//      	  int last_task = (bid + 1) * bundle_size;
          for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
                 tid++) {
              if (tid < tasks_packed) {
                /*Get an estimate for the max number of parts per cell in each bundle.
                 *  Used for determining the number of GPU CUDA blocks*/
                int count_i = fparti_fpartj_lparti_lpartj[tid].z
      															  - fparti_fpartj_lparti_lpartj[tid].x;
                parts_in_bundle_ci += count_i;
                max_parts_i = max(max_parts_i, count_i);
                int count_j = fparti_fpartj_lparti_lpartj[tid].w
      															  - fparti_fpartj_lparti_lpartj[tid].y;
                parts_in_bundle_cj += count_j;
                max_parts_j = max(max_parts_j, count_j);

//      		  last_task = tid;
              }
          }
          const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
          const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;
//////////////////////////////////
	      /* LAUNCH THE GPU KERNELS for ci & cj */
          // Setup 2d grid of GPU thread blocks for ci (number of tasks is
          // the y dimension and max_parts is the x dimension
          int numBlocks_y = 0;//tasks_left;
          int numBlocks_x = (bundle_n_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
          int bundle_part_0 = pack_vars->bundle_first_part[bid];
//          int bundle_first_task = pack_vars->bundle_first_task_list[bid];
//              fprintf(stderr, "bundle_part_0 %i bundle_first_task %i\n", bundle_part_0, bundle_first_task);

          /* Launch the kernel for ci using data for ci and cj */
          runner_dopair_branch_force_gpu_aos_f4(d_parts_send, d_parts_recv,
		        d_a, d_H, stream[bid], numBlocks_x, numBlocks_y, bundle_part_0, bundle_n_parts);

#ifdef CUDA_DEBUG
          cudaError_t cu_error = cudaPeekAtLastError(); // Get error code
	      if (cu_error != cudaSuccess) {
		    fprintf(stderr,
				"CUDA error with pair density kernel launch: %s cpuid id is: %i\n "
				"nbx %i nby %i max_parts_i %i max_parts_j %i\n",
				cudaGetErrorString(cu_error), r->cpuid, numBlocks_x, numBlocks_y, max_parts_i, max_parts_j);
		    exit(0);
	      }
#endif
      	}

    	for (int bid = 0; bid < nBundles_temp; bid++) {

          int max_parts_i = 0;
          int max_parts_j = 0;
          int parts_in_bundle_ci = 0;
          int parts_in_bundle_cj = 0;
          for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size;
               tid++) {
            if (tid < tasks_packed) {
              /*Get an estimate for the max number of parts per cell in each bundle.
               *  Used for determining the number of GPU CUDA blocks*/
              int count_i = fparti_fpartj_lparti_lpartj[tid].z
    															  - fparti_fpartj_lparti_lpartj[tid].x;
              parts_in_bundle_ci += count_i;
              max_parts_i = max(max_parts_i, count_i);
              int count_j = fparti_fpartj_lparti_lpartj[tid].w
    															  - fparti_fpartj_lparti_lpartj[tid].y;
              parts_in_bundle_cj += count_j;
              max_parts_j = max(max_parts_j, count_j);
            }
          }
          const int first_part_tmp_i = pack_vars->bundle_first_part[bid];
          const int bundle_n_parts = pack_vars->bundle_last_part[bid] - first_part_tmp_i;
///////////////////////////////////////////////////////////////////
        // Copy results back to CPU BUFFERS
      cudaMemcpyAsync(&parts_recv[first_part_tmp_i], &d_parts_recv[first_part_tmp_i],
    	  	  bundle_n_parts * sizeof(struct part_aos_f4_f_recv), cudaMemcpyDeviceToHost, stream[bid]);
      cudaEventRecord(pair_end[bid], stream[bid]);

#ifdef CUDA_DEBUG
      cudaError_t cu_error = cudaPeekAtLastError(); // cudaGetLastError();        //
										// Get error code
	  if (cu_error != cudaSuccess) {
		fprintf(stderr,
				"CUDA error with self density D2H memcpy: %s cpuid id is: %i\n ",
				cudaGetErrorString(cu_error), r->cpuid);
		error("Something's up with your cuda code");
	  }
#endif
	}/*End of looping over bundles to launch in streams*/

	/* Make sure all the kernels and copies back are finished */
//	cudaDeviceSynchronize();

    /*Time end of GPU work*/
	clock_gettime(CLOCK_REALTIME, &t1);
	*gpu_time += (t1.tv_sec - t0.tv_sec) +
			(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
    /*Time unpacking*/
//	clock_gettime(CLOCK_REALTIME, &t0);
	/* Now copy the data back from the CPU thread-local buffers to the cells */
	/* Pack length counter for use in unpacking */
	int pack_length_unpack=0;
	for (int bid = 0; bid < nBundles_temp; bid++){
	  clock_gettime(CLOCK_REALTIME, &t0);

//		cudaStreamSynchronize(stream[bid]);
	  cudaEventSynchronize(pair_end[bid]);

	  clock_gettime(CLOCK_REALTIME, &t1);
	  *gpu_time += (t1.tv_sec - t0.tv_sec) +
				(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;

	    /*Time unpacking*/
//	clock_gettime(CLOCK_REALTIME, &tp0);

	  for (int tid = bid * bundle_size; tid < (bid + 1) * bundle_size; tid++) {

		  if(tid < tasks_packed){

		    /*grab cell and task pointers*/
		    struct cell *cii = pack_vars->ci_list[tid];
		    struct cell *cjj = pack_vars->cj_list[tid];
		    struct task *tii = pack_vars->task_list[tid];

		  /*Let's lock ci*/
		    while(cell_locktree(cii)) {
			  ;  /* spin until we acquire the lock */
		    }
		  /*Let's lock cj*/
		    while(cell_locktree(cjj)) {
			  ;  /* spin until we acquire the lock */
		    }
		  /* Do the copy */
		    /*Time unpacking*/
		    clock_gettime(CLOCK_REALTIME, &tp0);
		    runner_do_ci_cj_gpu_unpack_neat_aos_f4_f(r, cii, cjj, parts_recv, 0,
				  &pack_length_unpack, tid, 2 * pack_vars->count_max_parts, e);

		    /* Record things for debugging */
		    cii->gpu_done_pair_f++;
		    cjj->gpu_done_pair_f++;

		    tii->gpu_done = 1;
		    /*Time end of unpacking*/
		    clock_gettime(CLOCK_REALTIME, &tp1);
		    *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
		    (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
		    /*schedule my dependencies (Only unpacks really)*/
		    enqueue_dependencies(s, tii);
		    /*Signal sleeping runners*/
		    signal_sleeping_runners(s, tii);

		    /* Release the locks */
		    cell_unlocktree(cii);
		    /* Release the locks */
		    cell_unlocktree(cjj);

		  }
	  }
	  /*Time end of unpacking*/
//	  clock_gettime(CLOCK_REALTIME, &tp1);
//	  *unpack_time += (tp1.tv_sec - tp0.tv_sec) +
//	  (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
//	  *packing_time += (tp1.tv_sec - tp0.tv_sec) +
//	  (tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
//		*packing_time += (tp1.tv_sec - tp0.tv_sec) +
//		(tp1.tv_nsec - tp0.tv_nsec) / 1000000000.0;
	}
	/* Zero counters for the next pack operations */
	pack_vars->count_parts = 0;
	pack_vars->tasks_packed = 0;
	/*Time end of unpacking*/
//	clock_gettime(CLOCK_REALTIME, &t1);
//	*packing_time += (t1.tv_sec - t0.tv_sec) +
//	(t1.tv_nsec - t0.tv_nsec) / 1000000000.0;
  } /*End of GPU work*/

