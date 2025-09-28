/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <sched.h>
#include <stdlib.h>
#include <string.h>
#ifdef SWIFT_DEBUG_THREADPOOL
#include <dlfcn.h>
#endif

/* This object's header. */
#include "threadpool.h"

/* Local headers. */
#include "atomic.h"
#include "clocks.h"
#include "error.h"
#include "minmax.h"

/* Keys for thread specific data. */
static pthread_key_t threadpool_tid;

/* Affinity mask shared by all threads, and if set. */
#ifdef HAVE_SETAFFINITY
static cpu_set_t thread_affinity;
static int thread_affinity_set = 0;
#endif

/* Local declarations. */
static void threadpool_apply_affinity_mask(void);

/* ---- Forward decls for queue extension (defined at end of file). ---- */
static int threadpool_queue_is_active(struct threadpool *tp);
static void threadpool_chomp_queue(struct threadpool *tp, int tid);
static void threadpool_queue_on_init(struct threadpool *tp);
static void threadpool_queue_on_clean(struct threadpool *tp);

#ifdef SWIFT_DEBUG_THREADPOOL
/**
 * @brief Store a log entry of the given chunk.
 */
static void threadpool_log(struct threadpool *tp, int tid, size_t chunk_size,
                           ticks tic, ticks toc) {
  struct mapper_log *log = &tp->logs[tid > 0 ? tid : 0];

  /* Check if we need to re-allocate the log buffer. */
  if (log->count == log->size) {
    log->size *= 2;
    struct mapper_log_entry *new_log;
    if ((new_log = (struct mapper_log_entry *)malloc(
             sizeof(struct mapper_log_entry) * log->size)) == NULL)
      error("Failed to re-allocate mapper log.");
    memcpy(new_log, log->log, sizeof(struct mapper_log_entry) * log->count);
    free(log->log);
    log->log = new_log;
  }

  /* Store the new entry. */
  struct mapper_log_entry *entry = &log->log[log->count];
  entry->tid = tid;
  entry->chunk_size = chunk_size;
  entry->tic = tic;
  entry->toc = toc;
  entry->map_function = tp->map_function;
  log->count++;
}

void threadpool_dump_log(struct threadpool *tp, const char *filename,
                         int reset) {

  /* Open the output file. */
  FILE *fd;
  if ((fd = fopen(filename, "w")) == NULL)
    error("Failed to create log file '%s'.", filename);

  /* Create a buffer of function names. */
  const int max_names = 100;
  struct name_entry {
    threadpool_map_function map_function;
    const char *name;
  };
  struct name_entry names[max_names];
  bzero(names, sizeof(struct name_entry) * max_names);

  /* Write a header. */
  fprintf(fd, "# map_function thread_id chunk_size tic toc");
  fprintf(fd, "# {'num_threads': %i, 'cpufreq': %lli}", tp->num_threads,
          clocks_get_cpufreq());

  /* Loop over the per-tid logs and dump them. */
  for (int k = 0; k < tp->num_threads; k++) {
    struct mapper_log *log = &tp->logs[k];

    /* Loop over the log entries and dump them. */
    for (int i = 0; i < log->count; i++) {

      struct mapper_log_entry *entry = &log->log[i];

      /* Look for the function pointer in the buffer. */
      int nid = 0;
      while (nid < max_names && names[nid].map_function != entry->map_function)
        nid++;

      /* If the name was not found, make a new entry. */
      if (nid == max_names) {
        for (int j = 1; j < max_names; j++) names[j - 1] = names[j];
        names[0].map_function = entry->map_function;
        Dl_info dl_info;
        dladdr(entry->map_function, &dl_info);
        names[0].name = dl_info.dli_sname;
        nid = 0;
      }

      /* Log a line to the file. */
      fprintf(fd,
              "%s %i %i %lli %lli
              ", names[nid].name, entry->tid,
              entry->chunk_size,
              entry->tic, entry->toc);
    }

    /* Clear the log if requested. */
    if (reset) log->count = 0;
  }

  /* Close the file. */
  fclose(fd);
}
#endif  // SWIFT_DEBUG_THREADPOOL

/**
 * @brief Runner main loop, get a chunk and call the mapper function.
 */
static void threadpool_chomp(struct threadpool *tp, int tid) {

  /* Store the thread ID as thread specific data. */
  int localtid = tid;
  pthread_setspecific(threadpool_tid, &localtid);

  /* Loop until we can't get a chunk. */
  while (1) {
    /* Compute the desired chunk size. */
    ptrdiff_t chunk_size;
    if (tp->map_data_chunk == threadpool_uniform_chunk_size) {
      chunk_size = ((tid + 1) * tp->map_data_size / tp->num_threads) -
                   (tid * tp->map_data_size / tp->num_threads);
    } else {
      chunk_size =
          (tp->map_data_size - tp->map_data_count) / (2 * tp->num_threads);
      if (chunk_size > tp->map_data_chunk) chunk_size = tp->map_data_chunk;
    }
    if (chunk_size < 1) chunk_size = 1;

    /* A chunk cannot exceed INT_MAX, as we use int elements in map_function. */
    if (chunk_size > INT_MAX) chunk_size = INT_MAX;

    /* Get a chunk and check its size. */
    size_t task_ind = atomic_add(&tp->map_data_count, chunk_size);
    if (task_ind >= tp->map_data_size) break;
    if (task_ind + chunk_size > tp->map_data_size)
      chunk_size = tp->map_data_size - task_ind;

/* Call the mapper function. */
#ifdef SWIFT_DEBUG_THREADPOOL
    ticks tic = getticks();
#endif

    tp->map_function((char *)tp->map_data + (tp->map_data_stride * task_ind),
                     chunk_size, tp->map_extra_data);

#ifdef SWIFT_DEBUG_THREADPOOL
    threadpool_log(tp, tid, chunk_size, tic, getticks());
#endif
  }
}

/**
 * @brief The thread start routine. Loops until told to exit.
 *
 * @param data the threadpool we are part of.
 */
static void *threadpool_runner(void *data) {

  /* Our threadpool. */
  struct threadpool *tp = (struct threadpool *)data;

  /* Our affinity, if set. */
  threadpool_apply_affinity_mask();

  /* Main loop. */
  while (1) {

    /* Let the controller know that this thread is waiting. */
    swift_barrier_wait(&tp->wait_barrier);

    /* Wait for the controller. */
    swift_barrier_wait(&tp->run_barrier);

    /* If no map function is specified, just die. We use this as a mechanism
       to shut down threads without leaving the barriers in an invalid state. */
    if (tp->map_function == NULL) pthread_exit(NULL);

    /* Do actual work (queue or chunk mode). */
    if (threadpool_queue_is_active(tp)) {
      threadpool_chomp_queue(tp, atomic_inc(&tp->num_threads_running));
    } else {
      threadpool_chomp(tp, atomic_inc(&tp->num_threads_running));
    }
  }
}

/**
 * @brief Initialises the #threadpool with a given number of threads.
 *
 * @param tp The #threadpool.
 * @param num_threads The number of threads.
 */
void threadpool_init(struct threadpool *tp, int num_threads) {

  /* Initialize the thread counters. */
  tp->num_threads = num_threads;

  /* Create thread local data areas. Only do this once for all threads. */
  pthread_key_create(&threadpool_tid, NULL);

  /* Store the main thread ID as thread specific data. */
  static int localtid = 0;
  pthread_setspecific(threadpool_tid, &localtid);

#ifdef SWIFT_DEBUG_THREADPOOL
  if ((tp->logs = (struct mapper_log *)malloc(sizeof(struct mapper_log) *
                                              num_threads)) == NULL)
    error("Failed to allocate mapper logs.");
  for (int k = 0; k < num_threads; k++) {
    tp->logs[k].size = threadpool_log_initial_size;
    tp->logs[k].count = 0;
    if ((tp->logs[k].log = (struct mapper_log_entry *)malloc(
             sizeof(struct mapper_log_entry) * tp->logs[k].size)) == NULL)
      error("Failed to allocate mapper log.");
  }
#endif

  /* If there is only a single thread, do nothing more as of here as
     we will just do work in the (blocked) calling thread. */
  if (num_threads == 1) return;

  /* Init the barriers. */
  if (swift_barrier_init(&tp->wait_barrier, NULL, num_threads) != 0 ||
      swift_barrier_init(&tp->run_barrier, NULL, num_threads) != 0)
    error("Failed to initialize barriers.");

  /* Set the task counter to zero. */
  tp->map_data_size = 0;
  tp->map_data_count = 0;
  tp->map_data_stride = 0;
  tp->map_data_chunk = 0;
  tp->map_function = NULL;

  /* Allocate the threads, one less than requested since the calling thread
     works as well. */
  if ((tp->threads = (pthread_t *)malloc(sizeof(pthread_t) *
                                         (num_threads - 1))) == NULL) {
    error("Failed to allocate thread array.");
  }

  /* Create and start the threads. */
  for (int k = 0; k < num_threads - 1; k++) {
    if (pthread_create(&tp->threads[k], NULL, &threadpool_runner, tp) != 0)
      error("Failed to create threadpool runner thread.");
  }

  /* Wait for all the threads to be up and running. */
  swift_barrier_wait(&tp->wait_barrier);

  /* Queue extension lazy init hook. */
  threadpool_queue_on_init(tp);
}

/**
 * @brief Map a function to an array of data in parallel using a #threadpool.
 *
 * The function @c map_function is called on each element of @c map_data
 * in parallel.
 *
 * @param tp The #threadpool on which to run.
 * @param map_function The function that will be applied to the map data.
 * @param map_data The data on which the mapping function will be called.
 * @param N Number of elements in @c map_data.
 * @param stride Size, in bytes, of each element of @c map_data.
 * @param chunk Number of map data elements to pass to the function at a time,
 *        or #threadpool_auto_chunk_size to choose the number dynamically
 *        depending on the number of threads and tasks (recommended), or
 *        #threadpool_uniform_chunk_size to spread the tasks evenly over the
 *        threads in one go.
 * @param extra_data Addtitional pointer that will be passed to the mapping
 *        function, may contain additional data.
 */
void threadpool_map(struct threadpool *tp, threadpool_map_function map_function,
                    void *map_data, size_t N, int stride, int chunk,
                    void *extra_data) {

#ifdef SWIFT_DEBUG_THREADPOOL
  ticks tic_total = getticks();
#endif

  /* If we just have a single thread, call the map function directly. */
  if (tp->num_threads == 1) {

    if (N <= INT_MAX) {
      map_function(map_data, N, extra_data);

#ifdef SWIFT_DEBUG_THREADPOOL
      tp->map_function = map_function;
      threadpool_log(tp, 0, N, tic_total, getticks());
#endif
    } else {

      /* N > INT_MAX, we need to do this in chunks as map_function only takes
       * an int. */
      size_t chunk_size = INT_MAX;
      size_t data_size = N;
      size_t data_count = 0;
      while (1) {

/* Call the mapper function. */
#ifdef SWIFT_DEBUG_THREADPOOL
        ticks tic = getticks();
#endif
        map_function((char *)map_data + (stride * data_count), chunk_size,
                     extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
        threadpool_log(tp, 0, chunk_size, tic, getticks());
#endif
        /* Get the next chunk and check its size. */
        data_count += chunk_size;
        if (data_count >= data_size) break;
        if (data_count + chunk_size > data_size)
          chunk_size = data_size - data_count;
      }
    }

    return;
  }

  /* Set the map data and signal the threads. */
  tp->map_data_stride = stride;
  tp->map_data_size = N;
  tp->map_data_count = 0;
  if (chunk == threadpool_auto_chunk_size) {
    tp->map_data_chunk =
        max((N / (tp->num_threads * threadpool_default_chunk_ratio)), 1U);
  } else if (chunk == threadpool_uniform_chunk_size) {
    tp->map_data_chunk = threadpool_uniform_chunk_size;
  } else {
    tp->map_data_chunk = chunk;
  }
  tp->map_function = map_function;
  tp->map_data = map_data;
  tp->map_extra_data = extra_data;
  tp->num_threads_running = 0;

  /* Wait for all the threads to be up and running. */
  swift_barrier_wait(&tp->run_barrier);

  /* Do some work while I'm at it. */
  threadpool_chomp(tp, tp->num_threads - 1);

  /* Wait for all threads to be done. */
  swift_barrier_wait(&tp->wait_barrier);

#ifdef SWIFT_DEBUG_THREADPOOL
  /* Log the total call time to thread id -1. */
  threadpool_log(tp, -1, N, tic_total, getticks());
#endif
}

/**
 * @brief Re-sets the log for this #threadpool.
 */
#ifdef SWIFT_DEBUG_THREADPOOL
void threadpool_reset_log(struct threadpool *tp) {
  for (int k = 0; k < tp->num_threads; k++) tp->logs[k].count = 0;
}
#endif

/**
 * @brief Frees up the memory allocated for this #threadpool.
 */
void threadpool_clean(struct threadpool *tp) {

  if (tp->num_threads > 1) {
    /* Destroy the runner threads by calling them with a NULL mapper function
     * and waiting for all the threads to terminate. This ensures that no
     * thread is still waiting at a barrier. */
    tp->map_function = NULL;
    swift_barrier_wait(&tp->run_barrier);
    for (int k = 0; k < tp->num_threads - 1; k++) {
      void *retval;
      pthread_join(tp->threads[k], &retval);
    }

    /* Release the barriers. */
    if (swift_barrier_destroy(&tp->wait_barrier) != 0 ||
        swift_barrier_destroy(&tp->run_barrier) != 0)
      error("Failed to destroy threadpool barriers.");

    /* Clean up memory. */
    free(tp->threads);
  }

#ifdef SWIFT_DEBUG_THREADPOOL
  for (int k = 0; k < tp->num_threads; k++) {
    free(tp->logs[k].log);
  }
  free(tp->logs);
#endif

  /* Queue extension cleanup. */
  threadpool_queue_on_clean(tp);
}

/**
 * @brief return the threadpool id of the current thread.
 */
int threadpool_gettid(void) {
  int *tid = (int *)pthread_getspecific(threadpool_tid);
  return *tid;
}

#ifdef HAVE_SETAFFINITY
/**
 * @brief set an affinity mask to be used for all threads.
 *
 * @param affinity the mask to use.
 */
void threadpool_set_affinity_mask(cpu_set_t *affinity) {
  memcpy(&thread_affinity, affinity, sizeof(cpu_set_t));
  thread_affinity_set = 1;
}
#endif

/**
 * @brief apply the affinity mask the current thread, if set.
 *
 */
static void threadpool_apply_affinity_mask(void) {
#ifdef HAVE_SETAFFINITY
  if (thread_affinity_set) {
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &thread_affinity);
  }
#endif
}

/* ========================================================================== */
/*                      Queue-based extension (end of file)                   */
/* ========================================================================== */

/* Minimal queue/stealing support to allow adding pointers to the map data
 * from within the mapper. No extra public task types. The queue stores small
 * tasks {ptr, size}. For user-added work, size==1. */

/* Task description. */
struct tpq_task {
  void *ptr;
  int size;
};

/* Per-thread queue. */
struct tpq_queue {
  pthread_mutex_t lock;
  struct tpq_task *buf;
  int cap;
  int head;
  int tail;
};

/* Per-pool state. */
struct tpq_state {
  struct threadpool *tp;
  struct tpq_queue *queues; /* [num_threads] */
  pthread_mutex_t sleep_lock;
  pthread_cond_t sleep_cv;
  int sleepers;
  int active;
  int tasks_in_flight;
  /* Cached map params. */
  threadpool_map_function map_function;
  void *map_extra_data;
};

/* Registry of pools. */
struct tpq_node {
  struct tpq_state *s;
  struct tpq_node *next;
};
static struct tpq_node *tpq_registry = NULL;

/* Helpers to find/create/destroy state. */
static struct tpq_state *threadpool_queue_get(struct threadpool *tp,
                                              int create) {
  struct tpq_node *n = tpq_registry;
  while (n) {
    if (n->s->tp == tp) return n->s;
    n = n->next;
  }
  if (!create) return NULL;
  struct tpq_state *s = (struct tpq_state *)malloc(sizeof(*s));
  if (s == NULL) error("Failed to allocate queue state.");
  s->tp = tp;
  s->queues =
      (struct tpq_queue *)calloc(tp->num_threads, sizeof(struct tpq_queue));
  if (s->queues == NULL) error("Failed to allocate queues.");
  for (int i = 0; i < tp->num_threads; i++) {
    pthread_mutex_init(&s->queues[i].lock, NULL);
    s->queues[i].buf = NULL;
    s->queues[i].cap = 0;
    s->queues[i].head = 0;
    s->queues[i].tail = 0;
  }
  pthread_mutex_init(&s->sleep_lock, NULL);
  pthread_cond_init(&s->sleep_cv, NULL);
  s->sleepers = 0;
  s->active = 0;
  &s->tasks_in_flight = 0;
  s->map_function = NULL;
  s->map_extra_data = NULL;
  struct tpq_node *nn = (struct tpq_node *)malloc(sizeof(*nn));
  if (nn == NULL) error("Failed to allocate queue registry node.");
  nn->s = s;
  nn->next = tpq_registry;
  tpq_registry = nn;
  return s;
}

static void threadpool_queue_destroy(struct threadpool *tp) {
  struct tpq_node **pp = &tpq_registry;
  while (*pp) {
    if ((*pp)->s->tp == tp) {
      struct tpq_node *dead = *pp;
      *pp = dead->next;
      struct tpq_state *s = dead->s;
      for (int i = 0; i < tp->num_threads; i++) {
        pthread_mutex_destroy(&s->queues[i].lock);
        free(s->queues[i].buf);
      }
      free(s->queues);
      pthread_mutex_destroy(&s->sleep_lock);
      pthread_cond_destroy(&s->sleep_cv);
      free(s);
      free(dead);
      return;
    }
    pp = &(*pp)->next;
  }
}

/* Queue buffer helpers. */
static void tpq_reserve(struct tpq_queue *q, int need) {
  if (q->cap >= need) return;
  int new_cap = q->cap ? q->cap : 64;
  while (new_cap < need) new_cap <<= 1;
  struct tpq_task *nb =
      (struct tpq_task *)malloc(sizeof(struct tpq_task) * new_cap);
  if (nb == NULL) error("Failed to allocate queue buffer.");
  const int n = q->tail - q->head;
  for (int i = 0; i < n; i++)
    nb[i] = q->cap ? q->buf[(q->head + i) & (q->cap - 1)]
                   : (struct tpq_task){NULL, 0};
  free(q->buf);
  q->buf = nb;
  q->cap = new_cap;
  q->head = 0;
  q->tail = n;
}

static void tpq_push_locked(struct tpq_queue *q, void *ptr, int size) {
  const int n = q->tail - q->head;
  if (q->cap - n < 1) tpq_reserve(q, q->cap ? 2 * q->cap : 64);
  const int idx = q->tail & (q->cap - 1);
  q->buf[idx].ptr = ptr;
  q->buf[idx].size = size;
  q->tail++;
}

static int tpq_owner_pop(struct tpq_queue *q, void **ptr, int *size) {
  int ok = 0;
  pthread_mutex_lock(&q->lock);
  if (q->tail - q->head > 0) {
    q->tail--;
    int idx = q->tail & (q->cap - 1);
    *ptr = q->buf[idx].ptr;
    *size = q->buf[idx].size;
    ok = 1;
  }
  pthread_mutex_unlock(&q->lock);
  return ok;
}

static int tpq_steal(struct tpq_queue *q, void **ptr, int *size) {
  int ok = 0;
  pthread_mutex_lock(&q->lock);
  if (q->tail - q->head > 0) {
    int idx = q->head & (q->cap - 1);
    *ptr = q->buf[idx].ptr;
    *size = q->buf[idx].size;
    q->head++;
    ok = 1;
  }
  pthread_mutex_unlock(&q->lock);
  return ok;
}

/* Hooks from init/clean. */
static void threadpool_queue_on_init(struct threadpool *tp) { (void)tp; }
static void threadpool_queue_on_clean(struct threadpool *tp) {
  threadpool_queue_destroy(tp);
}

/* Public: add N pointers as single-element tasks. */
void threadpool_queue_add(struct threadpool *tp, void **ptrs, int n_ptrs) {
  if (n_ptrs <= 0) return;
  struct tpq_state *s = threadpool_queue_get(tp, /*create*/ 1);
  int tid = threadpool_gettid();
  pthread_mutex_lock(&s->queues[tid].lock);
  for (int i = 0; i < n_ptrs; i++) tpq_push_locked(&s->queues[tid], ptrs[i], 1);
  pthread_mutex_unlock(&s->queues[tid].lock);
  atomic_add(&s->tasks_in_flight, n_ptrs);
  pthread_mutex_lock(&s->sleep_lock);
  if (s->sleepers > 0) pthread_cond_signal(&s->sleep_cv);
  pthread_mutex_unlock(&s->sleep_lock);
}

/* Public: queue-based mapper (same signature as normal, different entry). */
void threadpool_map_with_queue(struct threadpool *tp,
                               threadpool_map_function map_function,
                               void *map_data, size_t N, int stride, int chunk,
                               void *extra_data) {

#ifdef SWIFT_DEBUG_THREADPOOL
  ticks tic_total = getticks();
#endif

  /* Single-threaded path identical to threadpool_map(). */
  if (tp->num_threads == 1) {
    if (N <= INT_MAX) {
      map_function(map_data, N, extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
      tp->map_function = map_function;
      threadpool_log(tp, 0, N, tic_total, getticks());
#endif
    } else {
      size_t chunk_size = INT_MAX, data_count = 0;
      while (data_count < N) {
#ifdef SWIFT_DEBUG_THREADPOOL
        ticks tic = getticks();
#endif
        int c = (int)min(chunk_size, N - data_count);
        map_function((char *)map_data + (stride * data_count), c, extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
        threadpool_log(tp, 0, c, tic, getticks());
#endif
        data_count += (size_t)c;
      }
    }
    return;
  }

  /* Get/create state and set parameters. */
  struct tpq_state *s = threadpool_queue_get(tp, /*create*/ 1);
  s->map_function = map_function;
  s->map_extra_data = extra_data;

  /* Seed initial tasks (chunks) round-robin. */
  int chunk_size;
  if (chunk == threadpool_auto_chunk_size) {
    chunk_size =
        max((N / (tp->num_threads * threadpool_default_chunk_ratio)), 1U);
  } else if (chunk == threadpool_uniform_chunk_size) {
    chunk_size = (N + tp->num_threads - 1) / tp->num_threads;
  } else {
    chunk_size = (chunk > 0) ? chunk : 1;
  }

  size_t start = 0;
  int seed_tid = 0;
  int n_tasks = 0;
  while (start < N) {
    int c = (int)min((size_t)chunk_size, N - start);
    void *ptr = (char *)map_data + (stride * start);
    pthread_mutex_lock(&s->queues[seed_tid].lock);
    tpq_push_locked(&s->queues[seed_tid], ptr, c);
    pthread_mutex_unlock(&s->queues[seed_tid].lock);
    n_tasks++;
    start += (size_t)c;
    seed_tid = (seed_tid + 1) % tp->num_threads;
  }
  atomic_add(&s->tasks_in_flight, n_tasks);

  /* Activate queue mode and set standard fields for runner. */
  s->active = 1;
  tp->map_function = map_function;
  tp->map_data = map_data;
  tp->map_extra_data = extra_data;
  tp->num_threads_running = 0;

  /* Start and participate. */
  swift_barrier_wait(&tp->run_barrier);
  threadpool_chomp_queue(tp, tp->num_threads - 1);
  swift_barrier_wait(&tp->wait_barrier);

  s->active = 0;

#ifdef SWIFT_DEBUG_THREADPOOL
  threadpool_log(tp, -1, N, tic_total, getticks());
#endif
}

/* Active query used by runner. */
static int threadpool_queue_is_active(struct threadpool *tp) {
  struct tpq_state *s = threadpool_queue_get(tp, 0);
  return (s != NULL) && s->active;
}

/* Try get a task for tid: local pop then steal. */
static int threadpool_queue_take(struct tpq_state *s, int tid, void **ptr,
                                 int *size) {
  if (tpq_owner_pop(&s->queues[tid], ptr, size)) return 1;
  for (int off = 1; off < s->tp->num_threads; off++) {
    int vid = (tid + off) % s->tp->num_threads;
    if (tpq_steal(&s->queues[vid], ptr, size)) return 1;
  }
  return 0;
}

/* Queue-mode work loop. */
static void threadpool_chomp_queue(struct threadpool *tp, int tid) {

  /* Store the thread ID as thread specific data. */
  int localtid = tid;
  pthread_setspecific(threadpool_tid, &localtid);

  struct tpq_state *s = threadpool_queue_get(tp, 0);
  while (1) {
    void *ptr;
    int size;
    if (threadpool_queue_take(s, tid, &ptr, &size)) {
#ifdef SWIFT_DEBUG_THREADPOOL
      ticks tic = getticks();
#endif
      s->map_function(ptr, size, s->map_extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
      threadpool_log(tp, tid, size, tic, getticks());
#endif
      atomic_dec(&s->tasks_in_flight);
      continue;
    }

    /* Completion test: no tasks in flight and all queues empty. */
    if (s->tasks_in_flight == 0) {
      int empty = 1;
      for (int i = 0; i < tp->num_threads; i++) {
        pthread_mutex_lock(&s->queues[i].lock);
        empty &= (s->queues[i].tail - s->queues[i].head) == 0;
        pthread_mutex_unlock(&s->queues[i].lock);
        if (!empty) break;
      }
      if (empty) break;
    }

    /* Sleep until new work is signalled. */
    pthread_mutex_lock(&s->sleep_lock);
    s->sleepers++;
    pthread_cond_wait(&s->sleep_cv, &s->sleep_lock);
    s->sleepers--;
    pthread_mutex_unlock(&s->sleep_lock);
  }
}
