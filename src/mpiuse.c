/* This file is part of SWIFT.
 * Copyright (c) 2019 Peter W. Draper (p.w.draper@durham.ac.uk)
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
 *  @file mpiuse.c
 *  @brief file of routines to report about MPI tasks used in SWIFT.
 */
/* Config parameters. */
#include "../config.h"

#if defined(SWIFT_MPIUSE_REPORTS) && defined(WITH_MPI)

/* Standard includes. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

/* Local defines. */
#include "mpiuse.h"

/* Local includes. */
#include "atomic.h"
#include "clocks.h"
#include "engine.h"
#include "error.h"
#include "memuse_rnodes.h"

/* The initial size and increment of the log entries buffer. */
#define MPIUSE_INITLOG 1000000

/* A megabyte for conversions. */
#define MEGABYTE 1048576.0

/* Also recorded in logger. */
extern int engine_rank;
extern int engine_current_step;

/* Entry for logger of MPI send and recv requests in a step. */
struct mpiuse_log_entry {

  /* Type and subtype of MPI task. */
  int type;
  int subtype;

  /* Step of action. */
  int step;

  /* Whether an activation, send or recv, or if handoff completed. Not the
   * same as delivered, need to match across ranks to see that. */
  int activation;

  /* Memory of the request. */
  size_t size;

  /* Pointer to the request associated with the call. Needs to be
   * unique and match to the successful */
  union {
    void *ptr;
    uint8_t vptr[sizeof(uintptr_t)]; /* For rnode keys. */
  };

  /* Ticks at time of this action. */
  ticks tic;

  /* Time taken for handoff of this action. */
  ticks acttic;

  /* Whether request is still active, i.e. successful test not seen. */
  int active;

  /* Rank of otherside of communication. */
  int otherrank;

  /* The tag. */
  int tag;
};

/* The log of activations and handoffs. All volatile as accessed from threads
 * that use the value to synchronise. */
static struct mpiuse_log_entry *volatile mpiuse_log = NULL;
static volatile size_t mpiuse_log_size = 0;
static volatile size_t mpiuse_log_count = 0;
static volatile size_t mpiuse_log_done = 0;

/**
 * @brief reallocate the entries log if space is needed.
 */
static void mpiuse_log_reallocate(size_t ind) {

  if (ind == 0) {

    /* Need to perform initialization. Be generous. */
    if ((mpiuse_log = (struct mpiuse_log_entry *)malloc(
             sizeof(struct mpiuse_log_entry) * MPIUSE_INITLOG)) == NULL)
      error("Failed to allocate MPI use log.");

    /* Last action. */
    mpiuse_log_size = MPIUSE_INITLOG;

  } else {
    struct mpiuse_log_entry *new_log;
    if ((new_log = (struct mpiuse_log_entry *)malloc(
             sizeof(struct mpiuse_log_entry) *
             (mpiuse_log_size + MPIUSE_INITLOG))) == NULL)
      error("Failed to re-allocate MPI use log.");

    /* Wait for all writes to the old buffer to complete. */
    while (mpiuse_log_done < mpiuse_log_size)
      ;

    /* Copy to new buffer. */
    memcpy(new_log, mpiuse_log,
           sizeof(struct mpiuse_log_entry) * mpiuse_log_size);
    free(mpiuse_log);
    mpiuse_log = new_log;

    /* Last action, releases waiting threads. */
    atomic_add(&mpiuse_log_size, MPIUSE_INITLOG);
  }
}

/**
 * @brief Log an MPI request or handoff.
 *
 * @param type the task type (send or recv).
 * @param subtype the task subtype.
 * @param ptr pointer to the MPI request.
 * @param activation if not is a successful MPI_Test, not MPI_Isend or
 *        MPI_Irecv.
 * @param size the size in bytes of memory to be transfered or received.
 *             0 for a deactivation.
 * @param otherrank other rank associated with the transfer.
 * @param tag the MPI tag.
 */
void mpiuse_log_allocation(int type, int subtype, void *ptr, int activation,
                           size_t size, int otherrank, int tag) {

  size_t ind = atomic_inc(&mpiuse_log_count);

  /* If we are at the current size we need more space. */
  if (ind == mpiuse_log_size) mpiuse_log_reallocate(ind);

  /* Other threads wait for space. */
  while (ind > mpiuse_log_size)
    ;

  /* Record the log. */
  mpiuse_log[ind].step = engine_current_step;
  mpiuse_log[ind].type = type;
  mpiuse_log[ind].subtype = subtype;
  mpiuse_log[ind].activation = activation;
  mpiuse_log[ind].size = size;
  mpiuse_log[ind].ptr = ptr;
  mpiuse_log[ind].otherrank = otherrank;
  mpiuse_log[ind].tag = tag;
  mpiuse_log[ind].tic = getticks();
  mpiuse_log[ind].acttic = 0;
  mpiuse_log[ind].active = 1;
  atomic_inc(&mpiuse_log_done);
}

/**
 * @brief dump the log to a file and reset, if anything to dump.
 *
 * @param filename name of file for log dump.
 * @param stepticks the clock ticks at the start of step, if dumping a step,
 *                  otherwise some locally relative time that might help
 *                  synchronize across ranks.
 */
void mpiuse_log_dump(const char *filename, ticks stepticks) {

  /* Skip if nothing logged this step. */
  if (mpiuse_log_count == 0) return;

  // ticks tic = getticks();

  /* Create the radix tree root node. */
  struct memuse_rnode *memuse_rnode_root =
      (struct memuse_rnode *)calloc(1, sizeof(struct memuse_rnode));

  /* Stop any new logs from being processed while we are dumping. */
  size_t log_count = mpiuse_log_count;

  /* Open the output file. */
  FILE *fd;
  if ((fd = fopen(filename, "w")) == NULL) {
    message("Failed to create MPI use log file '%s', logs not dumped.",
            filename);
    return;
  }

  /* Write a header. */
  fprintf(fd,
          "# stic etic dtic step rank otherrank type itype subtype isubtype "
          "activation tag size sum\n");

  size_t mpiuse_current = 0;
  size_t mpiuse_max = 0;
  double mpiuse_sum = 0;
  size_t mpiuse_actcount = 0;
  for (size_t k = 0; k < log_count; k++) {

    /* Check if this address has already been recorded. */
    struct memuse_rnode *child = memuse_rnode_find_child(
        memuse_rnode_root, 0, mpiuse_log[k].vptr, sizeof(uintptr_t));

    if (child != NULL && child->value != -1) {

      /* Should be the handoff. Check that. */
      if (mpiuse_log[k].activation) {

        /* Used twice, this is an error, but just complain as not fatal. */
#if SWIFT_DEBUG_CHECKS
        message(
            "Used the same MPI request address twice "
            "(%s/%s: %d->%d: %zd/%d)",
            taskID_names[mpiuse_log[k].type],
            subtaskID_names[mpiuse_log[k].subtype], engine_rank,
            mpiuse_log[k].otherrank, mpiuse_log[k].size, mpiuse_log[k].tag);
#endif
        continue;
      }

      /* Free, update the missing fields, size of request is removed. */
      struct mpiuse_log_entry *oldlog = (struct mpiuse_log_entry *)child->value;
      mpiuse_log[k].size = -oldlog->size;
      mpiuse_log[k].otherrank = oldlog->otherrank;
      mpiuse_log[k].tag = oldlog->tag;

      /* Time taken to handoff. */
      mpiuse_log[k].acttic = mpiuse_log[k].tic - oldlog->tic;

      /* And deactivate this key. */
      child->value = -1;

      /* And mark this as handed off. */
      mpiuse_log[k].active = 0;
      oldlog->active = 0;

    } else if (child == NULL && mpiuse_log[k].activation) {

      /* Not found, so new send/recv which we store the log against the
       * address. */
      memuse_rnode_insert_child(memuse_rnode_root, 0, mpiuse_log[k].vptr,
                                sizeof(uintptr_t), (int64_t)&mpiuse_log[k]);

    } else if (child == NULL && !mpiuse_log[k].activation) {

      /* Unmatched handoff, not OK, but not fatal. */
#if SWIFT_DEBUG_CHECKS
      if (mpiuse_log[k].ptr != NULL) {
        message("Unmatched MPI_Test found: (%s/%s: %d->%d: %zd/%d)",
                taskID_names[mpiuse_log[k].type],
                subtaskID_names[mpiuse_log[k].subtype], engine_rank,
                mpiuse_log[k].otherrank, mpiuse_log[k].size, mpiuse_log[k].tag);
      }
#endif
      continue;
    } else if (mpiuse_log[k].activation) {

      /* Must be previously released request with the same address, so we
       * store. */
      memuse_rnode_insert_child(memuse_rnode_root, 0, mpiuse_log[k].vptr,
                                sizeof(uintptr_t), (int64_t)&mpiuse_log[k]);

    } else {
      message("Weird MPI log record found: (%s/%s: %d->%d: %zd/%d/%d/%p)",
              taskID_names[mpiuse_log[k].type],
              subtaskID_names[mpiuse_log[k].subtype], engine_rank,
              mpiuse_log[k].otherrank, mpiuse_log[k].size, mpiuse_log[k].tag,
              mpiuse_log[k].activation, mpiuse_log[k].ptr);
      continue;
    }

    /* Sum of memory in flight. */
    mpiuse_current += mpiuse_log[k].size;

    /* Gather for stats report. */
    if (mpiuse_log[k].activation) {
      if (mpiuse_log[k].size > mpiuse_max) mpiuse_max = mpiuse_log[k].size;
      mpiuse_sum += (double)mpiuse_log[k].size;
      mpiuse_actcount++;
    }

    /* And output. */
    fprintf(fd, "%lld %lld %lld %d %d %d %s %d %s %d %d %d %zd %zd\n",
            mpiuse_log[k].tic - stepticks,
            mpiuse_log[k].tic - clocks_start_ticks, mpiuse_log[k].acttic,
            mpiuse_log[k].step, engine_rank, mpiuse_log[k].otherrank,
            taskID_names[mpiuse_log[k].type], mpiuse_log[k].type,
            subtaskID_names[mpiuse_log[k].subtype], mpiuse_log[k].subtype,
            mpiuse_log[k].activation, mpiuse_log[k].tag, mpiuse_log[k].size,
            mpiuse_current);
  }

#ifdef MEMUSE_RNODE_DUMP
  /* Debug dump of tree. */
  // memuse_rnode_dump(0, memuse_rnode_root, 0);
#endif

  /* Write our statistics. */
  fprintf(fd, "##\n");
  fprintf(fd, "## Number of requests: %zd\n", mpiuse_actcount);
  fprintf(fd, "## Maximum request size: %.4f (MB)\n", mpiuse_max / MEGABYTE);
  fprintf(fd, "## Sum of all requests: %.4f (MB)\n", mpiuse_sum / MEGABYTE);
  fprintf(fd, "## Mean of all requests: %.4f (MB)\n",
          mpiuse_sum / (double)mpiuse_actcount / MEGABYTE);
  fprintf(fd, "##\n");

  /* Now check any still active logs, these are errors all should match. */
  if (mpiuse_current != 0) {
    message("Some MPI requests have not been completed");
    for (size_t k = 0; k < log_count; k++) {
      if (mpiuse_log[k].active)
        message("%s/%s: %d->%d: %zd/%d)", taskID_names[mpiuse_log[k].type],
                subtaskID_names[mpiuse_log[k].subtype], engine_rank,
                mpiuse_log[k].otherrank, mpiuse_log[k].size, mpiuse_log[k].tag);
    }
  }

  /* Finished with the rnodes. */
  memuse_rnode_cleanup(memuse_rnode_root);

  /* Clear the log. We expect this to clear step to step, unlike memory. */
  mpiuse_log_count = 0;
  mpiuse_log_done = 0;

  /* Close the file. */
  fflush(fd);
  fclose(fd);

  // message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
  //        clocks_getunit());
}

/**
 * @brief dump the log for using the given rank to generate a standard
 *        name for the output. Used when exiting in error.
 *
 * @param rank the rank exiting in error.
 */
void mpiuse_log_dump_error(int rank) {
  char filename[60];
  sprintf(filename, "mpiuse-error-report-rank%d.txt", rank);
  mpiuse_log_dump(filename, clocks_start_ticks);
}

#endif /* defined(SWIFT_MPIUSE_REPORTS) && defined(WITH_MPI) */
