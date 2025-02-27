/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* Standard includes. */
#include <errno.h>
#include <libgen.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local includes. */
#include "error.h"
#include "swift_lustre_api.h"

/* Lustre API */
#ifdef HAVE_LUSTREAPI
#include <lustre/lustre_user.h>
#include <lustre/lustreapi.h>
#endif

/* Number of OSTs to pre-allocate space for. */
#define PREALLOC (100)

/* Bytes in a TiB */
#define TiB (1024.0 * 1024.0 * 1024.0)

/* Bytes in a MiB */
#define MiB (1024.0 * 1024.0)

/**
 * @brief Allocate storage for a number of OSTs in a OST scan storage struct.
 *
 * Note does not reset count or fullcount. Zero these if you want an
 * empty struct.
 *
 * @param ost_infos pointer to the storage structure.
 * @param size number of OSTs to make space for.
 */
void swift_ost_store_alloc(struct swift_ost_store *ost_infos, int size) {
#ifdef HAVE_LUSTREAPI
  ost_infos->size = size;
  ost_infos->infos =
      (struct swift_ost_info *)malloc(sizeof(struct swift_ost_info) * size);
  if (ost_infos->infos == NULL)
    error("Failed to allocate space for an OST scan");
  memset(ost_infos->infos, 0, sizeof(struct swift_ost_info) * size);
#endif
}

/**
 * @brief Create a copy of an OST scan storage struct.
 *
 * @param ost_infos_src pointer to the storage structure to copy
 * @param ost_infos_dst pointer to a storage structure to populate with
 *                      the copy. Assumed to have no OST space so not
 *                      used or initialized.
 */
void swift_ost_store_copy(struct swift_ost_store *ost_infos_src,
                          struct swift_ost_store *ost_infos_dst) {
#ifdef HAVE_LUSTREAPI
  ost_infos_dst->size = ost_infos_src->fullcount; /* Used size. */
  ost_infos_dst->count = ost_infos_src->count;
  ost_infos_dst->fullcount = ost_infos_src->fullcount;
  ost_infos_dst->infos = (struct swift_ost_info *)malloc(
      sizeof(struct swift_ost_info) * ost_infos_dst->size);
  if (ost_infos_dst->infos == NULL)
    error("Failed to allocate space for an OST scan copy");
  memcpy(ost_infos_dst->infos, ost_infos_src->infos,
         sizeof(struct swift_ost_info) * ost_infos_dst->size);
#endif
}

/**
 * @brief Initialize an OST scan storage structure.
 *
 * @param ost_infos pointer to the storage structure.
 */
void swift_ost_store_init(struct swift_ost_store *ost_infos) {
#ifdef HAVE_LUSTREAPI
  swift_ost_store_alloc(ost_infos, PREALLOC);
  ost_infos->count = 0;
  ost_infos->fullcount = 0;
#endif
}

/**
 * @brief Release any storage associated with an OST scan storage structure.
 *
 * @param ost_infos pointer to the storage structure.
 */
void swift_ost_store_free(struct swift_ost_store *ost_infos) {
#ifdef HAVE_LUSTREAPI
  free(ost_infos->infos);
  ost_infos->infos = NULL;
  ost_infos->count = 0;
  ost_infos->fullcount = 0;
  ost_infos->size = 0;
#endif
}

/**
 * @brief Write about an OST storage structure to a given FILE.
 *
 * @param file FILE stream to write output to.
 * @param ost_infos pointer to the storage structure.
 */
void swift_ost_store_write(FILE *file, struct swift_ost_store *ost_infos) {
#ifdef HAVE_LUSTREAPI
  fprintf(file, "# %5s %21s %21s %21s\n", "Index", "Size (MiB)", "Used (MiB)",
          "Free (MiB)");
  size_t ssum = 0;
  size_t usum = 0;
  size_t smin = ost_infos->infos[0].size;
  size_t smax = 0;
  size_t umin = ost_infos->infos[0].used;
  size_t umax = 0;

  for (int i = 0; i < ost_infos->count; i++) {
    int msize = (int)(ost_infos->infos[i].size / MiB);
    int mused = (int)(ost_infos->infos[i].used / MiB);
    fprintf(file, "# %5d %21d %21d %21d\n", ost_infos->infos[i].index, msize,
            mused, msize - mused);

    ssum += ost_infos->infos[i].size;
    usum += ost_infos->infos[i].used;

    if (ost_infos->infos[i].size > smax) smax = ost_infos->infos[i].size;
    if (ost_infos->infos[i].size < smin) smin = ost_infos->infos[i].size;

    if (ost_infos->infos[i].used > umax) umax = ost_infos->infos[i].used;
    if (ost_infos->infos[i].used < umin) umin = ost_infos->infos[i].used;
  }
  if (ost_infos->count == ost_infos->fullcount) {
    /* Size is for used OSTs not all, so don't report as misleading. */
    fprintf(file,
            "# Filesystem size:%.2f TiB used:%.2f TiB free:%.2f TiB %.2f%%\n",
            ssum / TiB, usum / TiB, (ssum - usum) / TiB,
            100.0 * (double)(ssum - usum) / (double)ssum);
    fprintf(file, "# Min/max size: %.2f/%.2f TiB Min/max used: %.2f/%.2f TiB\n",
            smin / TiB, smax / TiB, umin / TiB, umax / TiB);
  } else {
    fprintf(file, "#\n");
  }
#endif
}

/**
 * @brief Print information about OST storage structure
 *
 * @param ost_infos pointer to the storage structure.
 * @param verbose if non zero additional information will be written
 *                to stdout.
 */
void swift_ost_store_print(struct swift_ost_store *ost_infos, int verbose) {
#ifdef HAVE_LUSTREAPI
  message("#  OSTs, using %d of %d", ost_infos->count, ost_infos->fullcount);
  if (verbose) swift_ost_store_write(stdout, ost_infos);
#endif
}

#ifdef HAVE_LUSTREAPI
/**
 * @brief Store information about an OST.
 *
 * @param ost_infos pointer to the storage structure.
 * @param index the index, zero based.
 * @param size the total size in bytes.
 * @param used the number of bytes used.
 */
static void swift_ost_store(struct swift_ost_store *ost_infos, int index,
                            size_t size, size_t used) {
  /* Add extra space if needed. Note not thread safe. */
  if (ost_infos->fullcount == ost_infos->size - 1) {
    size_t newsize = ost_infos->size + PREALLOC;
    struct swift_ost_info *newinfos = (struct swift_ost_info *)malloc(
        sizeof(struct swift_ost_info) * newsize);
    if (newinfos == NULL) error("Failed to allocate space for OST information");
    memset(newinfos, 0, sizeof(struct swift_ost_info) * newsize);
    memcpy(newinfos, ost_infos->infos,
           sizeof(struct swift_ost_info) * ost_infos->size);
    free(ost_infos->infos);
    ost_infos->infos = newinfos;
    ost_infos->size = newsize;
  }
  int count = ost_infos->count++;
  ost_infos->infos[count].index = index;
  ost_infos->infos[count].size = size;
  ost_infos->infos[count].used = used;
  ost_infos->fullcount = ost_infos->count;
}
#endif

/**
 * @brief Scan the OSTs associated with a lustre file system given a path.
 *
 * On exit the ost_infos struct will be populated with the
 * the number of OSTs found and details of the size and used bytes in each
 * OST.
 *
 * @param path a directory on the lustre file system, ideally the mount point.
 * @param ost_infos pointer to the storage structure.
 *
 * @return 0 on success, otherwise an error will have been reported to stdout.
 * If an error occurs the store will never be changed.
 */
int swift_ost_scan(const char *path, struct swift_ost_store *ost_infos) {

  int rc = 0;
#ifdef HAVE_LUSTREAPI
  char mntdir[PATH_MAX] = {0};
  char fsname[PATH_MAX] = {0};
  char cpath[PATH_MAX] = {0};

  /* Check this path exists. */
  if (!realpath(path, cpath)) {
    rc = errno;
    message("Not a filesystem path '%s': %s", path, strerror(rc));
  } else {

    /* Parse the path into the mount point and file system name. */
    if (llapi_search_mounts(cpath, 0, mntdir, fsname) == 0) {
      if (mntdir[0] != '\0') {
        struct obd_statfs stat_buf;
        struct obd_uuid uuid_buf;

        /* Loop while OSTs are located. */
        for (int index = 0;; index++) {
          memset(&stat_buf, 0, sizeof(struct obd_statfs));
          memset(&uuid_buf, 0, sizeof(struct obd_uuid));

          rc = llapi_obd_statfs(mntdir, LL_STATFS_LOV, index, &stat_buf,
                                &uuid_buf);
          rc = -rc;
          if (rc == ENODEV || rc == EAGAIN || rc == EINVAL || rc == EFAULT) {
            /* Nothing we can query here, so time to stop search. */
            break;
          }

          /* Inactive devices are empty. */
          if (rc == ENODATA) {
            swift_ost_store(ost_infos, index, 0, 0);
          } else {
            size_t used =
                (stat_buf.os_blocks - stat_buf.os_bfree) * stat_buf.os_bsize;
            size_t total = stat_buf.os_blocks * stat_buf.os_bsize;
            swift_ost_store(ost_infos, index, total, used);
          }
        }
        rc = 0;

      } else {
        message("No lustre mount point found for path: %s", path);
        rc = 1;
      }
    } else {
      message("Failed to locate a lustre mount point using path: %s", path);
      rc = 1;
    }
  }
#endif
  return rc;
}

#ifdef HAVE_LUSTREAPI
/** Comparison function for OST free space. */
static int ostcmp(const void *p1, const void *p2) {
  const struct swift_ost_info *i1 = (const struct swift_ost_info *)p1;
  const struct swift_ost_info *i2 = (const struct swift_ost_info *)p2;

  /* size_t ints so some care is needed to return an int. */
  size_t f1 = i1->size - i1->used;
  size_t f2 = i2->size - i2->used;
  if (f1 < f2) return 1;
  if (f1 > f2) return -1;
  return 0;
}
#endif

/**
 * @brief Sort the OSTs into decreasing free space culling those that do not
 * meet a free space threshold.
 *
 * @param ost_infos pointer to populated storage structure.
 * @param minfree the number of MiB that the OST should be capable of
 *                storing. Zero for no effect.
 */
void swift_ost_cull(struct swift_ost_store *ost_infos, int minfree) {
#ifdef HAVE_LUSTREAPI
  /* Sort by free space. */
  qsort(ost_infos->infos, ost_infos->count, sizeof(struct swift_ost_info),
        ostcmp);

  /* And cull if needed. */
  if (minfree > 0) {
    size_t bytesfree = minfree * (size_t)MiB;

    /* Always keep at least one! */
    for (int i = 1; i < ost_infos->count; i++) {
      struct swift_ost_info *curr = &ost_infos->infos[i];
      if ((curr->size - curr->used) < bytesfree) {

        /* Throw the rest away. Note fullcount now decoupled. */
        ost_infos->count = i;
      }
    }
  }
#endif
}

/**
 * @brief Get the next OST in an incrementing sequence.
 *
 * @param ost_infos pointer to populated storage structure.
 * @param arrayindex the last used array index, start with 0.
 *                   This will be wrapped as needed use as input for next
 *                   call.
 * @param count number of OSTs that will be used to stripe, that is the
 *              increment, usually 1. Only makes sense if the OST list is not
 *              culled as this implicitly assumes OSTs are in index order.
 * @return the selected OST index.
 */
int swift_ost_next(struct swift_ost_store *ost_infos, int *arrayindex,
                   int count) {
#ifdef HAVE_LUSTREAPI
  int index = (*arrayindex % ost_infos->count);
  *arrayindex = index + count;
  return ost_infos->infos[index].index;
#else
  return 0;
#endif
}

/**
 * @brief Remove an OST by index from the store.
 *
 * @param ost_infos pointer to populated storage structure.
 * @param index index of the OST to remove.
 */
void swift_ost_remove(struct swift_ost_store *ost_infos, int index) {

#ifdef HAVE_LUSTREAPI
  /* Find the array index. */
  int arrayindex = -1;
  for (int i = 0; i < ost_infos->fullcount; i++) {
    if (ost_infos->infos[i].index == index) {
      arrayindex = i;
      break;
    }
  }

  /* Do nothing if not found or we have the end array index. */
  if ((arrayindex != -1) && arrayindex != (ost_infos->fullcount - 1)) {

    /* Copy remaining infos down one place. Overlapping.. */
    memmove(&ost_infos->infos[arrayindex], &ost_infos->infos[arrayindex + 1],
            (ost_infos->fullcount - arrayindex - 1) *
                sizeof(struct swift_ost_info));
    if (arrayindex < ost_infos->count) ost_infos->count = ost_infos->count - 1;
    ost_infos->fullcount = ost_infos->fullcount - 1;

  } else if (arrayindex == ost_infos->fullcount - 1) {

    /* End array index, just adjust counts. */
    if (arrayindex < ost_infos->count) ost_infos->count = ost_infos->count - 1;
    ost_infos->fullcount = ost_infos->fullcount - 1;
  }
#endif
}

/**
 * @brief Create a file with a given OST index and number of OSTs to stripe.
 *
 * @param filename name of the file to create.
 * @param offset index of the first OST used with this file.
 * @param count number of OSTs to stripe this file over.
 * @param usedoffset the offset actually used by file.
 *
 * @return non-zero if there are problems creating the file.
 */
int swift_create_striped_file(const char *filename, int offset, int count,
                              int *usedoffset) {
  int rc = 0;

#ifdef HAVE_LUSTREAPI
  *usedoffset = offset;
  rc = llapi_file_create(filename, 0 /* Default block size */, offset, count,
                         LLAPI_LAYOUT_RAID0 /* Pattern default */);
  if (rc != 0) {
    rc = -rc;
    message("Cannot create file %s : %s", filename, strerror(rc));
  } else {

    /* Recover the file offset of first OST in case it is changed from
     * operational reasons. */
    /* Yuk, needs extra space for array os lov_user_ost_data. */
    size_t sizelum = sizeof(struct lov_user_md) +
                     LOV_MAX_STRIPE_COUNT * sizeof(struct lov_user_ost_data);
    struct lov_user_md *lum = (struct lov_user_md *)malloc(sizelum);

    rc = llapi_file_get_stripe(filename, lum);
    rc = -rc;
    if (rc == 0) {
      *usedoffset = lum->lmm_objects[0].l_ost_idx;
    } else {
      /* Shouldn't be fatal. */
      *usedoffset = offset;
    }
    free(lum);
  }
#endif
  return rc;
}

/**
 * @brief Scan for the available OSTs for a given file path.
 *
 * The OSTs will be sorted by free space on exit and may be further selected
 * to remove OSTs that are too full for use or cannot be written to. If too
 * many OSTs are rejected it will be considered to be a parameter error and
 * all OSTs, sorted by free space, will be returned.
 * We don't want to flood the OSTs with RPC calls so only one MPI rank
 * should make this call.
 *
 * @param ost_infos pointer to empty OST storage struct. Will contain the
 *                  selected free-space ordered OSTs found on exit. The OST
 *                  count will remain at zero if anything fails. Note this
 *                  will need to be freed as usual regardless.
 * @param filepath  path to a file on the lustre file system. Must not exist.
 *                  The containing directory must exist and be part of the
 *                  lustre file system.
 * @param minfree minimum free space to allow in MiB. -1 for use a guess based
 * on size of the current process, 0 to disable selection.
 * @param writetest whether to check if the OSTs are writable. If used
 *                  the path must be that of a non-existent file on the
 *                  file system that is writable by the process.
 * @param verbose if true information about the OSTs and the selections made
 *                will be output.
 *
 */
void swift_ost_select(struct swift_ost_store *ost_infos, const char *filepath,
                      int minfree, int writetest, int verbose) {

  /* Initialise the struct. */
  swift_ost_store_init(ost_infos);

  /* Get directory of filepath. */
  char *filepathc = strdup(filepath);
  char *dirp = dirname(filepathc);

  /* Scan for all OSTs. */
  int rc = swift_ost_scan(dirp, ost_infos);
  free(dirp);

  /* If does not succeed we do nothing, probably not a lustre mount. */
  if (rc == 0) {
    if (verbose) swift_ost_store_print(ost_infos, 1);

    /* Make a copy so we can undo any changes. */
    struct swift_ost_store ost_infos_full;
    swift_ost_store_copy(ost_infos, &ost_infos_full);

    /* Cull these so we do not use OSTs with too little free space. Also sorts
     * into most free space order. If given a value use that, otherwise we use
     * the resident set size of the process, dumps and restarts are always
     * smaller than that. */
    if (minfree != 0) {
      if (minfree < 0) {

        /* No guarantee this will work, hopefully will return 0 in those cases
         * and we do nothing. */
        long size, resident, shared, text, library, data, dirty;
        memuse_use(&size, &resident, &shared, &text, &data, &library, &dirty);

        /* KiB into MiB. */
        minfree = (int)(resident / 1024.0);
      }

      /* And cull and sort. */
      swift_ost_cull(ost_infos, minfree);
      if (verbose)
        message("Rejected %d OSTs using free space threshold %d (MiB)",
                ost_infos->fullcount - ost_infos->count, minfree);
    }

    if (writetest != 0) {
      /* Test writing to all OSTs and remove any that are not writable.  We do
       * this by creating our file on every OST and checking it was created on
       * it. */
      int usedindex = 0;
      int removed = 0;
      for (int i = ost_infos->count - 1; i >= 0; i--) {
        usedindex = ost_infos->infos[i].index;
        rc = swift_create_striped_file(filepath, ost_infos->infos[i].index, 1,
                                       &usedindex);

        if (rc != 0) {
          /* Failed so not likely to succeed next time. Probably file
           * exists, there is nothing we should do about that, the existing
           * stripe will be reused, along with the space of the existing file.
           */
          message("Failed testing file creation on OSTs, aborting test");
          break;
        }

        if (usedindex != ost_infos->infos[i].index) {
          /* Differing OST indices, so not what we asked for, bye. */
          swift_ost_remove(ost_infos, ost_infos->infos[i].index);
          removed++;
        }
        unlink(filepath);
      }
      if (verbose) message("Rejected %d OSTs as readonly", removed);
    }

    /* Safety first. If we have too few OSTs left after the above we will
     * make the choice to do nothing. */
    if ((ost_infos->fullcount * 0.25 > ost_infos->count) ||
        ost_infos->count < 2) {
      message("Too many OSTs have been rejected (%d of %d).",
              ost_infos->fullcount - ost_infos->count, ost_infos->fullcount);
      message("Assuming OST rejection is flawed and skipping.");
      swift_ost_store_copy(&ost_infos_full, ost_infos);

      /* Still good to use a sorted list. */
      swift_ost_cull(ost_infos, 0);
    }
    swift_ost_store_free(&ost_infos_full);
    if (verbose) swift_ost_store_print(ost_infos, 1);

  } else {

    /* If the scan failed we do nothing, this is probably not a lustre mount. */
    message("Lustre OST scan failed, is this a lustre mount?");
  }
}
