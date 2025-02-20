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
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local includes. */
#include "swift_lustre_api.h"

/* Lustre API */
#include <lustre/lustre_user.h>
#include <lustre/lustreapi.h>

/* Number of OSTs to pre-allocate space for. */
#define PREALLOC 100

/* Bytes in a TiB */
#define TiB (1024.0 * 1024.0 * 1024.0)

/**
 * @brief Initialize an OST scan storage structure.
 *
 * @param swift_ost_infos pointer to the storage structure.
 */
void swift_ost_store_init(struct swift_ost_store *ost_infos) {
  ost_infos->count = 0;
  ost_infos->fullcount = 0;
  ost_infos->size = PREALLOC;
  ost_infos->infos =
      (struct swift_ost_info *)malloc(sizeof(struct swift_ost_info) * PREALLOC);
  memset(ost_infos->infos, 0, sizeof(struct swift_ost_info) * PREALLOC);
}

/**
 * @brief Release any storage associated with an OST scan storage structure.
 *
 * @param ost_infos pointer to the storage structure.
 */
void swift_ost_store_free(struct swift_ost_store *ost_infos) {
  free(ost_infos->infos);
  ost_infos->infos = NULL;
  ost_infos->count = 0;
  ost_infos->fullcount = 0;
  ost_infos->size = 0;
}

/**
 * @brief Print an OST storage structure for debugging purposes.
 *
 * @param ost_infos pointer to the storage structure.
 */
void swift_ost_store_print(struct swift_ost_store *ost_infos) {
  printf("#  Listing of OSTs. Using %d of %d\n", ost_infos->count,
         ost_infos->fullcount);

  printf("%5s %21s %21s %21s\n", "Index", "Size (bytes)", "Used (bytes)",
         "Free (bytes)");
  size_t ssum = 0;
  size_t usum = 0;
  size_t smin = ost_infos->infos[0].size;
  size_t smax = 0;
  size_t umin = ost_infos->infos[0].used;
  size_t umax = 0;

  for (int i = 0; i < ost_infos->count; i++) {
    printf("%5d %21zd %21zd %21zd\n", ost_infos->infos[i].index,
           ost_infos->infos[i].size, ost_infos->infos[i].used,
           ost_infos->infos[i].size - ost_infos->infos[i].used);

    ssum += ost_infos->infos[i].size;
    usum += ost_infos->infos[i].used;

    if (ost_infos->infos[i].size > smax) smax = ost_infos->infos[i].size;
    if (ost_infos->infos[i].size < smin) smin = ost_infos->infos[i].size;

    if (ost_infos->infos[i].used > umax) umax = ost_infos->infos[i].used;
    if (ost_infos->infos[i].used < umin) umin = ost_infos->infos[i].used;
  }
  printf("# Filesystem size:%.2f TiB used:%.2f TiB free:%.2f TiB %.2f%%\n",
         ssum / TiB, usum / TiB, (ssum - usum) / TiB,
         100.0 * (double)(ssum - usum) / (double)ssum);
  printf("# Min/max size: %.2f/%.2f TiB Min/max used: %.2f/%.2f TiB\n",
         smin / TiB, smax / TiB, umin / TiB, umax / TiB);
}

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
 * @return 0 on success, otherwise an error will have been reported.
 */
int swift_ost_scan(const char *path, struct swift_ost_store *ost_infos) {

  char mntdir[PATH_MAX] = "";
  char fsname[PATH_MAX] = "";
  char cpath[PATH_MAX] = "";
  int rc = 0;

  /* Check this path exists. */
  if (!realpath(path, cpath)) {
    rc = -errno;
    fprintf(stderr, "Error: not a real path '%s': %s\n", path, strerror(-rc));
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
          if (rc == -ENODEV || rc == -EAGAIN || rc == -EINVAL ||
              rc == -EFAULT) {
            /* Nothing we can query here, so time to stop search. */
            break;
          }

          /* Inactive devices are empty. */
          if (rc == -ENODATA) {
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
        fprintf(stderr, "Error: no lustre mount point found for: %s\n", path);
        rc = 1;
      }
    } else {
      fprintf(stderr, "Error: failed to locate a lustre mount point for: %s\n",
              path);
      rc = 1;
    }
  }
  return rc;
}

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

/**
 * @brief Sort the OSTs into decreasing free space culling those that do not
 * meet a free space threshold.
 *
 * @param ost_infos pointer to populated storage structure.
 * @param minfree the number of bytes that the OST show be capable of
 *                storing. Zero for no effect.
 */
void swift_ost_cull(struct swift_ost_store *ost_infos, size_t minfree) {

  /* Sort by free space. */
  qsort(ost_infos->infos, ost_infos->count, sizeof(struct swift_ost_info),
        ostcmp);

  /* And cull if needed. */
  if (minfree > 0) {

    /* Always keep at least one! */
    for (int i = 1; i < ost_infos->count; i++) {
      struct swift_ost_info *curr = &ost_infos->infos[i];
      if ((curr->size - curr->used) < minfree) {

        /* Throw the rest away. Note fullcount now decoupled. */
        ost_infos->count = i;
      }
    }
  }
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
  int index = (*arrayindex % ost_infos->count);
  *arrayindex = index + count;
  return ost_infos->infos[index].index;
}

/**
 * @brief Remove an OST by index from the store.
 *
 * @param ost_infos pointer to populated storage structure.
 * @param index index of the OST to remove.
 */
void swift_ost_remove(struct swift_ost_store *ost_infos, int index) {

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

  *usedoffset = offset;
  int rc = llapi_file_create(filename, 0 /* Default block size */, offset,
                             count, LLAPI_LAYOUT_RAID0 /* Pattern default */);
  if (rc != 0) {
    fprintf(stderr, "Error: cannot create file %s : %s\n", filename,
            strerror(-rc));
  } else {

    /* Recover the file offset of first OST in case it is changed from
     * operational reasons. */
    /* Yuk, needs extra space for array os lov_user_ost_data. */
    size_t sizelum = sizeof(struct lov_user_md) +
                     LOV_MAX_STRIPE_COUNT * sizeof(struct lov_user_ost_data);
    struct lov_user_md *lum = (struct lov_user_md *)malloc(sizelum);

    rc = llapi_file_get_stripe(filename, lum);
    *usedoffset = lum->lmm_objects[0].l_ost_idx;
    free(lum);
  }
  return rc;
}
