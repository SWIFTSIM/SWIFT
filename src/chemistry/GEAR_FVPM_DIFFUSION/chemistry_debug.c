/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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

#if defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) || \
    defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)

#ifdef GEAR_FVPM_DIFF_DEBUG_PAIR_VISIT_COUNT

#include "chemistry_debug.h"
#include "error.h"
#include "part.h"

long long chemistry_fvpm_visit_acted = 0;
long long chemistry_fvpm_visit_skipped = 0;

/* Every visit that could cover a pair (SYM, ACT) plus every SKIP, keyed by
   particle address pair. Post-processed once at exit
   (chemistry_fvpm_visit_log_analyze()) into a dropped/double_applied
   verdict: a pair with zero covering records was dropped; two or more was
   double-applied. Not thread-safe. */
struct chemistry_fvpm_visit_record {
  const struct part *a; /* min(pi, pj) by address */
  const struct part *b; /* max(pi, pj) by address */
  enum chemistry_fvpm_visit_arm arm;
};

static struct chemistry_fvpm_visit_record *chemistry_fvpm_visit_log = NULL;
static size_t chemistry_fvpm_visit_log_count = 0;
static size_t chemistry_fvpm_visit_log_capacity = 0;

void chemistry_fvpm_visit_log_record(const struct part *pi,
                                     const struct part *pj,
                                     enum chemistry_fvpm_visit_arm arm) {
  if (chemistry_fvpm_visit_log_count == chemistry_fvpm_visit_log_capacity) {
    chemistry_fvpm_visit_log_capacity =
        chemistry_fvpm_visit_log_capacity
            ? chemistry_fvpm_visit_log_capacity * 2
            : 65536;
    chemistry_fvpm_visit_log = (struct chemistry_fvpm_visit_record *)realloc(
        chemistry_fvpm_visit_log,
        chemistry_fvpm_visit_log_capacity *
            sizeof(struct chemistry_fvpm_visit_record));
    if (chemistry_fvpm_visit_log == NULL)
      error("Failed to grow the chemistry_fvpm_visit_log diagnostic buffer.");
  }
  struct chemistry_fvpm_visit_record *r =
      &chemistry_fvpm_visit_log[chemistry_fvpm_visit_log_count++];
  if (pi < pj) {
    r->a = pi;
    r->b = pj;
  } else {
    r->a = pj;
    r->b = pi;
  }
  r->arm = arm;
}

static int chemistry_fvpm_visit_record_cmp(const void *x, const void *y) {
  const struct chemistry_fvpm_visit_record *rx =
      (const struct chemistry_fvpm_visit_record *)x;
  const struct chemistry_fvpm_visit_record *ry =
      (const struct chemistry_fvpm_visit_record *)y;
  if (rx->a != ry->a) return (rx->a < ry->a) ? -1 : 1;
  if (rx->b != ry->b) return (rx->b < ry->b) ? -1 : 1;
  return 0;
}

void chemistry_fvpm_visit_log_reset(void) {
  chemistry_fvpm_visit_log_count = 0;
}

void chemistry_fvpm_visit_log_analyze(void) {
  if (chemistry_fvpm_visit_log_count == 0) return;
  qsort(chemistry_fvpm_visit_log, chemistry_fvpm_visit_log_count,
        sizeof(struct chemistry_fvpm_visit_record),
        chemistry_fvpm_visit_record_cmp);
  long long distinct_pairs = 0, dropped = 0, double_applied = 0;
  size_t i = 0;
  while (i < chemistry_fvpm_visit_log_count) {
    size_t j = i;
    int covering = 0;
    while (j < chemistry_fvpm_visit_log_count &&
           chemistry_fvpm_visit_log[j].a == chemistry_fvpm_visit_log[i].a &&
           chemistry_fvpm_visit_log[j].b == chemistry_fvpm_visit_log[i].b) {
      if (chemistry_fvpm_visit_log[j].arm == CHEMISTRY_FVPM_VISIT_SYM ||
          chemistry_fvpm_visit_log[j].arm == CHEMISTRY_FVPM_VISIT_ACT)
        covering++;
      j++;
    }
    distinct_pairs++;
    if (covering == 0) dropped++;
    if (covering >= 2) double_applied++;
    i = j;
  }
  message(
      "GEAR_FVPM_DIFF_DEBUG_PAIR_VISIT_COUNT: per-pair log: distinct "
      "pairs=%lld dropped=%lld double_applied=%lld (both must be zero)",
      distinct_pairs, dropped, double_applied);
}

void chemistry_fvpm_print_visit_counts(void) {
  message(
      "GEAR_FVPM_DIFF_DEBUG_PAIR_VISIT_COUNT: mixed-band arm acted=%lld "
      "skipped=%lld",
      chemistry_fvpm_visit_acted, chemistry_fvpm_visit_skipped);
  chemistry_fvpm_visit_log_analyze();
}

#endif /* GEAR_FVPM_DIFF_DEBUG_PAIR_VISIT_COUNT */
#endif /* GEAR FVPM diffusion chemistry */
