/*******************************************************************************
 * This file is part of SWIFT.
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
#include <config.h>

/* Some standard headers. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "runner_radiation_feedback.h"
#include "swift.h"

/*! Number of random draws in the sampled checks. */
#define N_DRAWS 1000000

/**
 * @brief Next value of a splitmix64 stream, mapped to [0, 1).
 *
 * @param state The stream state.
 * @return A number in [0, 1).
 */
static double test_uniform(uint64_t *state) {
  uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z = z ^ (z >> 31);
  return (double)(z >> 11) * (1.0 / 9007199254740992.0);
}

/**
 * @brief Is a stored double finite? Read through its bits, which -ffast-math
 * cannot fold.
 *
 * @param x Pointer to the value.
 * @return 1 if finite.
 */
__attribute__((noinline)) static int test_is_finite(const double *x) {
  uint64_t bits;
  memcpy(&bits, x, sizeof(bits));
  return ((bits >> 52) & 0x7FFULL) != 0x7FFULL;
}

/**
 * @brief Print one check and return 1 if it failed.
 *
 * @param name What the check tests.
 * @param value The measured value.
 * @param bar The bar.
 * @param must_exceed 1 if the value must exceed the bar, 0 if it must not.
 * @return 1 if the check failed.
 */
static int report(const char *name, const double value, const double bar,
                  const int must_exceed) {
  const int ok =
      test_is_finite(&value) && (must_exceed ? value > bar : value <= bar);
  printf("%-4s %-72s %12.5g (bar %s %g)\n", ok ? "ok" : "FAIL", name, value,
         must_exceed ? ">" : "<=", bar);
  return !ok;
}

/**
 * @brief The claim decision of the claim loop before the rule existed.
 *
 * @param budget Photons left.
 * @param cost Candidate cost.
 * @param u The cone's draw.
 * @return 1 if claimed.
 */
static int old_claim(const double budget, const double cost, const float u) {
  if (cost <= budget) return 1;
  const float proba = budget / cost;
  return u <= proba;
}

/**
 * @brief Claimed share of holdable unaffordable candidates under a rule.
 *
 * @param rule The #hii_cone_claim_rule.
 * @param budget Photons left.
 * @param cost Candidate cost, above @p budget.
 * @param rate Recombination rate of the candidate.
 * @param cone_rate Photon rate of the cone.
 * @param fronts (return) Share of draws giving a front.
 * @return Share of draws giving a claim.
 */
static double claimed_share(const int rule, const double budget,
                            const double cost, const double rate,
                            const double cone_rate, double *fronts) {
  uint64_t state = 12345;
  long n_claimed = 0, n_front = 0;
  for (long i = 0; i < N_DRAWS; i++) {
    const float u = (float)test_uniform(&state);
    const enum hii_candidate_outcome o =
        runner_hii_candidate_outcome(rule, 0, budget, cost, rate, cone_rate, u);
    if (o == hii_candidate_claimed) n_claimed++;
    if (o == hii_candidate_front) n_front++;
  }
  *fronts = (double)n_front / N_DRAWS;
  return (double)n_claimed / N_DRAWS;
}

/**
 * @brief Checks of the claim decision.
 *
 * @return Number of failed checks.
 */
static int test_decision(void) {
  int fail = 0;
  uint64_t state = 777;

  /* The walk-on rule reproduces the old decision, draw for draw. */
  long differ = 0, passed_bad = 0;
  for (long i = 0; i < N_DRAWS; i++) {
    const double budget = pow(10.0, 60.0 * test_uniform(&state) - 10.0);
    const double cost = pow(10.0, 60.0 * test_uniform(&state) - 10.0);
    const float u = (float)test_uniform(&state);
    const enum hii_candidate_outcome o = runner_hii_candidate_outcome(
        hii_cone_claim_walk_on, 0, budget, cost, cost, budget, u);
    if ((o == hii_candidate_claimed) != old_claim(budget, cost, u)) differ++;
    if (o == hii_candidate_front) passed_bad++;
  }
  fail += report("walk-on rule: decisions differing from the old loop",
                 (double)differ, 0., 0);
  fail += report("walk-on rule: fronts returned (must be none)",
                 (double)passed_bad, 0., 0);

  /* Control: the stop rule with gas never holdable differs from the old
     decision. */
  long differ_stop = 0;
  for (long i = 0; i < N_DRAWS; i++) {
    const double budget = pow(10.0, 60.0 * test_uniform(&state) - 10.0);
    const double cost = pow(10.0, 60.0 * test_uniform(&state) - 10.0);
    const float u = (float)test_uniform(&state);
    const enum hii_candidate_outcome o = runner_hii_candidate_outcome(
        hii_cone_claim_stop, 0, budget, cost, 2.0, 1.0, u);
    if ((o == hii_candidate_claimed) != old_claim(budget, cost, u))
      differ_stop++;
  }
  fail += report("control: stop rule differs from the old loop (must)",
                 (double)differ_stop, 0., 1);

  /* Affordable gas is claimed by every rule, also at cost == budget. */
  long refused = 0;
  for (int rule = 0; rule < hii_cone_claim_rule_count; rule++) {
    for (long i = 0; i < 1000; i++) {
      const double budget = 1.0 + test_uniform(&state);
      const double cost = budget * (i == 0 ? 1.0 : test_uniform(&state));
      if (runner_hii_candidate_outcome(rule, 0, budget, cost, 1e30, 1.0,
                                       0.99f) != hii_candidate_claimed)
        refused++;
    }
  }
  fail +=
      report("affordable gas not claimed, all rules", (double)refused, 0., 0);

  /* The deterministic boundary claims every candidate. */
  long det = 0;
  for (int rule = 0; rule < hii_cone_claim_rule_count; rule++)
    if (runner_hii_candidate_outcome(rule, 1, 1.0, 10.0, 1e30, 1.0, 0.99f) !=
        hii_candidate_claimed)
      det++;
  fail += report("deterministic boundary: candidates not claimed", (double)det,
                 0., 0);

  /* The holdable gate: no roll on gas that costs more upkeep than the cone
     receives, even with the luckiest draw. */
  long gate_claims[hii_cone_claim_rule_count] = {0};
  for (int rule = 0; rule < hii_cone_claim_rule_count; rule++) {
    const double rate_ratio[3] = {1.0, 1.5, 1e4};
    for (int k = 0; k < 3; k++)
      if (runner_hii_candidate_outcome(rule, 0, 1.0, 1.0001, rate_ratio[k], 1.0,
                                       0.f) == hii_candidate_claimed)
        gate_claims[rule]++;
  }
  fail += report("stop rule: gas not holdable claimed with draw 0",
                 (double)gate_claims[1], 0., 0);
  fail += report("shadow rule: gas not holdable claimed with draw 0",
                 (double)gate_claims[2], 0., 0);
  fail += report("control: walk-on rule claims gas not holdable (must)",
                 (double)gate_claims[0], 0., 1);

  /* On holdable gas the roll keeps the expected spend at the budget: the
     claimed share is budget / cost, within 5 standard errors. */
  const double budget = 0.3, cost = 1.0;
  const double p = budget / cost;
  const double sigma = sqrt(p * (1. - p) / N_DRAWS);
  for (int rule = 0; rule < hii_cone_claim_rule_count; rule++) {
    double fronts = 0.;
    const double share = claimed_share(rule, budget, cost, 0.5, 1.0, &fronts);
    char name[128];
    snprintf(name, sizeof(name),
             "rule %d, holdable gas: |claimed share - budget/cost| / sigma",
             rule);
    fail += report(name, fabs(share - p) / sigma, 5., 0);
    snprintf(name, sizeof(name),
             "rule %d, holdable gas: |front share - expected| / sigma", rule);
    const double expect_front = rule == 0 ? 0. : 1. - p;
    fail += report(name, fabs(fronts - expect_front) / sigma, 5., 0);
  }
  /* Control: the same statistic against a wrong share must fail. */
  {
    double fronts = 0.;
    const double share = claimed_share(1, budget, cost, 0.5, 1.0, &fronts);
    fail += report("control: claimed share against 1.1 budget/cost (must)",
                   fabs(share - 1.1 * p) / sigma, 5., 1);
  }

  return fail;
}

/**
 * @brief Checks of the absorbed share at a front.
 *
 * @return Number of failed checks.
 */
static int test_absorbed(void) {
  int fail = 0;

  /* Blocking the 64 sub-cones one at a time leaves (64 - k) / 64 of the
     budget after k fronts, and exactly 0 after the last. */
  double budget = 3.7e58;
  const double b0 = budget;
  double worst = 0.;
  int n_open = hii_shadow_children;
  for (int k = 1; k <= hii_shadow_children; k++) {
    budget -= runner_hii_shadow_absorbed(budget, 1, n_open);
    n_open--;
    const double expect = b0 * (hii_shadow_children - k) / hii_shadow_children;
    if (k < hii_shadow_children)
      worst = fmax(worst, fabs(budget - expect) / expect);
  }
  fail += report("one sub-cone per front: worst relative error of the rest",
                 worst, 1e-12, 0);
  fail +=
      report("budget left after the last sub-cone closes", fabs(budget), 0., 0);
  fail += report("no sub-cone blocked: photons absorbed",
                 runner_hii_shadow_absorbed(5.0, 0, 10), 0., 0);
  fail += report("several sub-cones at once: |absorbed - 3/10 budget| / budget",
                 fabs(runner_hii_shadow_absorbed(5.0, 3, 10) - 1.5) / 5.0,
                 1e-15, 0);
  return fail;
}

#if defined(HAVE_CHEALPIX)
/**
 * @brief Draw a direction uniform on the sphere.
 *
 * @param state The stream state.
 * @param v (return) The unit direction.
 */
static void test_direction(uint64_t *state, double v[3]) {
  const double z = 2.0 * test_uniform(state) - 1.0;
  const double phi = 2.0 * M_PI * test_uniform(state);
  const double s = sqrt(1.0 - z * z);
  v[0] = s * cos(phi);
  v[1] = s * sin(phi);
  v[2] = z;
}

/**
 * @brief Checks of the sub-cone layout and of the blocking geometry.
 *
 * @return Number of failed checks.
 */
static int test_geometry(void) {
  int fail = 0;

  for (int nside = 1; nside <= 2; nside++) {
    const int n_pix = 12 * nside * nside;
    char name[128];

    /* Every sub-cone centre lies in its own cone. */
    long outside = 0;
    for (int p = 0; p < n_pix; p++) {
      struct hii_shadow_cone cone;
      runner_hii_shadow_cone_reset(&cone);
      runner_hii_shadow_fill_centres(&cone, nside, p);
      for (int k = 0; k < hii_shadow_children; k++) {
        const double v[3] = {cone.centre[k][0], cone.centre[k][1],
                             cone.centre[k][2]};
        int64_t ring = 0;
        vec2pix_ring64(nside, v, &ring);
        int pixel = -1, child = -1;
        runner_hii_shadow_locate(v, nside, &pixel, &child);
        if (ring != p || pixel != p || child != k) outside++;
      }
    }
    snprintf(name, sizeof(name),
             "nside %d: sub-cone centres outside their cone or index", nside);
    fail += report(name, (double)outside, 0., 0);

    /* The cone of the nested lookup is the cone the claim loop uses. */
    uint64_t state = 99 + nside;
    long mismatch = 0;
    for (long i = 0; i < N_DRAWS; i++) {
      double v[3];
      test_direction(&state, v);
      int pixel = -1, child = -1;
      runner_hii_shadow_locate(v, nside, &pixel, &child);
      const float dx[3] = {(float)-v[0], (float)-v[1], (float)-v[2]};
      if (pixel != runner_hii_get_pixel(dx, n_pix, NULL)) mismatch++;
    }
    snprintf(name, sizeof(name),
             "nside %d: nested cone differs from the claim-loop cone, 1e6 dirs",
             nside);
    fail += report(name, (double)mismatch, 0., 0);

    /* The same through the rotated frame of a star. */
    double R[9];
    runner_hii_cone_frame_from_id(4242, 0, R);
    const struct hii_cone_geometry geom = {
        .frame = R, .shadow = NULL, .nside = nside};
    long mismatch_rot = 0;
    for (long i = 0; i < N_DRAWS; i++) {
      double v[3];
      test_direction(&state, v);
      const float dx[3] = {(float)-v[0], (float)-v[1], (float)-v[2]};
      int pixel = -1, child = -1;
      runner_hii_shadow_direction(dx, &geom, &pixel, &child);
      if (pixel != runner_hii_get_pixel(dx, n_pix, R)) mismatch_rot++;
    }
    snprintf(name, sizeof(name),
             "nside %d, rotated frame: nested cone differs, 1e6 dirs", nside);
    fail += report(name, (double)mismatch_rot, 0., 0);
  }

  /* Blocking: a front covering angle a around a sub-cone centre blocks the
     sub-cones whose centres lie within a, counted independently from the
     angular coordinates of the centres. */
  const int nside = 1;
  long wrong = 0, total = 0;
  const double angles[5] = {0.0, 0.05, 0.2, 0.4, 0.8};
  for (int p = 0; p < 12; p++) {
    for (int k = 0; k < hii_shadow_children; k += 7) {
      for (int a = 0; a < 5; a++) {
        struct hii_shadow_cone cone;
        runner_hii_shadow_cone_reset(&cone);
        runner_hii_shadow_fill_centres(&cone, nside, p);
        int64_t parent = 0;
        ring2nest64(nside, p, &parent);
        double th0, ph0;
        pix2ang_nest64(8 * nside, parent * hii_shadow_children + k, &th0, &ph0);
        const double dir[3] = {sin(th0) * cos(ph0), sin(th0) * sin(ph0),
                               cos(th0)};
        const int m = runner_hii_shadow_block(&cone, dir, k, cos(angles[a]));
        /* Centres within 1e-6 of the edge may fall either side. */
        int n_in = 0, n_edge = 0;
        for (int j = 0; j < hii_shadow_children; j++) {
          double th, ph;
          pix2ang_nest64(8 * nside, parent * hii_shadow_children + j, &th, &ph);
          const double c =
              sin(th) * sin(th0) * cos(ph - ph0) + cos(th) * cos(th0);
          const double sep = acos(fmin(1.0, c));
          if (j == k || sep <= angles[a] * (1. - 1e-6))
            n_in++;
          else if (sep <= angles[a] * (1. + 1e-6) + 1e-9)
            n_edge++;
        }
        total++;
        if (m < n_in || m > n_in + n_edge) wrong++;
        if (cone.n_open != hii_shadow_children - m) wrong++;
      }
    }
  }
  fail += report("blocked sub-cone count differs from the independent count",
                 (double)wrong, 0., 0);

  /* Control: the independent count with twice the angle must disagree. */
  long wrong_ctl = 0;
  for (int k = 0; k < hii_shadow_children; k += 7) {
    struct hii_shadow_cone cone;
    runner_hii_shadow_cone_reset(&cone);
    runner_hii_shadow_fill_centres(&cone, 1, 4);
    const double dir[3] = {cone.centre[k][0], cone.centre[k][1],
                           cone.centre[k][2]};
    const int m = runner_hii_shadow_block(&cone, dir, k, cos(0.2));
    int n_wide = 0;
    for (int j = 0; j < hii_shadow_children; j++) {
      const double c = dir[0] * cone.centre[j][0] + dir[1] * cone.centre[j][1] +
                       dir[2] * cone.centre[j][2];
      if (j == k || c >= cos(0.4)) n_wide++;
    }
    if (m != n_wide) wrong_ctl++;
  }
  fail += report("control: counts at twice the angle disagree (must)",
                 (double)wrong_ctl, 0., 1);
  printf("     (%ld fronts tested)\n", total);

  /* A second front over the same sub-cones blocks nothing new. */
  {
    struct hii_shadow_cone cone;
    runner_hii_shadow_cone_reset(&cone);
    runner_hii_shadow_fill_centres(&cone, 1, 4);
    const double dir[3] = {cone.centre[10][0], cone.centre[10][1],
                           cone.centre[10][2]};
    const int m1 = runner_hii_shadow_block(&cone, dir, 10, cos(0.3));
    const int m2 = runner_hii_shadow_block(&cone, dir, 10, cos(0.3));
    fail += report("first front blocks sub-cones (must)", (double)m1, 0., 1);
    fail +=
        report("same front again: sub-cones newly blocked", (double)m2, 0., 0);
    const int m3 = runner_hii_shadow_block(&cone, dir, 10, -1.0);
    fail += report("front covering the sphere: sub-cones left open",
                   (double)cone.n_open, 0., 0);
    fail += report("front covering the sphere: |m1 + m3 - 64|",
                   fabs((double)(m1 + m3 - hii_shadow_children)), 0., 0);
  }

  return fail;
}

/**
 * @brief Checks of the holdable gate under the shadow mask: only the open
 * sub-cones of the candidate's footprint feed it.
 *
 * @return Number of failed checks.
 */
static int test_open_gate(void) {
  int fail = 0;
  const double cone_rate = 1.0;
  const int own = 27;

  struct hii_shadow_cone cone;
  runner_hii_shadow_cone_reset(&cone);
  runner_hii_shadow_fill_centres(&cone, 1, 4);
  const double dir[3] = {cone.centre[own][0], cone.centre[own][1],
                         cone.centre[own][2]};

  /* A far candidate covers only its own sub-cone. Its upkeep of half the
     cone's rate exceeds the 1/64 that reaches it: a front, whatever the
     draw. The whole-cone gate lets the luckiest draw claim it. */
  const uint64_t narrow = runner_hii_shadow_footprint(&cone, dir, own, 1.0);
  const double share_narrow = runner_hii_shadow_open_share(&cone, narrow);
  fail += report("narrow footprint: |open share - 1/64|",
                 fabs(share_narrow - 1.0 / hii_shadow_children), 0., 0);
  long narrow_claims = 0;
  for (int rule = 1; rule < hii_cone_claim_rule_count; rule++)
    if (runner_hii_candidate_outcome(rule, 0, 1.0, 1.0001, 0.5,
                                     share_narrow * cone_rate,
                                     0.f) != hii_candidate_front)
      narrow_claims++;
  fail += report("narrow footprint, upkeep 0.5: not a front with draw 0",
                 (double)narrow_claims, 0., 0);
  fail += report("control: whole-cone gate claims it with draw 0 (must)",
                 (double)(runner_hii_candidate_outcome(
                              hii_cone_claim_shadow_mask, 0, 1.0, 1.0001, 0.5,
                              cone_rate, 0.f) == hii_candidate_claimed),
                 0., 1);

  /* Partly blocked footprint: only its open sub-cones count, against an
     independent count over the centres. */
  const double wide = cos(0.3);
  const double side[3] = {cone.centre[own + 1][0], cone.centre[own + 1][1],
                          cone.centre[own + 1][2]};
  runner_hii_shadow_block(&cone, side, own + 1, cos(0.15));
  const uint64_t fp = runner_hii_shadow_footprint(&cone, dir, own, wide);
  int n_cover = 0, n_open_cover = 0;
  for (int j = 0; j < hii_shadow_children; j++) {
    const double c = dir[0] * cone.centre[j][0] + dir[1] * cone.centre[j][1] +
                     dir[2] * cone.centre[j][2];
    if (j != own && c < wide) continue;
    n_cover++;
    if (!runner_hii_shadow_is_blocked(&cone, j)) n_open_cover++;
  }
  const double share_fp = runner_hii_shadow_open_share(&cone, fp);
  fail += report("partly blocked footprint: own sub-cone blocked",
                 (double)runner_hii_shadow_is_blocked(&cone, own), 0., 0);
  fail += report("partly blocked footprint: sub-cones blocked (must)",
                 (double)(n_cover - n_open_cover), 0., 1);
  fail += report("partly blocked footprint: |open share - independent|",
                 fabs(share_fp * hii_shadow_children - n_open_cover), 0., 0);

  /* Upkeep between the open and the whole footprint share: a front. Below
     the open share: rolled, so the luckiest draw claims it. */
  const double upkeep_mid =
      0.5 * (n_open_cover + n_cover) / (double)hii_shadow_children;
  const double upkeep_low = 0.5 * share_fp;
  fail += report(
      "upkeep above the open share: claimed with draw 0",
      (double)(runner_hii_candidate_outcome(
                   hii_cone_claim_shadow_mask, 0, 1.0, 1.0001, upkeep_mid,
                   share_fp * cone_rate, 0.f) != hii_candidate_front),
      0., 0);
  fail += report(
      "upkeep below the open share: not claimed with draw 0",
      (double)(runner_hii_candidate_outcome(
                   hii_cone_claim_shadow_mask, 0, 1.0, 1.0001, upkeep_low,
                   share_fp * cone_rate, 0.f) != hii_candidate_claimed),
      0., 0);

  /* A footprint over the whole sphere sees the cone's open share. */
  const uint64_t all = runner_hii_shadow_footprint(&cone, dir, own, -1.0);
  fail += report("whole footprint: |open share - n_open / 64|",
                 fabs(runner_hii_shadow_open_share(&cone, all) -
                      (double)cone.n_open / hii_shadow_children),
                 0., 0);

  return fail;
}
#endif /* HAVE_CHEALPIX */

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  int fail = 0;
  fail += test_decision();
  fail += test_absorbed();
#if defined(HAVE_CHEALPIX)
  fail += test_geometry();
  fail += test_open_gate();
#else
  printf("skip geometry checks: HAVE_CHEALPIX is not defined\n");
#endif

  if (fail) printf("%d check(s) failed\n", fail);
  return fail ? 1 : 0;
}
