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

/* The checks of the photon budget need the GEAR feedback module with its HII
   radiation. */
#if defined(FEEDBACK_GEAR) && defined(GEAR_SUBGRID_RADIATION_HII)
#define TEST_HAS_HII_FEEDBACK
#endif

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
 * @brief The plain roll: claimed if affordable or on a won draw.
 *
 * @param budget Photons left.
 * @param cost Candidate cost.
 * @param u The cone's draw.
 * @return 1 if claimed.
 */
static int plain_roll_claim(const double budget, const double cost,
                            const float u) {
  if (cost <= budget) return 1;
  const float proba = budget / cost;
  return u <= proba;
}

/**
 * @brief The decision as it was when gas whose recombination rate reached the
 * cone's photon rate was never rolled for. Only the control checks use it.
 *
 * @param budget Photons left.
 * @param cost Candidate cost.
 * @param recombination_rate Recombination rate of the candidate.
 * @param photon_rate Photon rate of the cone.
 * @param u The cone's draw.
 * @return 1 if claimed.
 */
static int refusing_claim(const double budget, const double cost,
                          const double recombination_rate,
                          const double photon_rate, const float u) {
  if (cost <= budget) return 1;
  if (recombination_rate >= photon_rate) return 0;
  return u <= (float)(budget / cost);
}

/**
 * @brief Claimed share of unaffordable candidates over many draws.
 *
 * @param budget Photons left.
 * @param cost Candidate cost, above @p budget.
 * @param refuse 1 to use the refusing copy, with a gas recombination rate
 * above the cone's photon rate.
 * @param fronts (return) Share of draws giving a front.
 * @return Share of draws giving a claim.
 */
static double claimed_share(const double budget, const double cost,
                            const int refuse, double *fronts) {
  uint64_t state = 12345;
  long n_claimed = 0, n_front = 0;
  for (long i = 0; i < N_DRAWS; i++) {
    const float u = (float)test_uniform(&state);
    int claimed;
    if (refuse)
      claimed = refusing_claim(budget, cost, 8.0, 1.0, u);
    else
      claimed = runner_hii_candidate_outcome(0, budget, cost, u) ==
                hii_candidate_claimed;
    if (claimed)
      n_claimed++;
    else
      n_front++;
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

  /* A claim is decided by the plain roll, draw for draw; what is not claimed
     is a front. */
  long differ = 0, passed_bad = 0;
  for (long i = 0; i < N_DRAWS; i++) {
    const double budget = pow(10.0, 60.0 * test_uniform(&state) - 10.0);
    const double cost = pow(10.0, 60.0 * test_uniform(&state) - 10.0);
    const float u = (float)test_uniform(&state);
    const enum hii_candidate_outcome o =
        runner_hii_candidate_outcome(0, budget, cost, u);
    if ((o == hii_candidate_claimed) != plain_roll_claim(budget, cost, u))
      differ++;
    if (o == hii_candidate_passed) passed_bad++;
  }
  fail += report("claims differing from the plain roll", (double)differ, 0., 0);
  fail += report("candidates passed over (must be none)", (double)passed_bad,
                 0., 0);

  /* Affordable gas is claimed, also at cost == budget. */
  long refused = 0;
  for (long i = 0; i < 1000; i++) {
    const double budget = 1.0 + test_uniform(&state);
    const double cost = budget * (i == 0 ? 1.0 : test_uniform(&state));
    if (runner_hii_candidate_outcome(0, budget, cost, 0.99f) !=
        hii_candidate_claimed)
      refused++;
  }
  fail += report("affordable gas not claimed", (double)refused, 0., 0);

  /* The deterministic boundary claims every candidate. */
  fail += report("deterministic boundary: candidate not claimed",
                 (double)(runner_hii_candidate_outcome(1, 1.0, 10.0, 0.99f) !=
                          hii_candidate_claimed),
                 0., 0);

  /* A cone whose Stromgren mass is below one particle (budget 5 per cent of
     the cost, gas recombining faster than the cone emits) claims on a won
     roll: the claimed share is budget / cost, within 5 standard errors. */
  const double budget = 0.05, cost = 1.0;
  const double p = budget / cost;
  const double sigma = sqrt(p * (1. - p) / N_DRAWS);
  {
    double fronts = 0.;
    const double share = claimed_share(budget, cost, 0, &fronts);
    fail += report("sub-particle cone: |claimed share - budget/cost| / sigma",
                   fabs(share - p) / sigma, 5., 0);
    fail +=
        report("sub-particle cone: |front share - (1 - budget/cost)| / sigma",
               fabs(fronts - (1. - p)) / sigma, 5., 0);
  }

  /* Controls: a refusing copy and a wrong share must fail the same statistic.
   */
  {
    double fronts = 0.;
    double share = claimed_share(budget, cost, 1, &fronts);
    fail += report(
        "control, refusing copy: |claimed share - budget/cost| / sigma (must)",
        fabs(share - p) / sigma, 5., 1);
    share = claimed_share(budget, cost, 0, &fronts);
    fail += report("control: claimed share against 1.1 budget/cost (must)",
                   fabs(share - 1.1 * p) / sigma, 5., 1);
  }

  return fail;
}

/**
 * @brief The renewal decision as it was when a part-paid upkeep always
 * renewed. Only the control check uses it.
 *
 * @param paid Photons charged.
 * @param owed Photons owed.
 * @param u The draw.
 * @return 1, always.
 */
static int certain_renewal(const double paid, const double owed,
                           const double u) {
  (void)paid;
  (void)owed;
  (void)u;
  return 1;
}

/**
 * @brief Checks of the renewal of a held particle whose upkeep is part paid.
 *
 * @return Number of failed checks.
 */
static int test_renewal(void) {
  int fail = 0;
  const double owed = 1.0, paid = 0.13;
  const double sigma = sqrt(paid * (1. - paid) / N_DRAWS);
  uint64_t state = 4242;
  long n_ok = 0, n_certain = 0;
  for (long i = 0; i < N_DRAWS; i++) {
    const double u = test_uniform(&state);
    n_ok += runner_hii_upkeep_renewed(paid, owed, u);
    n_certain += certain_renewal(paid, owed, u);
  }
  fail += report("renewal: |renewed share - paid/owed| / sigma",
                 fabs((double)n_ok / N_DRAWS - paid) / sigma, 5., 0);
  fail += report(
      "control, certain renewal: |renewed share - paid/owed| / sigma (must)",
      fabs((double)n_certain / N_DRAWS - paid) / sigma, 5., 1);

  /* Fully paid or overpaid upkeep renews for every draw, also the largest. */
  const double draws[3] = {0., 0.5, 0.9999999999};
  long lost = 0;
  for (int k = 0; k < 3; k++) {
    lost += !runner_hii_upkeep_renewed(owed, owed, draws[k]);
    lost += !runner_hii_upkeep_renewed(1.5 * owed, owed, draws[k]);
  }
  fail += report("fully paid upkeep not renewed", (double)lost, 0., 0);
  fail +=
      report("control: part-paid upkeep with the largest draw lost (must)",
             (double)!runner_hii_upkeep_renewed(paid, owed, draws[2]), 0., 1);
  return fail;
}

#ifdef TEST_HAS_HII_FEEDBACK
/*! Number of checks in the no-debt sequence. */
#define N_NO_DEBT_CHECKS 6

/*! Names of the no-debt checks, in the order of the sequence. */
static const char *no_debt_names[N_NO_DEBT_CHECKS] = {
    "first pass: |budget - issued| / issued",
    "won roll: budget left (must be exactly 0)",
    "won roll: |photons taken - budget before| / issued",
    "spent pixel: photons taken by an upkeep (must be 0)",
    "pass after a won roll: |budget - issued| / issued",
    "negative budget before opening: |budget - issued| / issued"};

/**
 * @brief Open a budget as the code did when a won roll left a debt: the
 * photons issued plus the negative budget of the last pass.
 *
 * @param sp The star.
 * @param dt_back The time since the last pass.
 */
static void open_budget_with_debt(struct spart *sp, const double dt_back) {
  const double rate_now = sp->feedback_data.radiation.dot_N_ion_pix[0];
  const double rate_prev = sp->feedback_data.radiation.dot_N_ion_pix_prev[0];
  const double rate_used =
      rate_prev < 0. ? rate_now : 0.5 * (rate_prev + rate_now);
  const double debt = fmin(sp->feedback_data.radiation.N_ion_budget_pix[0], 0.);
  sp->feedback_data.radiation.N_ion_budget_pix[0] = debt + rate_used * dt_back;
  sp->feedback_data.radiation.dot_N_ion_pix_prev[0] = rate_now;
}

/**
 * @brief Consume a budget as the code did when a won roll left a debt: the
 * full charge, whatever is left.
 *
 * @param sp The star.
 * @param pixel The pixel.
 * @param cost The charge.
 * @return The charge.
 */
static double consume_budget_with_debt(struct spart *sp, const int pixel,
                                       const double cost) {
  sp->feedback_data.radiation.N_ion_budget_pix[pixel] -= cost;
  return cost;
}

/**
 * @brief Run the no-debt sequence on one pixel with the given budget functions.
 *
 * Each value is 0 when the code behaves as a budget that never goes below
 * zero.
 *
 * @param open Opens the budget of a pass.
 * @param consume Removes a charge from the budget.
 * @param v (return) One value per check.
 */
static void no_debt_sequence(void (*open)(struct spart *, double),
                             double (*consume)(struct spart *, int, double),
                             double v[N_NO_DEBT_CHECKS]) {
  static struct spart sp;
  memset(&sp, 0, sizeof(sp));
  sp.feedback_data.radiation.n_HII_pixels = 1;
  sp.feedback_data.radiation.dot_N_ion_pix[0] = 3.0e49;
  sp.feedback_data.radiation.dot_N_ion_pix_prev[0] = -1.0;
  const double dt = 0.02;
  const double issued = 3.0e49 * dt;

  open(&sp, dt);
  double budget = feedback_get_star_ionization_budget(&sp, 0);
  v[0] = fabs(budget - issued) / issued;

  /* A won roll on a candidate costing 1.7 times the budget. */
  const double taken = consume(&sp, 0, 1.7 * budget);
  budget = feedback_get_star_ionization_budget(&sp, 0);
  v[1] = budget == 0.0 ? 0. : 1.;
  v[2] = fabs(taken - issued) / issued;
  v[3] = consume(&sp, 0, 1.0e40);

  /* The next pass, same rate: its own photons and nothing else. */
  open(&sp, dt);
  budget = feedback_get_star_ionization_budget(&sp, 0);
  v[4] = fabs(budget - issued) / issued;

  /* A budget left negative by hand is ignored on opening. */
  sp.feedback_data.radiation.N_ion_budget_pix[0] = -0.5 * issued;
  open(&sp, dt);
  budget = feedback_get_star_ionization_budget(&sp, 0);
  v[5] = fabs(budget - issued) / issued;
}

/**
 * @brief Checks that a cone budget never goes below zero and that a pass
 * opens with its own photons only. The same sequence run on the debt-carrying
 * copies must fail every check except the first, which is a precondition.
 *
 * @return Number of failed checks.
 */
static int test_no_debt(void) {
  int fail = 0;
  double v[N_NO_DEBT_CHECKS], mutant[N_NO_DEBT_CHECKS];

  no_debt_sequence(radiation_open_ionizing_photon_budget,
                   radiation_consume_ionizing_photons, v);
  for (int k = 0; k < N_NO_DEBT_CHECKS; k++)
    fail += report(no_debt_names[k], v[k], 0., 0);

  no_debt_sequence(open_budget_with_debt, consume_budget_with_debt, mutant);
  char name[128];
  for (int k = 0; k < N_NO_DEBT_CHECKS; k++) {
    snprintf(name, sizeof(name), "control, debt copy: %s (%s)",
             no_debt_names[k], k > 0 ? "must fail" : "must hold");
    fail += report(name, mutant[k], 0., k > 0);
  }

  return fail;
}
#endif /* TEST_HAS_HII_FEEDBACK */

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  int fail = 0;
  fail += test_decision();
  fail += test_renewal();
#ifdef TEST_HAS_HII_FEEDBACK
  fail += test_no_debt();
#else
  printf("skip no-debt checks: GEAR HII feedback is not configured\n");
#endif

  if (fail) printf("%d check(s) failed\n", fail);
  return fail ? 1 : 0;
}
