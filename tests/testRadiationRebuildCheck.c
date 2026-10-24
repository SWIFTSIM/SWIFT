/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026.
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
#include <signal.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

/* Local headers. */
#include "runner_radiation_feedback.h"
#include "swift.h"

/* Not declared in cell.h (matches cell_set_super_hydro/_gravity, which are
 * likewise only called from within cell.c) -- declared here so the
 * orphaned-link regression test below can exercise it directly. */
void cell_set_super_radiation_subgrid(struct cell *c,
                                      struct cell *super_radiation);

/**
 * @brief Build a bare pair of cells with the fields relevant to the
 * radiation rebuild/split predicates set, everything else zeroed.
 */
static void setup_pair(struct cell *ci, struct cell *cj, float dmin,
                       float stars_h_max, float stars_h_hii_max,
                       float hydro_h_max, float stars_dx_max_part,
                       float hydro_dx_max_part) {
  bzero(ci, sizeof(struct cell));
  bzero(cj, sizeof(struct cell));

  ci->dmin = dmin;
  cj->dmin = dmin;
  ci->stars.h_max = stars_h_max;
  ci->stars.h_hii_max = stars_h_hii_max;
  cj->hydro.h_max = hydro_h_max;
  ci->stars.dx_max_part = stars_dx_max_part;
  cj->hydro.dx_max_part = hydro_dx_max_part;
}

/*! Stars:HII_radius_expansion_factor and Stars:HII_max_radius_expansion_tries
 * at their defaults. */
#define test_expansion_factor 1.1f
#define test_expansion_tries 5

/**
 * @brief Outcome of one modelled HII pass of a star that finds no gas.
 */
struct test_pass {

  /*! Limit that ended the radius growth. */
  enum hii_search_stop stop;

  /*! Last radius searched. */
  float radius;

  /*! h_hii_request the pass leaves on the star. */
  float h_request;
};

/**
 * @brief Run the radius growth of runner_dosub_stars_hii_ionization_feedback()
 * for a star that finds nothing at any radius.
 *
 * @param start_radius The shared start radius of the leaf.
 * @param margin cell_radiation_neighbour_margin() of the leaf.
 * @param max_search_radius Stars:HII_max_search_radius.
 * @param count_held Gas the star holds at the start of the pass.
 * @param with_request 0 to model the code before the request existed.
 */
static struct test_pass run_empty_pass(const float start_radius,
                                       const float margin,
                                       const float max_search_radius,
                                       const int count_held,
                                       const int with_request) {
  const float max_reachable = start_radius + margin;
  float radius = start_radius;
  int num_expansions = 0;
  int num_empty = 0;
  struct test_pass out = {hii_search_grow, 0.f, 0.f};
  while (1) {
    const int count_found = 0;
    out.stop = runner_hii_search_stop(radius, max_search_radius, max_reachable,
                                      num_expansions, test_expansion_tries,
                                      count_found, count_held, num_empty);
    if (out.stop == hii_search_stop_margin && with_request)
      out.h_request = runner_hii_requested_radius(radius, test_expansion_factor,
                                                  max_search_radius) /
                      (radiation_search_radius_factor * kernel_gamma);
    if (out.stop != hii_search_grow) break;
    if (count_found == 0 && count_held > 0)
      ++num_empty;
    else
      num_empty = 0;
    radius =
        min3(radius * test_expansion_factor, max_search_radius, max_reachable);
    ++num_expansions;
  }
  out.radius = radius;
  return out;
}

/**
 * @brief Build a radiation-level cell, one of its leaves and one linked
 * partner, wired as the radiation_in stencil of the HII search.
 */
static void setup_stencil(struct cell *level, struct cell *leaf,
                          struct cell *partner, struct task *t_self,
                          struct task *t_pair, struct link *l_self,
                          struct link *l_pair, const float level_dmin,
                          const float leaf_dmin, const float partner_dmin) {
  bzero(level, sizeof(struct cell));
  bzero(leaf, sizeof(struct cell));
  bzero(partner, sizeof(struct cell));
  bzero(t_self, sizeof(struct task));
  bzero(t_pair, sizeof(struct task));
  bzero(l_self, sizeof(struct link));
  bzero(l_pair, sizeof(struct link));
  level->dmin = level_dmin;
  leaf->dmin = leaf_dmin;
  partner->dmin = partner_dmin;
  level->stars.radiation_level = level;
  leaf->stars.radiation_level = level;
  t_self->type = task_type_self;
  t_self->ci = level;
  t_pair->type = task_type_pair;
  t_pair->ci = level;
  t_pair->cj = partner;
  l_self->t = t_self;
  l_self->next = l_pair;
  l_pair->t = t_pair;
  level->stars.radiation_in = l_self;
}

/**
 * @brief A star that only the neighbour cell margin stops asks for a wider
 * reach, and that reach reaches the rebuild check through the cell maxima.
 */
static void test_clamped_search_request(void) {

  const float to_h = 1.f / (radiation_search_radius_factor * kernel_gamma);
  struct cell level, leaf, partner;
  struct task t_self, t_pair;
  struct link l_self, l_pair;

  /* The margin is the thinnest of the leaf and its linked partners. */
  setup_stencil(&level, &leaf, &partner, &t_self, &t_pair, &l_self, &l_pair,
                1.f, 0.5f, 1.f);
  if (cell_radiation_neighbour_margin(&leaf) != 0.5f)
    error("Neighbour margin %g, expected the leaf width 0.5.",
          cell_radiation_neighbour_margin(&leaf));
  partner.dmin = 0.25f;
  if (cell_radiation_neighbour_margin(&leaf) != 0.25f)
    error("Neighbour margin %g, expected the finer partner width 0.25.",
          cell_radiation_neighbour_margin(&leaf));

  /* The margin is reached after k expansions, and binds only if k is below
   * the try count: start * f^(N - 1) >= start + margin. */
  float growth = 1.f;
  for (int k = 0; k < test_expansion_tries - 1; k++)
    growth *= test_expansion_factor;
  const float boundary = 1.f / (growth - 1.f);
  message("Margin binds only above a start radius of %.3f margins.", boundary);

  /* A lone starved star: its own h sets a start radius well below the
   * boundary, so the try count stops it and neither code requests. */
  {
    const float start = 0.2f * 0.5f;
    const struct test_pass old_pass = run_empty_pass(start, 0.5f, 10.f, 0, 0);
    const struct test_pass new_pass = run_empty_pass(start, 0.5f, 10.f, 0, 1);
    if (old_pass.stop != hii_search_stop_tries ||
        new_pass.stop != hii_search_stop_tries)
      error("Lone starved star: expected the try count to bind, got %d/%d.",
            old_pass.stop, new_pass.stop);
    if (new_pass.h_request != 0.f)
      error("Lone starved star requested %g although tries bound.",
            new_pass.h_request);
    message("Lone starved star: tries bind at radius %.4f, no request.",
            new_pass.radius);
  }

  /* A starved star sharing its leaf with a large region: the shared start
   * radius exceeds the boundary, and only the margin stops it. */
  struct spart region, starved;
  bzero(&region, sizeof(struct spart));
  bzero(&starved, sizeof(struct spart));
  region.h_hii = 1.2f * to_h;
  starved.h = 0.01f;

  for (int with_request = 0; with_request < 2; with_request++) {

    /* Radiation level just wider than the region: the first request must
     * fire. */
    setup_stencil(&level, &leaf, &partner, &t_self, &t_pair, &l_self, &l_pair,
                  1.5f, 0.5f, 1.5f);
    starved.h_hii_request = 0.f;
    int fired = 0;
    for (int pass = 0; pass < 20; pass++) {
      const float h_leaf = max(stars_get_hii_search_length(&region),
                               stars_get_hii_search_length(&starved));
      const float start = radiation_search_radius_factor * kernel_gamma *
                          max(h_leaf, starved.h);
      const struct test_pass out = run_empty_pass(
          start, cell_radiation_neighbour_margin(&leaf), 10.f, 0, with_request);
      if (out.stop != hii_search_stop_margin)
        error("Shared start radius: expected the margin to bind, got %d.",
              out.stop);
      starved.h_hii_request = out.h_request;
      level.stars.h_hii_max = max(stars_get_hii_search_length(&region),
                                  stars_get_hii_search_length(&starved));
      if (cell_need_rebuild_for_radiation_pair(&level, &partner)) {
        fired = pass + 1;
        break;
      }
    }
    if (with_request && fired != 1)
      error(
          "Shared start radius: the request fired after %d passes, "
          "expected 1.",
          fired);
    if (!with_request && fired)
      error(
          "Shared start radius: the code without a request fired a "
          "rebuild.");
    message(
        "Shared start radius, %s request: rebuild after %d passes "
        "(0 is never in 20).",
        with_request ? "with" : "without", fired);
  }

  /* Radiation level four cells wide: the request first grows the start
   * radius pass by pass, each request larger than the last, and fires once
   * it outgrows the wired neighbours. */
  {
    setup_stencil(&level, &leaf, &partner, &t_self, &t_pair, &l_self, &l_pair,
                  4.f, 0.5f, 4.f);
    starved.h_hii_request = 0.f;
    float last_request = 0.f;
    int fired = 0;
    for (int pass = 0; pass < 20; pass++) {
      const float h_leaf = max(stars_get_hii_search_length(&region),
                               stars_get_hii_search_length(&starved));
      const float start = radiation_search_radius_factor * kernel_gamma *
                          max(h_leaf, starved.h);
      const struct test_pass out = run_empty_pass(
          start, cell_radiation_neighbour_margin(&leaf), 10.f, 0, 1);
      if (!(out.h_request > last_request))
        error("Pass %d: request %g does not exceed the last one %g.", pass,
              out.h_request, last_request);
      last_request = out.h_request;
      starved.h_hii_request = out.h_request;
      level.stars.h_hii_max = max(stars_get_hii_search_length(&region),
                                  stars_get_hii_search_length(&starved));
      if (cell_need_rebuild_for_radiation_pair(&level, &partner)) {
        fired = pass + 1;
        break;
      }
    }
    if (fired < 2)
      error(
          "Wide radiation level: rebuild after %d passes, expected more "
          "than one.",
          fired);
    message(
        "Wide radiation level: rebuild after %d passes, requested radius "
        "%.3f.",
        fired, last_request / to_h);

    /* A tree whose level covers the request is silent: no second rebuild. */
    level.dmin = 2.f * last_request / to_h;
    partner.dmin = level.dmin;
    if (cell_need_rebuild_for_radiation_pair(&level, &partner))
      error("A level wide enough for the request still demands a rebuild.");
  }

  /* Not clamped: the ceiling, or the one-probe rule of a star holding gas,
   * stops these stars before the margin is the only limit. */
  {
    const struct test_pass ceiling = run_empty_pass(1.f, 0.5f, 1.2f, 0, 1);
    if (ceiling.stop != hii_search_stop_ceiling || ceiling.h_request != 0.f)
      error("Ceiling: stop %d, request %g.", ceiling.stop, ceiling.h_request);
    const struct test_pass probe = run_empty_pass(1.f, 0.05f, 10.f, 3, 1);
    if (probe.stop != hii_search_stop_one_probe || probe.h_request != 0.f)
      error("One probe: stop %d, request %g.", probe.stop, probe.h_request);
    message(
        "Not clamped: ceiling stop at %.3f, one-probe stop at %.3f, no "
        "request.",
        ceiling.radius, probe.radius);
  }

  /* The helper reads the larger of the region and the request. */
  starved.h_hii = 0.3f;
  starved.h_hii_request = 0.f;
  if (stars_get_hii_search_length(&starved) != 0.3f)
    error("Search length ignores h_hii.");
  starved.h_hii_request = 0.7f;
  if (stars_get_hii_search_length(&starved) != 0.7f)
    error("Search length ignores h_hii_request.");
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  struct cell ci, cj;

  /* Case 1: everything comfortably smaller than the cell -- no rebuild
   * needed for either the plain stars check or the radiation-aware one. */
  setup_pair(&ci, &cj, /*dmin=*/1.0f, /*stars_h_max=*/0.01f,
             /*stars_h_hii_max=*/0.01f, /*hydro_h_max=*/0.01f,
             /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
  if (cell_need_rebuild_for_stars_pair(&ci, &cj))
    error(
        "Unexpected rebuild flagged by cell_need_rebuild_for_stars_pair "
        "for a small, static configuration.");
  if (cell_need_rebuild_for_radiation_pair(&ci, &cj))
    error(
        "Unexpected rebuild flagged by cell_need_rebuild_for_radiation_pair "
        "for a small, static configuration.");

  /* Case 2: h_hii_max alone grown large enough to exceed the cell size,
   * while ordinary h_max/dx_max stay small. This is the scenario the two
   * checks must disagree on: the plain stars check has no h_hii term and
   * must stay silent, while the radiation-aware check must flag a
   * rebuild. */
  setup_pair(&ci, &cj, /*dmin=*/1.0f, /*stars_h_max=*/0.01f,
             /*stars_h_hii_max=*/2.0f, /*hydro_h_max=*/0.01f,
             /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
  if (cell_need_rebuild_for_stars_pair(&ci, &cj))
    error(
        "cell_need_rebuild_for_stars_pair must not react to h_hii_max "
        "growth (it has no h_hii term by design).");
  if (!cell_need_rebuild_for_radiation_pair(&ci, &cj))
    error(
        "cell_need_rebuild_for_radiation_pair failed to flag a rebuild "
        "after h_hii_max grew past the cell size.");

  /* Case 3: ordinary h_max grown large (not h_hii) -- both checks must
   * agree and flag a rebuild, since this is the term they share. */
  setup_pair(&ci, &cj, /*dmin=*/1.0f, /*stars_h_max=*/2.0f,
             /*stars_h_hii_max=*/0.0f, /*hydro_h_max=*/0.01f,
             /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
  if (!cell_need_rebuild_for_stars_pair(&ci, &cj))
    error(
        "cell_need_rebuild_for_stars_pair failed to flag a rebuild after "
        "stars.h_max grew past the cell size.");
  if (!cell_need_rebuild_for_radiation_pair(&ci, &cj))
    error(
        "cell_need_rebuild_for_radiation_pair failed to flag a rebuild "
        "after stars.h_max grew past the cell size.");

  /* Case 4: drift alone (dx_max_part) pushes the pair over the edge, with
   * h_max/h_hii_max both otherwise small. */
  setup_pair(&ci, &cj, /*dmin=*/1.0f, /*stars_h_max=*/0.01f,
             /*stars_h_hii_max=*/0.01f, /*hydro_h_max=*/0.01f,
             /*stars_dx_max_part=*/0.6f, /*hydro_dx_max_part=*/0.6f);
  if (!cell_need_rebuild_for_radiation_pair(&ci, &cj))
    error(
        "cell_need_rebuild_for_radiation_pair failed to flag a rebuild "
        "after dx_max_part drift pushed the pair past the cell size.");

  /* Case 5: h_hii_max sits in the 1.0x-radiation_search_radius_factor(x)
   * band -- inside the search radius the runner actually uses
   * (radiation_search_radius_factor * kernel_gamma * h_hii_max), but
   * outside the unfactored kernel_gamma * h_hii_max the check used before
   * the fold-in fix. Midpoint between 1.0 and the factor, so the case sits
   * comfortably inside the band regardless of the factor's exact value. */
  {
    const float mid_factor = 0.5f * (1.0f + radiation_search_radius_factor);
    const float dmin = 1.0f;
    const float h_hii_band = dmin / (mid_factor * kernel_gamma);
    setup_pair(&ci, &cj, dmin, /*stars_h_max=*/0.01f,
               /*stars_h_hii_max=*/h_hii_band, /*hydro_h_max=*/0.01f,
               /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
    if (kernel_gamma * h_hii_band >= dmin)
      error(
          "Test setup error: h_hii_band must sit BELOW the unfactored "
          "threshold (kernel_gamma * h_hii_max < dmin) for this case to "
          "discriminate the fold-in fix.");
    if (kernel_gamma * radiation_search_radius_factor * h_hii_band <= dmin)
      error(
          "Test setup error: h_hii_band must sit ABOVE the factored "
          "threshold (kernel_gamma * factor * h_hii_max > dmin) for this "
          "case to discriminate the fold-in fix.");
    if (!cell_need_rebuild_for_radiation_pair(&ci, &cj))
      error(
          "cell_need_rebuild_for_radiation_pair missed a rebuild for "
          "h_hii_max in the 1.0x-factor(x) search-radius band: the check "
          "must bound the SAME radius the runner searches with "
          "(radiation_search_radius_factor * kernel_gamma * h_hii_max), "
          "not the unfactored kernel_gamma * h_hii_max.");
  }

  /* Case 5b: unequal dmin (asymmetric radiation pair) -- generalizes the
   * "ci->dmin == cj->dmin" assumption baked into every other pair check in
   * this file. Both orientations below use the same reach (0.4), which
   * clears a dmin of 0.25 but not a dmin of 1.0. */
  {
    struct cell cell_fine, cell_coarse;
    bzero(&cell_fine, sizeof(struct cell));
    bzero(&cell_coarse, sizeof(struct cell));
    cell_fine.dmin = 0.25f;
    cell_coarse.dmin = 1.0f;

    /* Orientation A: the fine cell is ci (holds the star growth) and the
     * coarse cell is cj -- the discriminating case, where scaling against
     * cj->dmin alone (the pre-fix behaviour) would wrongly stay silent. */
    cell_fine.stars.dx_max_part = 0.4f;
    if (!cell_need_rebuild_for_radiation_pair(&cell_fine, &cell_coarse))
      error(
          "cell_need_rebuild_for_radiation_pair must scale against "
          "min(ci->dmin, cj->dmin), not cj->dmin alone -- missed a rebuild "
          "whose reach cleared the fine cell's dmin but not the coarse "
          "cell's.");

    /* Orientation B: same reach, but the fine cell is now cj and the
     * coarse cell (holding the star growth) is ci -- confirms min() picks
     * up the fine dmin regardless of which argument slot it occupies (the
     * pre-fix code would already have caught this orientation, since the
     * fine cell IS cj there -- a non-regression check, not a discriminator).
     */
    cell_fine.stars.dx_max_part = 0.0f;
    cell_coarse.stars.dx_max_part = 0.4f;
    if (!cell_need_rebuild_for_radiation_pair(&cell_coarse, &cell_fine))
      error(
          "cell_need_rebuild_for_radiation_pair failed when the fine cell "
          "(dmin=0.25) sits in the cj argument slot instead of ci.");
  }

  /* Case 5c: a split next to a cell that cannot split. The splitting cell's
   * progeny then pair with the coarse neighbour, whose reach the criterion
   * tests against the progeny width. Right after a rebuild (no drift) that
   * test must stay silent, or engine_rebuild() aborts. Here the coarse
   * cell's reach (0.278) fits its own width (0.5) but not the progeny width
   * (0.25), while the neighbour's own reach would let it split. */
  {
    struct cell coarse, parent, progeny;
    bzero(&coarse, sizeof(struct cell));
    bzero(&parent, sizeof(struct cell));
    bzero(&progeny, sizeof(struct cell));
    coarse.split = parent.split = 1;
    coarse.hydro.count = parent.hydro.count = progeny.hydro.count = 1000;
    coarse.dmin = parent.dmin = 0.5f;
    progeny.dmin = 0.25f;
    coarse.loc[0] = 0.;
    parent.loc[0] = 0.5;
    coarse.stars.h_hii_max = 0.119576f;
    parent.stars.h_hii_max = 0.0484152f;
    progeny.stars.h_hii_max = 0.0484152f;
    coarse.hydro.h_max = parent.hydro.h_max = progeny.hydro.h_max = 0.0347f;
    coarse.top = &coarse;
    parent.top = &parent;
    progeny.top = &parent;

    const int cdim[3] = {2, 1, 1};
    const float reach[2] = {coarse.stars.h_hii_max, parent.stars.h_hii_max};
    float reach_neighbours[2];
    space_radiation_neighbour_reach_for(reach, cdim, /*periodic=*/0,
                                        reach_neighbours);
    coarse.stars.h_reach_max_neighbours = reach_neighbours[0];
    parent.stars.h_reach_max_neighbours = reach_neighbours[1];

    if (cell_can_split_pair_radiation_subgrid_task(&coarse))
      error("Case 5c: the coarse cell's own reach must stop its split.");
    if (cell_can_split_pair_radiation_subgrid_task(&parent) &&
        cell_need_rebuild_for_radiation_pair(&coarse, &progeny))
      error(
          "Case 5c: a cell split next to a coarse neighbour whose reach "
          "exceeds the progeny width; the rebuild criterion fires with no "
          "drift, so no rebuild can satisfy it.");

    /* Control: with a small neighbour the same cell must still split. */
    const float reach_small[2] = {0.01f, parent.stars.h_hii_max};
    space_radiation_neighbour_reach_for(reach_small, cdim, /*periodic=*/0,
                                        reach_neighbours);
    parent.stars.h_reach_max_neighbours = reach_neighbours[1];
    if (!cell_can_split_pair_radiation_subgrid_task(&parent))
      error("Case 5c: a split was refused with no large reach nearby.");
  }

  /* Case 5d: a cell held coarse only by a big-reach star in a neighbouring
   * top cell, next to a cell that split two levels deeper. The coarse cell's
   * gas kernel support (0.35) fits its own width (1.0) but not the progeny
   * width (0.25). The HII search is centred on the star and never selects gas
   * by the gas's own smoothing length, so with no drift the rebuild check
   * must stay silent for this pair, in both directions. */
  {
    struct cell coarse, fine;
    bzero(&coarse, sizeof(struct cell));
    bzero(&fine, sizeof(struct cell));
    coarse.split = 1;
    coarse.hydro.count = fine.hydro.count = 1000;
    coarse.dmin = 1.0f;
    fine.dmin = 0.25f;
    coarse.hydro.h_max = 0.35f / kernel_gamma;
    coarse.stars.h_hii_max = 0.01f;
    fine.stars.h_hii_max = 0.01f;
    fine.hydro.h_max = 0.01f;
    coarse.top = &coarse;
    coarse.stars.h_reach_max_neighbours = 0.3f;

    if (!cell_can_split_pair_hydro_task(&coarse, /*r_cut=*/0.f))
      error("Case 5d: the coarse cell's own gas must allow a hydro split.");
    if (cell_can_split_pair_radiation_subgrid_task(&coarse))
      error("Case 5d: the neighbour reach must stop the coarse cell's split.");
    if (cell_need_rebuild_for_radiation_pair(&fine, &coarse) ||
        cell_need_rebuild_for_radiation_pair(&coarse, &fine))
      error(
          "Case 5d: the rebuild check fired on the coarse cell's gas "
          "smoothing length with no drift; no rebuild can satisfy it.");

    /* Control: gas drift beyond the progeny width must still fire. */
    fine.hydro.dx_max_part = 0.3f;
    if (!cell_need_rebuild_for_radiation_pair(&coarse, &fine))
      error("Case 5d: a gas drift beyond the progeny width was missed.");
  }

  /* Now check cell_can_split_pair/self_radiation_subgrid_task(): radiation
   * must NEVER split once the cell is at or below a hydro attach point
   * (cell_flag_at_or_below_hydro_attach set), for any combination of the
   * other fields (including small h_max/h_hii_max that would otherwise
   * geometrically allow splitting). This is the fix for the duplicate-unlock
   * crash at gas_mass=0.01: once a cell IS (or lies below) a hydro.super,
   * every one of its progeny would end up sharing that same hydro.super, so
   * splitting further would only produce sibling radiation tasks that all
   * wire dependencies to the identical hydro.super. See
   * cell_can_split_pair_radiation_subgrid_task()'s docstring in cell.h for
   * the full history. */
  struct cell c;
  const float dmins[] = {0.5f, 1.0f, 2.0f};
  const float h_values[] = {0.0f, 0.01f, 0.3f, 1.0f, 5.0f};
  for (int is = 0; is < 2; ++is) {
    for (int id = 0; id < 3; ++id) {
      for (int ih = 0; ih < 5; ++ih) {
        bzero(&c, sizeof(struct cell));
        c.split = is;
        c.dmin = dmins[id];
        c.hydro.count = 1; /* Non-gasless; see the dedicated block below. */
        c.hydro.h_max = h_values[ih];
        c.stars.h_max = h_values[ih];
        c.stars.h_hii_max = h_values[(ih + 1) % 5];
        c.sinks.h_max = h_values[ih];
        c.black_holes.h_max = h_values[ih];

        /* With the at-or-below-hydro-attach flag set, must always refuse
         * to split regardless of geometry. */
        cell_set_flag(&c, cell_flag_at_or_below_hydro_attach);
        if (cell_can_split_pair_radiation_subgrid_task(&c))
          error(
              "cell_can_split_pair_radiation_subgrid_task must never "
              "split once at-or-below a hydro attach point "
              "(split=%d dmin=%e h=%e h_hii=%e).",
              c.split, c.dmin, c.hydro.h_max, c.stars.h_hii_max);

        if (cell_can_split_self_radiation_subgrid_task(&c))
          error(
              "cell_can_split_self_radiation_subgrid_task must never "
              "split once at-or-below a hydro attach point "
              "(split=%d dmin=%e h=%e h_hii=%e).",
              c.split, c.dmin, c.hydro.h_max, c.stars.h_hii_max);

        /* With the flag cleared (strictly above any hydro attach point),
         * splitting must be governed purely by the geometric criterion --
         * identical to what the (now-removed) hydro.super-gated call would
         * have returned before this test replaced it. */
        cell_clear_flag(&c, cell_flag_at_or_below_hydro_attach);
        const int expect_geometric =
            c.split && c.hydro.count > 0 &&
            (space_stretch * kernel_gamma * c.hydro.h_max < 0.5f * c.dmin) &&
            (space_stretch * radiation_search_radius_factor * kernel_gamma *
                 c.stars.h_max <
             0.5f * c.dmin) &&
            (space_stretch * radiation_search_radius_factor * kernel_gamma *
                 c.stars.h_hii_max <
             0.5f * c.dmin) &&
            (space_stretch * kernel_gamma * c.sinks.h_max < 0.5f * c.dmin) &&
            (space_stretch * kernel_gamma * c.black_holes.h_max <
             0.5f * c.dmin);

        if (cell_can_split_pair_radiation_subgrid_task(&c) != expect_geometric)
          error(
              "cell_can_split_pair_radiation_subgrid_task disagreed with "
              "the geometric criterion while the attach flag is clear "
              "(split=%d dmin=%e h=%e h_hii=%e expected=%d).",
              c.split, c.dmin, c.hydro.h_max, c.stars.h_hii_max,
              expect_geometric);

        if (cell_can_split_self_radiation_subgrid_task(&c) != expect_geometric)
          error(
              "cell_can_split_self_radiation_subgrid_task disagreed with "
              "the geometric criterion while the attach flag is clear "
              "(split=%d dmin=%e h=%e h_hii=%e expected=%d).",
              c.split, c.dmin, c.hydro.h_max, c.stars.h_hii_max,
              expect_geometric);
      }
    }
  }

  /* A gasless cell (hydro.count == 0, e.g. a star-only branch with no gas
   * anywhere in its subtree) must never split, regardless of how permissive
   * the geometric terms otherwise are -- there is no gas at any depth to
   * gain resolution on, so descending further only wastes task-graph
   * churn. See cell_can_split_pair_radiation_subgrid_task()'s docstring. */
  bzero(&c, sizeof(struct cell));
  c.split = 1;
  c.dmin = 2.0f;
  c.hydro.count = 0;
  c.hydro.h_max = 0.0f;
  c.stars.h_max = 0.0f;
  c.stars.h_hii_max = 0.0f;
  c.sinks.h_max = 0.0f;
  c.black_holes.h_max = 0.0f;
  if (cell_can_split_pair_radiation_subgrid_task(&c))
    error(
        "cell_can_split_pair_radiation_subgrid_task must never split a "
        "gasless cell (hydro.count == 0), even with fully permissive "
        "geometry.");
  if (cell_can_split_self_radiation_subgrid_task(&c))
    error(
        "cell_can_split_self_radiation_subgrid_task must never split a "
        "gasless cell (hydro.count == 0), even with fully permissive "
        "geometry.");

  /* Confirm hydro's own split criteria stay decoupled from h_hii_max:
   * growing h_hii_max alone must NOT change cell_can_split_pair_hydro_task
   * -- hydro splitting is a pure hydro concept (see cell.h). */
  bzero(&c, sizeof(struct cell));
  c.split = 1;
  c.dmin = 1.0f;
  c.hydro.h_max = 0.01f;
  c.stars.h_max = 0.01f;
  c.stars.h_hii_max = 0.01f;
  c.sinks.h_max = 0.01f;
  c.black_holes.h_max = 0.01f;
  if (!cell_can_split_pair_hydro_task(&c, /*r_cut=*/0.f))
    error(
        "Expected cell_can_split_pair_hydro_task to allow splitting for "
        "a small configuration.");
  c.stars.h_hii_max = 2.0f;
  if (!cell_can_split_pair_hydro_task(&c, /*r_cut=*/0.f))
    error(
        "cell_can_split_pair_hydro_task must stay decoupled from "
        "h_hii_max -- it should still allow splitting after h_hii_max "
        "alone grew past the cell size.");

  /* Regression test for cell_set_super_radiation_subgrid()'s "orphaned
   * radiation_in link" SWIFT_DEBUG_CHECKS assertion (cell.c). Invariant
   * discovered 2026-07-12: the flat radiation_in walk in
   * runner_dosub_stars_hii_ionization_feedback only ever reads
   * radiation_level's OWN link list -- it never descends into
   * radiation_level's children to pick up links attached deeper. So every
   * cell holding a radiation_in link must BE its own radiation_level.
   * This can be violated when a pair task with one neighbour stops
   * splitting early (its link lands on a coarse ancestor, e.g. because
   * that neighbour's hydro.super sits at a shallower depth) while a
   * DIFFERENT, matching-depth neighbour's pair (or the cell's own
   * self-splitting) attaches its link to a deeper descendant: the
   * "topmost link found first wins, propagated to all descendants" rule
   * in cell_set_super_radiation_subgrid then makes the deeper link
   * permanently invisible to the flat walk -- a silent, directional
   * search-completeness bug, not a crash. */
  {
    struct cell top, child;
    struct link top_link, child_link;
    bzero(&top, sizeof(struct cell));
    bzero(&child, sizeof(struct cell));
    bzero(&top_link, sizeof(struct link));
    bzero(&child_link, sizeof(struct link));

    /* Case A: no orphaning -- only the top cell holds a link (the common,
     * correct case: a pair/self task that never split below the top
     * cell). radiation_level must land on top for both cells; no error. */
    top.split = 1;
    top.progeny[0] = &child;
    top.stars.radiation_in = &top_link;
    child.stars.radiation_in = NULL;
    cell_set_super_radiation_subgrid(&top, NULL);
    if (top.stars.radiation_level != &top)
      error(
          "cell_set_super_radiation_subgrid: expected radiation_level == "
          "top for the unsplit case.");
    if (child.stars.radiation_level != &top)
      error(
          "cell_set_super_radiation_subgrid: expected the child to "
          "inherit radiation_level == top.");

    /* Case B: orphaning -- BOTH top and child hold their own link (top's
     * simulates an unsplit, mismatched-neighbour pair; child's simulates a
     * properly-split sibling/neighbour pair or self-task). radiation_level
     * is pulled to top (topmost wins), silently orphaning child's link.
     * Must be caught by the SWIFT_DEBUG_CHECKS assertion -- run in a
     * forked child so the expected abort doesn't kill the test. If this
     * test binary was built without SWIFT_DEBUG_CHECKS, the assertion
     * doesn't exist and this case is skipped (nothing to catch it with). */
#ifdef SWIFT_DEBUG_CHECKS
    child.stars.radiation_in = &child_link;

    const pid_t pid = fork();
    if (pid == 0) {
      /* Child process: silence stderr (the error() message is expected
       * output, not a test failure), then trigger the assertion. */
      if (freopen("/dev/null", "w", stderr) == NULL) _exit(43);
      cell_set_super_radiation_subgrid(&top, NULL);
      /* Reached only if the assertion did NOT fire -- signal failure with
       * a distinguishable exit code (real error() exits with status 1). */
      _exit(42);
    } else if (pid > 0) {
      int status;
      waitpid(pid, &status, 0);
      /* error() aborts via swift_abort(), which is exit(1) normally but
       * abort() (SIGABRT) under SWIFT_DEVELOP_MODE -- accept either. */
      const int exited_with_error =
          WIFEXITED(status) && WEXITSTATUS(status) == 1;
      const int aborted = WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
      if (!exited_with_error && !aborted)
        error(
            "cell_set_super_radiation_subgrid failed to detect an "
            "orphaned radiation_in link (WIFEXITED=%d WEXITSTATUS=%d "
            "WIFSIGNALED=%d WTERMSIG=%d) -- the SWIFT_DEBUG_CHECKS "
            "assertion did not fire as expected.",
            WIFEXITED(status), WIFEXITED(status) ? WEXITSTATUS(status) : -1,
            WIFSIGNALED(status), WIFSIGNALED(status) ? WTERMSIG(status) : -1);
    } else {
      error("fork() failed in the orphaned-link regression test.");
    }
#endif
  }

  /* Regression tests for the periodic cdim==3 rebuild-floor fix: at that
   * grid, the 27-cell radiation stencil wires every top-level pair, so
   * cell_unskip.c's call site skips cell_need_rebuild_for_radiation_pair
   * for a pair of UNSPLIT top-level cells once
   * space_radiation_top_stencil_covers_box() is true (see its docstring in
   * space.h and the guard in cell_unskip_radiation_tasks()). That predicate
   * is called directly here (not threaded through
   * cell_need_rebuild_for_radiation_pair, which stays cell.h-only and
   * space/engine-agnostic), so it is exercised as the real function, not a
   * re-implementation that could silently drift from space_regrid.c. */

  /* Case 6: space_radiation_top_stencil_covers_box() itself -- true only
   * for a periodic box at exactly cdim == {3,3,3}. */
  {
    struct space test_space;
    const struct {
      int periodic;
      int cdim[3];
      int expected;
    } flag_cases[] = {
        {1, {3, 3, 3}, 1}, {1, {4, 3, 3}, 0}, {1, {3, 4, 3}, 0},
        {1, {3, 3, 4}, 0}, {1, {2, 3, 3}, 0}, {0, {3, 3, 3}, 0},
    };
    for (size_t i = 0; i < sizeof(flag_cases) / sizeof(flag_cases[0]); ++i) {
      bzero(&test_space, sizeof(struct space));
      test_space.periodic = flag_cases[i].periodic;
      test_space.cdim[0] = flag_cases[i].cdim[0];
      test_space.cdim[1] = flag_cases[i].cdim[1];
      test_space.cdim[2] = flag_cases[i].cdim[2];
      const int computed = space_radiation_top_stencil_covers_box(&test_space);
      if (computed != flag_cases[i].expected)
        error(
            "space_radiation_top_stencil_covers_box mismatch for case %zu "
            "(periodic=%d cdim=[%d %d %d]): got %d, expected %d.",
            i, flag_cases[i].periodic, flag_cases[i].cdim[0],
            flag_cases[i].cdim[1], flag_cases[i].cdim[2], computed,
            flag_cases[i].expected);
    }
  }

  /* Case 7: the cell_unskip.c guard's full condition, mirrored here since
   * exercising the real call site needs a fully-wired engine/scheduler/task
   * graph. A pair whose h_hii_max violates the rebuild criterion must have
   * its rebuild demand suppressed ONLY when BOTH the box-covers flag is set
   * AND both cells are unsplit top-level cells (ci->top == ci); a split
   * sub-cell pair (ci->top != ci) must still demand a rebuild even with the
   * flag set, since its own smaller dmin came from a geometric split
   * decision the top-level stencil argument says nothing about. */
  {
    const float dmin = 1.0f;
    setup_pair(&ci, &cj, dmin, /*stars_h_max=*/0.01f,
               /*stars_h_hii_max=*/2.0f, /*hydro_h_max=*/0.01f,
               /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
    if (!cell_need_rebuild_for_radiation_pair(&ci, &cj))
      error(
          "Test setup error: this pair must violate the radiation rebuild "
          "criterion on its own for the cases below to mean anything.");

    struct space test_space;
    bzero(&test_space, sizeof(struct space));
    test_space.periodic = 1;
    test_space.cdim[0] = test_space.cdim[1] = test_space.cdim[2] = 3;

    for (int box_covers = 0; box_covers <= 1; ++box_covers) {
      for (int top_level = 0; top_level <= 1; ++top_level) {
        /* box_covers toggles between the case-6-verified true/false space;
         * top_level toggles whether ci/cj present as unsplit top cells. */
        struct space *s = box_covers ? &test_space : NULL;
        ci.top = top_level ? &ci : NULL;
        cj.top = top_level ? &cj : NULL;

        const int guard_skips = s != NULL &&
                                space_radiation_top_stencil_covers_box(s) &&
                                ci.top == &ci && cj.top == &cj;
        const int would_rebuild =
            !guard_skips && (cell_need_rebuild_for_radiation_pair(&ci, &cj) ||
                             cell_need_rebuild_for_radiation_pair(&cj, &ci));

        const int expect_skip = box_covers && top_level;
        if (expect_skip && would_rebuild)
          error(
              "cell_unskip_radiation_tasks's guard must not demand a "
              "rebuild for a criterion-violating pair of unsplit top-level "
              "cells once the box-covers flag is set.");
        if (!expect_skip && !would_rebuild)
          error(
              "cell_unskip_radiation_tasks's guard must still demand a "
              "rebuild for a criterion-violating pair when box_covers=%d "
              "top_level=%d (must only skip when BOTH hold).",
              box_covers, top_level);
      }
    }
  }

  /* Case 8: the criterion and the top-level grid must agree. If a rebuild
   * criterion demands a rebuild and the unskip guard does not exempt the
   * pair, the rebuild that follows MUST produce a different grid (or,
   * once no grid can satisfy it, space_regrid()'s own "too few cells"
   * abort must catch it) -- otherwise unskip demands a rebuild again
   * immediately and engine_rebuild() aborts with "engine_unskip failed
   * after a rebuild!". Two independent criteria read a star's smoothing
   * lengths, so this case covers both: cell_need_rebuild_for_radiation_pair()
   * bounds a star's reach by max(h_hii_max, h_max), factored and exempted
   * once the radiation stencil covers the box (8a sweeps h_hii, 8b sweeps
   * h_max, h_max held/held-at-0 respectively);
   * cell_need_rebuild_for_stars_pair() separately bounds it by the unfactored,
   * unexempted stars.h_max alone (8c).
   *
   * Geometry is TaskGraphPair's (examples/SubgridTests/StellarFeedback/
   * HIIRegions/TaskGraphPair/params.yml): a periodic box of 0.16 with
   * Scheduler:max_top_level_cells 4 and the default min_top_level_cells 3,
   * so the grid sits at cdim = {4,4,4} with dmin = 0.04. cell_min and
   * cell_max_width mirror space_init()'s derivation (src/space.c:1271 and
   * :1297); the grid sizing itself is NOT mirrored -- "what the next
   * regrid would build" calls the real space_regrid_star_radiation_term_for(),
   * space_regrid_star_h_max_no_hii_term_for() and
   * space_regrid_cell_width_for(), the exact functions space_regrid() itself
   * calls at every one of its three accumulation sites (the
   * local_cells_with_particles_top loop, the cells_top loop, and the
   * flat-particle-array fallback), so this case's ROUTING FORMULAS cannot
   * drift from space_regrid.c's the way an in-test re-derivation would.
   * This does NOT cover wiring: the test supplies its own arguments to
   * these functions directly, so a future edit that stops calling either
   * function from one of the three sites (reverting to an inline
   * re-implementation there instead) would not be caught here -- only a
   * direct read of space_regrid.c's three call sites catches that.
   *
   * All three sweeps deliberately run well past
   * space_regrid_radiation_cap_for() (the coarsening threshold below which
   * the radiation-pair demand cannot be satisfied), so 8a/8b also exercise
   * the capped branch (radiation term == cap, grid at min_top_level_cells),
   * not only the transitional window around first onset; 8c's h_max_no_hii
   * route is uncapped, so its sweep instead crosses into the range where
   * no grid can satisfy the demand and the "too few cells" abort is the
   * correct outcome (see 8c's own comment). */
  {
    const double dim = 0.16;
    const int maxtcells = 4;
    const int mintcells = 3;
    const double tol = fmax(1.0 - 1.0 / (double)(maxtcells * maxtcells), 0.99);
    const double cell_min = tol * dim / (double)maxtcells;
    const double cell_max_width = dim / (double)mintcells;
    const float h_max_floor = (float)cell_min / kernel_gamma / space_stretch;
    const double dmin = dim / (double)maxtcells;

    struct space test_space;
    bzero(&test_space, sizeof(struct space));
    test_space.periodic = 1;
    test_space.cdim[0] = test_space.cdim[1] = test_space.cdim[2] = maxtcells;

    /* A star-free grid must be untouched by anything the radiation term
     * does: with no h_hii and no star anywhere the width stays at the
     * cell_min floor and the grid stays at max_top_level_cells. Without
     * this, scaling the radiation term while its accumulator still starts
     * at the floor value would coarsen every run in the tree, radiation or
     * not, and no example-level test would attribute it. */
    {
      const float star_radiation_term = space_regrid_star_radiation_term_for(
          /*h_hii_max=*/0.f, /*h_max=*/0.f, cell_max_width);
      const double width = space_regrid_cell_width_for(
          h_max_floor, star_radiation_term, cell_min, cell_max_width);
      if (!(width > 0.))
        error("Star-free regrid width is not a positive number (%g).", width);
      const int cd = (int)floor(dim / width);
      if (cd != maxtcells)
        error(
            "A star-free grid (h_hii = h_max = 0) must stay at "
            "Scheduler:max_top_level_cells (%d), got cdim = %d. The "
            "radiation term must not contribute when no star is present.",
            maxtcells, cd);
    }

    /* Case 8a: h_hii-driven, h_max held small (0.01, well below the sweep
     * range) so max(h_hii_max, h_max) == h_hii_max throughout. */
    {
      int fired = 0;
      for (int i = 0; i <= 290; ++i) {
        const float h_hii = 0.016f + 0.0001f * (float)i;

        setup_pair(&ci, &cj, (float)dmin, /*stars_h_max=*/0.01f,
                   /*stars_h_hii_max=*/h_hii, /*hydro_h_max=*/0.01f,
                   /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
        ci.top = &ci;
        cj.top = &cj;

        const int exempt =
            space_radiation_top_stencil_covers_box(&test_space) &&
            ci.top == &ci;
        if (!cell_need_rebuild_for_radiation_pair(&ci, &cj) || exempt) continue;
        ++fired;

        const float star_radiation_term = space_regrid_star_radiation_term_for(
            h_hii, /*h_max=*/0.01f, cell_max_width);
        const double width = space_regrid_cell_width_for(
            h_max_floor, star_radiation_term, cell_min, cell_max_width);
        if (!(width > 0.))
          error("Regrid width is not a positive number (%g) at h_hii = %g.",
                width, h_hii);
        const int cd = (int)floor(dim / width);

        if (cd == maxtcells)
          error(
              "Case 8a: cell_need_rebuild_for_radiation_pair demands a "
              "rebuild at h_hii = %g (dmin = %g, unskip guard does not "
              "exempt the pair), but the regrid rebuilds the identical "
              "cdim = %d grid (width %g, radiation term %g). The demand can "
              "never be satisfied, so engine_unskip() re-raises it and "
              "engine_rebuild() aborts.",
              h_hii, dmin, cd, width, star_radiation_term);
      }
      if (fired == 0)
        error(
            "Case 8a swept h_hii across [0.016, 0.045] and no sample ever "
            "made cell_need_rebuild_for_radiation_pair fire a non-exempt "
            "demand; the sweep range no longer exercises the invariant it "
            "is meant to check.");
    }

    /* Case 8b: h_max-driven, h_hii held at 0 so max(h_hii_max, h_max) ==
     * h_max throughout -- a star with no active HII region, resting on its
     * own ordinary smoothing length in a sparse top-level cell. */
    {
      int fired = 0;
      for (int i = 0; i <= 290; ++i) {
        const float stars_h_max = 0.016f + 0.0001f * (float)i;

        setup_pair(&ci, &cj, (float)dmin, /*stars_h_max=*/stars_h_max,
                   /*stars_h_hii_max=*/0.f, /*hydro_h_max=*/0.01f,
                   /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
        ci.top = &ci;
        cj.top = &cj;

        const int exempt =
            space_radiation_top_stencil_covers_box(&test_space) &&
            ci.top == &ci;
        if (!cell_need_rebuild_for_radiation_pair(&ci, &cj) || exempt) continue;
        ++fired;

        const float star_radiation_term = space_regrid_star_radiation_term_for(
            /*h_hii_max=*/0.f, stars_h_max, cell_max_width);
        const double width = space_regrid_cell_width_for(
            h_max_floor, star_radiation_term, cell_min, cell_max_width);
        if (!(width > 0.))
          error("Regrid width is not a positive number (%g) at h_max = %g.",
                width, stars_h_max);
        const int cd = (int)floor(dim / width);

        if (cd == maxtcells)
          error(
              "Case 8b: cell_need_rebuild_for_radiation_pair demands a "
              "rebuild at stars.h_max = %g with h_hii = 0 (dmin = %g, "
              "unskip guard does not exempt the pair), but the regrid "
              "rebuilds the identical cdim = %d grid (width %g, radiation "
              "term %g). A star with no HII region can hit the same "
              "unsatisfiable-rebuild deadlock as an HII-driven star unless "
              "stars.h_max is folded into the same factored/capped term as "
              "h_hii.",
              stars_h_max, dmin, cd, width, star_radiation_term);
      }
      if (fired == 0)
        error(
            "Case 8b swept stars.h_max across [0.016, 0.045] and no sample "
            "ever made cell_need_rebuild_for_radiation_pair fire a "
            "non-exempt demand; the sweep range no longer exercises the "
            "invariant it is meant to check.");
    }

    /* Case 8c: the SAME invariant as 8a/8b, but for
     * cell_need_rebuild_for_stars_pair() -- unfactored, uncapped, and
     * unexempted, unlike the radiation-pair criterion above. It reads
     * ci->stars.h_max directly (cell_unskip.c:2593-2594 calls it
     * unconditionally on every pair task), so stars.h_max must feed
     * h_max_no_hii via space_regrid_star_h_max_no_hii_term_for() -- not
     * only h_max_radiation via space_regrid_star_radiation_term_for() --
     * or this criterion can demand a rebuild space_regrid() never
     * satisfies and never aborts on either, since h_max_no_hii never sees
     * the growth. h_hii is irrelevant to this criterion, so it is held at
     * 0 throughout. */
    {
      int fired = 0;
      for (int i = 0; i <= 290; ++i) {
        const float stars_h_max = 0.016f + 0.0001f * (float)i;

        setup_pair(&ci, &cj, (float)dmin, /*stars_h_max=*/stars_h_max,
                   /*stars_h_hii_max=*/0.f, /*hydro_h_max=*/0.01f,
                   /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);

        if (!cell_need_rebuild_for_stars_pair(&ci, &cj)) continue;
        ++fired;

        /* What the next regrid would accumulate: stars.h_max folded into
         * BOTH terms, exactly as space_regrid()'s three sites now do. */
        const float h_max_no_hii_sample =
            space_regrid_star_h_max_no_hii_term_for(h_max_floor, stars_h_max);
        const float star_radiation_term = space_regrid_star_radiation_term_for(
            /*h_hii_max=*/0.f, stars_h_max, cell_max_width);
        const double width = space_regrid_cell_width_for(
            h_max_no_hii_sample, star_radiation_term, cell_min, cell_max_width);
        if (!(width > 0.))
          error("Regrid width is not a positive number (%g) at h_max = %g.",
                width, stars_h_max);
        const int cd = (int)floor(dim / width);

        /* The precise invariant: unlike 8a/8b (where "cdim changed" and
         * "the criterion is satisfied" coincide by construction of the
         * radiation cap), here they do NOT -- the radiation-term route
         * alone can still coarsen cdim away from maxtcells while
         * saturating below this criterion's own unfactored bound, which
         * scales differently (no radiation_search_radius_factor). So
         * check the bound directly: a fresh check against the produced
         * width must NOT fire again. Beyond cell_max_width (where no grid
         * can satisfy it), space_regrid()'s own "too few cells" abort
         * (h_max_no_hii's search radius > cell_max_width) must be the one
         * to catch it instead. */
        if (kernel_gamma * (double)stars_h_max <= cell_max_width) {
          if (!(kernel_gamma * (double)stars_h_max <= width))
            error(
                "Case 8c: cell_need_rebuild_for_stars_pair demands a "
                "rebuild at stars.h_max = %g (dmin = %g, no exemption "
                "exists for this criterion), but the regrid's own width "
                "(%g, cdim = %d) still does not satisfy the criterion's "
                "unfactored bound (kernel_gamma * h_max = %g). "
                "stars.h_max is not reaching h_max_no_hii, so the demand "
                "can never be satisfied and engine_rebuild() aborts with "
                "\"engine_unskip failed after a rebuild!\".",
                stars_h_max, dmin, width, cd,
                kernel_gamma * (double)stars_h_max);
        } else {
          /* Unsatisfiable by any grid (dmin is already pinned at
           * cell_max_width / mintcells): space_regrid()'s "too few cells"
           * abort requires h_max_no_hii's OWN search radius -- not just
           * the capped star_radiation_term -- to exceed cell_max_width.
           * If stars.h_max never reaches h_max_no_hii, this never
           * triggers either, and the unsatisfiable demand becomes
           * invisible: cdim stays pinned (via the radiation-term route
           * alone), no abort ever fires, and
           * cell_need_rebuild_for_stars_pair keeps demanding a rebuild
           * forever. */
          if (!(h_max_no_hii_sample * kernel_gamma * space_stretch >
                cell_max_width))
            error(
                "Case 8c: at stars.h_max = %g the pair criterion's own "
                "unfactored bound needs a search radius beyond "
                "cell_max_width (%g), which no grid can satisfy -- but "
                "h_max_no_hii (%g) does not reflect that growth, so "
                "space_regrid()'s \"too few cells\" abort will never fire "
                "either. The demand is unsatisfiable AND invisible: "
                "engine_rebuild() loops instead of aborting loudly.",
                stars_h_max, cell_max_width, h_max_no_hii_sample);
        }
      }
      if (fired == 0)
        error(
            "Case 8c swept stars.h_max across [0.016, 0.045] and no sample "
            "ever made cell_need_rebuild_for_stars_pair fire; the sweep "
            "range no longer exercises the invariant it is meant to "
            "check.");
    }
  }

  /* cell_boxes_touch_under_shift(): the geometric facing test that replaces
   * cell_split_pairs/space_getsid for asymmetric radiation pair descent
   * (scheduler_splittasks.c). Only loc/width matter -- everything else is
   * left zeroed. */
  {
    struct cell a, b;

    /* Case 1: equal-size cells sharing a face (abutting, not overlapping)
     * -- must count as touching. */
    bzero(&a, sizeof(struct cell));
    bzero(&b, sizeof(struct cell));
    a.loc[0] = 0.0;
    a.loc[1] = 0.0;
    a.loc[2] = 0.0;
    a.width[0] = a.width[1] = a.width[2] = 1.0;
    b.loc[0] = 1.0;
    b.loc[1] = 0.0;
    b.loc[2] = 0.0;
    b.width[0] = b.width[1] = b.width[2] = 1.0;
    double zero_shift[3] = {0.0, 0.0, 0.0};
    if (!cell_boxes_touch_under_shift(&a, &b, zero_shift))
      error("cell_boxes_touch_under_shift missed a face-abutting pair.");

    /* Case 2: a one-cell gap -- must NOT count as touching. */
    b.loc[0] = 2.0;
    if (cell_boxes_touch_under_shift(&a, &b, zero_shift))
      error(
          "cell_boxes_touch_under_shift flagged a pair separated by a "
          "full cell width as touching.");

    /* Case 3: periodic wrap -- b sits far away in real coordinates but its
     * periodic image (via shift) abuts a. Must only touch WITH the shift
     * applied, not without it. */
    b.loc[0] = 9.0;
    double wrap_shift[3] = {-10.0, 0.0, 0.0};
    if (!cell_boxes_touch_under_shift(&a, &b, wrap_shift))
      error(
          "cell_boxes_touch_under_shift missed a pair touching only "
          "through the periodic image.");
    if (cell_boxes_touch_under_shift(&a, &b, zero_shift))
      error(
          "cell_boxes_touch_under_shift touched a pair separated in real "
          "space without the periodic shift that brings them together.");

    /* Case 4: unequal-size facing progeny -- the actual asymmetric-pair
     * scenario. Parent P = [0,2] splits into progeny p_near = [0,1] and
     * p_far = [1,2]^3... using x only for clarity, y/z stay full-width so
     * only the x-axis discriminates. Partner Q = [2,4]^3 (unsplit, twice
     * p's width). p_far shares P's face with Q and must touch; p_near
     * does not. */
    struct cell p_near, p_far, Q;
    bzero(&p_near, sizeof(struct cell));
    bzero(&p_far, sizeof(struct cell));
    bzero(&Q, sizeof(struct cell));
    p_near.loc[0] = 0.0;
    p_far.loc[0] = 1.0;
    p_near.width[0] = p_far.width[0] = 1.0;
    p_near.width[1] = p_far.width[1] = 2.0;
    p_near.width[2] = p_far.width[2] = 2.0;
    Q.loc[0] = 2.0;
    Q.width[0] = Q.width[1] = Q.width[2] = 2.0;
    if (!cell_boxes_touch_under_shift(&p_far, &Q, zero_shift))
      error(
          "cell_boxes_touch_under_shift missed the fine progeny that "
          "actually faces the coarse partner (the asymmetric-pair case "
          "cell_split_pairs/space_getsid cannot handle).");
    if (cell_boxes_touch_under_shift(&p_near, &Q, zero_shift))
      error(
          "cell_boxes_touch_under_shift touched a fine progeny that does "
          "NOT face the coarse partner -- would create a spurious pair "
          "task.");

    /* Case 5: corner-only touch between unequal-size cells -- all three
     * axes abut exactly, no axis overlaps. Perturbing one axis by a full
     * cell width must break it. */
    struct cell corner_a, corner_b;
    bzero(&corner_a, sizeof(struct cell));
    bzero(&corner_b, sizeof(struct cell));
    corner_a.width[0] = corner_a.width[1] = corner_a.width[2] = 1.0;
    corner_b.width[0] = corner_b.width[1] = corner_b.width[2] = 2.0;
    corner_b.loc[0] = corner_b.loc[1] = corner_b.loc[2] = 1.0;
    if (!cell_boxes_touch_under_shift(&corner_a, &corner_b, zero_shift))
      error("cell_boxes_touch_under_shift missed a corner-only touch.");
    corner_b.loc[2] = 3.5; /* Gap on z alone -- no longer touching. */
    if (cell_boxes_touch_under_shift(&corner_a, &corner_b, zero_shift))
      error(
          "cell_boxes_touch_under_shift touched a pair with a real gap on "
          "one axis.");
  }

  test_clamped_search_request();

  return 0;
}
