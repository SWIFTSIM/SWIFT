/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Unit test of the gas-sink overlap test of sink formation.
 *
 * It checks:
 * - sink_prepare_part_sink_formation_sink_criteria(): overlap inside and
 *   outside the reach, the periodic wrap, dead sinks and gas that cannot form
 *   a sink.
 * - cell_can_recurse_in_pair_aperture_task() and
 *   cell_can_recurse_in_pair_sink_aperture_task(): the margin for the
 *   particle movement in the pair recursion.
 *
 * The overlap test only exists for the GEAR sink model. */

/* Config parameters. */
#include <config.h>

/* Standard headers. */
#include <string.h>

/* Local headers. */
#include "swift.h"

#ifdef SINK_GEAR

/* Run the overlap test for a gas particle at pos_p and a sink at pos_s. The
 * gas particle radius is r_cut and the sink radius is r_cut too. */
static int overlaps(const double pos_p[3], const double pos_s[3],
                    const double dim[3], const float r_cut,
                    const int can_form_sink, const double sink_age) {

  struct cosmology cosmo;
  memset(&cosmo, 0, sizeof(cosmo));
  cosmo.a = 1.;

  struct sink_props props;
  memset(&props, 0, sizeof(props));
  props.use_fixed_r_cut = 1;
  props.cut_off_radius = r_cut;
  props.age_threshold_unlimited = 100.;

  struct part p;
  struct xpart xp;
  struct sink s;
  memset(&p, 0, sizeof(p));
  memset(&xp, 0, sizeof(xp));
  memset(&s, 0, sizeof(s));

  for (int k = 0; k < 3; k++) {
    p.x[k] = pos_p[k];
    s.x[k] = pos_s[k];
  }
  p.sink_data.can_form_sink = can_form_sink;
  p.sink_data.is_overlapping_sink = 0;

  /* The accretion radius of the sink is h * kernel_gamma. */
  s.h = r_cut / kernel_gamma;

  /* Non-cosmological run: the age is time - birth time. */
  const double time = 1000.;
  s.birth_data.time = time - sink_age;

  sink_prepare_part_sink_formation_sink_criteria(
      /*e=*/NULL, &p, &xp, &s, /*with_cosmology=*/0, &cosmo, &props, time,
      /*r_acc_p=*/r_cut * cosmo.a, dim);

  return p.sink_data.is_overlapping_sink;
}

static void test_overlap_criteria(void) {

  const float r_cut = 0.02f;
  const double box[3] = {1., 1., 1.};
  const double no_box[3] = {0., 0., 0.};
  const double gas[3] = {0.5, 0.5, 0.5};

  /* Reach is 2 * r_cut = 0.04. Inside and outside. */
  const double sink_in[3] = {0.5 + 0.039, 0.5, 0.5};
  const double sink_out[3] = {0.5 + 0.041, 0.5, 0.5};
  if (!overlaps(gas, sink_in, no_box, r_cut, 1, 0.))
    error("A sink inside 2 * r_cut was not found.");
  if (overlaps(gas, sink_out, no_box, r_cut, 1, 0.))
    error("A sink outside 2 * r_cut was found.");

  /* Periodic wrap: the two particles are 0.01 apart across the box edge. */
  const double gas_edge[3] = {0.005, 0.5, 0.5};
  const double sink_edge[3] = {0.995, 0.5, 0.5};
  if (!overlaps(gas_edge, sink_edge, box, r_cut, 1, 0.))
    error("The periodic wrap did not find the sink across the box edge.");
  if (overlaps(gas_edge, sink_edge, no_box, r_cut, 1, 0.))
    error("A sink across the box edge was found in a non-periodic box.");

  /* Wrap along the other axes and with a far sink. */
  const double gas_y[3] = {0.5, 0.995, 0.5};
  const double sink_y[3] = {0.5, 0.005, 0.5};
  if (!overlaps(gas_y, sink_y, box, r_cut, 1, 0.))
    error("The periodic wrap along y did not find the sink.");
  const double sink_far[3] = {0.5, 0.5, 0.5};
  const double gas_far[3] = {0.005, 0.5, 0.5};
  if (overlaps(gas_far, sink_far, box, r_cut, 1, 0.))
    error("A far sink was found in a periodic box.");

  /* A dead sink (older than age_threshold_unlimited = 100) is ignored. */
  if (overlaps(gas, sink_in, no_box, r_cut, 1, 101.))
    error("A dead sink was found.");
  if (!overlaps(gas, sink_in, no_box, r_cut, 1, 99.))
    error("A sink younger than the age threshold was not found.");

  /* Gas that cannot form a sink is left alone. */
  if (overlaps(gas, sink_in, no_box, r_cut, 0, 0.))
    error("The overlap flag was set for gas that cannot form a sink.");
}

static struct cell make_bare_cell(const double dmin, const float dx_gas,
                                  const float dx_sink) {
  struct cell c;
  memset(&c, 0, sizeof(struct cell));
  c.dmin = dmin;
  c.hydro.dx_max_part_old = dx_gas;
  c.sinks.dx_max_part_old = dx_sink;
  return c;
}

static void test_recursion_margin(void) {

  const double dmin = 1.;
  const float r_cut = 0.2f;

  /* Gas-gas: r_cut + dx_i + dx_j < 0.5 * dmin. The margin is 0.3. */
  {
    struct cell ci = make_bare_cell(dmin, 0.10f, 0.f);
    struct cell cj = make_bare_cell(dmin, 0.15f, 0.f);
    if (!cell_can_recurse_in_pair_aperture_task(&ci, &cj, r_cut))
      error("Gas-gas: recursion refused below the limit.");
    cj = make_bare_cell(dmin, 0.25f, 0.f);
    if (cell_can_recurse_in_pair_aperture_task(&ci, &cj, r_cut))
      error("Gas-gas: recursion allowed above the limit.");
    if (cell_can_recurse_in_pair_aperture_task(&cj, &ci, r_cut))
      error("Gas-gas: recursion allowed above the limit (swapped).");
  }

  /* Gas-sink: 2 * r_cut + the worse of the two directions < 0.5 * dmin.
   * The margin is 0.1. */
  {
    /* Direction 1: gas of ci with sinks of cj = 0.03 + 0.04 = 0.07.
     * Direction 2: gas of cj with sinks of ci = 0.02 + 0.03 = 0.05. */
    struct cell ci = make_bare_cell(dmin, 0.03f, 0.03f);
    struct cell cj = make_bare_cell(dmin, 0.02f, 0.04f);
    if (!cell_can_recurse_in_pair_sink_aperture_task(&ci, &cj, r_cut))
      error("Gas-sink: recursion refused below the limit.");

    /* Direction 2 becomes the worse one: 0.09 + 0.03 = 0.12 > 0.1. */
    cj = make_bare_cell(dmin, 0.09f, 0.04f);
    if (cell_can_recurse_in_pair_sink_aperture_task(&ci, &cj, r_cut))
      error("Gas-sink: recursion allowed above the limit.");
    if (cell_can_recurse_in_pair_sink_aperture_task(&cj, &ci, r_cut))
      error("Gas-sink: recursion allowed above the limit (swapped).");

    /* No movement: 2 * r_cut = 0.4 < 0.5. And a radius that is too big. */
    ci = make_bare_cell(dmin, 0.f, 0.f);
    cj = make_bare_cell(dmin, 0.f, 0.f);
    if (!cell_can_recurse_in_pair_sink_aperture_task(&ci, &cj, r_cut))
      error("Gas-sink: recursion refused without movement.");
    if (cell_can_recurse_in_pair_sink_aperture_task(&ci, &cj, 0.26f))
      error("Gas-sink: recursion allowed with a radius above 0.25.");
  }
}

int main(int argc, char *argv[]) {

  test_overlap_criteria();
  test_recursion_margin();

  message("All sink overlap criteria and recursion margin checks passed.");
  return 0;
}

#else

int main(int argc, char *argv[]) {
  message("Skipped: this test needs the GEAR sink model.");
  return 0;
}

#endif /* SINK_GEAR */
