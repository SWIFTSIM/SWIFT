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
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

/* HDF5. */
#include <hdf5.h>

/* Local headers. */
#include "swift.h"

/* Not declared in radiation.h: unlike the other five radiation_read_*_array()
 * entry points, this one stays file-local to radiation_table_io.c (only this
 * test needs it from outside); declared here so the test below can exercise
 * it directly, matching testRadiationRebuildCheck.c's own
 * cell_set_super_radiation_subgrid() forward declaration. */
void radiation_read_grid_metadata(hid_t group_id,
                                  struct radiation_grid_metadata *grid);

/**
 * @brief Build a minimal, synthetic 2D ("M,Z") Data/Radiation group with
 * every attribute radiation_read_grid_metadata() requires, plus the two new
 * generic edge_policy_q_h_below/above and
 * edge_policy_mean_excess_energy_below/above attributes
 * pychem's write_h5_table_v2() always writes alongside the existing
 * per-variant ones.
 *
 * @param file_id Open HDF5 file id to create the group in.
 * @param source_value The group's "source" attribute value: an
 * arbitrary/unrecognised string (standing in for a future pychem mode such
 * as "parsec_spectral"), to prove radiation_read_grid_metadata() never
 * dispatches on it.
 * @param omit_attr If not NULL, this one attribute is skipped entirely
 * (used to build the negative, missing-required-attribute cases). Ignored
 * for a name this function does not itself write.
 * @return The newly created, still-open "Radiation" group id.
 */
static hid_t build_radiation_group(hid_t file_id, const char *source_value,
                                   const char *omit_attr) {

  const hid_t grp =
      H5Gcreate(file_id, "/Radiation", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (grp < 0) error("Failed to create the test 'Radiation' group.");

#define WRITE_STR(name, value)                             \
  if (omit_attr == NULL || strcmp(omit_attr, name) != 0) { \
    io_write_attribute_s(grp, name, value);                \
  }

  WRITE_STR("dimensionality", "M,Z");
  io_write_attribute_f(grp, "m0", -1.0f);
  io_write_attribute_f(grp, "dm", 0.2f);
  io_write_attribute_i(grp, "nm", 5);
  io_write_attribute_i(grp, "nz", 3);

  const float metallicity[3] = {0.001f, 0.005f, 0.02f};
  const hsize_t dim = 3;
  const hid_t h_space = H5Screate_simple(1, &dim, NULL);
  if (h_space < 0) error("Failed to create the test 'Metallicity' dataspace.");
  const hid_t h_dset = H5Dcreate(grp, "Metallicity", H5T_NATIVE_FLOAT, h_space,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_dset < 0) error("Failed to create the test 'Metallicity' dataset.");
  if (H5Dwrite(h_dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               metallicity) < 0)
    error("Failed to write the test 'Metallicity' dataset.");
  H5Dclose(h_dset);
  H5Sclose(h_space);

  /* below="zero"/above="constant" so the positive-case assertions below can
   * tell the luminosity/Q_H/mean-excess-energy edge policies apart from one
   * another and from a boundary_condition_const default. */
  WRITE_STR("edge_policy_luminosity_below", "zero");
  WRITE_STR("edge_policy_luminosity_above", "constant");
  WRITE_STR("edge_policy_q_h_below", "zero");
  WRITE_STR("edge_policy_q_h_above", "constant");
  WRITE_STR("edge_policy_mean_excess_energy_below", "constant");
  WRITE_STR("edge_policy_mean_excess_energy_above", "constant");

  WRITE_STR("source", source_value);

#undef WRITE_STR

  return grp;
}

/**
 * @brief Confirm radiation_read_grid_metadata() still error()s when one of
 * the two new generic required attributes is missing -- mandatory-attribute
 * enforcement must survive the removal of the "source"-keyed dispatch, just
 * no longer keyed on "source" itself.
 *
 * Runs the read in a forked child (mirroring testRadiationRebuildCheck.c's
 * own orphaned-link regression test) since error() aborts the process.
 *
 * @param omit_attr Name of the required attribute to omit from the
 * synthetic file.
 */
static void run_negative_case(const char *omit_attr) {

  char filename[256];
  snprintf(filename, sizeof(filename), "testRadiationGridMetadata_neg_%s.hdf5",
           omit_attr);
  unlink(filename);

  const hid_t file_id =
      H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id < 0) error("Failed to create test HDF5 file '%s'.", filename);
  const hid_t grp =
      build_radiation_group(file_id, "parsec_spectral", omit_attr);
  H5Gclose(grp);
  H5Fclose(file_id);

  const pid_t pid = fork();
  if (pid == 0) {
    /* Child process: silence stderr (the error() message is expected
     * output, not a test failure), then trigger it. */
    if (freopen("/dev/null", "w", stderr) == NULL) _exit(43);
    const hid_t read_file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    const hid_t read_grp = H5Gopen(read_file_id, "/Radiation", H5P_DEFAULT);
    struct radiation_grid_metadata grid;
    memset(&grid, 0, sizeof(grid));
    radiation_read_grid_metadata(read_grp, &grid);
    /* Reached only if the missing-attribute error did NOT fire -- signal
     * failure with a distinguishable exit code (real error() exits 1). */
    _exit(42);
  } else if (pid > 0) {
    int status;
    waitpid(pid, &status, 0);
    /* error() aborts via swift_abort(), which is exit(1) normally but
     * abort() (SIGABRT) under SWIFT_DEVELOP_MODE -- accept either. */
    const int exited_with_error = WIFEXITED(status) && WEXITSTATUS(status) == 1;
    const int aborted = WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
    if (!exited_with_error && !aborted)
      error(
          "radiation_read_grid_metadata failed to error() on a missing "
          "'%s' attribute (WIFEXITED=%d WEXITSTATUS=%d WIFSIGNALED=%d "
          "WTERMSIG=%d) -- mandatory-attribute enforcement must survive "
          "the removal of the 'source'-keyed dispatch.",
          omit_attr, WIFEXITED(status),
          WIFEXITED(status) ? WEXITSTATUS(status) : -1, WIFSIGNALED(status),
          WIFSIGNALED(status) ? WTERMSIG(status) : -1);
  } else {
    error("fork() failed in the missing-attribute regression test.");
  }

  unlink(filename);
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  /* Positive case: a full, valid 2D group carrying the two new generic
   * edge_policy_* attributes and an arbitrary/unknown "source" string (a
   * stand-in for a future pychem mode like "parsec_spectral"). Must load
   * successfully and must never dispatch behaviour on "source". */
  {
    const char *filename = "testRadiationGridMetadata.hdf5";
    unlink(filename);
    const hid_t file_id =
        H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) error("Failed to create test HDF5 file '%s'.", filename);
    const hid_t grp =
        build_radiation_group(file_id, "totally_bogus_source_mode", NULL);

    struct radiation_grid_metadata grid;
    memset(&grid, 0, sizeof(grid));
    radiation_read_grid_metadata(grp, &grid);

    if (!grid.is_2d)
      error("Expected a 2D ('M,Z') grid to be recognised as such.");
    if (grid.edge_policy_luminosity != boundary_condition_zero_const)
      error(
          "edge_policy_luminosity mismatch: expected "
          "boundary_condition_zero_const, got %d.",
          grid.edge_policy_luminosity);
    if (grid.edge_policy_q_h != boundary_condition_zero_const)
      error(
          "edge_policy_q_h mismatch: expected boundary_condition_zero_const "
          "(from the generic edge_policy_q_h_below/above attributes), got "
          "%d. An unrecognised 'source' value must not prevent this from "
          "loading correctly.",
          grid.edge_policy_q_h);
    if (grid.edge_policy_dot_e_excess != boundary_condition_const)
      error(
          "edge_policy_dot_e_excess mismatch: expected "
          "boundary_condition_const (from the generic "
          "edge_policy_mean_excess_energy_below/above attributes), got %d.",
          grid.edge_policy_dot_e_excess);

    free(grid.metallicity);
    H5Gclose(grp);
    H5Fclose(file_id);
    unlink(filename);
  }

  /* Negative cases: each of the two new generic required attributes must
   * still be enforced -- just no longer keyed on "source". io_read_attribute
   * ()/radiation_read_string_attribute() already error() loudly on a missing
   * attribute, so this only needs to confirm that behaviour still fires for
   * these two specific names now that the "source" strcmp/error dispatch is
   * gone. */
  run_negative_case("edge_policy_q_h_below");
  run_negative_case("edge_policy_mean_excess_energy_below");

  return 0;
}
