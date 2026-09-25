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
#ifndef SWIFT_FEEDBACK_IO_GEAR_H
#define SWIFT_FEEDBACK_IO_GEAR_H

#include "feedback.h"
#include "io_properties.h"

/**
 * @brief Snapshot converter for #IsIonizedFlags, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_is_ionized(const struct engine *e,
                                           const struct part *p,
                                           const struct xpart *xp, char *ret) {
  ret[0] = feedback_is_part_tagged_as_ionized(p, xp);
}

/**
 * @brief Snapshot converter for #HIIStarIDs, see #feedback_write_particles.
 */
INLINE static void convert_part_HII_star_id(const struct engine *e,
                                            const struct part *p,
                                            const struct xpart *xp,
                                            long long *ret) {
  ret[0] = feedback_get_part_ionized_star_id(p, xp);
}

/**
 * @brief Snapshot converter for #PESpecificEnergy, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_u_PE(const struct engine *e,
                                     const struct part *p,
                                     const struct xpart *xp, float *ret) {
  ret[0] = feedback_get_part_u_PE(p);
}

/**
 * @brief Snapshot converter for #LWSpecificEnergy, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_u_LW(const struct engine *e,
                                     const struct part *p,
                                     const struct xpart *xp, float *ret) {
  ret[0] = feedback_get_part_u_LW(p);
}

/**
 * @brief Snapshot converter for #LWPhotonSpecificEnergies, see
 * #feedback_write_particles. Distinct from the RT schemes' PhotonEnergies
 * (src/rt/GEAR/rt_io.h, src/rt/SPHM1RT/rt_io.h): those are raw per-group
 * energies, not mass-specific and not band-prefixed.
 */
INLINE static void convert_part_u_LW_PHOTON(const struct engine *e,
                                            const struct part *p,
                                            const struct xpart *xp,
                                            float *ret) {
  ret[0] = feedback_get_part_u_LW_PHOTON(p);
}

/**
 * @brief Snapshot converter for #PEArtificialDissipationCoefficients, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_dissipation_alpha_PE(const struct engine *e,
                                                     const struct part *p,
                                                     const struct xpart *xp,
                                                     float *ret) {
  ret[0] = feedback_get_part_dissipation_alpha_PE(p);
}

/**
 * @brief Snapshot converter for #LWArtificialDissipationCoefficients, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_dissipation_alpha_LW(const struct engine *e,
                                                     const struct part *p,
                                                     const struct xpart *xp,
                                                     float *ret) {
  ret[0] = feedback_get_part_dissipation_alpha_LW(p);
}

/**
 * @brief Snapshot converter for #LWPhotonArtificialDissipationCoefficients,
 * see #feedback_write_particles. Decorative: #ISRF_MOMENT_LW_PHOTON shares
 * #ISRF_OPERATOR_LW with #ISRF_MOMENT_LW, so this reads the identical
 * operator field as #convert_part_dissipation_alpha_LW.
 */
INLINE static void convert_part_dissipation_alpha_LW_PHOTON(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  ret[0] = feedback_get_part_dissipation_alpha_LW_PHOTON(p);
}

/**
 * @brief Snapshot converter for #PESpecificFluxDivergences, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_div_specific_flux_PE(const struct engine *e,
                                                     const struct part *p,
                                                     const struct xpart *xp,
                                                     float *ret) {
  ret[0] = feedback_get_part_div_specific_flux_PE(p);
}

/**
 * @brief Snapshot converter for #LWSpecificFluxDivergences, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_div_specific_flux_LW(const struct engine *e,
                                                     const struct part *p,
                                                     const struct xpart *xp,
                                                     float *ret) {
  ret[0] = feedback_get_part_div_specific_flux_LW(p);
}

/**
 * @brief Snapshot converter for #LWPhotonSpecificFluxDivergences, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_div_specific_flux_LW_PHOTON(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  ret[0] = feedback_get_part_div_specific_flux_LW_PHOTON(p);
}

/**
 * @brief Snapshot converter for #PESpecificFluxes, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_specific_flux_PE(const struct engine *e,
                                                 const struct part *p,
                                                 const struct xpart *xp,
                                                 float *ret) {
  feedback_get_part_specific_flux_PE(p, ret);
}

/**
 * @brief Snapshot converter for #LWSpecificFluxes, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_specific_flux_LW(const struct engine *e,
                                                 const struct part *p,
                                                 const struct xpart *xp,
                                                 float *ret) {
  feedback_get_part_specific_flux_LW(p, ret);
}

/**
 * @brief Snapshot converter for #LWPhotonSpecificFluxes, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_specific_flux_LW_PHOTON(const struct engine *e,
                                                        const struct part *p,
                                                        const struct xpart *xp,
                                                        float *ret) {
  feedback_get_part_specific_flux_LW_PHOTON(p, ret);
}

/**
 * @brief Snapshot converter for #PEMinimumSpecificEnergies, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_u_min_since_snapshot_PE(const struct engine *e,
                                                        const struct part *p,
                                                        const struct xpart *xp,
                                                        float *ret) {
  ret[0] = feedback_get_part_u_min_since_snapshot_PE(p, e);
}

/**
 * @brief Snapshot converter for #LWMinimumSpecificEnergies, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_u_min_since_snapshot_LW(const struct engine *e,
                                                        const struct part *p,
                                                        const struct xpart *xp,
                                                        float *ret) {
  ret[0] = feedback_get_part_u_min_since_snapshot_LW(p, e);
}

/**
 * @brief Snapshot converter for #LWPhotonMinimumSpecificEnergies, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_u_min_since_snapshot_LW_PHOTON(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  ret[0] = feedback_get_part_u_min_since_snapshot_LW_PHOTON(p, e);
}

/**
 * @brief Snapshot converter for #PECumulativeInjectedSpecificEnergies, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_cumulative_injected_PE(const struct engine *e,
                                                       const struct part *p,
                                                       const struct xpart *xp,
                                                       float *ret) {
  ret[0] = feedback_get_part_cumulative_injected_PE(p);
}

/**
 * @brief Snapshot converter for #LWCumulativeInjectedSpecificEnergies, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_cumulative_injected_LW(const struct engine *e,
                                                       const struct part *p,
                                                       const struct xpart *xp,
                                                       float *ret) {
  ret[0] = feedback_get_part_cumulative_injected_LW(p);
}

/**
 * @brief Snapshot converter for #LWPhotonCumulativeInjectedSpecificEnergies,
 * see #feedback_write_particles.
 */
INLINE static void convert_part_cumulative_injected_LW_PHOTON(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  ret[0] = feedback_get_part_cumulative_injected_LW_PHOTON(p);
}

/**
 * @brief Snapshot converter for #PECumulativeAbsorbedSpecificEnergies, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_cumulative_absorbed_PE(const struct engine *e,
                                                       const struct part *p,
                                                       const struct xpart *xp,
                                                       float *ret) {
  ret[0] = feedback_get_part_cumulative_absorbed_PE(p);
}

/**
 * @brief Snapshot converter for #LWCumulativeAbsorbedSpecificEnergies, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_cumulative_absorbed_LW(const struct engine *e,
                                                       const struct part *p,
                                                       const struct xpart *xp,
                                                       float *ret) {
  ret[0] = feedback_get_part_cumulative_absorbed_LW(p);
}

/**
 * @brief Snapshot converter for #LWPhotonCumulativeAbsorbedSpecificEnergies,
 * see #feedback_write_particles.
 */
INLINE static void convert_part_cumulative_absorbed_LW_PHOTON(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  ret[0] = feedback_get_part_cumulative_absorbed_LW_PHOTON(p);
}

/**
 * @brief Snapshot converter for #HyperbolicPropagationSpeeds, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_c_hyp(const struct engine *e,
                                      const struct part *p,
                                      const struct xpart *xp, float *ret) {
  ret[0] = feedback_get_part_c_hyp(p);
}

/**
 * @brief Specifies which particle fields to read from the ICs.
 *
 * Reads "PESpecificEnergy"/"LWSpecificEnergy" (singular: the codebase's
 * own IC-read convention, e.g. hydro's "Density"/"SmoothingLength", vs. the
 * plural "Densities"/"SmoothingLengths" used for the matching snapshot
 * output) as OPTIONAL fields into #feedback_isrf_moment_data.u.
 * The snapshot output for the same quantity
 * (#convert_part_u_PE/convert_part_u_LW, this file)
 * deliberately uses the plural "PESpecificEnergies"/"LWSpecificEnergies"
 * instead, so a snapshot cannot be fed back in as an IC unmodified. This
 * is a validation/testing tool, not a normal production IC input, and is
 * not meant to make that round-trip easy. It lets a test set an arbitrary,
 * analytically-known initial LW/PE field shape (a value range across
 * particles, a pulse, a step) and watch only Grackle's chemistry, or only
 * the LW/PE propagation PDE, evolve it, decoupled from the star and
 * injection machinery. It only makes physical sense in a run with no star,
 * though this reader does not itself enforce that.
 *
 * An IC without these fields is unaffected: #radiation_first_init_part no
 * longer zeroes #feedback_isrf_moment_data.u (and seeds
 * #feedback_isrf_moment_data.u_prev from it, not from 0.f, so a supplied value
 * also survives the very first propagation update when
 * `GEARFeedback:ISRF_propagation` is on) so that a supplied value survives
 * first-init, but every #part is bzero'd before this read runs
 * (single_io.c/parallel_io.c/serial_io.c), so a missing field still leaves
 * exactly 0.f, matching pre-existing behaviour.
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int feedback_read_particles(struct part *parts,
                                          struct io_props *list) {

  /* Both are PHYSICAL, mass-specific quantities: no scale-factor exponent
     beyond what UNIT_CONV_ENERGY_PER_UNIT_MASS implies, and an input field
     carries no a-exponent slot at all, so an IC value is taken verbatim. */

  list[0] = io_make_input_field("PESpecificEnergy", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                                feedback_data.isrf_moment[ISRF_MOMENT_PE].u);
  list[1] = io_make_input_field("LWSpecificEnergy", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                                feedback_data.isrf_moment[ISRF_MOMENT_LW].u);

  /* LWPhotonSpecificEnergy = 0 together with a nonzero LWSpecificEnergy is
     not representable after first-init: radiation_first_init_part()
     overwrites it with the LW value, for either sign, since a seeded LW
     field with no attribution is a field at the reference photon energy by
     definition. An IC author who wants no photon moment carried must also
     zero LWSpecificEnergy. */
  list[2] =
      io_make_input_field("LWPhotonSpecificEnergy", FLOAT, 1, OPTIONAL,
                          UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                          feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].u);

  return 3;
}

/**
 * @brief Specifies which particle fields to write to a dataset.
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology switched on?
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int feedback_write_particles(
    const struct part *parts, const struct xpart *xparts, struct io_props *list,
    const int with_cosmology) {

  int num = 24;

  list[0] = io_make_output_field_convert_part(
      "IsIonizedFlags", CHAR, 1, UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
      convert_part_is_ionized,
      "Were the particles flagged as ionized by HII ionzation subgrid model?");

  list[1] = io_make_output_field_convert_part(
      "HIIStarIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
      convert_part_HII_star_id,
      "Star particle IDs that ionized these gas particles due to HII ionzation "
      "subgrid model?");

  list[2] = io_make_output_field_convert_part(
      "PESpecificEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f,
      parts, xparts, convert_part_u_PE,
      "Local specific PE-band (6-11.2 eV) interstellar radiation field. "
      "Physical, mass-specific: no scale-factor exponent of its own.");

  list[3] = io_make_output_field_convert_part(
      "LWSpecificEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f,
      parts, xparts, convert_part_u_LW,
      "Local specific Lyman-Werner-band (11.2-13.6 eV) interstellar "
      "radiation field.");

  list[4] = io_make_output_field_convert_part(
      "PEArtificialDissipationCoefficients", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
      parts, xparts, convert_part_dissipation_alpha_PE,
      "Negativity-triggered artificial-dissipation coefficient of the "
      "PE-band hyperbolic propagation, in "
      "[0, max(ISRF_dissipation_alpha_max, ISRF_dissipation_alpha_floor)]. "
      "Only meaningful when ISRF_propagation is on.");

  list[5] = io_make_output_field_convert_part(
      "LWArtificialDissipationCoefficients", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
      parts, xparts, convert_part_dissipation_alpha_LW,
      "Same as PEArtificialDissipationCoefficients, Lyman-Werner band.");

  list[6] = io_make_output_field_convert_part(
      "PESpecificFluxDivergences", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS_PER_TIME, 0.f, parts, xparts,
      convert_part_div_specific_flux_PE,
      "`(1/rho) div(rho F)` accumulator of the PE-band hyperbolic "
      "propagation, accumulated in the force loop from the step's relaxed "
      "flux. Physical, like the "
      "specific energy it is a rate of change of. Only meaningful when "
      "ISRF_propagation is on.");

  list[7] = io_make_output_field_convert_part(
      "LWSpecificFluxDivergences", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS_PER_TIME, 0.f, parts, xparts,
      convert_part_div_specific_flux_LW,
      "Same as PESpecificFluxDivergences, Lyman-Werner band.");

  list[8] = io_make_output_field_convert_part(
      "PESpecificFluxes", FLOAT, 3, UNIT_CONV_ENERGY_PER_UNIT_MASS_VELOCITY,
      0.f, parts, xparts, convert_part_specific_flux_PE,
      "Tracked specific flux moment of the PE-band hyperbolic propagation, "
      "mass-specific like PESpecificEnergies. Physical: no scale-factor "
      "exponent of its own. Only meaningful when ISRF_propagation is on.");

  list[9] = io_make_output_field_convert_part(
      "LWSpecificFluxes", FLOAT, 3, UNIT_CONV_ENERGY_PER_UNIT_MASS_VELOCITY,
      0.f, parts, xparts, convert_part_specific_flux_LW,
      "Same as PESpecificFluxes, Lyman-Werner band.");

  list[10] = io_make_output_field_convert_part(
      "PEMinimumSpecificEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS,
      0.f, parts, xparts, convert_part_u_min_since_snapshot_PE,
      "Most negative PESpecificEnergies value the propagation update wrote "
      "since the previous snapshot, 0 if none was negative. The number of "
      "nonzero entries is the count of particles that undershot. The "
      "interval is since the last increment of engine.snapshot_output_count, "
      "which also happens when a FOF seeding catalogue is dumped "
      "(FOF:dump_catalogue_when_seeding), not only at a real snapshot. "
      "Always 0 unless the code is configured with --enable-debugging-checks.");

  list[11] = io_make_output_field_convert_part(
      "LWMinimumSpecificEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS,
      0.f, parts, xparts, convert_part_u_min_since_snapshot_LW,
      "Same as PEMinimumSpecificEnergies, Lyman-Werner band.");

  list[12] = io_make_output_field_convert_part(
      "PECumulativeInjectedSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_cumulative_injected_PE,
      "Cumulative mass-specific PE-band dose this particle has drawn from "
      "the dose reservoir since first init, rescaled by c_hyp/c exactly as "
      "PESpecificEnergies' own update rescales it, but not relaxed by the "
      "per-step phi factor: the raw amount attempted every step, summed. "
      "Energy-conservation diagnostic (with "
      "PECumulativeAbsorbedSpecificEnergies and PESpecificEnergies); "
      "always 0 unless the code is configured with "
      "--enable-debugging-checks.");

  list[13] = io_make_output_field_convert_part(
      "LWCumulativeInjectedSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_cumulative_injected_LW,
      "Same as PECumulativeInjectedSpecificEnergies, Lyman-Werner band.");

  list[14] = io_make_output_field_convert_part(
      "PECumulativeAbsorbedSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_cumulative_absorbed_PE,
      "Cumulative mass-specific PE-band energy PESpecificEnergies' own "
      "exact-relaxation update has attributed to decay (dust absorption "
      "and the cosmological redshift term) plus the fraction of each "
      "step's source and transport terms that never reached the field "
      "because the step was optically thick, summed since first init. "
      "Energy-conservation diagnostic: PESpecificEnergies plus this field "
      "minus PECumulativeInjectedSpecificEnergies isolates the transport "
      "and artificial-dissipation residual the closed-form split does not "
      "attribute to either term, which the SPH divergence's kernel-sum "
      "identity and the dissipation's pairwise antisymmetry drive towards "
      "0 when summed over every particle. Always 0 unless the code is "
      "configured with --enable-debugging-checks.");

  list[15] = io_make_output_field_convert_part(
      "LWCumulativeAbsorbedSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_cumulative_absorbed_LW,
      "Same as PECumulativeAbsorbedSpecificEnergies, Lyman-Werner band.");

  list[16] = io_make_output_field_convert_part(
      "HyperbolicPropagationSpeeds", FLOAT, 1, UNIT_CONV_SPEED, 0.f, parts,
      xparts, convert_part_c_hyp,
      "Kernel-local hyperbolic propagation speed the band updates and the "
      "pairwise transport operators ran with, shared by both bands. "
      "Physical: built from the physical smoothing length and a physical "
      "timestep, so no scale-factor exponent of its own. The conserved "
      "ledger of the consistent-variable-c schemes is `sum m u / c_hyp` "
      "rather than `sum m u`, which is what this field makes measurable "
      "from a snapshot. Only meaningful when ISRF_propagation is on.");

  list[17] = io_make_output_field_convert_part(
      "LWPhotonSpecificEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f,
      parts, xparts, convert_part_u_LW_PHOTON,
      "Lyman-Werner-band photon-number moment, energy-equivalent at a fixed "
      "reference photon energy: see LWSpecificEnergies. Distinct from the "
      "RT schemes' PhotonEnergies, which are raw per-group energies, not "
      "mass-specific and not band-prefixed.");

  list[18] = io_make_output_field_convert_part(
      "LWPhotonArtificialDissipationCoefficients", FLOAT, 1, UNIT_CONV_NO_UNITS,
      0.f, parts, xparts, convert_part_dissipation_alpha_LW_PHOTON,
      "Same as LWArtificialDissipationCoefficients: the photon-number "
      "moment shares the Lyman-Werner operator, so these two fields always "
      "read numerically identical.");

  list[19] = io_make_output_field_convert_part(
      "LWPhotonSpecificFluxDivergences", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS_PER_TIME, 0.f, parts, xparts,
      convert_part_div_specific_flux_LW_PHOTON,
      "Same as LWSpecificFluxDivergences, Lyman-Werner-band photon-number "
      "moment.");

  list[20] = io_make_output_field_convert_part(
      "LWPhotonSpecificFluxes", FLOAT, 3,
      UNIT_CONV_ENERGY_PER_UNIT_MASS_VELOCITY, 0.f, parts, xparts,
      convert_part_specific_flux_LW_PHOTON,
      "Same as LWSpecificFluxes, Lyman-Werner-band photon-number moment. "
      "Distinct from the RT schemes' PhotonFluxes, which are raw per-group "
      "fluxes, not mass-specific and not band-prefixed.");

  list[21] = io_make_output_field_convert_part(
      "LWPhotonMinimumSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_u_min_since_snapshot_LW_PHOTON,
      "Same as LWMinimumSpecificEnergies, Lyman-Werner-band photon-number "
      "moment. Always 0 unless the code is configured with "
      "--enable-debugging-checks.");

  list[22] = io_make_output_field_convert_part(
      "LWPhotonCumulativeInjectedSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_cumulative_injected_LW_PHOTON,
      "Same as LWCumulativeInjectedSpecificEnergies, Lyman-Werner-band "
      "photon-number moment. Always 0 unless the code is configured with "
      "--enable-debugging-checks.");

  list[23] = io_make_output_field_convert_part(
      "LWPhotonCumulativeAbsorbedSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_cumulative_absorbed_LW_PHOTON,
      "Same as LWCumulativeAbsorbedSpecificEnergies, Lyman-Werner-band "
      "photon-number moment. Always 0 unless the code is configured with "
      "--enable-debugging-checks.");

  return num;
}

/**
 * @brief Specifies which star particle fields to write to a dataset.
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology switched on?
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int feedback_write_sparticles(
    const struct spart *sparts, struct io_props *list,
    const int with_cosmology) {

  int num = 2;

  list[0] = io_make_output_field(
      "FinalHIIRegionRadii", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, sparts,
      feedback_data.radiation.final_HII_radius,
      "Co-moving HII region radius of the star particles before they die or "
      "were not eligible to form HII regions anymore. Same algorithm's "
      "bookkeeping caveat as the live HIIRegionRadii it is retired from.");

  list[1] = io_make_output_field(
      "FinalHIIRegionMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      feedback_data.radiation.final_HII_mass,
      "Ionized gas mass of the star particles' HII region before they die or "
      "were not eligible to form HII regions anymore. Same algorithm's "
      "bookkeeping caveat as the live HIIRegionMasses it is retired from.");

  return num;
}

/**
 * @brief Specifies which star particle fields to read from the ICs or a
 * restart file. Empty by design, not by omission: no feedback-owned #spart
 * field is currently read from either source (see #stars_read_particles).
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int feedback_read_sparticles(struct spart *sparts,
                                           struct io_props *list) {
  return 0;
}

#endif /* SWIFT_FEEDBACK_IO_GEAR_H */
