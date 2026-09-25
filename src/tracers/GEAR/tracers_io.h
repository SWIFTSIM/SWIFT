/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_TRACERS_NONE_IO_H
#define SWIFT_TRACERS_NONE_IO_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "feedback.h"
#include "io_properties.h"
#include "tracers.h"

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of tracers to the file.
 *
 * @param h_grp The HDF5 group in which to write
 * @param tracers The #tracers_function_data
 */
__attribute__((always_inline)) INLINE static void tracers_write_flavour(
    hid_t h_grp) {

  io_write_attribute_s(h_grp, "Tracers", "GEAR");
}
#endif

INLINE static void convert_sink_averaged_SFR(const struct engine *e,
                                             const struct sink *sink,
                                             float *ret) {

  for (int i = 0; i < num_snapshot_triggers_sink; ++i) {
    if (e->snapshot_recording_triggers_started_sink[i]) {
      ret[i] = sink->tracers_data.averaged_SFR[i] /
               e->snapshot_recording_triggers_sink[i];
    } else {
      ret[i] = 0.f;
    }
  }
}

INLINE static void convert_bpart_averaged_accretion_rate(const struct engine *e,
                                                         const struct bpart *bp,
                                                         float *ret) {

  for (int i = 0; i < num_snapshot_triggers_bpart; ++i) {
    if (e->snapshot_recording_triggers_started_bpart[i]) {
      ret[i] = bp->tracers_data.averaged_accretion_rate[i] /
               e->snapshot_recording_triggers_bpart[i];
    } else {
      ret[i] = 0.f;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (ret[i] < 0.f)
      error(
          "Negative averaged accretion rate for black hole id=%lld "
          "trigger=%d value=%e",
          bp->id, i, ret[i]);
#endif
  }
}

INLINE static void convert_sink_averaged_accretion_rate(const struct engine *e,
                                                        const struct sink *sink,
                                                        float *ret) {

  for (int i = 0; i < num_snapshot_triggers_sink; ++i) {
    if (e->snapshot_recording_triggers_started_sink[i]) {
      ret[i] = sink->tracers_data.averaged_accretion_rate[i] /
               e->snapshot_recording_triggers_sink[i];
    } else {
      ret[i] = 0.f;
    }
  }
}

/**
 * @brief Snapshot converter for #IsIonizedFlags, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_is_ionized(const struct engine *e,
                                           const struct part *p,
                                           const struct xpart *xp, char *ret) {
  ret[0] = feedback_is_part_tagged_as_ionized(p, xp);
}

/**
 * @brief Snapshot converter for #HIIStarIDs, see #tracers_write_particles.
 */
INLINE static void convert_part_HII_star_id(const struct engine *e,
                                            const struct part *p,
                                            const struct xpart *xp,
                                            long long *ret) {
  ret[0] = feedback_get_part_ionized_star_id(p, xp);
}

/**
 * @brief Snapshot converter for #PESpecificEnergy, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_u_PE(const struct engine *e,
                                     const struct part *p,
                                     const struct xpart *xp, float *ret) {
  ret[0] = feedback_get_part_u_PE(p);
}

/**
 * @brief Snapshot converter for #LWSpecificEnergy, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_u_LW(const struct engine *e,
                                     const struct part *p,
                                     const struct xpart *xp, float *ret) {
  ret[0] = feedback_get_part_u_LW(p);
}

/**
 * @brief Snapshot converter for #LWPhotonSpecificEnergies, see
 * #tracers_write_particles. Distinct from the RT schemes' PhotonEnergies
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
 * #tracers_write_particles.
 */
INLINE static void convert_part_dissipation_alpha_PE(const struct engine *e,
                                                     const struct part *p,
                                                     const struct xpart *xp,
                                                     float *ret) {
  ret[0] = feedback_get_part_dissipation_alpha_PE(p);
}

/**
 * @brief Snapshot converter for #LWArtificialDissipationCoefficients, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_dissipation_alpha_LW(const struct engine *e,
                                                     const struct part *p,
                                                     const struct xpart *xp,
                                                     float *ret) {
  ret[0] = feedback_get_part_dissipation_alpha_LW(p);
}

/**
 * @brief Snapshot converter for #LWPhotonArtificialDissipationCoefficients,
 * see #tracers_write_particles. Decorative: #ISRF_MOMENT_LW_PHOTON shares
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
 * #tracers_write_particles.
 */
INLINE static void convert_part_div_specific_flux_PE(const struct engine *e,
                                                     const struct part *p,
                                                     const struct xpart *xp,
                                                     float *ret) {
  ret[0] = feedback_get_part_div_specific_flux_PE(p);
}

/**
 * @brief Snapshot converter for #LWSpecificFluxDivergences, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_div_specific_flux_LW(const struct engine *e,
                                                     const struct part *p,
                                                     const struct xpart *xp,
                                                     float *ret) {
  ret[0] = feedback_get_part_div_specific_flux_LW(p);
}

/**
 * @brief Snapshot converter for #LWPhotonSpecificFluxDivergences, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_div_specific_flux_LW_PHOTON(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  ret[0] = feedback_get_part_div_specific_flux_LW_PHOTON(p);
}

/**
 * @brief Snapshot converter for #PESpecificFluxes, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_specific_flux_PE(const struct engine *e,
                                                 const struct part *p,
                                                 const struct xpart *xp,
                                                 float *ret) {
  feedback_get_part_specific_flux_PE(p, ret);
}

/**
 * @brief Snapshot converter for #LWSpecificFluxes, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_specific_flux_LW(const struct engine *e,
                                                 const struct part *p,
                                                 const struct xpart *xp,
                                                 float *ret) {
  feedback_get_part_specific_flux_LW(p, ret);
}

/**
 * @brief Snapshot converter for #LWPhotonSpecificFluxes, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_specific_flux_LW_PHOTON(const struct engine *e,
                                                        const struct part *p,
                                                        const struct xpart *xp,
                                                        float *ret) {
  feedback_get_part_specific_flux_LW_PHOTON(p, ret);
}

/**
 * @brief Snapshot converter for #PEMinimumSpecificEnergies, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_u_min_since_snapshot_PE(const struct engine *e,
                                                        const struct part *p,
                                                        const struct xpart *xp,
                                                        float *ret) {
  ret[0] = feedback_get_part_u_min_since_snapshot_PE(p, e);
}

/**
 * @brief Snapshot converter for #LWMinimumSpecificEnergies, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_u_min_since_snapshot_LW(const struct engine *e,
                                                        const struct part *p,
                                                        const struct xpart *xp,
                                                        float *ret) {
  ret[0] = feedback_get_part_u_min_since_snapshot_LW(p, e);
}

/**
 * @brief Snapshot converter for #LWPhotonMinimumSpecificEnergies, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_u_min_since_snapshot_LW_PHOTON(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  ret[0] = feedback_get_part_u_min_since_snapshot_LW_PHOTON(p, e);
}

/**
 * @brief Snapshot converter for #PECumulativeInjectedSpecificEnergies, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_cumulative_injected_PE(const struct engine *e,
                                                       const struct part *p,
                                                       const struct xpart *xp,
                                                       float *ret) {
  ret[0] = feedback_get_part_cumulative_injected_PE(p);
}

/**
 * @brief Snapshot converter for #LWCumulativeInjectedSpecificEnergies, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_cumulative_injected_LW(const struct engine *e,
                                                       const struct part *p,
                                                       const struct xpart *xp,
                                                       float *ret) {
  ret[0] = feedback_get_part_cumulative_injected_LW(p);
}

/**
 * @brief Snapshot converter for #LWPhotonCumulativeInjectedSpecificEnergies,
 * see #tracers_write_particles.
 */
INLINE static void convert_part_cumulative_injected_LW_PHOTON(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  ret[0] = feedback_get_part_cumulative_injected_LW_PHOTON(p);
}

/**
 * @brief Snapshot converter for #PECumulativeAbsorbedSpecificEnergies, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_cumulative_absorbed_PE(const struct engine *e,
                                                       const struct part *p,
                                                       const struct xpart *xp,
                                                       float *ret) {
  ret[0] = feedback_get_part_cumulative_absorbed_PE(p);
}

/**
 * @brief Snapshot converter for #LWCumulativeAbsorbedSpecificEnergies, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_cumulative_absorbed_LW(const struct engine *e,
                                                       const struct part *p,
                                                       const struct xpart *xp,
                                                       float *ret) {
  ret[0] = feedback_get_part_cumulative_absorbed_LW(p);
}

/**
 * @brief Snapshot converter for #LWPhotonCumulativeAbsorbedSpecificEnergies,
 * see #tracers_write_particles.
 */
INLINE static void convert_part_cumulative_absorbed_LW_PHOTON(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  ret[0] = feedback_get_part_cumulative_absorbed_LW_PHOTON(p);
}

/**
 * @brief Snapshot converter for #HyperbolicPropagationSpeeds, see
 * #tracers_write_particles.
 */
INLINE static void convert_part_c_hyp(const struct engine *e,
                                      const struct part *p,
                                      const struct xpart *xp, float *ret) {
  ret[0] = feedback_get_part_c_hyp(p);
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int tracers_write_particles(
    const struct part *parts, const struct xpart *xparts, struct io_props *list,
    const int with_cosmology) {

  int num = 32;

  /* The tag core (is_ionized/star_id) lives on struct part's feedback_data,
     not tracers_xpart_data, so read it through the feedback-model dispatch
     wrapper: this must compile under any --with-feedback choice paired
     with --with-tracers=GEAR, not only --with-feedback=GEAR. */
  list[0] = io_make_output_field_convert_part(
      "IsIonizedFlags", CHAR, 1, UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
      convert_part_is_ionized,
      "Were the particles flagged as ionized by HII ionzation subgrid model?");

  list[1] = io_make_output_field_convert_part(
      "HIIStarIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
      convert_part_HII_star_id,
      "Star particle IDs that ionized these gas particles due to HII ionzation "
      "subgrid model?");

  list[2] = io_make_physical_output_field(
      "CumulativeMomentumFromSupernovae", FLOAT, 1, UNIT_CONV_MOMENTUM, 0.f,
      xparts, tracers_data.feedback_cumulative.momentum_supernovae,
      /*can convert to comoving=*/0,
      "Cumulative |delta_p| per event from supernovae over this particle's "
      "lifetime (scalar sum, not vector: isotropic kicks would else "
      "cancel).");

  list[3] = io_make_physical_output_field(
      "CumulativeMomentumFromWinds", FLOAT, 1, UNIT_CONV_MOMENTUM, 0.f, xparts,
      tracers_data.feedback_cumulative.momentum_winds,
      /*can convert to comoving=*/0,
      "Same convention as CumulativeMomentumFromSupernovae, for stellar "
      "winds.");

  list[4] = io_make_physical_output_field(
      "CumulativeMomentumFromRadiationPressure", FLOAT, 1, UNIT_CONV_MOMENTUM,
      0.f, xparts, tracers_data.feedback_cumulative.momentum_radiation,
      /*can convert to comoving=*/0,
      "Same convention as CumulativeMomentumFromSupernovae, for radiation "
      "pressure.");

  list[5] = io_make_physical_output_field(
      "CumulativeEnergyFromSupernovae", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, xparts,
      tracers_data.feedback_cumulative.energy_supernovae,
      /*can convert to comoving=*/0,
      "Cumulative specific internal energy received from supernovae over "
      "this particle's lifetime.");

  list[6] = io_make_physical_output_field(
      "CumulativeEnergyFromWinds", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS,
      0.f, xparts, tracers_data.feedback_cumulative.energy_winds,
      /*can convert to comoving=*/0,
      "Cumulative specific internal energy received from stellar winds. "
      "A conservation residual, not strictly positive: can go negative "
      "when the gas was already moving towards the star before the kick.");

  list[7] = io_make_physical_output_field(
      "MaxKickVelocityFromSupernovae", FLOAT, 1, UNIT_CONV_SPEED, 0.f, xparts,
      tracers_data.feedback_cumulative.max_kick_velocity_supernovae,
      /*can convert to comoving=*/0,
      "Largest single-event kick velocity this particle received from "
      "supernovae (outflow diagnostic).");

  list[8] = io_make_physical_output_field(
      "MaxKickVelocityFromWinds", FLOAT, 1, UNIT_CONV_SPEED, 0.f, xparts,
      tracers_data.feedback_cumulative.max_kick_velocity_winds,
      /*can convert to comoving=*/0,
      "Same convention as MaxKickVelocityFromSupernovae, for stellar "
      "winds.");

  list[9] = io_make_physical_output_field(
      "MaxKickVelocityFromRadiationPressure", FLOAT, 1, UNIT_CONV_SPEED, 0.f,
      xparts, tracers_data.feedback_cumulative.max_kick_velocity_radiation,
      /*can convert to comoving=*/0,
      "Same convention as MaxKickVelocityFromSupernovae, for radiation "
      "pressure.");

  /* Same feedback-model dispatch reasoning as IsIonizedFlags above: must
     compile under any --with-feedback choice paired with
     --with-tracers=GEAR. */
  list[10] = io_make_output_field_convert_part(
      "PESpecificEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f,
      parts, xparts, convert_part_u_PE,
      "Local specific PE-band (6-11.2 eV) interstellar radiation field. "
      "Physical, mass-specific: no scale-factor exponent of its own.");

  list[11] = io_make_output_field_convert_part(
      "LWSpecificEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f,
      parts, xparts, convert_part_u_LW,
      "Local specific Lyman-Werner-band (11.2-13.6 eV) interstellar "
      "radiation field.");

  list[12] = io_make_output_field_convert_part(
      "PEArtificialDissipationCoefficients", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
      parts, xparts, convert_part_dissipation_alpha_PE,
      "Negativity-triggered artificial-dissipation coefficient of the "
      "PE-band hyperbolic propagation, in "
      "[0, max(ISRF_dissipation_alpha_max, ISRF_dissipation_alpha_floor)]. "
      "Only meaningful when ISRF_propagation is on.");

  list[13] = io_make_output_field_convert_part(
      "LWArtificialDissipationCoefficients", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
      parts, xparts, convert_part_dissipation_alpha_LW,
      "Same as PEArtificialDissipationCoefficients, Lyman-Werner band.");

  list[14] = io_make_output_field_convert_part(
      "PESpecificFluxDivergences", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS_PER_TIME, 0.f, parts, xparts,
      convert_part_div_specific_flux_PE,
      "`(1/rho) div(rho F)` accumulator of the PE-band hyperbolic "
      "propagation, accumulated in the force loop from the step's relaxed "
      "flux. Physical, like the "
      "specific energy it is a rate of change of. Only meaningful when "
      "ISRF_propagation is on.");

  list[15] = io_make_output_field_convert_part(
      "LWSpecificFluxDivergences", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS_PER_TIME, 0.f, parts, xparts,
      convert_part_div_specific_flux_LW,
      "Same as PESpecificFluxDivergences, Lyman-Werner band.");

  list[16] = io_make_output_field_convert_part(
      "PESpecificFluxes", FLOAT, 3, UNIT_CONV_ENERGY_PER_UNIT_MASS_VELOCITY,
      0.f, parts, xparts, convert_part_specific_flux_PE,
      "Tracked specific flux moment of the PE-band hyperbolic propagation, "
      "mass-specific like PESpecificEnergies. Physical: no scale-factor "
      "exponent of its own. Only meaningful when ISRF_propagation is on.");

  list[17] = io_make_output_field_convert_part(
      "LWSpecificFluxes", FLOAT, 3, UNIT_CONV_ENERGY_PER_UNIT_MASS_VELOCITY,
      0.f, parts, xparts, convert_part_specific_flux_LW,
      "Same as PESpecificFluxes, Lyman-Werner band.");

  list[18] = io_make_output_field_convert_part(
      "PEMinimumSpecificEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS,
      0.f, parts, xparts, convert_part_u_min_since_snapshot_PE,
      "Most negative PESpecificEnergies value the propagation update wrote "
      "since the previous snapshot, 0 if none was negative. The number of "
      "nonzero entries is the count of particles that undershot. The "
      "interval is since the last increment of engine.snapshot_output_count, "
      "which also happens when a FOF seeding catalogue is dumped "
      "(FOF:dump_catalogue_when_seeding), not only at a real snapshot. "
      "Always 0 unless the code is configured with --enable-debugging-checks.");

  list[19] = io_make_output_field_convert_part(
      "LWMinimumSpecificEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS,
      0.f, parts, xparts, convert_part_u_min_since_snapshot_LW,
      "Same as PEMinimumSpecificEnergies, Lyman-Werner band.");

  list[20] = io_make_output_field_convert_part(
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

  list[21] = io_make_output_field_convert_part(
      "LWCumulativeInjectedSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_cumulative_injected_LW,
      "Same as PECumulativeInjectedSpecificEnergies, Lyman-Werner band.");

  list[22] = io_make_output_field_convert_part(
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

  list[23] = io_make_output_field_convert_part(
      "LWCumulativeAbsorbedSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_cumulative_absorbed_LW,
      "Same as PECumulativeAbsorbedSpecificEnergies, Lyman-Werner band.");

  list[24] = io_make_output_field_convert_part(
      "HyperbolicPropagationSpeeds", FLOAT, 1, UNIT_CONV_SPEED, 0.f, parts,
      xparts, convert_part_c_hyp,
      "Kernel-local hyperbolic propagation speed the band updates and the "
      "pairwise transport operators ran with, shared by both bands. "
      "Physical: built from the physical smoothing length and a physical "
      "timestep, so no scale-factor exponent of its own. The conserved "
      "ledger of the consistent-variable-c schemes is `sum m u / c_hyp` "
      "rather than `sum m u`, which is what this field makes measurable "
      "from a snapshot. Only meaningful when ISRF_propagation is on.");

  list[25] = io_make_output_field_convert_part(
      "LWPhotonSpecificEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f,
      parts, xparts, convert_part_u_LW_PHOTON,
      "Lyman-Werner-band photon-number moment, energy-equivalent at a fixed "
      "reference photon energy: see LWSpecificEnergies. Distinct from the "
      "RT schemes' PhotonEnergies, which are raw per-group energies, not "
      "mass-specific and not band-prefixed.");

  list[26] = io_make_output_field_convert_part(
      "LWPhotonArtificialDissipationCoefficients", FLOAT, 1, UNIT_CONV_NO_UNITS,
      0.f, parts, xparts, convert_part_dissipation_alpha_LW_PHOTON,
      "Same as LWArtificialDissipationCoefficients: the photon-number "
      "moment shares the Lyman-Werner operator, so these two fields always "
      "read numerically identical.");

  list[27] = io_make_output_field_convert_part(
      "LWPhotonSpecificFluxDivergences", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS_PER_TIME, 0.f, parts, xparts,
      convert_part_div_specific_flux_LW_PHOTON,
      "Same as LWSpecificFluxDivergences, Lyman-Werner-band photon-number "
      "moment.");

  list[28] = io_make_output_field_convert_part(
      "LWPhotonSpecificFluxes", FLOAT, 3,
      UNIT_CONV_ENERGY_PER_UNIT_MASS_VELOCITY, 0.f, parts, xparts,
      convert_part_specific_flux_LW_PHOTON,
      "Same as LWSpecificFluxes, Lyman-Werner-band photon-number moment. "
      "Distinct from the RT schemes' PhotonFluxes, which are raw per-group "
      "fluxes, not mass-specific and not band-prefixed.");

  list[29] = io_make_output_field_convert_part(
      "LWPhotonMinimumSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_u_min_since_snapshot_LW_PHOTON,
      "Same as LWMinimumSpecificEnergies, Lyman-Werner-band photon-number "
      "moment. Always 0 unless the code is configured with "
      "--enable-debugging-checks.");

  list[30] = io_make_output_field_convert_part(
      "LWPhotonCumulativeInjectedSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_cumulative_injected_LW_PHOTON,
      "Same as LWCumulativeInjectedSpecificEnergies, Lyman-Werner-band "
      "photon-number moment. Always 0 unless the code is configured with "
      "--enable-debugging-checks.");

  list[31] = io_make_output_field_convert_part(
      "LWPhotonCumulativeAbsorbedSpecificEnergies", FLOAT, 1,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_cumulative_absorbed_LW_PHOTON,
      "Same as LWCumulativeAbsorbedSpecificEnergies, Lyman-Werner-band "
      "photon-number moment. Always 0 unless the code is configured with "
      "--enable-debugging-checks.");

  return num;
}

__attribute__((always_inline)) INLINE static int tracers_write_sparticles(
    const struct spart *sparts, struct io_props *list,
    const int with_cosmology) {

  int num = 11;

  list[0] = io_make_output_field(
      "FinalHIIRegionRadii", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, sparts,
      tracers_data.final_HII_radius,
      "Co-moving HII region radius of the star particles before they die or "
      "were not eligible to form HII regions anymore. Same algorithm's "
      "bookkeeping caveat as the live HIIRegionRadii it is retired from.");

  list[1] = io_make_output_field(
      "FinalHIIRegionMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      tracers_data.final_HII_mass,
      "Ionized gas mass of the star particles' HII region before they die or "
      "were not eligible to form HII regions anymore. Same algorithm's "
      "bookkeeping caveat as the live HIIRegionMasses it is retired from.");

  list[2] = io_make_output_field(
      "NumberOfSNIIEvents", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      tracers_data.snii_events.n_events,
      "Number of SNII events this star produced over its lifetime so far "
      "(fractional for a continuously-sampled population particle; always "
      "0 or 1 for a discrete star).");

  list[3] = io_make_physical_output_field(
      "DensityAtLastSNIIEvent", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts,
      tracers_data.snii_events.density_at_last_event,
      /*can convert to comoving=*/0,
      "Gas density at the star's location at its most recent SNII event. "
      "0 if it has never had one.");

  if (with_cosmology) {
    list[4] = io_make_physical_output_field(
        "ScaleFactorAtLastSNIIEvent", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
        tracers_data.snii_events.last_event_scale_factor,
        /*can convert to comoving=*/0,
        "Scale-factor at this star's most recent SNII event. 0 if it has "
        "never had one.");
  } else {
    list[4] = io_make_output_field(
        "TimeAtLastSNIIEvent", FLOAT, 1, UNIT_CONV_TIME, 0.f, sparts,
        tracers_data.snii_events.last_event_time,
        "Simulation time at this star's most recent SNII event. 0 if it "
        "has never had one.");
  }

  list[5] = io_make_output_field(
      "NumberOfSNIaEvents", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      tracers_data.snia_events.n_events,
      "Number of SNIa events this star produced over its lifetime so far. "
      "Always 0 for a discrete (single_star) particle: SNIa is a "
      "population-level channel in this model.");

  list[6] = io_make_physical_output_field(
      "DensityAtLastSNIaEvent", FLOAT, 1, UNIT_CONV_DENSITY, 0.f, sparts,
      tracers_data.snia_events.density_at_last_event,
      /*can convert to comoving=*/0,
      "Same as DensityAtLastSNIIEvent, for the SNIa channel.");

  if (with_cosmology) {
    list[7] = io_make_physical_output_field(
        "ScaleFactorAtLastSNIaEvent", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
        tracers_data.snia_events.last_event_scale_factor,
        /*can convert to comoving=*/0,
        "Same as ScaleFactorAtLastSNIIEvent, for the SNIa channel.");
  } else {
    list[7] = io_make_output_field(
        "TimeAtLastSNIaEvent", FLOAT, 1, UNIT_CONV_TIME, 0.f, sparts,
        tracers_data.snia_events.last_event_time,
        "Same as TimeAtLastSNIIEvent, for the SNIa channel.");
  }

  list[8] = io_make_output_field(
      "CumulativeWindsMassEjected", DOUBLE, 1, UNIT_CONV_MASS, 0.f, sparts,
      tracers_data.winds.mass_ejected,
      "Mass this star ejected by stellar winds over its lifetime so far, "
      "counted once per injection step from the star's own budget, but only "
      "for a step whose preceding step left the star with gas neighbours; "
      "the star's mass budget is spent regardless, so this undercounts the "
      "star's true total wind mass loss whenever a step had none.");

  list[9] = io_make_physical_output_field(
      "CumulativeWindsEnergyEjected", DOUBLE, 1, UNIT_CONV_ENERGY, 0.f, sparts,
      tracers_data.winds.energy_ejected,
      /*can convert to comoving=*/0,
      "Total (not specific) stellar-wind energy this star ejected over its "
      "lifetime so far, after the winds efficiency factor.");

  list[10] = io_make_physical_output_field(
      "CumulativeWindsMomentumEjected", DOUBLE, 1, UNIT_CONV_MOMENTUM, 0.f,
      sparts, tracers_data.winds.momentum_ejected,
      /*can convert to comoving=*/0,
      "Sum over injection steps of the wind momentum budget sqrt(2 m_ej E_ej), "
      "in the star's rest frame (scalar sum, not vector: isotropic ejecta "
      "would else cancel). Excludes the m_ej v_star term that the gas "
      "CumulativeMomentumFromWinds includes, which is why the received/ejected "
      "ratio is not 1 for a moving star.");

  return num;
}

__attribute__((always_inline)) INLINE static int tracers_write_bparticles(
    const struct bpart *bparts, struct io_props *list,
    const int with_cosmology) {

  list[0] = io_make_output_field_convert_bpart(
      "AveragedAccretionRates", FLOAT, num_snapshot_triggers_bpart,
      UNIT_CONV_MASS_PER_UNIT_TIME, 0.f, bparts,
      convert_bpart_averaged_accretion_rate,
      "Accretion rates of the black holes averaged over the period set by "
      "the first N snapshot triggers");

  return 1;
}

__attribute__((always_inline)) INLINE static int tracers_write_sinkparticles(
    const struct sink *sinks, struct io_props *list, const int with_cosmology) {

  list[0] = io_make_output_field_convert_sink(
      "AveragedAccretionRates", FLOAT, num_snapshot_triggers_sink,
      UNIT_CONV_MASS_PER_UNIT_TIME, 0.f, sinks,
      convert_sink_averaged_accretion_rate,
      "Accretion rates of the sinks averaged over the period set by the "
      "first N snapshot triggers");

  list[1] = io_make_output_field_convert_sink(
      "AveragedStarFormationRates", FLOAT, num_snapshot_triggers_sink,
      UNIT_CONV_SFR, 0.f, sinks, convert_sink_averaged_SFR,
      "Star formation rates of the particles averaged over the period set by "
      "the first N snapshot triggers");

  return 2;
}

#endif /* SWIFT_TRACERS_NONE_IO_H */
