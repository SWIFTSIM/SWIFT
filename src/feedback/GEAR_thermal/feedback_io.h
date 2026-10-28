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

#ifdef GEAR_SUBGRID_RADIATION_HII
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

#endif /* GEAR_SUBGRID_RADIATION_HII */

#ifdef GEAR_SUBGRID_RADIATION_ISRF
/**
 * @brief Energy waiting in the arrival bins of one band of a slot.
 *
 * @param t The band of the slot.
 */
INLINE static double packet_bins_energy(
    const struct feedback_isrf_packet_slot_band_data *t) {
  double sum = 0.;
  for (int i = 0; i < RADIATION_PACKET_BINS; i++) sum += t->bin_energy[i];
  return sum;
}

/**
 * @brief Photon number waiting in the arrival bins of one band of a slot.
 *
 * @param t The band of the slot.
 */
INLINE static double packet_bins_photons(
    const struct feedback_isrf_packet_slot_band_data *t) {
  double sum = 0.;
  for (int i = 0; i < RADIATION_PACKET_BINS; i++) sum += t->bin_photons[i];
  return sum;
}

/**
 * @brief Snapshot converter of the packet inbox energy per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_inbox_energy(const struct engine *e,
                                                    const struct part *p,
                                                    const struct xpart *xp,
                                                    float *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    double sum = 0.;
    for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
      sum += packet_bins_energy(&p->feedback_data.isrf_packet.slot[s].band[b]);
    ret[b] = (float)sum;
  }
}

/**
 * @brief Snapshot converter of the packet inbox photons per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_inbox_photons(const struct engine *e,
                                                     const struct part *p,
                                                     const struct xpart *xp,
                                                     float *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    double sum = 0.;
    for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
      sum += packet_bins_photons(&p->feedback_data.isrf_packet.slot[s].band[b]);
    ret[b] = (float)sum;
  }
}

/**
 * @brief Snapshot converter of the packet reservoir energy per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_reservoir_energy(const struct engine *e,
                                                        const struct part *p,
                                                        const struct xpart *xp,
                                                        float *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    double sum = 0.;
    for (int s = 0; s < RADIATION_PACKET_SLOTS; s++) {
      const struct feedback_isrf_packet_slot_band_data *t =
          &p->feedback_data.isrf_packet.slot[s].band[b];
      if (e->feedback_props->ISRF_packet_light_speed !=
          isrf_packet_light_speed_unlimited) {
        for (int c = 0; c < RADIATION_PACKET_FAR_CELLS; c++)
          sum += t->far_energy[c];
      } else {
        sum += t->reservoir_energy;
      }
    }
    ret[b] = (float)sum;
  }
}

/**
 * @brief Snapshot converter of the packet reservoir photons per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_reservoir_photons(const struct engine *e,
                                                         const struct part *p,
                                                         const struct xpart *xp,
                                                         float *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    double sum = 0.;
    for (int s = 0; s < RADIATION_PACKET_SLOTS; s++) {
      const struct feedback_isrf_packet_slot_band_data *t =
          &p->feedback_data.isrf_packet.slot[s].band[b];
      if (e->feedback_props->ISRF_packet_light_speed !=
          isrf_packet_light_speed_unlimited) {
        for (int c = 0; c < RADIATION_PACKET_FAR_CELLS; c++)
          sum += t->far_photons[c];
      } else {
        sum += t->reservoir_photons;
      }
    }
    ret[b] = (float)sum;
  }
}

/**
 * @brief Snapshot converter of the packet direction of each reception slot,
 * see #feedback_write_particles.
 */
INLINE static void convert_part_packet_slot_directions(const struct engine *e,
                                                       const struct part *p,
                                                       const struct xpart *xp,
                                                       float *ret) {
  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
    for (int k = 0; k < 3; k++)
      ret[3 * s + k] = p->feedback_data.isrf_packet.slot[s].direction[k];
}

/**
 * @brief Snapshot converter of the origin distance of each reception slot,
 * see #feedback_write_particles.
 */
INLINE static void convert_part_packet_slot_origin_distances(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
    ret[s] = p->feedback_data.isrf_packet.slot[s].origin_distance;
}

/**
 * @brief Snapshot converter of the inbox energy of each reception slot and
 * band, see #feedback_write_particles.
 */
INLINE static void convert_part_packet_slot_inbox_energy(const struct engine *e,
                                                         const struct part *p,
                                                         const struct xpart *xp,
                                                         float *ret) {
  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
      ret[ISRF_OPERATOR_COUNT * s + b] =
          packet_bins_energy(&p->feedback_data.isrf_packet.slot[s].band[b]);
}

/**
 * @brief Snapshot converter of the direction of the inbox light of each
 * reception slot, see #feedback_write_particles.
 */
INLINE static void convert_part_packet_slot_inbox_directions(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++) {
    const struct feedback_isrf_packet_slot_data *sl =
        &p->feedback_data.isrf_packet.slot[s];
    double energy = 0., norm = 0.;
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
      energy += packet_bins_energy(&sl->band[b]);
    for (int k = 0; k < 3; k++)
      norm += sl->inbox_origin_offset[k] * sl->inbox_origin_offset[k];
    norm = sqrt(norm);
    for (int k = 0; k < 3; k++)
      ret[3 * s + k] = (energy > 0. && norm > 0.)
                           ? (float)(-sl->inbox_origin_offset[k] / norm)
                           : 0.f;
  }
}

/**
 * @brief Snapshot converter of the packet energy per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_energy(const struct engine *e,
                                              const struct part *p,
                                              const struct xpart *xp,
                                              float *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    double sum = 0.;
    for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
      sum += p->feedback_data.isrf_packet.slot[s].band[b].energy;
    ret[b] = (float)sum;
  }
}

/**
 * @brief Snapshot converter of the packet photons per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_photons(const struct engine *e,
                                               const struct part *p,
                                               const struct xpart *xp,
                                               float *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    double sum = 0.;
    for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
      sum += p->feedback_data.isrf_packet.slot[s].band[b].photons;
    ret[b] = (float)sum;
  }
}

/**
 * @brief Snapshot converter of the packet absorbed energy per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_absorbed_energy(const struct engine *e,
                                                       const struct part *p,
                                                       const struct xpart *xp,
                                                       double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = p->feedback_data.isrf_packet.band[b].absorbed_energy;
}

/**
 * @brief Snapshot converter of the packet absorbed photons per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_absorbed_photons(const struct engine *e,
                                                        const struct part *p,
                                                        const struct xpart *xp,
                                                        double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = p->feedback_data.isrf_packet.band[b].absorbed_photons;
}

/**
 * @brief Snapshot converter of the packet escaped energy per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_escaped_energy(const struct engine *e,
                                                      const struct part *p,
                                                      const struct xpart *xp,
                                                      double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = p->feedback_data.isrf_packet.band[b].escaped_energy;
}

/**
 * @brief Snapshot converter of the packet escaped photons per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_escaped_photons(const struct engine *e,
                                                       const struct part *p,
                                                       const struct xpart *xp,
                                                       double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = p->feedback_data.isrf_packet.band[b].escaped_photons;
}

/**
 * @brief Snapshot converter of the packet redshifted energy per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_redshifted_energy(const struct engine *e,
                                                         const struct part *p,
                                                         const struct xpart *xp,
                                                         double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = p->feedback_data.isrf_packet.band[b].redshifted_energy;
}

/**
 * @brief Snapshot converter of the packet redshifted photons per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_redshifted_photons(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = p->feedback_data.isrf_packet.band[b].redshifted_photons;
}

/**
 * @brief Snapshot converter of the packet field energy per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_field_energy(const struct engine *e,
                                                    const struct part *p,
                                                    const struct xpart *xp,
                                                    float *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = p->feedback_data.isrf_packet.band[b].field_energy;
}

/**
 * @brief Snapshot converter of the packet field photons per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_field_photons(const struct engine *e,
                                                     const struct part *p,
                                                     const struct xpart *xp,
                                                     float *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = p->feedback_data.isrf_packet.band[b].field_photons;
}

/**
 * @brief Snapshot converter of the packet ngb mean field energy per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_ngb_mean_field_energy(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = p->feedback_data.isrf_packet.band[b].ngb_mean_field_energy;
}

/**
 * @brief Snapshot converter of the packet ngb mean field photons per band, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_packet_ngb_mean_field_photons(
    const struct engine *e, const struct part *p, const struct xpart *xp,
    float *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = p->feedback_data.isrf_packet.band[b].ngb_mean_field_photons;
}

/**
 * @brief Snapshot converter of a star's emitted packet energy per band, see
 * #feedback_write_sparticles.
 */
INLINE static void convert_spart_packet_emitted_energy(const struct engine *e,
                                                       const struct spart *sp,
                                                       double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = sp->feedback_data.radiation.packet.emitted_energy[b];
}

/**
 * @brief Snapshot converter of a star's emitted packet photon number per
 * band, see #feedback_write_sparticles.
 */
INLINE static void convert_spart_packet_emitted_photons(const struct engine *e,
                                                        const struct spart *sp,
                                                        double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = sp->feedback_data.radiation.packet.emitted_photons[b];
}

/**
 * @brief Snapshot converter of a star's escaped packet energy per band, see
 * #feedback_write_sparticles.
 */
INLINE static void convert_spart_packet_escaped_energy(const struct engine *e,
                                                       const struct spart *sp,
                                                       double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = sp->feedback_data.radiation.packet.escaped_energy[b];
}

/**
 * @brief Snapshot converter of a star's escaped packet photon number per
 * band, see #feedback_write_sparticles.
 */
INLINE static void convert_spart_packet_escaped_photons(const struct engine *e,
                                                        const struct spart *sp,
                                                        double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = sp->feedback_data.radiation.packet.escaped_photons[b];
}

/**
 * @brief Snapshot converter of the energy in flight in the removed gas per
 * band, see #feedback_write_sparticles.
 */
INLINE static void convert_spart_packet_removed_energy(const struct engine *e,
                                                       const struct spart *sp,
                                                       double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = sp->feedback_data.radiation.packet.removed.removed_energy[b];
}

/**
 * @brief Snapshot converter of the photon number in flight in the removed gas
 * per band, see #feedback_write_sparticles.
 */
INLINE static void convert_spart_packet_removed_photons(const struct engine *e,
                                                        const struct spart *sp,
                                                        double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = sp->feedback_data.radiation.packet.removed.removed_photons[b];
}

/**
 * @brief Snapshot converter of the absorbed and escaped energy of the removed
 * gas per band, see #feedback_write_sparticles.
 */
INLINE static void convert_spart_packet_inherited_energy(const struct engine *e,
                                                         const struct spart *sp,
                                                         double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = sp->feedback_data.radiation.packet.removed.inherited_energy[b];
}

/**
 * @brief Snapshot converter of the absorbed and escaped photon number of the
 * removed gas per band, see #feedback_write_sparticles.
 */
INLINE static void convert_spart_packet_inherited_photons(
    const struct engine *e, const struct spart *sp, double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = sp->feedback_data.radiation.packet.removed.inherited_photons[b];
}

/**
 * @brief Snapshot converter of the energy in flight in the removed gas per
 * band, see #feedback_write_sinkparticles.
 */
INLINE static void convert_sink_packet_removed_energy(const struct engine *e,
                                                      const struct sink *s,
                                                      double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = s->feedback_data.isrf_packet.removed_energy[b];
}

/**
 * @brief Snapshot converter of the photon number in flight in the removed gas
 * per band, see #feedback_write_sinkparticles.
 */
INLINE static void convert_sink_packet_removed_photons(const struct engine *e,
                                                       const struct sink *s,
                                                       double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = s->feedback_data.isrf_packet.removed_photons[b];
}

/**
 * @brief Snapshot converter of the absorbed and escaped energy of the removed
 * gas per band, see #feedback_write_sinkparticles.
 */
INLINE static void convert_sink_packet_inherited_energy(const struct engine *e,
                                                        const struct sink *s,
                                                        double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = s->feedback_data.isrf_packet.inherited_energy[b];
}

/**
 * @brief Snapshot converter of the absorbed and escaped photon number of the
 * removed gas per band, see #feedback_write_sinkparticles.
 */
INLINE static void convert_sink_packet_inherited_photons(const struct engine *e,
                                                         const struct sink *s,
                                                         double *ret) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    ret[b] = s->feedback_data.isrf_packet.inherited_photons[b];
}

/**
 * @brief Snapshot converter for #PESpecificEnergy, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_u_PE(const struct engine *e,
                                     const struct part *p,
                                     const struct xpart *xp, double *ret) {
  ret[0] = feedback_get_part_u_PE(p);
}

/**
 * @brief Snapshot converter for #LWSpecificEnergy, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_u_LW(const struct engine *e,
                                     const struct part *p,
                                     const struct xpart *xp, double *ret) {
  ret[0] = feedback_get_part_u_LW(p);
}

/**
 * @brief Snapshot converter for #LWPhotonSpecificEnergies, see
 * #feedback_write_particles.
 */
INLINE static void convert_part_u_LW_PHOTON(const struct engine *e,
                                            const struct part *p,
                                            const struct xpart *xp,
                                            double *ret) {
  ret[0] = feedback_get_part_u_LW_PHOTON(p);
}

#endif /* GEAR_SUBGRID_RADIATION_ISRF */

/**
 * @brief Specifies which particle fields to read from the ICs.
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int feedback_read_particles(struct part *parts,
                                          struct io_props *list) {

#ifdef GEAR_SUBGRID_RADIATION_ISRF
  /* Physical, mass-specific quantities: an IC value is taken verbatim. */

  list[0] = io_make_input_field("PESpecificEnergy", DOUBLE, 1, OPTIONAL,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                                feedback_data.isrf_moment[ISRF_MOMENT_PE].u);
  list[1] = io_make_input_field("LWSpecificEnergy", DOUBLE, 1, OPTIONAL,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                                feedback_data.isrf_moment[ISRF_MOMENT_LW].u);

  /* A zero value beside a nonzero LWSpecificEnergy is overwritten at first
     init, see radiation_first_init_part(). */
  list[2] =
      io_make_input_field("LWPhotonSpecificEnergy", DOUBLE, 1, OPTIONAL,
                          UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                          feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].u);

  return 3;
#else
  return 0;
#endif
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

  int num = 0;

#ifdef GEAR_SUBGRID_RADIATION_HII
  list[num++] = io_make_output_field_convert_part(
      "IsIonizedFlags", CHAR, 1, UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
      convert_part_is_ionized,
      "Were the particles flagged as ionized by HII ionzation subgrid model?");

  list[num++] = io_make_output_field_convert_part(
      "HIIStarIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
      convert_part_HII_star_id,
      "Star particle IDs that ionized these gas particles due to HII ionzation "
      "subgrid model?");
#endif

#ifdef GEAR_SUBGRID_RADIATION_ISRF
  list[num++] = io_make_output_field_convert_part(
      "PESpecificEnergies", DOUBLE, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f,
      parts, xparts, convert_part_u_PE,
      "Local specific PE-band (6-11.2 eV) interstellar radiation field. "
      "Physical, mass-specific: no scale-factor exponent of its own. "
      "DOUBLE.");

  list[num++] = io_make_output_field_convert_part(
      "LWSpecificEnergies", DOUBLE, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f,
      parts, xparts, convert_part_u_LW,
      "Local specific Lyman-Werner-band (11.2-13.6 eV) interstellar "
      "radiation field. DOUBLE: see PESpecificEnergies.");

  list[num++] = io_make_output_field_convert_part(
      "LWPhotonSpecificEnergies", DOUBLE, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS,
      0.f, parts, xparts, convert_part_u_LW_PHOTON,
      "Lyman-Werner-band photon-number moment, energy-equivalent at a fixed "
      "reference photon energy: see LWSpecificEnergies. Distinct from the "
      "RT schemes' PhotonEnergies, which are raw per-group energies, not "
      "mass-specific and not band-prefixed. DOUBLE: see LWSpecificEnergies.");

  list[num++] = io_make_output_field_convert_part(
      "PacketInboxEnergies", FLOAT, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY, 0.f,
      parts, xparts, convert_part_packet_inbox_energy,
      "Packet energy received since the particle's last step, sent at its "
      "next step. Two components, PE (6-11.2 eV) then LW (11.2-13.6 eV). "
      "Zero unless GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketInboxPhotonNumbers", FLOAT, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, parts, xparts, convert_part_packet_inbox_photons,
      "Photon number of the inbox times the band's lower-edge photon "
      "energy, so in energy units: divide by 6.0 eV (PE) or 11.2 eV (LW) "
      "for a count. Two components, PE (6-11.2 eV) then LW (11.2-13.6 "
      "eV). Zero unless GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketReservoirEnergies", FLOAT, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, parts, xparts, convert_part_packet_reservoir_energy,
      "Packet energy that waits beyond the arrival bins of the inbox: with "
      "a speed limit, light that arrives more than three steps ahead; "
      "without one, energy from senders on longer steps, released over "
      "their step. Two components, PE (6-11.2 eV) then LW (11.2-13.6 eV). "
      "Zero unless GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketReservoirPhotonNumbers", FLOAT, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY, 0.f, parts, xparts,
      convert_part_packet_reservoir_photons,
      "Photon number of the reservoir, in the units of "
      "PacketInboxPhotonNumbers. Two components, PE (6-11.2 eV) then LW "
      "(11.2-13.6 eV). Zero unless GEARFeedback:ISRF_propagation is "
      "on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketHeldEnergies", FLOAT, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY, 0.f,
      parts, xparts, convert_part_packet_energy,
      "Packet energy held for sending this step; zero between steps. Two "
      "components, PE (6-11.2 eV) then LW (11.2-13.6 eV). Zero unless "
      "GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketHeldPhotonNumbers", FLOAT, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, parts, xparts, convert_part_packet_photons,
      "Photon number held, in the units of PacketInboxPhotonNumbers. Two "
      "components, PE (6-11.2 eV) then LW (11.2-13.6 eV). Zero unless "
      "GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketAbsorbedEnergies", DOUBLE, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, parts, xparts, convert_part_packet_absorbed_energy,
      "Packet energy absorbed by dust in the particle since the start of "
      "the run. Two components, PE (6-11.2 eV) then LW (11.2-13.6 eV). "
      "Zero unless GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketAbsorbedPhotonNumbers", DOUBLE, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY, 0.f, parts, xparts,
      convert_part_packet_absorbed_photons,
      "Photon number absorbed, in the units of PacketInboxPhotonNumbers. "
      "Two components, PE (6-11.2 eV) then LW (11.2-13.6 eV). Zero unless "
      "GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketEscapedEnergies", DOUBLE, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, parts, xparts, convert_part_packet_escaped_energy,
      "Packet energy that left the particle with no downstream neighbour, "
      "since the start of the run. Two components, PE (6-11.2 eV) then LW "
      "(11.2-13.6 eV). Zero unless GEARFeedback:ISRF_propagation is "
      "on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketEscapedPhotonNumbers", DOUBLE, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY, 0.f, parts, xparts, convert_part_packet_escaped_photons,
      "Photon number escaped, in the units of PacketInboxPhotonNumbers. "
      "Two components, PE (6-11.2 eV) then LW (11.2-13.6 eV). Zero unless "
      "GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketRedshiftedEnergies", DOUBLE, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, parts, xparts, convert_part_packet_redshifted_energy,
      "Packet energy removed from the band by the expansion in the shares "
      "booked into the particle since the start of the run (redshift and "
      "band-edge loss), minus the energy gained from the other band. Two "
      "components, PE (6-11.2 eV) then LW (11.2-13.6 eV). Zero in a static "
      "box and unless GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketRedshiftedPhotonNumbers", DOUBLE, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY, 0.f, parts, xparts,
      convert_part_packet_redshifted_photons,
      "Photon number removed from the band by the expansion, in the units of "
      "PacketInboxPhotonNumbers, net of the gain from the other band. Two "
      "components, PE (6-11.2 eV) then LW (11.2-13.6 eV). Zero in a static "
      "box and unless GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketSpecificEnergies", FLOAT, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_packet_field_energy,
      "Own field of the packet transport, (E / m) l_par / (c dt), the specific "
      "energy of the light that crosses the particle in one step. Two "
      "components, PE then LW. Zero unless GEARFeedback:ISRF_propagation is "
      "on. Its running mean over "
      "GEARFeedback:ISRF_packet_read_mean_time_Myr when that is above 0.");

  list[num++] = io_make_output_field_convert_part(
      "PacketSpecificPhotonNumbers", FLOAT, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_packet_field_photons,
      "Own field of the photon number, per unit mass, times the band's "
      "lower-edge photon energy (6.0 eV PE, 11.2 eV LW).");

  list[num++] = io_make_output_field_convert_part(
      "PacketKernelMeanSpecificEnergies", FLOAT, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_packet_ngb_mean_field_energy,
      "Mass- and kernel-weighted mean of PacketSpecificEnergies over the "
      "particle's support, itself included, at its last step.");

  list[num++] = io_make_output_field_convert_part(
      "PacketKernelMeanSpecificPhotonNumbers", FLOAT, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY_PER_UNIT_MASS, 0.f, parts, xparts,
      convert_part_packet_ngb_mean_field_photons,
      "Kernel mean of PacketSpecificPhotonNumbers, as "
      "PacketKernelMeanSpecificEnergies.");

  list[num++] = io_make_output_field_convert_part(
      "PacketSlotDirections", FLOAT, 3 * RADIATION_PACKET_SLOTS,
      UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
      convert_part_packet_slot_directions,
      "Unit direction of the packet of each of the " RADIATION_PACKET_SLOT_TEXT
      " Zero for a slot never merged, and unless GEARFeedback:ISRF_propagation "
      "is "
      "on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketSlotOriginDistances", FLOAT, RADIATION_PACKET_SLOTS,
      UNIT_CONV_LENGTH, 1.f, parts, xparts,
      convert_part_packet_slot_origin_distances,
      "Co-moving distance from the packet's origin to the particle, per "
      "reception slot. Zero unless GEARFeedback:ISRF_propagation is "
      "on.");

  list[num++] = io_make_output_field_convert_part(
      "PacketSlotInboxEnergies", FLOAT,
      RADIATION_PACKET_SLOTS * ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY, 0.f,
      parts, xparts, convert_part_packet_slot_inbox_energy,
      "Packet energy in the inbox of each reception slot: slot 0 first, "
      "then PE (6-11.2 eV) and LW (11.2-13.6 eV). Sums over the slots to "
      "PacketInboxEnergies.");

  list[num++] = io_make_output_field_convert_part(
      "PacketSlotInboxDirections", FLOAT, 3 * RADIATION_PACKET_SLOTS,
      UNIT_CONV_NO_UNITS, 0.f, parts, xparts,
      convert_part_packet_slot_inbox_directions,
      "Unit direction from the energy-weighted mean origin of the inbox of "
      "each reception slot to the particle: slot 0 first, then x, y, z. Zero "
      "for an empty inbox.");
#endif

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

  int num = 0;

#ifdef GEAR_SUBGRID_RADIATION_HII
  list[num++] = io_make_output_field(
      "FinalHIIRegionRadii", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, sparts,
      feedback_data.radiation.final_HII_radius,
      "Co-moving HII region radius of the star particles before they die or "
      "were not eligible to form HII regions anymore. Same algorithm's "
      "bookkeeping caveat as the live HIIRegionRadii it is retired from.");

  list[num++] = io_make_output_field(
      "FinalHIIRegionMasses", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      feedback_data.radiation.final_HII_mass,
      "Ionized gas mass of the star particles' HII region before they die or "
      "were not eligible to form HII regions anymore. Same algorithm's "
      "bookkeeping caveat as the live HIIRegionMasses it is retired from.");
#endif

#ifdef GEAR_SUBGRID_RADIATION_ISRF
  list[num++] = io_make_output_field_convert_spart(
      "PacketEmittedEnergies", DOUBLE, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, sparts, convert_spart_packet_emitted_energy,
      "Packet energy emitted by the star since its birth. Two components, PE "
      "then LW. Zero unless GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_spart(
      "PacketEmittedPhotonNumbers", DOUBLE, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY, 0.f, sparts, convert_spart_packet_emitted_photons,
      "Photon number emitted by the star, times the band's lower-edge photon "
      "energy (6.0 eV PE, 11.2 eV LW).");

  list[num++] = io_make_output_field_convert_spart(
      "PacketEscapedEnergies", DOUBLE, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, sparts, convert_spart_packet_escaped_energy,
      "Packet energy of the star's empty emission cones since its birth. "
      "Two components, PE then LW. Zero unless GEARFeedback:ISRF_propagation "
      "is on.");

  list[num++] = io_make_output_field_convert_spart(
      "PacketEscapedPhotonNumbers", DOUBLE, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY, 0.f, sparts, convert_spart_packet_escaped_photons,
      "Photon number of the star's empty emission cones, times the band's "
      "lower-edge photon energy (6.0 eV PE, 11.2 eV LW).");

  list[num++] = io_make_output_field_convert_spart(
      "PacketRemovedEnergies", DOUBLE, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, sparts, convert_spart_packet_removed_energy,
      "Packet energy in flight in the gas particle the star was converted "
      "from, lost with it. Two components, PE then LW.");

  list[num++] = io_make_output_field_convert_spart(
      "PacketRemovedPhotonNumbers", DOUBLE, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY, 0.f, sparts, convert_spart_packet_removed_photons,
      "Photon number in flight in the gas particle the star was converted "
      "from, in the units of PacketEmittedPhotonNumbers.");

  list[num++] = io_make_output_field_convert_spart(
      "PacketInheritedEnergies", DOUBLE, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, sparts, convert_spart_packet_inherited_energy,
      "Gas PacketAbsorbedEnergies plus PacketEscapedEnergies plus "
      "PacketRedshiftedEnergies of the gas particle the star was converted "
      "from. Two components, PE then LW.");

  list[num++] = io_make_output_field_convert_spart(
      "PacketInheritedPhotonNumbers", DOUBLE, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY, 0.f, sparts, convert_spart_packet_inherited_photons,
      "Absorbed, escaped and redshifted photon number of the gas particle the "
      "star was converted from, in the units of PacketEmittedPhotonNumbers.");
#endif

  return num;
}

/**
 * @brief Specifies which sink particle fields to write to a dataset.
 *
 * @param sinks The sink particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology switched on?
 *
 * @return Returns the number of fields to write.
 */
INLINE static int feedback_write_sinkparticles(const struct sink *sinks,
                                               struct io_props *list,
                                               const int with_cosmology) {

  int num = 0;

#ifdef GEAR_SUBGRID_RADIATION_ISRF
  list[num++] = io_make_output_field_convert_sink(
      "PacketRemovedEnergies", DOUBLE, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, sinks, convert_sink_packet_removed_energy,
      "Packet energy in flight in the gas the sink formed from or swallowed, "
      "and in the sinks it merged with, lost with it. Two components, PE then "
      "LW. Zero unless GEARFeedback:ISRF_propagation is on.");

  list[num++] = io_make_output_field_convert_sink(
      "PacketRemovedPhotonNumbers", DOUBLE, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY, 0.f, sinks, convert_sink_packet_removed_photons,
      "Photon number in flight in the gas the sink formed from or swallowed, "
      "times the band's lower-edge photon energy (6.0 eV PE, 11.2 eV LW).");

  list[num++] = io_make_output_field_convert_sink(
      "PacketInheritedEnergies", DOUBLE, ISRF_OPERATOR_COUNT, UNIT_CONV_ENERGY,
      0.f, sinks, convert_sink_packet_inherited_energy,
      "Gas PacketAbsorbedEnergies plus PacketEscapedEnergies plus "
      "PacketRedshiftedEnergies of the gas the sink formed from or swallowed. "
      "Two components, PE then LW.");

  list[num++] = io_make_output_field_convert_sink(
      "PacketInheritedPhotonNumbers", DOUBLE, ISRF_OPERATOR_COUNT,
      UNIT_CONV_ENERGY, 0.f, sinks, convert_sink_packet_inherited_photons,
      "Absorbed, escaped and redshifted photon number of the gas the sink "
      "formed from or swallowed, in the units of PacketRemovedPhotonNumbers.");
#endif

  return num;
}

/**
 * @brief Specifies which star particle fields to read from the ICs or a
 * restart file. Empty: no feedback-owned #spart field is read.
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
