/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_GEAR_RADIATION_PART_H
#define SWIFT_FEEDBACK_GEAR_RADIATION_PART_H

/**
 * @file src/feedback/GEAR/radiation_part.h
 * @brief Accessors of the subgrid radiation fields of a gas particle and the
 * records of its HII state exchanged over MPI, shared by the GEAR feedback
 * modules. Include this file after the module's #feedback_part_data.
 */

#ifdef GEAR_SUBGRID_RADIATION_ISRF
/*! Bit of #feedback_part_data.ISRF_H2_line_absorption: the H2 lines absorb
    the transported LW light on its path. */
#define RADIATION_ISRF_H2_LINES 1

/*! Bit of #feedback_part_data.ISRF_H2_line_absorption: the stars give their
    LW light their own H2 dissociation factor. */
#define RADIATION_ISRF_LW_POPULATION 2

/**
 * @brief Do the H2 lines absorb the LW light received by this particle?
 *
 * @param fd The feedback data of the particle.
 */
__attribute__((always_inline)) static inline int feedback_part_has_H2_lines(
    const struct feedback_part_data *fd) {
  return (fd->ISRF_H2_line_absorption & RADIATION_ISRF_H2_LINES) != 0;
}

/**
 * @brief Does this particle carry the H2-dissociating part of its LW light?
 *
 * @param fd The feedback data of the particle.
 */
__attribute__((always_inline)) static inline int
feedback_part_tracks_LW_dissociating(const struct feedback_part_data *fd) {
  return fd->ISRF_H2_line_absorption != 0;
}
#endif

#if defined(WITH_MPI) && defined(GEAR_SUBGRID_RADIATION_HII)

/**
 * @brief HII claim sent back by the rank that ran the star's search to the
 * rank that owns the gas (task_subtype_part_hii_tag).
 *
 * One entry per particle of the cell. The owner merges only the entries
 * stamped by #claimed_this_pass.
 */
struct hii_tag_report {

  /*! Id of the star that claimed the particle. */
  long long star_id;

  /*! Simulation time until which the particle stays tagged. */
  double end_time;

  /*! Photons charged for this claim, counted as forfeited if it loses. */
  double cost;

  /*! Squared distance to the claiming star at claim time. The owner may hold
      no copy of that star, so it is sent, not recomputed. */
  float r2;

  /*! Ionization tag of the claim. */
  char is_ionized;

  /*! Set if this entry is a claim of this pass; others are not merged. */
  char claimed_this_pass;
};

/**
 * @brief Post-cooling HII state sent by the gas owner to every rank that
 * holds a foreign copy (task_subtype_part_hii_state).
 *
 * The search on the other rank then sees the owner's tag and eligibility
 * caches of this step, not the pre-cooling snapshot of the xv/rho channel.
 * The unpack writes only these fields.
 */
struct hii_state_update {

  /*! Id of the star holding the tag. */
  long long star_id;

  /*! Simulation time until which the particle stays tagged. */
  double end_time;

  /*! Photons debited for the current claim. */
  double cost;

  /*! Eligibility temperature cache, see #feedback_part_data. */
  float T_eligibility;

  /*! Mean molecular weight cache, see #feedback_part_data. */
  float mu_eligibility;

  /*! Neutral hydrogen mass fraction cache, see #feedback_part_data. */
  float neutral_H_frac;

  /*! Squared distance to the star holding the tag, at claim time. */
  float r2;

  /*! Ionization tag. */
  char is_ionized;
};

#else

/*! Empty payloads: without HII under MPI these subtypes carry nothing. */
struct hii_tag_report {};
struct hii_state_update {};

#endif

#endif /* SWIFT_FEEDBACK_GEAR_RADIATION_PART_H */
