/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_GEAR_STAR_FORMATION_STRUCT_H
#define SWIFT_GEAR_STAR_FORMATION_STRUCT_H

/**
 * @brief Star-formation-related properties stored in the extended particle
 * data.
 */
struct star_formation_xpart_data {
	/*!star formation rate (mass/(time*volume))*/
	double SFR; 
    /*! Estimation of local turbulence (squared) */
    float sigma2;
	};

/* Starformation struct */
struct star_formation {
	/*! Njeans number. We request that the SPH mass resolution is Njeans times smaller than the Jeans mass*/
	int Njeans;
	/*!star formation efficency. If the particle can create a star, it happens with this probablity*/
	double star_formation_rate;
	/*!do we include local turbulence*/
	int with_sigma;
	/*!Maximum temperature to allow star formation* (should be about 10'000 or 30'000 K*/
	double Max_temperature;
	/*!some other useful elements*/
	const struct hydro_props* restrict hydro_props;
	/*!units*/
	const struct unit_system* restrict us;
	/*! physical constants*/
	const struct phys_const* phys_const;
	};

#endif /* SWIFT_GEAR_STAR_FORMATION_STRUCT_H */
