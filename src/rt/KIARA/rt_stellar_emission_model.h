/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_STELLAR_EMISSION_MODEL_KIARA_H
#define SWIFT_RT_STELLAR_EMISSION_MODEL_KIARA_H

/**
 * @file src/rt/KIARA/rt_stellar_emission_model.h
 * @brief Main header file for the KIARA M1 closure radiative transfer scheme
 * stellar radiation emission models.
 */

enum rt_stellar_emission_models {
  rt_stellar_emission_model_none = 0,
  rt_stellar_emission_model_const,
  rt_stellar_emission_model_IlievTest,
  rt_stellar_emission_model_BPASS,
  rt_stellar_emission_model_count
};

/**
 * @brief Compute the energy emitted from a star during the time step dt.
 * This is for the constant emission rate model.
 *
 * @param emission_this_step (return) the emitted radiation energy of a star
 * during the time interval dt
 * @param const_stellar_emission_rates the constant emission rates used in this
 * run
 * @param dt time step size (in internal units)
 */
__attribute__((always_inline)) INLINE static void
rt_get_emission_this_step_const(
    double emission_this_step[RT_NGROUPS],
    const double const_stellar_emission_rates[RT_NGROUPS], double dt) {

  /* The read-in constant stellar emisison rates are in units of L_sol.
   * But they have been read in assuming they are in cgs. Convert this
   * only now to proper internal units to avoid float overflows. We only
   * store the energy that is to be distributed from this spart to its
   * neighbours in this step in internal units.*/
  const double solar_luminosity = 3.828e33; /* erg/s */
  for (int g = 0; g < RT_NGROUPS; g++) {
    const double emission_rate_internal_units =
        const_stellar_emission_rates[g] * solar_luminosity;
    emission_this_step[g] = emission_rate_internal_units * dt;
  }
}

/**
 * @brief Compute the energy emitted from a star during the time step dt.
 * This is for the Iliev+2006 Test 4.
 *
 * @param emission_this_step (return) the emitted radiation energy of a star
 * during the time interval dt
 * @param M the star mass (in internal units)
 * @param dt time step size (in internal units)
 * @param photon_number_integral Integrated photon numbers over frequency
 * interval
 * @param average_photon_energy average photon energy in each frequency bin, in
 * erg
 * @param phys_const struct holding physical constants
 * @param internal_units units struct containing internal units
 */
__attribute__((always_inline)) INLINE static void
rt_get_emission_this_step_IlievTest(
    double emission_this_step[RT_NGROUPS], float M, const double dt,
    const double photon_number_integral[RT_NGROUPS],
    const double average_photon_energy[RT_NGROUPS],
    const struct phys_const* phys_const,
    const struct unit_system* internal_units) {

  /* Note that this model uses the halo mass to determine the luminosity
   * of a source. I'm cheating the system here by storing the required halo
   * mass as the star mass. This is only ok because the test is supposed to
   * run with all dynamics and gravity turned off. */

  const double Omega_b = 0.043;
  const double Omega_0 = 0.27;
  const double m_p = phys_const->const_proton_mass;
  const double t_s = 3e6 * phys_const->const_year;
  const double f_gamma = 250.;
  const double Ndot_gamma = (f_gamma * M * Omega_b) / (Omega_0 * m_p * t_s);

  double Nsum = 0.;
  for (int g = 0; g < RT_NGROUPS; g++) Nsum += photon_number_integral[g];

  if (Nsum <= 0.) error("No photons in spectrum...???");

  const double energy_units =
      units_cgs_conversion_factor(internal_units, UNIT_CONV_ENERGY);
  for (int g = 0; g < RT_NGROUPS; g++) {
    const double fi = photon_number_integral[g] / Nsum;
    const double Ndot_i = fi * Ndot_gamma;
    /* average photon densities are in cgs! */
    const double Edot_i = average_photon_energy[g] * Ndot_i / energy_units;
    emission_this_step[g] = Edot_i * dt;
  }
}

/*
 * @brief compute the index of a lower bound bin for a non-uniform spaced 1D array_x.
 * It takes an array_x that indicate the positions, its size and a value x for which
 * we wish to find the lower bound. It returns only the index of the lower bound bin.
 *
 *
 * @param array_x Values of gridpoints
 * @param size    Length of array_x
 * @param x       Value within range of array_x for which we want to find the lower bound
 */
static int interpolate_1D_non_uniform_LowerBound(const double* array_x,
						const int size,
                                                const double x) {

  // Handle below bounds explicitly
    if (x <= array_x[0])
        return 0;

    // Handle above bounds explicitly
    if (x >= array_x[size - 1])
        return size - 2;

    // Find the lower bin: array_x[i] <= x < array_x[i+1]
    for (int i = 0; i < size - 1; ++i) {
        if (x < array_x[i + 1])
            return i;
    }

    // Should not reach here
    return size - 2;
}

/*
 * @brief compute an interpolated value from a 2d table (array_y) whose values are
 * given on a non-uniformly spaced grids. The first dimension uses the outer bounds
 * of the grid when given values below/above their min/max. The second dimension does the
 * same for its max but gives an error when below their min.
 *
 * @param array_x1 Values of gridpoints over first dimension
 * @param array_x2 Values of gridpoints over second dimension
 * @param array_y  Values of a 2D table that we interpolate, lengths must match size1 and size2
 * @param size1    Length of array_x1
 * @param size2    Length of array_x2
 * @param x1       Value within range of array_x1 to interpolate to 
 * here is specifically metallicity
 * @param x2       Value within range of array_x2 to interpolate to
 * here is specifically stellar age
 */

static double interpolate_2D_non_uniform(const double* array_x1,const double* array_x2,
                                                double** array_y,
                                                const int size1, const int size2,
                                                double x1, double x2) {

// First we check the bounds and if x1 and x2 are within it
// return 0 if Time is above star age threshold, give error when below lowest bin value (should always be 0 though if properly binned)
// check for star age
/* TODO: add debug macro for the error check. */

if (array_y == NULL) {
    printf("ERROR: array_y is NULL!\n");
    return 0; // Handle error properly
}

if (x2 > array_x2[size2-1]){
    //printf("\n");
    //printf("x2 Has exceeded the upper bound =  %f", x2);
    //printf("\n");
    //printf("Setting it to: %f",array_x2[size2-1]);
    x2 = array_x2[size2-1];

}
if (x2 < array_x2[0]){ // MAKE THIS A REAL ERROR
    //printf("\n");
    //printf("x2 =  %f", x2);
    //printf("\n");
    //printf("x2 IS TOO LOW!!!");
    error("Stellar age %e is below 0!", x2);
    //return 0;  // MAKE THIS A REAL ERROR
}

//check for metallicity
if (x1 < array_x1[0]){
    //printf("\n");
    //printf("x1 =  %f", x1);

    x1 = array_x1[0];
    //printf("\n");
    //printf("x1 =  %f", x1);
}

if (x1 > array_x1[size1-1]){
    x1 = array_x1[size1-1];
}
// Now we finished bound checks, consider doing this seperate?

// First we want to find the 4 indices for the corner values.
// We get the lower index for the Metals
const int Ind1 = interpolate_1D_non_uniform_LowerBound(array_x1,size1,x1);
// Now we get the lower index for the Times
const int Ind2 = interpolate_1D_non_uniform_LowerBound(array_x2,size2,x2);

if (Ind1 >= size1 - 1 || Ind2 >= size2 - 1) {
    printf("Index out of bounds: Ind1=%d, Ind2=%d (size1=%d, size2=%d)\n", Ind1, Ind2, size1, size2);
    printf("x1 = %e, x2 = %e\n", x1, x2);
    exit(1); // or return a safe default
}

if (Ind1 < 0 || Ind2 < 0) {
    error("Interpolation index out of bounds! Ind1=%d Ind2=%d", Ind1, Ind2);
}

double dx1 = array_x1[Ind1+1] - array_x1[Ind1];
double dx2 = array_x2[Ind2+1] - array_x2[Ind2];

// Handle very small differences (avoid division by zero)
if (fabs(dx1) < 1e-10 || fabs(dx2) < 1e-10) {
	printf("Warning: Small step sizes detected: dx1 = %g, dx2 = %g\n", dx1, dx2);
        exit(1);  // Or return a default value if this is recoverable
}

// We know the values we need are at [Ind1,Ind2], [Ind1 + 1,Ind2], [Ind1,Ind2 + 1] and [Ind1 + 1,Ind2 + 1]
// We calculate the offsets over the first dimension
const double offset1 =
      (array_x1[Ind1+1] - x1) / dx1;
// And now the second dimension
const double offset2 =
      (array_x2[Ind2+1] - x2) / dx2;

// First Interpolation:
// Interpolate [Ind1,Ind2] --- [Ind1+1,Ind2]   and interpolate [Ind1,Ind2+1] --- [Ind1+1,Ind2+1] with offset1 to make 2 new values.

if (array_y[Ind1] == NULL || array_y[Ind1+1] == NULL) {
    printf("ERROR: array_y[%d] or array_y[%d] is NULL!\n", Ind1, Ind1+1);
    exit(1);
}

const double Interpol1 = offset1 * array_y[Ind1][Ind2] + (1. - offset1) * array_y[Ind1+1][Ind2];
const double Interpol2 = offset1 * array_y[Ind1][Ind2+1] + (1. - offset1) * array_y[Ind1+1][Ind2+1];

const double result = offset2 * Interpol1 + (1. - offset2) * Interpol2;

if (isnan(result) || isinf(result)) {
    printf("Interpolation result is NaN or Inf\n");
    printf("Metallicity (x1) = %g, age (x2) = %g, offsets: %g, %g\n", x1, x2, offset1, offset2);
    printf("array_y values: %g %g %g %g\n",
        array_y[Ind1][Ind2], array_y[Ind1+1][Ind2],
        array_y[Ind1][Ind2+1], array_y[Ind1+1][Ind2+1]);
    exit(1);
}

return result;
}

/**
 * @brief Compute the energy emitted from a star during the time step dt.
 * This is using BPASS model with three photon bins.
 *
 * @param emission_this_step (return) the emitted radiation energy of a star
 * during the time interval dt
 * @param M the star mass (in internal units)
 * @param Metallicity the stellar metallicity (in internal units)
 * @param star_age_begin_of_step the star age at the beginning of this timestep
 * (in internal units)
 * @param star_age the star age at the end of this timestep(in internal units)
 * @param ionizing_tables BPASS table that stores the total emitted photon number 
 * at different stellar age
 * @param average_photon_energy average photon energy in each frequency bin, in
 * erg
 * @param phys_const struct holding physical constants
 * @param internal_units units struct containing internal units
 */
__attribute__((always_inline)) INLINE static void
rt_get_emission_this_step_BPASS(
    double emission_this_step[RT_NGROUPS], float M,float Metallicity, 
    double star_age_begin_of_step, double star_age,
    double ***ionizing_tables,
    const double average_photon_energy[RT_NGROUPS],
    const struct phys_const* phys_const,
    const struct unit_system* internal_units,
    const double f_esc) {

  /* Get the array for the table range.
   * TODO: Hard code in now, need to replace it. */ 
  static const double Metallicities[] = {1e-5, 1e-4, 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.1, 0.14, 0.2, 0.3, 0.4};
  int size_Metals = sizeof(Metallicities) / sizeof(Metallicities[0]);
  static const double age_100Myr[] = {
      0.0,        // 0
      1.0,        // 10^0
      1.258925,   // 10^0.1
      1.584893,   // 10^0.2
      1.995262,   // 10^0.3
      2.511886,   // 10^0.4
      3.162278,   // 10^0.5
      3.981072,   // 10^0.6
      5.011872,   // 10^0.7
      6.309573,   // 10^0.8
      7.943282,   // 10^0.9
      10.0,       // 10^1.0
      12.589254,  // 10^1.1
      15.848932,  // 10^1.2
      19.952623,  // 10^1.3
      25.118864,  // 10^1.4
      31.622777,  // 10^1.5
      39.810717,  // 10^1.6
      50.118723,  // 10^1.7
      63.095734,  // 10^1.8
      79.432823,  // 10^1.9
      100.0       // 10^2.0
  };
  int size_Times = sizeof(age_100Myr) / sizeof(age_100Myr[0]);
  const double normalized_mass = 1; //units of solar mass
  
  /* Convert some quantities. */
  const double time_to_Myr = units_cgs_conversion_factor(internal_units, UNIT_CONV_TIME) /
          (365.25f * 24.f * 60.f * 60.f * 1e6f);
  const double energy_units =
      units_cgs_conversion_factor(internal_units, UNIT_CONV_ENERGY);

  /* Calculate internal mass to solar mass conversion factor */
  const double Msun_cgs = phys_const->const_solar_mass *
                          units_cgs_conversion_factor(internal_units, UNIT_CONV_MASS);
  const double unit_mass_cgs = units_cgs_conversion_factor(internal_units, UNIT_CONV_MASS);
  const double mass_to_solar_mass = unit_mass_cgs / Msun_cgs;

  /* Get the converted quantities. */
  const double star_age_before_Myr = star_age_begin_of_step * time_to_Myr;
  const double star_age_now_Myr = star_age * time_to_Myr;
  const double star_mass_Msolar = M * mass_to_solar_mass;
  const double M_star_fraction = star_mass_Msolar / normalized_mass;

  if (Metallicity < Metallicities[0])
    Metallicity = Metallicities[0];
  else if (Metallicity > Metallicities[size_Metals - 1])
    Metallicity = Metallicities[size_Metals - 1];
  
  /* TODO: add debug check macro for error check. */
  if (star_age_before_Myr < 0 || isnan(star_age_before_Myr)) {
    error("Invalid star_age_before_Myr = %e", star_age_before_Myr);
  }

  for (int g = 0; g < RT_NGROUPS; g++) {
    if (star_age_before_Myr > 100.) {
	/* If the stellar age before this timestep is above 100Myr, we set the emission to zero*/
    	emission_this_step[g] = 0.;
    } else {
	    if (!ionizing_tables[g]) error("ionizing_tables[%d] is NULL", g);
	    for (int i = 0; i < size_Metals; i++) {
		    if (!ionizing_tables[g][i]) error("ionizing_tables[%d][%d] is NULL", g, i);
	    }
	    double N_total_before = interpolate_2D_non_uniform(Metallicities, age_100Myr,
                                                ionizing_tables[g],
                                                size_Metals, size_Times,
                                                Metallicity, star_age_before_Myr);
	    double N_total_now = interpolate_2D_non_uniform(Metallicities, age_100Myr,
                                                ionizing_tables[g],
                                                size_Metals, size_Times,
                                                Metallicity, star_age_now_Myr);
	    double N_emission_this_step = N_total_now - N_total_before;
	    
	    /* average photon densities are in cgs! */
	    const double E_g = f_esc * average_photon_energy[g] * N_emission_this_step * M_star_fraction / energy_units;
	    emission_this_step[g] = E_g;
	    
	    if (E_g < 0) {
		    error("Negative Photons??, N_total_before: %e, N_total_now: %e,star_age_before_Myr: %e, star_age_now_Myr: %e",N_total_before, N_total_now, star_age_before_Myr, star_age_now_Myr);
	    }

    //message("energy this step: %e, N_total_before: %e, N_total_now: %e, average_photon_energy[g]: %e, star_age_before_Myr: %e, star_age_now_Myr: %e, M_star_fraction: %e, N_emission_this_step:%e, Metalicities: %e", E_g, N_total_before, N_total_now, average_photon_energy[g], star_age_before_Myr, star_age_now_Myr, M_star_fraction, N_emission_this_step, Metallicity);
    }
  }
}

#endif /* SWIFT_RT_STELLAR_EMISSION_MODEL_KIARA_H */
