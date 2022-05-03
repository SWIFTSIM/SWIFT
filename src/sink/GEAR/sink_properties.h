/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_DEFAULT_SINK_PROPERTIES_H
#define SWIFT_DEFAULT_SINK_PROPERTIES_H

/**
 * @brief Properties of sink in the Default model.
 */
struct sink_props {

  /*! Cut off radius */
  float cut_off_radius;

  /*! Maximal gas temperature for forming a star. */
  float maximal_temperature;

  /*! Minimal gas density for forming a star. */
  float density_threashold;
  
  /*! Mass of the stellar particle representing the low mass stars 
   * (continuous IMF sampling). */
  float stellar_particle_mass;  
  
  /*! Minimal mass of stars represented by discrete particles */
  float minimal_discrete_mass;  

  /*! Size of the calibration sample used to determine the probabilities
   * to form stellar particles with mass stellar_particle_mass */
  int size_of_calibration_sample;  
  
  
};

/**
 * @brief Initialise the sink properties from the parameter file.
 *
 * @param sp The #sink_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
INLINE static void sink_props_init(struct sink_props *sp,
                                   const struct phys_const *phys_const,
                                   const struct unit_system *us,
                                   struct swift_params *params,
                                   const struct cosmology *cosmo) {

  sp->cut_off_radius =
      parser_get_param_float(params, "GEARSink:cut_off_radius");

  sp->maximal_temperature =
      parser_get_param_float(params, "GEARSink:maximal_temperature");

  sp->density_threashold =
      parser_get_param_float(params, "GEARSink:density_threashold");
      
  sp->stellar_particle_mass =
      parser_get_param_float(params, "GEARSink:stellar_particle_mass");      
      
  sp->minimal_discrete_mass =
      parser_get_param_float(params, "GEARSink:minimal_discrete_mass");   
      
  sp->size_of_calibration_sample =
      parser_get_param_int(params, "GEARSink:size_of_calibration_sample");        
      
      
      

  /* Apply unit change */
  sp->maximal_temperature /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  sp->density_threashold /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

  sp->stellar_particle_mass*=phys_const->const_solar_mass;

  message("maximal_temperature    = %g", sp->maximal_temperature);
  message("density_threashold     = %g", sp->density_threashold);
  message("stellar_particle_mass  = %g", sp->stellar_particle_mass);
}

/**
 * @brief Write a sink_props struct to the given FILE as a stream of
 * bytes.
 *
 * @param props the sink properties struct
 * @param stream the file stream
 */
INLINE static void sink_struct_dump(const struct sink_props *props,
                                    FILE *stream) {
  restart_write_blocks((void *)props, sizeof(struct sink_props), 1, stream,
                       "sink props", "Sink props");
}

/**
 * @brief Restore a sink_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the sink properties struct
 * @param stream the file stream
 */
INLINE static void sink_struct_restore(const struct sink_props *props,
                                       FILE *stream) {
  restart_read_blocks((void *)props, sizeof(struct sink_props), 1, stream, NULL,
                      "Sink props");
}

#endif /* SWIFT_DEFAULT_SINK_PROPERTIES_H */
