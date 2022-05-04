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

/* Standard header */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>

/* Local header */
#include <random.h>

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


struct f_zero_params {
    double Mc;
    double Md; 
    double mmin;
    double mmax;
    double exp;
    int n;
    gsl_rng * r;
};


INLINE static double f_zero(double Pc, void *params) {

  struct f_zero_params *p
    = (struct f_zero_params *) params;

  double Mc   = p->Mc;
  double Md   = p->Md; 
  double mmin = p->mmin;
  double mmax = p->mmax;
  double exp  = p->exp;
  int   n    = p->n;
  gsl_rng *r = p->r;
  

  /* mass of the continuous part of the IMF */ 
  double Mcs = 0;  
  
  /* mass of the discrete part of the IMF */ 
  double Mds = 0;
  
  for (int i = 0; i < n; i++) {
    
    double xc = gsl_rng_uniform(r);
    
    /* discrete part */
    if (xc < Pc)
      Mcs = Mcs + Mc;
      
    /* continuous part */  
    else {
      double m = random_sample_power_law(mmin, mmax, exp, gsl_rng_uniform(r));
      Mds = Mds + m;
    }
  }
    
  return (double) (Mcs/Mds) - (Mc/Md);
};




/**
 * @brief Initialise the probabilities to get a stellar mass (continuous
 * sampling of the IMF)
 *
 * @param sp The #sink_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param cosmo The cosmological model.
 */
INLINE static void sink_props_init_probabilities(struct sink_props *sp) {

  /* random generator */
  gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlux);
  
  /* set the seed */
  gsl_rng_set(r,0);
    
  
  struct f_zero_params params;
  params.Mc      = 45.02843206491189;
  params.Md      = 4.971567793485943; 
  params.mmin    = 8;
  params.mmax    = 30;
  params.exp     = -1.3;
  params.n       = 100000;
  params.r       = r;
  
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  
  F.function = &f_zero;
  F.params = &params;
  
  T = gsl_root_fsolver_bisection;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, 0, 0.99);  
  
  int status;
  int iter = 0;
  int max_iter = 100;
  
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      double root = gsl_root_fsolver_root (s);
      double x_lo = gsl_root_fsolver_x_lower (s);
      double x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-4);
        
      message("%5d %.7f", iter, root);
    }
  while (status == GSL_CONTINUE && iter < max_iter);





  /* free the solver */
  gsl_root_fsolver_free (s);
  

  /* free the random generator */
  gsl_rng_free (r);

  exit(-1);

}


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
  
  
  sink_props_init_probabilities(sp);
  

  message("maximal_temperature         = %g", sp->maximal_temperature);
  message("density_threashold          = %g", sp->density_threashold);
  message("stellar_particle_mass       = %g", sp->stellar_particle_mass);
  message("minimal_discrete_mass       = %g", sp->minimal_discrete_mass);
  message("size_of_calibration_sample  = %d", sp->size_of_calibration_sample);


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
