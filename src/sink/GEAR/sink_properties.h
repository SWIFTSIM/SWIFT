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
#include <feedback_properties.h>


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
  
  /*! Size of the calibration sample used to determine the probabilities
   * to form stellar particles with mass stellar_particle_mass */
  int size_of_calibration_sample;    
  
  /*! Mass of the stellar particle representing the low mass stars 
   * (continuous IMF sampling). */
  float stellar_particle_mass;  
  
  /*! Minimal mass of stars represented by discrete particles */
  float minimal_discrete_mass;  

  /*! Mass of the stellar particle representing the low mass stars 
   * (continuous IMF sampling). First stars */
  float stellar_particle_mass_first_stars;  
  
  /*! Minimal mass of stars represented by discrete particles. 
   * First stars. */
  float minimal_discrete_mass_first_stars;  
  
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
INLINE static void sink_props_init_probabilities(struct sink_props *sp, struct initial_mass_function *imf,
  const struct phys_const *phys_const,int first_stars) {

  /* random number generator */
  gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlux);
  /* set the seed */
  gsl_rng_set(r,0);
  

  /* get the slope of the last IMF portion */
  float exp      = imf->exp[imf->n_parts-1];
  
  /* get the IMF mass limits (all in Msol) */
  float mass_min = imf->mass_min;
  float mass_max = imf->mass_max;
  
  float minimal_discrete_mass;
  float stellar_particle_mass;
  
  if (!first_stars) {
    minimal_discrete_mass = sp->minimal_discrete_mass;
    stellar_particle_mass = sp->stellar_particle_mass/phys_const->const_solar_mass;
  } else {
    minimal_discrete_mass = sp->minimal_discrete_mass_first_stars;
    stellar_particle_mass = sp->stellar_particle_mass_first_stars/phys_const->const_solar_mass;  
  }  
  
  
  /* sanity check */
  if (minimal_discrete_mass < imf->mass_limits[imf->n_parts-1])  
    error("minimal_discrete_mass (=%8.3f) cannot be smaller than the mass limit (=%8.3f) of the last IMF segment,",
    minimal_discrete_mass, imf->mass_limits[imf->n_parts-1]);
    
  /* Compute the IMF mass below the minimal IMF discrete mass (continuous part) */
  double Mtot, Md, Mc;
  Mc = initial_mass_function_get_imf_mass_fraction(imf,mass_min,minimal_discrete_mass);
    
  
  if (Mc > 0) {
    Mtot = stellar_particle_mass/Mc;
    Md   = Mtot-stellar_particle_mass;
    Mc   = stellar_particle_mass;
  } else {
    Mtot = stellar_particle_mass;
    Md   = Mtot;
    Mc   = 0;
  }  
  
  message("Mass of the continuous part : %g",Mc);
  message("Mass of the discrete   part : %g",Md);
  message("Total IMF mass              : %g",Mtot);
  
  /* if no continous part, return */
  if (Mc == 0){
    imf->sink_Pc = 0;
    imf->sink_stellar_particle_mass = 0;
    message("probability of the continuous part : %g",0.); 
    message("probability of the discrete   part : %g",1.);
    return;
  }  
  
  
  struct f_zero_params params;
  
  params.Mc      = Mc;
  params.Md      = Md; 
  
  params.mmin    = minimal_discrete_mass;
  params.mmax    = mass_max;
  params.exp     = exp;
  params.n       = sp->size_of_calibration_sample;
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
  double root;
  
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      root = gsl_root_fsolver_root (s);
      double x_lo = gsl_root_fsolver_x_lower (s);
      double x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-4);
        
      message("%5d %.7f", iter, root);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  /* the final probabilities */
  double Pc = root;
  double Pd = 1 - Pc;
  imf->sink_Pc = Pc;
  imf->sink_stellar_particle_mass = Mc;
  
  message("probability of the continuous part : %g",Pc); 
  message("probability of the discrete   part : %g",Pd);
  

  /* free the solver */
  gsl_root_fsolver_free (s);
  

  /* free the random generator */
  gsl_rng_free (r);

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
INLINE static void sink_props_init(struct sink_props *sp, struct feedback_props *fp,
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
      
  sp->size_of_calibration_sample =
      parser_get_param_int(params, "GEARSink:size_of_calibration_sample");         
      
  sp->stellar_particle_mass =
      parser_get_param_float(params, "GEARSink:stellar_particle_mass");      
      
  sp->minimal_discrete_mass =
      parser_get_param_float(params, "GEARSink:minimal_discrete_mass");   
      
  sp->stellar_particle_mass_first_stars =
      parser_get_param_float(params, "GEARSink:stellar_particle_mass_first_stars");      
      
  sp->minimal_discrete_mass_first_stars =
      parser_get_param_float(params, "GEARSink:minimal_discrete_mass_first_stars");      
      
      
  /* Apply unit change */
  sp->maximal_temperature /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  sp->density_threashold /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

  sp->stellar_particle_mass*=phys_const->const_solar_mass;
  sp->stellar_particle_mass_first_stars*=phys_const->const_solar_mass;
  
  
  /* here, we need to differenciate between the stellar models */
  struct initial_mass_function *imf;
  struct stellar_model* sm;
  
  sm = &fp->stellar_model;
  imf = &sm->imf;
  
  
  /* Initialize for the stellar models (PopII) */
  sink_props_init_probabilities(sp, imf, phys_const, 0);
    
  /* Now initialize the first stars. */
  if (fp->metallicity_max_first_stars != -1) {
     sm = &fp->stellar_model_first_stars;
     imf = &sm->imf;
     sink_props_init_probabilities(sp, imf, phys_const, 1);  
  }
  
  message("maximal_temperature               = %g", sp->maximal_temperature);
  message("density_threashold                = %g", sp->density_threashold);
  message("size_of_calibration_sample        = %d", sp->size_of_calibration_sample);
  
  message("stellar_particle_mass             = %g", sp->stellar_particle_mass);
  message("minimal_discrete_mass             = %g", sp->minimal_discrete_mass);
  
  message("stellar_particle_mass_first_stars = %g", sp->stellar_particle_mass_first_stars);
  message("minimal_discrete_mass_first_stars = %g", sp->minimal_discrete_mass_first_stars);
 
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
