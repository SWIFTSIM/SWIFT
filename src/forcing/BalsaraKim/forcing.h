#ifndef SWIFT_FORCING_BALSARAKIM_H
#define SWIFT_FORCING_BALSARAKIM_H

/* Config parameters */
#include <config.h>

/* Standard includes. */
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* Local includes. */
#include "error.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "units.h"
#include "hydro.h"
#include "periodic.h"
#include "timeline.h"
#include "timestep_sync_part.h"

enum mechanism {
  SET_U       = 0,
  SET_CONST_U = 1,
  SET_V       = 2,
  SET_CONST_V = 3
};

struct forcing_terms {
    /* Amount of SN injections */
    size_t size;

    /* injection energy */
    double E_inj;

    /* injection radius */
    double r_inj;

    /* injection volume */
    double V_inj;

    /* injection times */
    double *times;

    /* injection coordinates */
    double *x_SN;
    double *y_SN;
    double *z_SN;

    /* next event injection index */
    int t_index;

    /* injection mechanism */
    enum mechanism injection_model;

    /* injection model: set_const_u, specific energy injection */ 
    double u_inj;

    /* injection model: set_const_v, velocity kick */
    double vel_inj;

    /* keep track of times time-condition was valid */
    int counter;

    /* keep track if the final injection has happened already */
    int final_injection;

    /* amplification factor for the dedner field */
    float dedner_amp;

    /* min and max timestep, for timestep computation purposes */
    float dt_min;
    float dt_max;

};

/**
 * @brief Computes the forcing terms.
 *
 * @param time The current time.
 * @param terms The properties of the forcing terms.
 * @param s The #space we act on.
 * @param phys_const The physical constants in internal units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_terms_apply(
    const double time, struct forcing_terms* terms, const struct space* s,
    const struct phys_const* phys_const, struct part* p, struct xpart* xp) {

  const int t_index = terms->t_index;

  if ((time >= terms->times[t_index]) &&
      (terms->final_injection == 0)) {
    terms->counter++;

    float SN_loc[3];
    SN_loc[0] = terms->x_SN[t_index];
    SN_loc[1] = terms->y_SN[t_index];
    SN_loc[2] = terms->z_SN[t_index];
    
    const int periodic = s->periodic;
    const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

    double dx = SN_loc[0] - p->x[0];
    double dy = SN_loc[1] - p->x[1];
    double dz = SN_loc[2] - p->x[2];

    if (periodic) {
          dx = nearest(dx, dim[0]);
          dy = nearest(dy, dim[1]);
          dz = nearest(dz, dim[2]);
        }
    
    double distance = sqrtf( dx * dx + dy * dy + dz * dz);

    if (distance <= terms->r_inj) {
      message("applying injection to particle %lld at t: %f", p->id, 
        time);

      /* store old specific energy and velocity */
      const double u_old = xp->u_full;
      const double v_old = sqrtf(xp->v_full[0]*xp->v_full[0] + 
                                 xp->v_full[1]*xp->v_full[1] +
                                 xp->v_full[2]*xp->v_full[2]);

      double u_new;
      double v_new;

      /* inject energy according to specified model */
      enum mechanism injection_model = terms->injection_model;

      switch (injection_model) {

        case SET_U:
        
          /* compute new specific energy */
          u_new = terms->E_inj / (p->rho * terms->V_inj);

          /* set the specific energy */
          xp->u_full = u_new;

          /* store injected energy */
          xp->forcing_data.forcing_injected_energy += (u_new - u_old) * p->mass;

          break;
        
        case SET_CONST_U:

          /* get new specific energy */
          u_new = terms->u_inj;

          /* set the specific energy */
          xp->u_full = u_new;

          /* store injected energy */
          xp->forcing_data.forcing_injected_energy += (u_new - u_old) * p->mass;
          break;

        case SET_V:

          /* compute new velocity */
          v_new = sqrtf(2 * terms->E_inj / (p->rho * terms->V_inj));
        
          /* set the velocity */
          xp->v_full[0] = (dx / distance) * v_new;
          xp->v_full[1] = (dy / distance) * v_new;
          xp->v_full[2] = (dz / distance) * v_new;

          /* store injected energy */
          xp->forcing_data.forcing_injected_energy += 
                p->mass * (v_new * v_new - v_old * v_old) / 2;
          break;

        case SET_CONST_V:

          /* get new velocity */
          v_new = terms->vel_inj;
        
          /* set the velocity */
          xp->v_full[0] = (dx / distance) * v_new;
          xp->v_full[1] = (dy / distance) * v_new;
          xp->v_full[2] = (dz / distance) * v_new;

          /* store injected energy */
          xp->forcing_data.forcing_injected_energy += 
                p->mass * (v_new * v_new - v_old * v_old) / 2;
          break;
        
        default:

          error("no injection model specified");
          
      }

      /* increase dedner field */
      xp->mhd_data.psi_over_ch_full *= terms->dedner_amp;
      
      timestep_sync_part(p);
    }
  }
}

/**
 * @brief Computes the time-step condition due to the forcing terms.
 *
 * Nothing to do here. --> Return FLT_MAX.
 *
 * @param time The current time.
 * @param terms The properties of the forcing terms.
 * @param phys_const The physical constants in internal units.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static float forcing_terms_timestep(
    const double time, const struct forcing_terms *terms, 
    const struct space *s, const struct phys_const *phys_const, 
    const struct part *p, const struct xpart *xp) {
 
  const int t_index = terms->t_index;

  /* Don't do anything if there is no injection happening in the near future */
  if (((terms->times[t_index] - time) > terms->dt_max) ||
      (terms->final_injection == 1)) {
    return FLT_MAX;
  }
  
  /* Get injection coordinates of the next injection */
  float SN_loc[3];
  SN_loc[0] = terms->x_SN[t_index];
  SN_loc[1] = terms->y_SN[t_index];
  SN_loc[2] = terms->z_SN[t_index];
  
  /* Compute distance to center injection volume */
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  
  double dx = SN_loc[0] - p->x[0];
  double dy = SN_loc[1] - p->x[1];
  double dz = SN_loc[2] - p->x[2];

  if (periodic) {
    dx = nearest(dx, dim[0]);
    dy = nearest(dy, dim[1]);
    dz = nearest(dz, dim[2]);
  }

  /* Distance to the injection volume, is allowed to be negative */
  double distance = sqrtf( dx * dx + dy * dy + dz * dz ) - terms->r_inj;

  /* Compute an upper limit estimate of the signal velocity */
  float v_abs2 = xp->v_full[0] * xp->v_full[0] +
		             xp->v_full[1] * xp->v_full[1] +
		             xp->v_full[2] * xp->v_full[2];

  float v_sig = sqrtf( p->viscosity.v_sig * p->viscosity.v_sig + v_abs2 );
		   
  /* Could the particle be in the injection volume in dt_max time*/
  if (distance / v_sig < terms->dt_max) {

    /* Compute timestep such that particle would be active */
    float new_dt_forcing;
    new_dt_forcing = fmaxf(terms->times[t_index] - time, terms->dt_min);

    return new_dt_forcing;
  }
  else {
    return FLT_MAX;
  }  
}

/**  
 * @brief updates the forcing terms
 *
 * increases the current supernova index after one has happend 
 * 
 * @param terms The #forcing_terms properties of the run
 * @param time_old The previous system time
 */
void forcing_update(struct forcing_terms *terms, const double time_old);

/**
 * @brief Prints the properties of the forcing terms to stdout.
 *
 * @param terms The #forcing_terms properties of the run.
 */
static INLINE void forcing_terms_print(const struct forcing_terms* terms) {
  /* Print energy injection mechanism */
  enum mechanism injection_model = terms->injection_model;

  switch (injection_model) {

    case SET_U:
    
      message("Balsara-Kim density dependent specific energy injection with E_inj: %f",
        terms->E_inj);
      break;
    
    case SET_CONST_U:

      message("Balsara-Kim specific energy injection u_inj: %f", terms->u_inj);
      break;

    case SET_V:

      message("Balsara-Kim density dependent velocity kick with E_inj: %f", 
        terms->E_inj);
      break;

    case SET_CONST_V:

      message("Balsara-Kim constant velocity kick v_inj: %f", terms->vel_inj);
      break;
    
    default:

      error("no injection model specified");
      
  }

  message("Injection radius r_inj: %f", terms->r_inj);
  
  size_t i;
  for (i = 0; i < terms->size; i++) {
    message("SN at t: %f, at [x,y,z]: [%.2f,%.2f,%.2f]", terms->times[i],
    terms->x_SN[i], terms->y_SN[i], terms->z_SN[i]);
  }
}

/**
 * @brief Initialises the forcing term properties
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param s The #space object.
 * @param terms The forcing term properties to initialize
 */
static INLINE void forcing_terms_init(struct swift_params* parameter_file,
                                      const struct phys_const* phys_const,
                                      const struct unit_system* us,
                                      const struct space* s,
                                      struct forcing_terms* terms) {

  /* Read the filename containing the SN times & coordinates */
  char coords_filename[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(parameter_file,"BalsaraKimForcing:coords", 
                          coords_filename);

  /* 
   * Read the file containing the SN times & coordinates
   * Store the read values in terms
   */
  
  /* Open file */
  FILE *file = fopen(coords_filename, "r");
  if (file == NULL) error("Error opening file '%s'", coords_filename);

  /* Count number of lines */
  size_t len = 0;
  char *line = NULL;
  size_t nber_line = 0;
  while (getline(&line, &len, file) != -1) nber_line++;

  terms->size = nber_line - 1; /* Do not count header */

  /* Return to start of file and initialize time array */
  fseek(file, 0, SEEK_SET);
  terms->times = (double *)malloc(sizeof(double) * terms->size);
  terms->x_SN = (double *)malloc(sizeof(double) * terms->size);
  terms->y_SN = (double *)malloc(sizeof(double) * terms->size);
  terms->z_SN = (double *)malloc(sizeof(double) * terms->size);

  if ((!terms->times) || (!terms->z_SN)) {
    error(
        "Unable to malloc output_list. "
        "Try reducing the number of lines in %s",
        coords_filename);
  }

  /* Read header */
  if (getline(&line, &len, file) == -1)
    error("Unable to read header in file '%s'", coords_filename);

  /* Remove end of line character */
  line[strcspn(line, "\n")] = 0;

  /* Read file */
  size_t ind = 0;
  int read_successfully = 0;
  while (getline(&line, &len, file) != -1) {

    double *time = &terms->times[ind];
    double *x_SN = &terms->x_SN[ind];
    double *y_SN = &terms->y_SN[ind];
    double *z_SN = &terms->z_SN[ind];

    /* Write data to output_list */
    read_successfully = sscanf(line, "%lf,%lf,%lf,%lf", time,
                               x_SN, y_SN, z_SN) == 4;

    if (!read_successfully) {
      error(
          "Tried parsing output_list but found '%s' with illegal "
          "characters in file '%s'.",
          line, coords_filename);
    }
    
    /* convert the times [Myr] and coordinates [L_box] to IU */
    terms->times[ind] *= 1.e6 * phys_const->const_year;
    terms->x_SN[ind]  *= s->dim[0];
    terms->y_SN[ind]  *= s->dim[1];
    terms->z_SN[ind]  *= s->dim[2];
    
    ind += 1;
  }

  /* Cleanup */
  free(line);

  if (ind != terms->size)
    error("Did not read the correct number of output times.");

  /* Check that the list is in monotonic order */
  for (size_t i = 1; i < terms->size; ++i) {

    if (terms->times[i] <= terms->times[i - 1])
      error("Output list not having monotonically increasing ages.");
  }

  fclose(file);

  /* Read the injection radius, defaults to 5 pc*/
  double r_inj = parser_get_opt_param_double(parameter_file,
        "BalsaraKimForcing:r_inj", 5);

  /* Calculate & store injection volume, saves calculations */
  double V_inj = (4/3) * M_PI * pow(r_inj, 3);

  /* What energy injection scheme do we use? */
  char injection_model[20];
  parser_get_param_string(parameter_file, "BalsaraKimForcing:inj_model",
        injection_model);
  if (strcmp(injection_model, "set_u") == 0) {
    terms->injection_model = SET_U;
  }
  else if (strcmp(injection_model, "set_const_u") == 0) {
    terms->injection_model = SET_CONST_U;
  }
  else if (strcmp(injection_model, "set_v") == 0) {
    terms->injection_model = SET_V;
  }
  else if (strcmp(injection_model, "set_const_v") == 0) {
    terms->injection_model = SET_CONST_V;
  }
  else {
    error("unknown injection model '%s'", injection_model);
  }
  
  /* Read injection parameters*/
  double E_inj = parser_get_opt_param_double(parameter_file,
        "BalsaraKimForcing:E_inj", 1e51);
  double u_inj = parser_get_opt_param_double(parameter_file,
	      "BalsaraKimForcing:u_inj", 2.8263e15);
  double vel_inj = parser_get_opt_param_double(parameter_file,
        "BalsaraKimForcing:v_inj", 200) * 1.e5;

  /* Read dedner handling */
  terms->dedner_amp = parser_get_opt_param_float(parameter_file,
         "BalsaraKimForcing:psi_amp_factor", 1.);

  /* Read min and max timestep to store here as well */
  terms->dt_min = parser_get_param_float(parameter_file,
          "TimeIntegration:dt_min");
  terms->dt_max = parser_get_param_float(parameter_file,
          "TimeIntegration:dt_max");

  /* convert everything to internal units */
  double parsec_ui = phys_const->const_parsec;
  double erg_cgs = units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
  double v_cgs = us->UnitLength_in_cgs / us->UnitTime_in_cgs;

  terms->r_inj = r_inj * parsec_ui;
  terms->V_inj = V_inj * pow(parsec_ui, 3);
  terms->E_inj = E_inj / erg_cgs;
  terms->u_inj = u_inj / pow(v_cgs, 2);
  terms->vel_inj = vel_inj / v_cgs;

  /* initialise some indices and counters */
  terms->t_index = 0;
  terms->counter = 0;
  terms->final_injection = 0;
}
#endif /* SWIFT_FORCING_BALSARAKIM_H */
