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
    
    /* do we set a constant specific energy as injection */
    int set_const_u;

    /* if yes, to what value */ 
    double u_inj;

    /* do we kinetically kick particles, based on density */
    int set_v;

    /* do we kinetically kick particles, with constant v */
    int set_const_v;

    /* if yes, what velocity do we give particles */
    double vel_inj;

    /* keep track of times time-condition was valid */
    int counter;

    /* keep track if the final injection has happened already */
    int final_injection;
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
        /*
         * Since we are looking at one particle, we don't know
         * the amount of particles affected. Thus approximate
         * N_part * m = rho * V_inj
         *
	      */
      if (terms->set_const_u == 1) {
        xp->u_full = terms->u_inj; 
      }
      else if (terms->set_v == 1) {
        double v = sqrtf(2 * terms->E_inj / (p->rho * terms->V_inj));
        
        xp->v_full[0] = (dx / distance) * v;
        xp->v_full[1] = (dy / distance) * v;
        xp->v_full[2] = (dz / distance) * v;
      }
      else if (terms->set_const_v == 1) {
        xp->v_full[0] = (dx / distance) * terms->vel_inj;
        xp->v_full[1] = (dy / distance) * terms->vel_inj;
        xp->v_full[2] = (dz / distance) * terms->vel_inj;
      }
      else {
        xp->u_full = terms->E_inj / (p->rho * terms->V_inj);
      }
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
    double time, const struct forcing_terms* terms,
    const struct phys_const* phys_const, const struct part* p,
    const struct xpart* xp) {

  return FLT_MAX;
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
  if (terms->set_const_u == 1) {
    message("Balsara-Kim specific energy injection u_inj: %f", terms->u_inj);
  }
  else if (terms->set_const_v == 1) {
    message("Balsara-Kim constant velocity kick v_inj: %f", terms->vel_inj);
  }
  else if (terms->set_v == 1) {
    message("Balsara-Kim density dependent velocity kick with E_inj: %f", 
      terms->E_inj);
  }
  else {
    message("Balsara-Kim density dependent specific energy injection with E_inj: %f",
      terms->E_inj);
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
    
  /* Read the injection energy, defaults to 1e51 erg*/
  double E_inj = parser_get_opt_param_double(parameter_file,
        "BalsaraKimForcing:E_inj", 1e51);

  /* What energy injection scheme do we use? */
  terms->set_const_u = parser_get_opt_param_int(parameter_file,
	"BalsaraKimForcing:set_const_u", 0);
  double u_inj = parser_get_opt_param_double(parameter_file,
	"BalsaraKimForcing:u_inj", 2.8263e15);

  terms->set_const_v = parser_get_opt_param_int(parameter_file,
	"BalsaraKimForcing:set_const_v", 0);
  double vel_inj = parser_get_opt_param_double(parameter_file,
  "BalsaraKimForcing:v_inj", 200) * 1.e5;

  terms->set_v = parser_get_opt_param_int(parameter_file,
	"BalsaraKimForcing:set_v", 0);

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
