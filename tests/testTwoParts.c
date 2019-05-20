/*******************************************************************************
 *  * This file is part of SWIFT.
 *   * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *    *
 *     * This program is free software: you can redistribute it and/or modify
 *      * it under the terms of the GNU Lesser General Public License as published
 *       * by the Free Software Foundation, either version 3 of the License, or
 *        * (at your option) any later version.
 *         *
 *          * This program is distributed in the hope that it will be useful,
 *           * but WITHOUT ANY WARRANTY; without even the implied warranty of
 *            * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *             * GNU General Public License for more details.
 *              *
 *               * You should have received a copy of the GNU Lesser General Public License
 *                * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *                 *
 *                  ******************************************************************************/
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local includes */
#include "swift.h"


int main(int argc, char *argv[]){


  struct cell c;
  int i;

  float dt = 0.01;
  struct part *parts = NULL;
  parts = (struct part *) malloc(2*sizeof(struct part));
  bzero(parts, 2*sizeof(struct part));
  struct xpart *xparts = NULL;
  xparts = (struct xpart *)malloc(2*sizeof(struct xpart));
  bzero(xparts, 2*sizeof(struct xpart));

  parts[0].x[0] = 0.35;
  parts[0].x[1] = 0.0;
  parts[0].x[2] = 0.0;

  parts[1].x[0] = 0.65;
  parts[1].x[1] = 0.0;
  parts[1].x[2] = 0.0;

  parts[0].v[0] = 0.01;
  parts[0].v[1] = 0.;
  parts[0].v[2] = 0.;
  parts[0].v_minus1[0] = 0.01;
  parts[0].v_minus1[1] = 0.;
  parts[0].v_minus1[2] = 0.;

  parts[1].v[0] = -0.01;
  parts[1].v[1] = 0.;
  parts[1].v[2] = 0.;
  parts[1].v_minus1[0] = -0.01;
  parts[1].v_minus1[1] = 0.;
  parts[1].v_minus1[2] = 0.;

  parts[0].a_hydro[0] = 0.0;
  parts[0].a_hydro[1] = 0.0;
  parts[0].a_hydro[2] = 0.0;
  parts[1].a_hydro[0] = 0.0;
  parts[1].a_hydro[1] = 0.0;
  parts[1].a_hydro[2] = 0.0;

  xparts[0].v_full[0] = 0.01;
  xparts[0].v_full[1] = 0.;
  xparts[0].v_full[2] = 0.;

  xparts[1].v_full[0] = -0.01;
  xparts[1].v_full[1] = 0.;
  xparts[1].v_full[2] = 0.;

  parts[0].h = 0.13;
  parts[1].h = 0.13;
  parts[0].rho = 1000.0;
  parts[1].rho = 1000.0;
  parts[0].rho_t_minus1 = 1000.0;
  parts[1].rho_t_minus1 = 1000.0;
  parts[0].mass = 10.0;
  parts[1].mass = 10.0;

  c.hydro.parts = parts;
  c.hydro.xparts = xparts;
  c.hydro.count = 2;
  c.split = 0;

  eos.soundspeed = 2.21;
  eos.soundspeed_squared = 2.210*2.210;
  eos.density_reference = 1000.0;
  eos_print(&eos);
  struct engine eng;

  eng.time = 0.;
  eng.time_begin = 0.;
  eng.time_end = 1.;
  eng.dt_min = 0.01;
  eng.dt_max = 0.01;

  for(i = 0; i < 2048; i++) {
    eng.time_old = eng.time;
    eng.time += dt;
    
    hydro_reset_acceleration(&parts[0]);
    hydro_reset_acceleration(&parts[1]);

    float dx[3];
    dx[0] = parts[0].x[0] - parts[1].x[0];
    dx[1] = parts[0].x[1] - parts[1].x[1];
    dx[2] = parts[0].x[2] - parts[1].x[2];
    float r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    runner_iact_force(r2, dx, parts[0].h, parts[1].h, &parts[0], &parts[1], 0.0, 0.0);
    if(i%50 == 0){
      for(int j = 0; j < c.hydro.count; j++){
        struct part *p = &parts[j];
        float temp = p->rho;
        p->rho = p->rho + p->drho_dt*dt;
        p->rho_t_minus1 = temp;

        p->x[0] = p->x[0] + p->v[0]*dt + 0.5*p->a_hydro[0]*dt*dt;

        temp = p->v[0];
        p->v[0] = p->v[0] + p->a_hydro[0]*dt;
        p->v_minus1[0] = temp;

        temp = p->v[1];
        p->v[1] = p->v[1] + p->a_hydro[1]*dt;
        p->v_minus1[1] = temp;

        temp = p->v[2];
        p->v[2] = p->v[2] + p->a_hydro[2]*dt;
        p->v_minus1[2] = temp;
        p->pressure = pressure_from_density(p->rho);

      } 
    }else{ 
      for(int j = 0; j < c.hydro.count; j++){
        struct part *p = &parts[j];
        float temp = p->rho;
        p->rho = p->rho_t_minus1 + 2*dt*p->drho_dt;
        p->rho_t_minus1 = temp;
        p->x[0] = p->x[0] + p->v[0]*dt + 0.5*p->a_hydro[0]*dt*dt;
        p->x[1] = p->x[1] + p->v[1]*dt + 0.5*p->a_hydro[1]*dt*dt;
        p->x[2] = p->x[2] + p->v[2]*dt + 0.5*p->a_hydro[2]*dt*dt;
        temp = p->v[0];
        p->v[0] = p->v_minus1[0] + 2*p->a_hydro[0]*dt;
        p->v_minus1[0] = temp;
  
        temp = p->v[1];
        p->v[1] = p->v_minus1[1] + 2*p->a_hydro[1]*dt;
        p->v_minus1[1] = temp;
  
        temp = p->v[2];
        p->v[2] = p->v_minus1[2] + 2*p->a_hydro[2]*dt;
        p->v_minus1[2] = temp;
        p->pressure = pressure_from_density(p->rho);
      }
    }
    printf("parts 0 v = %f, parts 1 v = %f\n", parts[0].v[0], parts[1].v[0]);
  }

  printf("Testing Kernel\n");
  float wi_dx, wi;
  kernel_deval(0.8, &wi, &wi_dx);
  printf("Kernel from 0.8 is %f %f\n", wi, wi_dx);

  free(parts);
  free(xparts);

  




}
