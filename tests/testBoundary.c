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

  float dt = 1e-3;
  struct part *parts = NULL;
  parts = (struct part *) malloc(7*sizeof(struct part));
  bzero(parts, 7*sizeof(struct part));
  struct xpart *xparts = NULL;
  xparts = (struct xpart *)malloc(7*sizeof(struct xpart));
  bzero(xparts, 7*sizeof(struct xpart));

  for(i = 0; i <= 2; i++){
    parts[i].x[0] = 0.4 + (0.1 * (float)i);
    parts[i].x[1] = 0.2;
    parts[i].x[2] = 0.5;

    parts[i+3].x[0] = 0.4 + (0.1 * (float)i);
    parts[i+3].x[1] = 0.3;
    parts[i+3].x[2] = 0.5;
  }
  
  parts[6].x[0] = 0.5;
  parts[6].x[1] = 1.35;
  parts[6].x[2] = 0.5;



  for(i = 0; i < 7; i++){
    parts[i].v[0] = 0.00;
    parts[i].v[1] = 0.;
    parts[i].v[2] = 0.;
    parts[i].v_minus1[0] = 0.0;
    parts[i].v_minus1[1] = 0.;
    parts[i].v_minus1[2] = 0.;
    if(i < 6){
      parts[i].is_boundary = 1;
      parts[i].a_constant[0] = 0.0;
      parts[i].a_constant[1] = 0.0;
      parts[i].a_constant[2] = 0.0;
    }else{
      parts[i].is_boundary = 0;
      parts[i].a_constant[0] = 0.0;
      parts[i].a_constant[1] = -9.81;
      parts[i].a_constant[2] = 0.0;
    }
    parts[i].a_hydro[0] = 0.0;
    parts[i].a_hydro[1] = 0.0;
    parts[i].a_hydro[2] = 0.0;
    parts[i].h = 0.13;
    parts[i].rho = 1000.0;
    parts[i].rho_t_minus1 = 1000.0;
    parts[i].mass = 10.0;
    xparts[i].v_full[0] = 0.00;
    xparts[i].v_full[1] = 0.;
    xparts[i].v_full[2] = 0.;
  }



  c.hydro.parts = parts;
  c.hydro.xparts = xparts;
  c.hydro.count = 7;
  c.split = 0;

  eos.soundspeed = 221.0;
  eos.soundspeed_squared = 221.0*221.0;
  eos.density_reference = 1000.0;
  eos_print(&eos);
  struct engine eng;

  eng.time = 0.;
  eng.time_begin = 0.;
  eng.time_end = 2.;
  eng.dt_min = 0.0001;
  eng.dt_max = 0.0001;

  for(i = 0; i < 20000048 && eng.time < eng.time_end; i++) {
    eng.time_old = eng.time;
    eng.time += dt;
    for(int j = 0; j < c.hydro.count; j++){
      hydro_reset_acceleration(&parts[j]);
    }
    for(int j = 0; j < c.hydro.count; j++){ 

      for(int k = j+1; k < c.hydro.count; k++){
        float dx[3];
        dx[0] = parts[j].x[0] - parts[k].x[0];
        dx[1] = parts[j].x[1] - parts[k].x[1];
        dx[2] = parts[j].x[2] - parts[k].x[2];
        float r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        runner_iact_force(r2, dx, parts[j].h, parts[k].h, &parts[j], &parts[k], 0.0, 0.0);
      }
    }

    for(int j = 0; j < c.hydro.count; j++){
/*      if(!parts[j].is_boundary){
         printf("a_hydro[1] = %f\n", parts[j].a_hydro[1]);
      }*/
      parts[j].a_hydro[0] += parts[j].a_constant[0];
      parts[j].a_hydro[1] += parts[j].a_constant[1];
      parts[j].a_hydro[2] += parts[j].a_constant[2];
      if(parts[j].is_boundary){
        parts[j].a_hydro[0] = 0.0;
        parts[j].a_hydro[1] = 0.0;
        parts[j].a_hydro[2] = 0.0;
      }
    }

    if(i%50 == 0){
      for(int j = 0; j < c.hydro.count; j++){
        struct part *p = &parts[j];
        float temp = p->rho;
        p->rho = p->rho + p->drho_dt*dt;
        p->rho_t_minus1 = temp;

        p->x[0] = p->x[0] + p->v[0]*dt + 0.5*p->a_hydro[0]*dt*dt;
        p->x[1] = p->x[1] + p->v[1]*dt + 0.5*p->a_hydro[1]*dt*dt;
        p->x[2] = p->x[2] + p->v[2]*dt + 0.5*p->a_hydro[2]*dt*dt;

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
    printf("%f, %f, %f\n", parts[6].v[1], parts[6].x[1], parts[6].rho);
  }


  printf("Final position: [%f %f %f]\n", parts[6].x[0], parts[6].x[1], parts[6].x[2]);
  free(parts);
  free(xparts);

  




}
