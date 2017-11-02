/***********************************************************************
/
/ Grackle c wrapper
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/
#ifndef SWIFT_COOLING_GRACKLE_WRAPPER_H
#define SWIFT_COOLING_GRACKLE_WRAPPER_H

#include "config.h"
#include "error.h"
#include <grackle.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



int wrap_init_cooling(char* CloudyTable,int UVbackground, double udensity, double ulength, double utime, 
                       int grackle_chemistry);
		       
		       
		       
//int wrap_update_UVbackground_rates(double auni);

int wrap_set_UVbackground_On();
int wrap_set_UVbackground_Off();

int wrap_get_cooling_time(double rho, double u, double Z, double a_now, double *coolingtime);

int wrap_do_cooling(double density, double *energy, double dtime, double Z, double a_now);




#endif /* SWIFT_COOLING_GRACKLE_WRAPPER_H */
