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

#include <chemistry_data.h>
#include <grackle.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "error.h"

int wrap_init_cooling(char *CloudyTable, int UVbackground, double udensity,
                      double ulength, double utime, int grackle_chemistry);

int wrap_init_cooling(char *CloudyTable, int UVbackground, double udensity,
                      double ulength, double utime, int grackle_chemistry);

int wrap_set_UVbackground_on();

int wrap_set_UVbackground_off();

int wrap_get_cooling_time(double rho, double u, double Z, double a_now,
                          double *coolingtime);

int wrap_do_cooling(double density, double *energy, double dtime, double Z,
                    double a_now);

void grackle_print_data();

void cloudy_print_data(const cloudy_data c, const int print_mmw);

#endif /* SWIFT_COOLING_GRACKLE_WRAPPER_H */
