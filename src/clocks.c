/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/**
 *  @file clocks.c
 *  @brief support for measuring intervals in seconds. Use cycle.h
 *  or timers.h for relative times.
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <stdio.h>
#include <unistd.h>

/* Local headers. */
#include "clocks.h"

/* 0.25 of a second in nanoseconds. */
#define SLEEPTIME 250000000

/* The CPU frequency used to convert ticks to seconds. */
static unsigned long long cpufreq = 0;

/* Local prototypes. */
static void estimate_cpufreq();

/**
 * @brief Get the current time.
 *
 * @param time the current time.
 */
void clocks_gettime(struct clockstime *time) {

#ifdef HAVE_CLOCK_GETTIME
  clock_gettime(CLOCK_REALTIME, &time->time);
#else
  time->time = getticks();
#endif
}

/**
 * @brief Get difference in milli-seconds between two times.
 *
 * @param start the start time.
 * @param end the end time.
 *
 * @return the difference in milli-secinds.
 */
double clocks_diff(struct clockstime *start, struct clockstime *end)
{
#ifdef HAVE_CLOCK_GETTIME
  struct timespec temp;
  if ((end->time.tv_nsec - start->time.tv_nsec) < 0) {
    temp.tv_sec = end->time.tv_sec - start->time.tv_sec - 1;
    temp.tv_nsec = 1000000000 + end->time.tv_nsec - start->time.tv_nsec;
  } else {
    temp.tv_sec = end->time.tv_sec - start->time.tv_sec;
    temp.tv_nsec = end->time.tv_nsec - start->time.tv_nsec;
  }
  return (double)temp.tv_sec * 1000.0 + (double)temp.tv_nsec * 1.0E-6;
#else
  return elapsed(end->time, start-time) / clocks_get_cpufreq() * 1000;
#endif

}

/**
 * @brief Set the CPU frequency.
 *
 * This function should be called at least once to set the CPU frequency.
 * To use the builtin estimation techniques give a value of 0.
 *
 * @param freq the CPU frequency in Hz or 0 to estimate one.
 */
void clocks_set_cpufreq(unsigned long long freq) {
  if ( freq > 0 ) {
    cpufreq = freq;
  } else {
    estimate_cpufreq();
  }
}

/**
 * @brief Get the CPU frequency in Hz.
 *
 * @result the CPU frequency.
 */
unsigned long long clocks_get_cpufreq() {

  if (cpufreq > 0)
    return cpufreq;

  /* It not already set estimate it. */
  estimate_cpufreq();
  return cpufreq;

}

/**
 * @brief Estimate the CPU frequency in Hz.
 *
 * If already set return the CPU frequency, then estimate the CPU frequency.
 *
 * The technique is either use a clock timed nanosleep (this was the best
 * method on i7), to read the value from the cpuinfo_max_freq
 * file (probably a overestimate), to use the macro value CPU_TPS or
 * finally just use a value of 1.
 */
static void estimate_cpufreq() {

#ifdef HAVE_CLOCK_GETTIME
  /* Try to time a nanosleep() in ticks. */
  struct clockstime time1;
  struct clockstime time2;

  struct timespec sleep;
  sleep.tv_sec = 0;
  sleep.tv_nsec = SLEEPTIME;

  clocks_gettime(&time1);
  ticks tic = getticks();

  /* Could do some calculation, but constant_tsc should protect us. */
  nanosleep(&sleep, NULL);

  clocks_gettime(&time2);
  ticks toc = getticks();
  double realsleep = clocks_diff(&time1, &time2);

  cpufreq = (signed long long) (double)(toc - tic) * 1.0/realsleep * 1000.0;
#endif

  /* Look for the system value, if available. Tends to be too large. */
#ifdef __linux__
  if (cpufreq == 0) {
    FILE *file = fopen("/sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq",
                       "r");
    if (file != NULL) {
      unsigned long long maxfreq;
      if (fscanf(file, "%llu", &maxfreq) == 1) {
        cpufreq = maxfreq * 1000;
      }
      fclose(file);
    }
  }
#endif

  /* Nearly final attempt */
#ifdef CPU_TPS
  if (cpufreq == 0)
    cpufreq = CPU_TPS;
#endif

  /* If all fails just report ticks in any times. */
  if (cpufreq == 0)
    cpufreq = 1;
}

/**
 * @brief Return the difference between two ticks in seconds.
 *
 * Only an approximation as based on how well we have estimated the
 * rtc frequency. Should be good for machines that support constant_rtc
 * and clock_gettime().
 *
 * @param tic a number of ticks returned by the cycle.h getticks() function.
 * @param toc a number of ticks returned by the cycle.h getticks() function.
 *
 * @result the absolute difference in approximated seconds.
 */
double clocks_diff_ticks(ticks tic, ticks toc)
{
    return clocks_from_ticks(tic - toc);
}

/**
 * @brief Convert a number of ticks into seconds.
 *
 * Only an approximation as based on how well we have estimated the
 * rtc frequency. Should be good for machines that support constant_rtc
 * and clock_gettime().
 *
 * @param tics a number of ticks returned by the cycle.h getticks() function.
 *
 * @result the approximated seconds.
 */
double clocks_from_ticks(ticks tics)
{
    return ((double)tics / (double)clocks_get_cpufreq() * 1000.0);
}
