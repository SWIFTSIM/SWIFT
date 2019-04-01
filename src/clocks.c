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
 *  @brief support for measuring intervals in milli seconds, when that
 *  is possible, otherwise ticks.
 *
 *  Use cycle.h or timers.h for relative times.
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <limits.h>
#include <stdio.h>
#include <unistd.h>

/* Local headers. */
#include "clocks.h"

/* 0.25 of a second in nanoseconds. */
#define SLEEPTIME 250000000

/* The CPU frequency used to convert ticks to seconds. */
static unsigned long long clocks_cpufreq = 0;

/* Ticks when the CPU frequency was initialised, this marks the start of
 * time. */
ticks clocks_start_ticks = 0;

/* The units of any returned times. */
static const char *clocks_units[] = {"ms", "~ms"};
static int clocks_units_index = 0;
static double clocks_units_scale = 1000.0;

/* Local prototypes. */
static void clocks_estimate_cpufreq(void);

/**
 * @brief Get the current time.
 *
 * @param time the current time.
 */
void clocks_gettime(struct clocks_time *time) {

#ifdef HAVE_CLOCK_GETTIME
  clock_gettime(CLOCK_REALTIME, &time->time);
#else
  time->time = getticks();
#endif
}

/**
 * @brief Get difference in between two times.
 *
 * @param start the start time.
 * @param end the end time.
 *
 * @return the difference.
 */
double clocks_diff(struct clocks_time *start, struct clocks_time *end) {
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
  return elapsed(end->time, start->time) / clocks_get_cpufreq() *
         clocks_units_scale;
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
  if (freq > 0) {
    clocks_cpufreq = freq;
  } else {
    clocks_estimate_cpufreq();
  }
  clocks_start_ticks = getticks();
}

/**
 * @brief Get the CPU frequency in Hz.
 *
 * @result the CPU frequency.
 */
unsigned long long clocks_get_cpufreq(void) {

  if (clocks_cpufreq > 0) return clocks_cpufreq;

  /* It not already set estimate it. */
  clocks_estimate_cpufreq();
  return clocks_cpufreq;
}

/**
 * @brief Estimate the CPU frequency in Hz.
 *
 * If already set return the CPU frequency, then estimate the CPU frequency.
 *
 * The technique is either use a clock timed nanosleep (this was the best
 * method on i7), to read the value from the cpuinfo_max_freq
 * file (probably a overestimate) or finally just use a value of 1 with
 * time units of ticks.
 */
static void clocks_estimate_cpufreq(void) {

#ifdef HAVE_CLOCK_GETTIME
  /* Try to time a nanosleep() in ticks. */
  struct clocks_time time1;
  struct clocks_time time2;

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

  clocks_cpufreq =
      (signed long long)(double)(toc - tic) * 1.0 / realsleep * 1000.0;
  clocks_units_index = 0;
  clocks_units_scale = 1000.0;
#endif

/* Look for the system value, if available. Tends to be too large. */
#ifdef __linux__
  if (clocks_cpufreq == 0) {
    FILE *file =
        fopen("/sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq", "r");
    if (file != NULL) {
      unsigned long long maxfreq;
      if (fscanf(file, "%llu", &maxfreq) == 1) {
        clocks_cpufreq = maxfreq * 1000;
        clocks_units_index = 0;
        clocks_units_scale = 1000.0;
      }
      fclose(file);
    }
  }
#endif

  /* If all fails just report ticks for a notional 2.6GHz machine. */
  if (clocks_cpufreq == 0) {
    unsigned long long maxfreq = 2600000;
    clocks_cpufreq = maxfreq * 1000;
    clocks_units_index = 1;
    clocks_units_scale = 1000.0;
  }
}

/**
 * @brief Return the difference between two ticks.
 *
 * Only an approximation as based on how well we have estimated the
 * rtc frequency. Should be good for machines that support constant_rtc
 * and clock_gettime().
 *
 * @param tic a number of ticks returned by the cycle.h getticks() function.
 * @param toc a number of ticks returned by the cycle.h getticks() function.
 *
 * @result the difference.
 */
double clocks_diff_ticks(ticks tic, ticks toc) {
  return clocks_from_ticks(tic - toc);
}

/**
 * @brief Convert a number of ticks into milli seconds, if possible.
 *
 * Only an approximation as based on how well we have estimated the
 * rtc frequency. Should be good for machines that support constant_rtc
 * and clock_gettime(), and reasonable for most Linux machines, otherwise
 * ticks will just be returned. See clocks_getunit() for the actual units.
 *
 * @param tics a number of ticks returned by the cycle.h getticks() function.
 *
 * @result the milli seconds, if possible.
 */
double clocks_from_ticks(ticks tics) {
  return ((double)tics / (double)clocks_get_cpufreq() * clocks_units_scale);
}

/**
 * @brief Convert a number of milli seconds into ticks, if possible.
 *
 * Only an approximation as based on how well we have estimated the
 * rtc frequency. Should be good for machines that support constant_rtc
 * and clock_gettime(), and reasonable for most Linux machines, otherwise
 * a guess will just be returned. See clocks_getunit() for the actual units.
 *
 * @param ms a number of "milliseconds" to convert to ticks
 *
 * @result the number of ticks, if possible.
 */
ticks clocks_to_ticks(double ms) {
  return (ticks)(ms * (double)clocks_get_cpufreq() / clocks_units_scale);
}

/**
 * @brief return the time units.
 *
 * Normally "ms" for milliseconds, but can be "ticks" when no conversion
 * factor is available.
 *
 * @result the current time units.
 */
const char *clocks_getunit(void) { return clocks_units[clocks_units_index]; }

/**
 * @brief returns the time since the start of the execution in seconds
 *
 * Need to call clocks_set_cpufreq() to mark the start of execution.
 *
 * The time is return in the format [sssss.s].
 *
 * @result the time since the start of the execution
 */
const char *clocks_get_timesincestart(void) {

  static char buffer[40];

  sprintf(buffer, "[%07.1f]",
          clocks_diff_ticks(getticks(), clocks_start_ticks) / 1000.0);

  return buffer;
}

/**
 * Returns the wall-clock time since the start of execution in hours.
 *
 * Need to call clocks_set_cpufreq() to mark the start of execution.
 *
 * @result the time since the start of the execution
 */
double clocks_get_hours_since_start(void) {
  return clocks_diff_ticks(getticks(), clocks_start_ticks) / (3600. * 1000.0);
}

/**
 * @brief return the cpu time used.
 *
 * Uses the times(2) function to access the user cpu times and returns the sum
 * of these for the process tree, i.e. current process plus "waited-for"
 * children. This may be pthread implementation specific as to what that
 * exactly means. Note we do not include the system time as that includes
 * spin times and we don't want to give credit for that.
 *
 * @result cpu time used in sysconf(_SC_CLK_TCK) ticks, usually 100/s not our
 *         usual ticks.
 */
double clocks_get_cputime_used(void) {

  struct tms tmstic;
  times(&tmstic);
  return (double)(tmstic.tms_utime + tmstic.tms_cutime);
}

/**
 * @brief Return an integer based on the current time.
 *
 * Normally this will be the remainder of the current number of nanoseconds
 * so not very dissimilar in the most significant figures unless the time
 * between calls is greater than INT_MAX nanoseconds. For faster calls use
 * fewer figures, if that matters.
 *
 * @result an integer.
 */
int clocks_random_seed(void) {
#ifdef HAVE_CLOCK_GETTIME
  struct timespec timespec;
  clock_gettime(CLOCK_REALTIME, &timespec);
  return (timespec.tv_nsec % INT_MAX);
#else
  return (getticks() % INT_MAX);
#endif
}
