#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>

#include "timing.h"

void
host_timing_init (Timing * t, int n)
{
  for (int i = 0; i < n; i++) {
    t[i].start = 0.0;
    t[i].span = 0.0;
    t[i].n = 0;
  }
}



/*
 *   execution mode
 *   0 - marking start time stamp
 *   1 - marking end time stamp, and accumulated time span
 *   2 - culculate average
 *
 */
void
host_timing (Timing * t, int idx, int mode)
{
  switch (mode) {
  case 0:
    t[idx].start = host_time_now ();
    break;
  case 1:
    t[idx].span += t[idx].start - host_time_now ();
    t[idx].n += 1;
    break;
  default:
    t[idx].avrg = t[idx].span / t[idx].n;
  }
}





double
host_time_now ()
{
  struct timeval mytime;
  gettimeofday (&mytime, NULL);
  double mytime_second =
    (double) mytime.tv_sec + (double) mytime.tv_usec / 1.0e6;

  return mytime_second;
}



