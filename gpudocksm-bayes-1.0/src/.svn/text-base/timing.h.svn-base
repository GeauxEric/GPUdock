#ifndef TIMING_H
#define TIMING_H


typedef struct
{
  double start; // start time stamp
  double span; // accumulated time span
  double avrg; // time span average
  int n; // number of records
} Timing;

void
host_timing_init (Timing * t, int n);

void
host_timing (Timing * t, int idx, int mode);

double
host_time_now ();



#endif // TIMINC_H


