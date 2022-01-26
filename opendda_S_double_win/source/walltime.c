/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Timing function, microsecond accurate

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "walltime.h"

#ifndef __USE_BSD
	struct timezone
	{
		int tz_minuteswest; /* Minutes west of GMT.  */
		int tz_dsttime;  /* Nonzero if DST is ever in effect.  */
	};
#endif

#if BUILD_PLATFORM == WINDOWS_BUILD
double walltime(void){

   double time;
   LARGE_INTEGER time_val, time_freq;

   QueryPerformanceFrequency(&time_freq);
   QueryPerformanceCounter(&time_val);

   time=((double)time_val.LowPart/(double)time_freq.LowPart)*1e6;

   return(time);
}
#else
double walltime(void){

   double time;
   struct timeval tv;
   struct timezone tz;

   gettimeofday(&tv,&tz);

   time=(double)tv.tv_sec+(double)tv.tv_usec*1e-6;

   return(time);
}
#endif
