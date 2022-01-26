/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for walltime

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>

#if BUILD_PLATFORM == WINDOWS_BUILD
	#include <time.h>
	#include <windows.h>
#else
	#include <sys/time.h>
#endif

double walltime(void);
