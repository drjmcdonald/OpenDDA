/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for fcomplex_mlbicgstab_ss_OMP

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "fcomplex_type.h"
#include "fcomplex_math_OMP.h"
#include "definitions.h"
#include "degree_radian_conversion.h"
#include "dftmatvec_OMP.h"
#include "print_details.h"
#include "sfmt_rng.h" /* SIMD oriented Fast Mersenne Twister RNG */
#include "walltime.h"

void fcomplex_mlbicgstab_ss_OMP(void);
