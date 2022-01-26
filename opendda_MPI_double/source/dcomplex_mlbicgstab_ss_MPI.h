/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for dcomplex_mlbicgstab_ss_MPI

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "dcomplex_type.h"
#include "dcomplex_math_MPI.h"
#include "definitions.h"
#include "degree_radian_conversion.h"
#include "dftmatvec_MPI.h"
#include "print_details.h"
#include "sfmt_rng.h" /* SIMD oriented Fast Mersenne Twister RNG */

void dcomplex_mlbicgstab_ss_MPI(void);
