/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for dcomplex_qmr_MPI

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <mpi.h>
#include <math.h>
#include "dcomplex_type.h"
#include "dcomplex_math_MPI.h"
#include "definitions.h"
#include "degree_radian_conversion.h"
#include "dftmatvec_MPI.h"
#include "print_details.h"

void dcomplex_qmr_MPI(void);
