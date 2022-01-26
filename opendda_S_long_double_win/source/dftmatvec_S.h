/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for dftmatvec_S

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <string.h>
#include <fftw3.h>
#include "ldcomplex_type.h"
#include "ldcomplex_math_S.h"
#include "definitions.h"
#include "walltime.h"

void dftmatvec_S(ldcomplex *vector);
