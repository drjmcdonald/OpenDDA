/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for fcomplex_math_S

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "fcomplex_type.h"
#include "definitions.h"

float fcomplex_abs(fcomplex in);
float fcomplex_abs2(fcomplex in);
fcomplex fcomplex_add(fcomplex in0,fcomplex in1);
fcomplex fcomplex_sub(fcomplex in0,fcomplex in1);
fcomplex fcomplex_mul(fcomplex in0,fcomplex in1);
fcomplex fcomplex_div(fcomplex in0,fcomplex in1);
fcomplex fcomplex_scale(fcomplex in,float scale);
fcomplex fcomplex_no_r_scale(fcomplex in,float scale);
float fcomplex_2vectornorm(fcomplex *in);
fcomplex fcomplex_exp(fcomplex in);
fcomplex fcomplex_no_r_exp(fcomplex in);
fcomplex fcomplex_conj(fcomplex in);
fcomplex fcomplex_sqrt(fcomplex in);
