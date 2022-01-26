/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for dcomplex_math_OMP

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include "dcomplex_type.h"
#include "definitions.h"

double dcomplex_abs(dcomplex in);
double dcomplex_abs2(dcomplex in);
dcomplex dcomplex_add(dcomplex in0,dcomplex in1);
dcomplex dcomplex_sub(dcomplex in0,dcomplex in1);
dcomplex dcomplex_mul(dcomplex in0,dcomplex in1);
dcomplex dcomplex_div(dcomplex in0,dcomplex in1);
dcomplex dcomplex_scale(dcomplex in,double scale);
dcomplex dcomplex_no_r_scale(dcomplex in,double scale);
double dcomplex_2vectornorm(dcomplex *in);
dcomplex dcomplex_exp(dcomplex in);
dcomplex dcomplex_no_r_exp(dcomplex in);
dcomplex dcomplex_conj(dcomplex in);
dcomplex dcomplex_sqrt(dcomplex in);
