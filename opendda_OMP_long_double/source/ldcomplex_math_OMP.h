/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for ldcomplex_math_OMP

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include "ldcomplex_type.h"
#include "definitions.h"

long double ldcomplex_abs(ldcomplex in);
long double ldcomplex_abs2(ldcomplex in);
ldcomplex ldcomplex_add(ldcomplex in0,ldcomplex in1);
ldcomplex ldcomplex_sub(ldcomplex in0,ldcomplex in1);
ldcomplex ldcomplex_mul(ldcomplex in0,ldcomplex in1);
ldcomplex ldcomplex_div(ldcomplex in0,ldcomplex in1);
ldcomplex ldcomplex_scale(ldcomplex in,long double scale);
ldcomplex ldcomplex_no_r_scale(ldcomplex in,long double scale);
long double ldcomplex_2vectornorm(ldcomplex *in);
ldcomplex ldcomplex_exp(ldcomplex in);
ldcomplex ldcomplex_no_r_exp(ldcomplex in);
ldcomplex ldcomplex_conj(ldcomplex in);
ldcomplex ldcomplex_sqrt(ldcomplex in);
