/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for interpolation

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <float.h>
#include "fcomplex_type.h"
#include "definitions.h"
#include "print_details.h"

void interpolation(int degree,int size,float *xdata,fcomplex *ydata,int points,float *xinterp,fcomplex *yinterp,int part);
float power(float x,int n);
