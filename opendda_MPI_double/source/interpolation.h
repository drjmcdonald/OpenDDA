/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for interpolation

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <float.h>
#include "dcomplex_type.h"
#include "definitions.h"
#include "print_details.h"

void interpolation(int degree,int size,double *xdata,dcomplex *ydata,int points,double *xinterp,dcomplex *yinterp,int part);
double power(double x,int n);
