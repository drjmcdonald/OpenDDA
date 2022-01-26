/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for rotation

   Copyright (C) 2006 James Mc Donald
   Computational Astrophysics Laboratory
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <math.h>
#include "ldcomplex_type.h"
#include "ldcomplex_math_MPI.h"
#include "definitions.h"
#include "degree_radian_conversion.h"

void scattering_vectors(long double phi,long double theta);
void rotation_euler(long double phi,long double costheta,long double psi);
