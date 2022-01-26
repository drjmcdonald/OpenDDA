/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for rotation

   Copyright (C) 2006 James Mc Donald
   Computational Astrophysics Laboratory
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <math.h>
#include "fcomplex_type.h"
#include "fcomplex_math_S.h"
#include "definitions.h"
#include "degree_radian_conversion.h"

void scattering_vectors(float phi,float theta);
void rotation_euler(float phi,float costheta,float psi);
