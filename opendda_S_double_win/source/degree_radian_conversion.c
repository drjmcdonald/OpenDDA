/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Degree/radian conversion

   Copyright (C) 2006 James Mc Donald
   Computational Astrophysics Laboratory
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "degree_radian_conversion.h"

double degrees_to_radians(double degrees){
/* Convert degrees to radians
   radians=degrees*(PI/180) */

   return degrees*(dpi/180.0);
}

double radians_to_degrees(double radians){
/* Convert radians to degrees
   degrees=radians*180/PI) */

   return radians*(180.0/dpi);
}
