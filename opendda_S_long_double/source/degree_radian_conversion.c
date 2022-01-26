/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Degree/radian conversion

   Copyright (C) 2006 James Mc Donald
   Computational Astrophysics Laboratory
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "degree_radian_conversion.h"

long double degrees_to_radians(long double degrees){
/* Convert degrees to radians
   radians=degrees*(PI/180) */

   return degrees*(ldpi/180.0L);
}

long double radians_to_degrees(long double radians){
/* Convert radians to degrees
   degrees=radians*180/PI) */

   return radians*(180.0L/ldpi);
}
