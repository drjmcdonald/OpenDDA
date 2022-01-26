/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Set K, J and P using the shape and the construction parameters

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "set_domain.h"

void set_domain(void){

   if(target.from_file==0){ /* Built-in target construction */   
      if((strcmp(target.shape,"ellipsoid")==0)||(strcmp(target.shape,"cuboid")==0)){
         target.K=target.construction_parameters[0]; /* K */
         target.J=target.construction_parameters[1]; /* J */
         target.P=target.construction_parameters[2]; /* P */
      }
   }
   else{ /* Custom target being read from file */
      target.K=target.construction_parameters[0]; /* K */
      target.J=target.construction_parameters[1]; /* J */
      target.P=target.construction_parameters[2]; /* P */
   }
}
