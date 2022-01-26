/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Build the target

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "build_target.h"

void build_target(void){

   if(target.from_file==0){ /* Built-in target construction */
      if(strcmp(target.shape,"ellipsoid")==0){
         target_ellipsoid();
      }
      else if(strcmp(target.shape,"cuboid")==0){
         target_cuboid();
      }
   }
   else{
      if(parallel.myrank==0){ /* Restrict to master */   
         target_from_file();
      }
      broadcast_target_from_file();
   }
}
