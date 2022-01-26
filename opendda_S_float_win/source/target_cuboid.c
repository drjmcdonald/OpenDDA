/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Build a homogeneous & isotropic cuboidal target

   i,j,k,index0*: Loop and array index control variables

   Nd: Nd is the number of lattice sites actually occupied by dipoles

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "target_cuboid.h"

void target_cuboid(void){

   long int index0i;

   target.Nd=0; /* Initialise Nd */
   for(index0i=0;index0i<target.N;index0i++){
      target.occupied[(int)((double)index0i/32.0)]|=(1<<(index0i%32)); /* Set occupied flag */
      /* Homogeneous & isotropic cuboid. Set composition data array target.composition[0..2]=0 already */
      target.Nd++; /* Increment Nd */
   }
   if(target.output_xyz_vtk==1){
      /* target xyz coordinates are output to a vtk file for subsequent 3D visualisation */
      target_xyz_vtk();
   }
}
