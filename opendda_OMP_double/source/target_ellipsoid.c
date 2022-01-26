/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Build a homogeneous & isotropic ellipsoidal target

   Equation of a general ellipsoid

   (x-x0)^2/a^2 + (y-y0)^2/b^2 + (k-k0)^2/c^2 = 1

   a,b,c are the semi axis lengths (radii)

   (x-x0)^2/(A/2)^2 + (y-y0)^2/(B/2)^2 + (k-k0)^2/(C^2)^2 = 1

   A,B,C are the major axis lengths (diameters)

   (x-x0)^2/A^2 + (y-y0)^2/B^2 + (k-k0)^2/C^2 = 0.25

   (x-x0)^2/K^2 + (y-y0)^2/J^2 + (k-k0)^2/P^2 = 0.25

   i0,j0,k0: Origin
   i,j,k,index0*: Loop and array index control variables
   ii,jj,kk: Spatial coordinate variables
             ii=(x-x0)
             jj=(y-y0)
             kk=(z-z0)
   Ksqd,
   Jsqd,
   Psqd: Ksqd=K^{2},Jsqd=J^{2},Psqd=P^{2}

   isqrd,
   jsqrd,
   ksqrd: isqrd=(x-x0)^{2}/K^{2}
          jsqrd=(y-y0)^{2}/J^{2}
          ksqrd=(z-z0)^{2}/P^{2}
   Nd: Nd is the number of lattice sites actually occupied by dipoles

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "target_ellipsoid.h"

void target_ellipsoid(void){

   int i,j,k;
   long int index0k,index0j,index0i,temp;
   double ii,jj,kk,i0,j0,k0,isqrd,jsqrd,ksqrd,Ksqd,Jsqd,Psqd;

   Ksqd=(double)(target.K*target.K);
   Jsqd=(double)(target.J*target.J);
   Psqd=(double)(target.P*target.P);

/* Origin */
   i0=((double)target.K-1.0)/2.0;
   j0=((double)target.J-1.0)/2.0;
   k0=((double)target.P-1.0)/2.0;

   temp=0.0;   
#pragma omp parallel for private(index0k,kk,ksqrd,j,index0j,jj,jsqrd,i,ii,isqrd,index0i) reduction(+:temp) schedule(dynamic)
   for(k=0;k<target.P;k++){
      index0k=k*target.M;
      kk=(double)k-k0;
      ksqrd=(kk*kk)/Psqd;
      if(ksqrd<0.25){
         for(j=0;j<target.J;j++){
            index0j=index0k+j*target.K;
            jj=(double)j-j0;
            jsqrd=(jj*jj)/Jsqd;
            if((jsqrd+ksqrd)<0.25){
               for(i=0;i<target.K;i++){
                  ii=(double)i-i0;
                  isqrd=(ii*ii)/Ksqd;
                  if((isqrd+jsqrd+ksqrd)<0.25){
                     index0i=index0j+i;
                     #pragma omp critical
                     target.occupied[(int)((double)index0i/32.0)]|=(1<<(index0i%32)); /* Set occupied flag */
                     /* Homogeneous & isotropic ellipsoid. Set composition data array target.composition[0..2]=0 already */
                     temp++; /* Increment Nd */
                  }
               }
            }
         }
      }
   }
   target.Nd=temp;
   if(target.output_xyz_vtk==1){
      /* target xyz coordinates are output to a vtk file for subsequent 3D visualisation */
      target_xyz_vtk();
   }
}
