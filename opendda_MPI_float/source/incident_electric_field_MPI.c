/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Compute the electric field at the lattice sites

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "incident_electric_field_MPI.h"

void incident_electric_field_MPI(void){

   int i,j,k,vc;
   long int index0i,index0j,index0k,element;
   float wave_vector_TF[3],kr0,kr1,kr2;
   fcomplex current_polarisation_TF[3],eikr,ikr;

   if(polarisation_state==0){ /* Incident polarisation e0 */ 
      current_polarisation_TF[0]=incident_TF.polarisation[0]; /* e0x in the target frame */
      current_polarisation_TF[1]=incident_TF.polarisation[1]; /* e0y in the target frame */
      current_polarisation_TF[2]=incident_TF.polarisation[2]; /* e0z in the target frame */
   }
   else if(polarisation_state==1){ /* Orthonormal polarisation e1 */
      current_polarisation_TF[0]=incident_TF.orthonormal[0]; /* e1x in the target frame */
      current_polarisation_TF[1]=incident_TF.orthonormal[1]; /* e1y in the target frame */
      current_polarisation_TF[2]=incident_TF.orthonormal[2]; /* e1z in the target frame */
   }
   else{
      print_error("Incident electric field calculation error","Incident polarisation error","The current polarisation state MUST be 0 OR 1");
   }

   /* Wave vector in the target frame */
   wave_vector_TF[0]=wavenumber*incident_TF.n_inc[0]*target.dipole_spacing;
   wave_vector_TF[1]=wavenumber*incident_TF.n_inc[1]*target.dipole_spacing;
   wave_vector_TF[2]=wavenumber*incident_TF.n_inc[2]*target.dipole_spacing;

   /*   Note: Dot product k.r=kx*rx+ky*ry+kz*rz */
   for(k=1;k<=target.P;k++){
      if((k-1)%parallel.np==parallel.myrank){ /* Plane k is on proc (k%np) */
         parallel.plane=((int)((double)(k-1)/parallel.dnp)); /* Plane k is plane ((int)(k/np)) on that proc */
         index0k=(k-1)*target.M;
         element=target.populated[parallel.plane]; /* Start index for plane k for the stored non-zero lattice sites */
         kr2=wave_vector_TF[2]*((float)k); /* kz*rz */
         for(j=1;j<=target.J;j++){
            index0j=index0k+(j-1)*target.K;
            kr1=kr2+wave_vector_TF[1]*((float)j); /* ky*ry+kz*rz */
            for(i=1;i<=target.K;i++){
               index0i=index0j+(i-1);
               kr0=kr1+wave_vector_TF[0]*((float)i); /* k.r=kx*rx+ky*ry+kz*rz */
               ikr.dat[0]=0.0F;ikr.dat[1]=kr0; /* ikr */
               eikr.dat[0]=cosf(ikr.dat[1]);eikr.dat[1]=sinf(ikr.dat[1]); /* e^{ikr} */
               /* Complex exponential where real part of 'in' is zero i.e. in=a+ib=0.0+ib
               out=exp(in)=exp(a)*[cos(b)+isin(b)]=exp(0.0)*[cos(b)+isin(b)]=1.0*[cos(b)+isin(b)]=cos(b)+isin(b) */
               /* incident_E_field=eikr*current_polarisation_TF[vc] */
               if(target.occupied[(int)((double)index0i/32.0)]&(1<<(index0i%32))){ /* Lattice site is occupied */
                  for(vc=0;vc<3;vc++){
                     incident_E_field[element+vc*parallel.alloc_vector]=fcomplex_mul(eikr,current_polarisation_TF[vc]);
                  }
                  element++;
               }
            }
         }
      }
   }   
}
