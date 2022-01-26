/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Calculates the various cross sections

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "cross_sections_MPI.h"

void cross_sections_MPI(void){

   long int index0i;
   float k3,temp,temp1,result;
   fcomplex current_polarisation_TF[3],E2;

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
      print_error("Cross section calculation error","Incident polarisation error","The current polarisation state MUST be 0 OR 1");
   }

   /* |E|^{2}=sum_{i}[current_polarisation_TF[i]*conj(current_polarisation_TF[i])]
      Should not be necessary in practice as the polarisations are normalised i.e.
      |E|^{2}=1 */
   E2=fcomplex_add(fcomplex_mul(current_polarisation_TF[0],fcomplex_conj(current_polarisation_TF[0])),
                  fcomplex_add(fcomplex_mul(current_polarisation_TF[1],fcomplex_conj(current_polarisation_TF[1])),
                  fcomplex_mul(current_polarisation_TF[2],fcomplex_conj(current_polarisation_TF[2]))));

/* C_{ext} Extinction cross section
   cross_section.extinction=
   [(4*pi*k)/(|E|^{2})]*sum_{j}IM[conj(incident_E_field_{j}*dipole_polarisation_{j})] */
   temp=0.0F;result=0.0F;
   for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
      temp+=(incident_E_field[index0i].dat[0]*dipole_polarisation[index0i].dat[1]-
            incident_E_field[index0i].dat[1]*dipole_polarisation[index0i].dat[0]);
   }
   MPI_Reduce(&temp,&result,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      cross_section.extinction=(4.0F*fpi*wavenumber*result)/E2.dat[0];
   }

/* C_{abs} Absorption cross section
   cross_section.absorption=
   [(4*pi*k)/(|E|^{2})]*sum_{j}{IM[dipole_polarisation_{j}*(1/polarisability_{j})*
   conj(dipole_polarisation_{j})]-(2/3)k^{3}|dipole_polarisation_{j}|^{2}} */
   k3=(2.0F/3.0F)*wavenumber*wavenumber*wavenumber; /* (2/3)k^{3} */
   temp=0.0F;result=0.0F;
   for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
      temp1=fcomplex_abs2(dipole_polarisation[index0i]);
      temp+=(-interaction_matrix_diagonal[index0i].dat[1]*temp1-(k3*temp1));
   }
   MPI_Reduce(&temp,&result,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      cross_section.absorption=((4.0F*fpi*wavenumber)/E2.dat[0])*result;
      cross_section.scattering=cross_section.extinction-cross_section.absorption;
   }
}
