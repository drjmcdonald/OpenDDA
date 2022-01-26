/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Calculates the various cross sections

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "cross_sections_OMP.h"

void cross_sections_OMP(void){

   long int index0i;
   double k3,temp,temp1;
   dcomplex current_polarisation_TF[3],E2;

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
   E2=dcomplex_add(dcomplex_mul(current_polarisation_TF[0],dcomplex_conj(current_polarisation_TF[0])),
                  dcomplex_add(dcomplex_mul(current_polarisation_TF[1],dcomplex_conj(current_polarisation_TF[1])),
                  dcomplex_mul(current_polarisation_TF[2],dcomplex_conj(current_polarisation_TF[2]))));

/* C_{ext} Extinction cross section
   cross_section.extinction=
   [(4*pi*k)/(|E|^{2})]*sum_{j}IM[conj(incident_E_field_{j}*dipole_polarisation_{j})] */
   temp=0.0;
#pragma omp parallel for reduction(+:temp) schedule(dynamic,target.M)
   for(index0i=0;index0i<target.Nd3;index0i++){
      temp+=incident_E_field[index0i].dat[0]*dipole_polarisation[index0i].dat[1]-
            incident_E_field[index0i].dat[1]*dipole_polarisation[index0i].dat[0];
   }
   cross_section.extinction=(4.0*dpi*wavenumber*temp)/E2.dat[0];

/* C_{abs} Absorption cross section
   cross_section.absorption=
   [(4*pi*k)/(|E|^{2})]*sum_{j}{IM[dipole_polarisation_{j}*(1/polarisability_{j})*
   conj(dipole_polarisation_{j})]-(2/3)k^{3}|dipole_polarisation_{j}|^{2}} */
   k3=(2.0/3.0)*wavenumber*wavenumber*wavenumber; /* (2/3)k^{3} */
   temp=0.0;
#pragma omp parallel for private(temp1) reduction(+:temp) schedule(dynamic,target.M)
   for(index0i=0;index0i<target.Nd3;index0i++){
      temp1=dcomplex_abs2(dipole_polarisation[index0i]);
      temp+=(-interaction_matrix_diagonal[index0i].dat[1]*temp1-(k3*temp1));
   }
   cross_section.absorption=((4.0*dpi*wavenumber)/E2.dat[0])*temp;
   cross_section.scattering=cross_section.extinction-cross_section.absorption;
}
