/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Set the initial guess for the iterative solver

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "set_initial_guess.h"

void set_initial_guess(void){

   long int index0i;

   if(iterative.initial_guess==0){ /* Set the initial guess to x=0 */
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=zero;
      }
   }
   else if(iterative.initial_guess==1){ /* Set the initial guess to x=1 */
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=one;
      }
   }
   else if(iterative.initial_guess==2){ /* Set the initial guess to x=b i.e. incident electric-field */
      memcpy(dipole_polarisation,incident_E_field,(size_t)target.Nd3*sizeof(ldcomplex));
   }
   else if(iterative.initial_guess==3){ /* Set the initial guess to the matrix diagonal i.e. 1/polarisability */
      memcpy(dipole_polarisation,interaction_matrix_diagonal,(size_t)target.Nd3*sizeof(ldcomplex));
   }
}
