/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Jacobi-Point Preconditioning initialisation

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "preconditioner_OMP.h"

void Point_Jacobi_preconditioner_OMP(void){
/* Point-Jacobi Preconditioning i.e.
   Divide by the diagonal of the interation matrix

   This function simply inverts the diagonal so that
   the actual preconditioning is implemented via
   multiplication rather than division */

   long int index0i;

#pragma omp parallel for schedule(dynamic,target.M)
   for(index0i=0;index0i<target.Nd3;index0i++){
      point_jacobi[index0i]=fcomplex_div(one,interaction_matrix_diagonal[index0i]);
   }
}
