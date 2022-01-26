/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Call the appropriate iterative scheme

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "iterative_solution.h"

void iterative_solution(void){

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Setup the iterative.preconditioner */
   if(iterative.precond==1){
      /* Point-Jacobi Preconditioning */
      Point_Jacobi_preconditioner_OMP();
   }

   if(timing.enabled){ /* Increment the number of times the iterative kernel has been executed */
      timing.iterative_kernel_count++;
   }
      
   if(strcmp(iterative.scheme,"bicg")==0){
      ldcomplex_bicg_OMP();
   }
   else if(strcmp(iterative.scheme,"bicg_sym")==0){
      ldcomplex_bicg_sym_OMP();
   }
   else if(strcmp(iterative.scheme,"bicgstab")==0){
      ldcomplex_bicgstab_OMP();
   }
   else if(strcmp(iterative.scheme,"cg")==0){
      ldcomplex_cg_OMP();
   }
   else if(strcmp(iterative.scheme,"cgs")==0){
      ldcomplex_cgs_OMP();
   }
   else if(strcmp(iterative.scheme,"mlbicgstab")==0){
      ldcomplex_mlbicgstab_OMP();
   }
   else if(strcmp(iterative.scheme,"mlbicgstab_orig")==0){
      ldcomplex_mlbicgstab_orig_OMP();
   }
   else if(strcmp(iterative.scheme,"mlbicgstab_ss")==0){
      ldcomplex_mlbicgstab_ss_OMP();
   }
   else if(strcmp(iterative.scheme,"qmr")==0){
      ldcomplex_qmr_OMP();
   }
   else if(strcmp(iterative.scheme,"qmr_sym")==0){
      ldcomplex_qmr_sym_OMP();
   }
   else if(strcmp(iterative.scheme,"rbicgstab")==0){
      ldcomplex_rbicgstab_OMP();
   }
   else if(strcmp(iterative.scheme,"tfqmr")==0){
      ldcomplex_tfqmr_OMP();
   }
}
