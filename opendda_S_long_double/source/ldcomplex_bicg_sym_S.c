/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   BiConjugate-Gradients for symmetric systems: BiCG_sym

   Conjugate gradient-type methods for linear systems with
   complex symmetric coefficient matrices, Freund, R. W.,
   SIAM Journal on Scientific and Statistical Computing,
   Special issue on iterative methods in numerical linear
   algebra, Vol. 13, No. 1, Pages 425-448, 1992.

   vc,index0i,array*: Loop and array index control variables
   iteration: Iteration number
   iterative.precond: 0: None
                      1: Point-Jacobi Preconditioning
   r[3],p[3]: Local vectors for iterative scheme
   delta,rho,rhok,
   alpha,beta: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "ldcomplex_bicg_sym_S.h"

void ldcomplex_bicg_sym_S(void){

   int vc,iteration=0;
   long int index0i,array0,array1,array2;
   ldcomplex rho,rhok,alpha,beta,delta;
   long double ratio,norm0,norm1,terminate;
   double start,finish;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=walltime();}

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 1. r=incident_E_field-interaction_matrix*dipole_polarisation */
   if(iterative.initial_guess>0){
      dftmatvec_S(dipole_polarisation);
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            r[array2]=ldcomplex_sub(incident_E_field[array2],zero_padded_vector[array1+index0i]);
         }
      }
   }
   else{ /* r=incident_E_field */
      memcpy(r,incident_E_field,(size_t)target.Nd3*sizeof(ldcomplex));
   }

   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* r=M^{-1}r */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
         for(index0i=0;index0i<target.Nd3;index0i++){
            r[index0i]=ldcomplex_mul(r[index0i],point_jacobi[index0i]);
         }
      }
   }

   norm0=ldcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0L; 
   };

   print_section_title(output.iter,"Iterative convergence details");
   fprintf(output.iter,"\n  Wavelength [micromemters]: %.*Lg\n",LDBLP,wavelength.value);
   fprintf(output.iter,"  Effective radius [micromemters]: %.*Lg\n",LDBLP,radius.value);
   fprintf(output.iter,"  Size parameter: %.*Lg\n",LDBLP,size_parameter);
   fprintf(output.iter,"  Dipole spacing [micromemters]: %.*Lg\n",LDBLP,target.dipole_spacing);
   fprintf(output.iter,"  Dipoles per wavelength: %.*Lg\n",LDBLP,target.dipoles_per_wavelength);
   fprintf(output.iter,"  Euler phi: %.*Lg\n",LDBLP,euler_phi.value);
   fprintf(output.iter,"  Euler theta: %.*Lg\n",LDBLP,radians_to_degrees(acosl(euler_theta.value)));
   fprintf(output.iter,"  Euler psi: %.*Lg\n",LDBLP,euler_psi.value);
   fprintf(output.iter,"  Polarisation state: %d\n\n",polarisation_state);

   terminate=iterative.tolerance;

   while(iteration<iterative.maximum){ /* 11. Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. rhok=r^{T}r */
      rhok=zero;
      for(index0i=0;index0i<target.Nd3;index0i++){
         rhok=ldcomplex_add(rhok,ldcomplex_mul(r[index0i],r[index0i]));
      }

      if(iteration==0){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. p=r */
         memcpy(p,r,(size_t)target.Nd3*sizeof(ldcomplex));
      }
      else{
         /* Check for iterative.breakdown */
         if((ldcomplex_abs(rho))<iterative.breakdown){
            print_error("Iterative scheme error","BiCG_sym() iterative.breakdown","|rho|<iterative.breakdown"); 
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. beta=rhok/rho */
         beta=ldcomplex_div(rhok,rho);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. p=r+beta*p */
         for(index0i=0;index0i<target.Nd3;index0i++){
            p[index0i]=ldcomplex_add(r[index0i],ldcomplex_mul(beta,p[index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. interaction_matrix*p */
/* interaction_matrix*p is stored in zero_padded_vector */
      dftmatvec_S(p);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*p=M^{-1}interaction_matrix*p */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=vc*target.Nv;
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=ldcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. alpha=p^{T}interaction_matrix*p */
/* interaction_matrix*p is stored in zero_padded_vector */
      alpha=zero;
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            alpha=ldcomplex_add(alpha,ldcomplex_mul(p[array0+index0i],zero_padded_vector[array1+index0i]));
         }
      }
      /* Check for iterative.breakdown */
      if((ldcomplex_abs(alpha))<iterative.breakdown){
         print_error("Iterative scheme error","BiCG_sym() iterative.breakdown","|alpha|<iterative.breakdown"); 
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. delta=rhok/alpha */
      delta=ldcomplex_div(rhok,alpha);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. dipole_polarisation=dipole_polarisation+delta*p */
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=ldcomplex_add(dipole_polarisation[index0i],ldcomplex_mul(delta,p[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. r=r-delta*interaction_matrix*p */
/* interaction_matrix*p is stored in zero_padded_vector */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            r[array2]=ldcomplex_sub(r[array2],ldcomplex_mul(delta,zero_padded_vector[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. Check stopping criterion */
      norm1=ldcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      fprintf(output.iter,"  At iteration %d the convergence ratio is %.*Lg and the residual norm is %.*Lg\n",iteration,LDBLP,ratio,LDBLP,norm1);
      if(ratio<(terminate+LDBL_EPSILON)){
         break;
      }
      rho=rhok;
      iteration++;
   }
   if(timing.enabled){ /* Finalise timing */
      finish=walltime();
      fprintf(output.iter,"\n  Time for convergence: %.*gs\n\n",DBLP,finish-start);
      timing.iterative_solution+=finish-start;
   }
}
