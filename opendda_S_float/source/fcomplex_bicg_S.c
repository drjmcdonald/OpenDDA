/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   BiConjugate-Gradients: BiCG

   Conjugate gradient methods for indefinite systems, Fletcher, R.,
   Proceedings of the Dundee Biennial Conference on Numerical Analysis,
   Ed. Watson, G.~A., Pages 73-89, 1975.

   Adapted from:
   Templates for the Solution of Linear Systems: Building Blocks for
   Iterative Methods, Barrett, R. and Berry, M. and Chan, T.F. and Demmel,
   J. and Donato, J. and Dongarra, J. and Eijkhout, V. and Pozo, R. and
   Romine, C. and Van der Vorst, H, SIAM, Philadelphia, 2nd Ed., 1994. 
   http://www.netlib.org/linalg/html_templates/Templates.html

   vc,index0i,array*: Loop and array index control variables
   iteration: Iteration number
   iterative.precond: 0: None
                      1: Point-Jacobi Preconditioning
   r[3],rtilde[3],p[3],
   ptilde[3],q[3]: Local vectors for iterative scheme
   rho,rhok,alpha,beta: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "fcomplex_bicg_S.h"

void fcomplex_bicg_S(void){

   int vc,iteration=0;
   long int index0i,array0,array1,array2;
   fcomplex rho,rhok,alpha,beta;
   float ratio,norm0,norm1,terminate;
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
            r[array2]=fcomplex_sub(incident_E_field[array2],zero_padded_vector[array1+index0i]);
         }
      }
   }
   else{ /* r=incident_E_field */
      memcpy(r,incident_E_field,(size_t)target.Nd3*sizeof(fcomplex));
   }      

   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* r=M^{-1}r */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
         for(index0i=0;index0i<target.Nd3;index0i++){
            r[index0i]=fcomplex_mul(r[index0i],point_jacobi[index0i]);
         }
      }
   }

   norm0=fcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0F; 
   };

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. rtilde=r */
   memcpy(rtilde,r,(size_t)target.Nd3*sizeof(fcomplex));

   print_section_title(output.iter,"Iterative convergence details");
   fprintf(output.iter,"\n  Wavelength [micromemters]: %.*g\n",FLTP,wavelength.value);
   fprintf(output.iter,"  Effective radius [micromemters]: %.*g\n",FLTP,radius.value);
   fprintf(output.iter,"  Size parameter: %.*g\n",FLTP,size_parameter);
   fprintf(output.iter,"  Dipole spacing [micromemters]: %.*g\n",FLTP,target.dipole_spacing);
   fprintf(output.iter,"  Dipoles per wavelength: %.*g\n",FLTP,target.dipoles_per_wavelength);
   fprintf(output.iter,"  Euler phi: %.*g\n",FLTP,euler_phi.value);
   fprintf(output.iter,"  Euler theta: %.*g\n",FLTP,radians_to_degrees(acosf(euler_theta.value)));
   fprintf(output.iter,"  Euler psi: %.*g\n",FLTP,euler_psi.value);
   fprintf(output.iter,"  Polarisation state: %d\n\n",polarisation_state);

   terminate=iterative.tolerance;

   while(iteration<iterative.maximum){ /* 17. Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. rhok=r^{T}rtilde */
      rhok=zero;
      for(index0i=0;index0i<target.Nd3;index0i++){
         rhok=fcomplex_add(rhok,fcomplex_mul(r[index0i],rtilde[index0i]));
      }

      if(iteration==0){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. p=r */
         memcpy(p,r,(size_t)target.Nd3*sizeof(fcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. ptilde=rtilde */
         memcpy(ptilde,rtilde,(size_t)target.Nd3*sizeof(fcomplex));
      }
      else{
         /* Check for iterative.breakdown */
         if((fcomplex_abs(rho))<iterative.breakdown){
            print_error("Iterative scheme error","BiCG() iterative.breakdown","|rho|<iterative.breakdown"); 
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. beta=rhok/rho */
         beta=fcomplex_div(rhok,rho);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. p=r+beta*p */
         for(index0i=0;index0i<target.Nd3;index0i++){
            p[index0i]=fcomplex_add(r[index0i],fcomplex_mul(beta,p[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. ptilde=rtilde+beta*ptilde */
         for(index0i=0;index0i<target.Nd3;index0i++){
            ptilde[index0i]=fcomplex_add(rtilde[index0i],fcomplex_mul(beta,ptilde[index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. q=interaction_matrix*p */
      dftmatvec_S(p);
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            q[array0+index0i]=zero_padded_vector[array1+index0i];
         }
      }

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* q=M^{-1}q */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(index0i=0;index0i<target.Nd3;index0i++){
               q[index0i]=fcomplex_mul(q[index0i],point_jacobi[index0i]);
            }
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. qtilde=interaction_matrix^{T}ptilde */
/* qtilde is stored in zero_padded_vector */
      dftmatvec_S(ptilde);

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix^{T}ptilde=M^{-1}interaction_matrix^{T}ptilde */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=vc*target.Nv;
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=fcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. alpha=rhok/(ptilde^{T}q) */
      alpha=zero;
      for(index0i=0;index0i<target.Nd3;index0i++){
         alpha=fcomplex_add(alpha,fcomplex_mul(ptilde[index0i],q[index0i]));
      }
      /* Check for iterative.breakdown */
      if((fcomplex_abs(alpha))<iterative.breakdown){
         print_error("Iterative scheme error","BiCG() iterative.breakdown","|alpha|<iterative.breakdown"); 
      }
      alpha=fcomplex_div(rhok,alpha);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. dipole_polarisation=dipole_polarisation+alpha*p */
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=fcomplex_add(dipole_polarisation[index0i],fcomplex_mul(alpha,p[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. r=r-alpha*q */
      for(index0i=0;index0i<target.Nd3;index0i++){
         r[index0i]=fcomplex_sub(r[index0i],fcomplex_mul(alpha,q[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. rtilde=rtilde-alpha*qtilde */
/* qtilde is stored in zero_padded_vector */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            rtilde[array2]=fcomplex_sub(rtilde[array2],fcomplex_mul(alpha,zero_padded_vector[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. Check stopping criterion */
      norm1=fcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      fprintf(output.iter,"  At iteration %d the convergence ratio is %.*g and the residual norm is %.*g\n",iteration,FLTP,ratio,FLTP,norm1);
      if(ratio<(terminate+FLT_EPSILON)){
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
