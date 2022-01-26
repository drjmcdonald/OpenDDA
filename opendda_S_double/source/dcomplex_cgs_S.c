/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   Conjugate-Gradients Squared: CGS

   CGS, a fast Lanczos-type solver for nonsymmetric linear systems, Sonneveld, P.,
   SIAM J. Sci. Stat. Comput., Vol. 10, No. 1, Pages 36-52, 1989.

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
   u[3],q[3]: Local vectors for iterative scheme
   rho,rhok,alpha,beta: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dcomplex_cgs_S.h"

void dcomplex_cgs_S(void){

   int vc,iteration=0;
   long int index0i,array0,array1,array2;
   dcomplex rho,rhok,alpha,beta;
   double ratio,norm0,norm1,terminate,start,finish;

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
            r[array2]=dcomplex_sub(incident_E_field[array2],zero_padded_vector[array1+index0i]);
         }
      }
   }
   else{ /* r=incident_E_field */
      memcpy(r,incident_E_field,(size_t)target.Nd3*sizeof(dcomplex));
   }

   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* r=M^{-1}r */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
         for(index0i=0;index0i<target.Nd3;index0i++){
            r[index0i]=dcomplex_mul(r[index0i],point_jacobi[index0i]);
         }
      }
   }

   norm0=dcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0; 
   };

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. rtilde=r */
   memcpy(rtilde,r,(size_t)target.Nd3*sizeof(dcomplex));

   print_section_title(output.iter,"Iterative convergence details");
   fprintf(output.iter,"\n  Wavelength [micromemters]: %.*g\n",DBLP,wavelength.value);
   fprintf(output.iter,"  Effective radius [micromemters]: %.*g\n",DBLP,radius.value);
   fprintf(output.iter,"  Size parameter: %.*g\n",DBLP,size_parameter);
   fprintf(output.iter,"  Dipole spacing [micromemters]: %.*g\n",DBLP,target.dipole_spacing);
   fprintf(output.iter,"  Dipoles per wavelength: %.*g\n",DBLP,target.dipoles_per_wavelength);
   fprintf(output.iter,"  Euler phi: %.*g\n",DBLP,euler_phi.value);
   fprintf(output.iter,"  Euler theta: %.*g\n",DBLP,radians_to_degrees(acos(euler_theta.value)));
   fprintf(output.iter,"  Euler psi: %.*g\n",DBLP,euler_psi.value);
   fprintf(output.iter,"  Polarisation state: %d\n\n",polarisation_state);

   terminate=iterative.tolerance;

   while(iteration<iterative.maximum){ /* 17. Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. rhok=conj(rtilde^{T})r */
      rhok=zero;
      for(index0i=0;index0i<target.Nd3;index0i++){
         rhok=dcomplex_add(rhok,dcomplex_mul(dcomplex_conj(rtilde[index0i]),r[index0i]));
      }
      if(iteration==0){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. u=r */
         memcpy(u,r,(size_t)target.Nd3*sizeof(dcomplex));
   
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. p=u */
         memcpy(p,u,(size_t)target.Nd3*sizeof(dcomplex));
      }
      else{
         /* Check for iterative.breakdown */
         if((dcomplex_abs(rho))<iterative.breakdown){
            print_error("Iterative scheme error","CGS() iterative.breakdown","|rho|<iterative.breakdown"); 
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. beta=rhok/rho */
         beta=dcomplex_div(rhok,rho);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. u=r+beta*q */
         for(index0i=0;index0i<target.Nd3;index0i++){
            u[index0i]=dcomplex_add(r[index0i],dcomplex_mul(beta,q[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. p=u+beta(q+beta*p) */
         for(index0i=0;index0i<target.Nd3;index0i++){
            p[index0i]=dcomplex_add(u[index0i],dcomplex_mul(beta,dcomplex_add(q[index0i],dcomplex_mul(beta,p[index0i]))));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. vhat=interaction_matrix*p */
/* vhat is stored in zero_padded_vector */
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
                  zero_padded_vector[array2]=dcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. alpha=rhok/(conj(rtilde^{T})vhat) */
/* vhat is stored in zero_padded_vector */
      alpha=zero;
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            alpha=dcomplex_add(alpha,dcomplex_mul(dcomplex_conj(rtilde[array0+index0i]),zero_padded_vector[array1+index0i]));
         }
      }
      /* Check for iterative.breakdown */
      if((dcomplex_abs(alpha))<iterative.breakdown){
         print_error("Iterative scheme error","CGS() iterative.breakdown","|alpha|<iterative.breakdown");
      }
      alpha=dcomplex_div(rhok,alpha);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. q=u-alpha*vhat */
/* vhat is stored in zero_padded_vector */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            q[array2]=dcomplex_sub(u[array2],dcomplex_mul(alpha,zero_padded_vector[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. u=u+q */
      for(index0i=0;index0i<target.Nd3;index0i++){
         u[index0i]=dcomplex_add(u[index0i],q[index0i]);
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. dipole_polarisation=dipole_polarisation+alpha*u */
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=dcomplex_add(dipole_polarisation[index0i],dcomplex_mul(alpha,u[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. qhat=interaction_matrix*u */
/* qhat is stored in zero_padded_vector */
      dftmatvec_S(u);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*u=M^{-1}interaction_matrix*u */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=vc*target.Nv;
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=dcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. r=r-alpha*qhat */
/* qhat is stored in zero_padded_vector */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            r[array2]=dcomplex_sub(r[array2],dcomplex_mul(alpha,zero_padded_vector[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. Check stopping criterion */
      norm1=dcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      fprintf(output.iter,"  At iteration %d the convergence ratio is %.*g and the residual norm is %.*g\n",iteration,DBLP,ratio,DBLP,norm1);
      if(ratio<(terminate+DBL_EPSILON)){
         break;
      }
      rho=rhok;
      iteration++;
   }
   if(timing.enabled){ /* Finalise timing */
      finish=walltime();
      ratio=finish-start;
      fprintf(output.iter,"\n  Time for convergence: %.*gs\n\n",DBLP,ratio);
      timing.iterative_solution+=ratio;
   }
}
