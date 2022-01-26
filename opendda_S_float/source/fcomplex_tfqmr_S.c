/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   Transpose-free quasi-minimal residual: TFQMR

   A transpose-free quasi-minimal residual algorithm for
   non-Hermitian linear systems, Freund, R., SIAM J. Sci. Comput.,
   Vol. 14, No. 2, Pages 470-482, 1993. Algorithm 5.1, Pages 477-478.

   Adapted from:
   PIM 2.2: The parallel iterative methods package for systems of linear
   equations. User's Guide (Fortran 77 version)}, da Cunha, R.D. and
   Hopkins, T., 1997.

   m,vc,index0i,array*: Loop and array index control variables
   iteration: Iteration number
   iterative.precond: 0: None
                      1: Point-Jacobi Preconditioning
   r[3],rtilde[3],omega[3],
   v[3],d[3],ya[3],yb[3]: Local vectors for iterative scheme
   eta,eta0,rho,rhok,sigma,
   alpha,beta,temp,tau,theta,
   tau0,theta0,c,kappa: Local variables for iterative scheme
   error,ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "fcomplex_tfqmr_S.h"

void fcomplex_tfqmr_S(void){

   int m,vc,iteration=0,exitflag=0;
   long int index0i,array0,array1,array2;
   fcomplex eta={{0.0F,0.0F}},eta0,rho,rhok,sigma,alpha,beta,temp,*ptr;
   float ratio,norm0,norm1,tau,theta=0.0F,tau0,theta0,c,kappa,error;
   float terminate;
   double start,finish;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=walltime();}

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 1. r=incident_E_field-interaction_matrix(dipole_polarisation) */
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
   ratio=1.0F;
   error=iterative.tolerance*norm0;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. omega=yb=r */
   memcpy(omega,r,(size_t)target.Nd3*sizeof(fcomplex));
   memcpy(yb,r,(size_t)target.Nd3*sizeof(fcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. v=interaction_matrix*yb */
   dftmatvec_S(yb);
   for(vc=0;vc<3;vc++){
      array0=vc*target.Nd;
      array1=vc*target.Nv;
      for(index0i=0;index0i<target.Nd;index0i++){
         v[array0+index0i]=zero_padded_vector[array1+index0i];
      }
   }

   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* v=M^{-1}v */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
         for(index0i=0;index0i<target.Nd3;index0i++){
            v[index0i]=fcomplex_mul(v[index0i],point_jacobi[index0i]);
         }
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. d=0 Set to zero using calloc in memory allocation */ 

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. tau=||r|| */
   tau=fcomplex_2vectornorm(r);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. theta=0.0,eta=0.0+0.0i Set at variable declaration stage */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. rtilde=r */
   memcpy(rtilde,r,(size_t)target.Nd3*sizeof(fcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. rhok=rtilde^{T}r */
   rhok=zero;
   for(index0i=0;index0i<target.Nd3;index0i++){
      rhok=fcomplex_add(rhok,fcomplex_mul(rtilde[index0i],r[index0i]));
   }

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

   while(iteration<iterative.maximum){ /* Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. sigma=rtilde^{T}v */
      sigma=zero;
      for(index0i=0;index0i<target.Nd3;index0i++){
         sigma=fcomplex_add(sigma,fcomplex_mul(rtilde[index0i],v[index0i]));
      }
      /* Check for iterative.breakdown */
      if(fcomplex_abs(sigma)<iterative.breakdown){
         print_error("Iterative scheme error","TFQMR() iterative.breakdown","|sigma|<iterative.breakdown");
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. alpha=rhok/sigma */
      alpha=fcomplex_div(rhok,sigma);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. ya=yb,yb=ya-alpha*v */
      for(index0i=0;index0i<target.Nd3;index0i++){
         ya[index0i]=fcomplex_sub(yb[index0i],fcomplex_mul(alpha,v[index0i]));
      }
      /* ya=yb,yb=ya-alpha*v */
      ptr=ya;ya=yb;yb=ptr;

      for(m=2*iteration-1;m<=2*iteration;m++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. omega=omega-alpha*interaction_matrix*ya */
         dftmatvec_S(ya);
         if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*ya=M^{-1}interaction_matrix*ya */
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
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=vc*target.Nv;
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=array0+index0i;
               omega[array2]=fcomplex_sub(omega[array2],fcomplex_mul(alpha,zero_padded_vector[array1+index0i]));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. theta=(||omega||)/(tau) */
         theta0=theta;
         tau0=tau;
         /* Check for iterative.breakdown */
         if(tau0<iterative.breakdown){
            print_error("Iterative scheme error","TFQMR() iterative.breakdown","tau0<iterative.breakdown");
         }
         theta=fcomplex_2vectornorm(omega)/tau0;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. c=(1)/(sqrt(1+theta^{2})) */
         c=1.0F/(sqrtf(1.0F+theta*theta));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. tau=tau0*theta*c */
         tau=tau0*theta*c;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. eta=c^{2}*alpha */
         eta0=eta;
         eta=fcomplex_scale(alpha,c*c);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. d=ya+(theta0^{2}*eta0/alpha)*d */
         temp=fcomplex_scale(fcomplex_div(eta0,alpha),theta0*theta0);
         for(index0i=0;index0i<target.Nd3;index0i++){
            d[index0i]=fcomplex_add(ya[index0i],fcomplex_mul(temp,d[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. dipole_polarisation=dipole_polarisation+eta*d */
         for(index0i=0;index0i<target.Nd3;index0i++){
            dipole_polarisation[index0i]=fcomplex_add(dipole_polarisation[index0i],fcomplex_mul(eta,d[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. kappa=tau*sqrt(m+1) */
         kappa=tau*sqrtf((float)m+1.0F);
         
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 20. Check stopping criterion */
         if(kappa<error){
/* r=incident_E_field-interaction_matrix*dipole_polarisation)*/
            dftmatvec_S(dipole_polarisation);
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=vc*target.Nv;
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array0+index0i;
                  r[array2]=fcomplex_sub(incident_E_field[array2],zero_padded_vector[array1+index0i]);
               }
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

            norm1=fcomplex_2vectornorm(r);
            ratio=norm1/norm0;
            fprintf(output.iter,"  At iteration %d the convergence ratio is %.*g and the residual norm is %.*g\n",iteration,FLTP,ratio,FLTP,norm1);
            if(ratio<(terminate+FLT_EPSILON)){
               exitflag=1;
               break;
            }
         }

/* ya=yb */
         if(m==(2*iteration-1)){
            ptr=ya;ya=yb;yb=ptr;
         }
      } /* End m for loop */
      if(exitflag){
         break;
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 21. rhok=rtilde^{T}omega */
      rho=rhok;
      rhok=zero;
      for(index0i=0;index0i<target.Nd3;index0i++){
         rhok=fcomplex_add(rhok,fcomplex_mul(rtilde[index0i],omega[index0i]));
      }
      /* Check for iterative.breakdown */
      if(fcomplex_abs(rho)<iterative.breakdown){
         print_error("Iterative scheme error","TFQMR() iterative.breakdown","|rho|<iterative.breakdown"); 
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 22. beta=rhok/rho   */
      beta=fcomplex_div(rhok,rho);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 23. yb=omega+beta*ya */
      for(index0i=0;index0i<target.Nd3;index0i++){
         yb[index0i]=fcomplex_add(omega[index0i],fcomplex_mul(beta,ya[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 24. v=interaction_matrix*yb+beta*(interaction_matrix*ya+beta*v) */
      dftmatvec_S(ya);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*ya=M^{-1}interaction_matrix*ya */
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
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array1+index0i;
            zero_padded_vector[array2]=fcomplex_add(zero_padded_vector[array2],fcomplex_mul(beta,v[array0+index0i]));
         }
      }
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            v[array0+index0i]=fcomplex_mul(beta,zero_padded_vector[array1+index0i]);
         }
      }
      dftmatvec_S(yb);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*yb=M^{-1}interaction_matrix*yb */
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
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            v[array2]=fcomplex_add(v[array2],zero_padded_vector[array1+index0i]);
         }
      }
      iteration++;
   } /* End while loop */
   if(timing.enabled){ /* Finalise timing */
      finish=walltime();
      fprintf(output.iter,"\n  Time for convergence: %.*gs\n\n",DBLP,finish-start);
      timing.iterative_solution+=finish-start;
   }
}
