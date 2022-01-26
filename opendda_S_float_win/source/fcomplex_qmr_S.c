/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   Quasi-minimal residual with coupled two-term recurrences: QMR

   A quasi-minimal residual method for non-Hermitian
   linear systems, Freund, R. and Nachtigal, N.,
   Numer. Math., Vol. 60, Pages 315-339, 1991.

   An implementation of the QMR method based on coupled
   two-term recurrences, Freund, R. and Nachtigal, N.,
   SIAM J. Sci. Comput., Vol. 15, No. 2, Pages 313-337, 1992.

   A parallel version of the unsymmetric Lanczos algorithm
   and its application to QMR, B\"ucker, H.M. and Sauren, M.,
   Central Institute for Applied Mathematics, Research Centre
   J\"ulich, Germany, Tech. Report No. = KFA-ZAM-IB-9606, 1996.

   Adapted from:
   PIM 2.2: The parallel iterative methods package for systems of linear
   equations. User's Guide (Fortran 77 version)}, da Cunha, R.D. and
   Hopkins, T., 1997.

   vc,index0i,array*: Loop and array index control variables
   iteration: Iteration number
   iterative.precond: 0: None
                      1: Point-Jacobi Preconditioning
   r[3],d[3],s[3],p[3],q[3],
   omegatilde[3],vtilde[3]: Local vectors for iterative scheme
   tempdc,kappa,kappa0,rho,
   rho0,epsilon,mu,tau,tau0,
   gamma0,gamma1,xi,xi0,
   lambda,theta,denom: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "fcomplex_qmr_S.h"

void fcomplex_qmr_S(void){

   int vc,iteration=0;
   long int index0i,array0,array1,array2;
   fcomplex tempdc={{0.0F,0.0F}},kappa={{-1.0F,0.0F}},kappa0,rho,rho0;
   fcomplex epsilon,mu={{0.0F,0.0F}},tau,tau0;
   float ratio,norm0,norm1,gamma0,gamma1,xi,xi0,lambda=1.0F;
   float theta=-1.0F,denom,tempd,terminate;
   double start,finish;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=walltime();}

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 1. lambda=1,kappa=theta=-1 Set at variable declaration stage */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. r=incident_E_field-interaction_matrix(dipole_polarisation),omegatilde=vtilde=r */
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

   /* omegatilde=vtilde=r */
   memcpy(omegatilde,r,(size_t)target.Nd3*sizeof(fcomplex));
   memcpy(vtilde,r,(size_t)target.Nd3*sizeof(fcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. p=q=d=s=0 Set to zero using calloc in memory allocation */ 

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. gamma0=||vtilde||
      xi=||omegatilde||
      rho=[omegatilde]^{T}vtilde
      epsilon=(interaction_matrix^{T}omegatilde)^{T}vtilde
      mu=0 set at variable declaration stage */
   epsilon=zero;
/*   gamma0=fcomplex_2vectornorm(vtilde);
   xi=fcomplex_2vectornorm(omegatilde); */
   gamma0=norm0;xi=norm0;
   rho=zero;
   for(index0i=0;index0i<target.Nd3;index0i++){
      rho=fcomplex_add(rho,fcomplex_mul(omegatilde[index0i],vtilde[index0i]));
   }
   dftmatvec_S(omegatilde);
   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*omegatilde=M^{-1}interaction_matrix*omegatilde */
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
         epsilon=fcomplex_add(epsilon,fcomplex_mul(zero_padded_vector[array1+index0i],vtilde[array0+index0i]));
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. tau=epsilon/rho    */
   /* Check for iterative.breakdown */
   if(fcomplex_abs(rho)<iterative.breakdown){
      print_error("Iterative scheme error","QMR() iterative.breakdown","|rho|<iterative.breakdown"); 
   }
   tau=fcomplex_div(epsilon,rho);

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
/* 6. p=vtilde/gamma0-mu*p */
      /* Check for iterative.breakdown */
      if(gamma0<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","gamma0<iterative.breakdown"); 
      }
      tempd=1.0F/gamma0;
      for(index0i=0;index0i<target.Nd3;index0i++){
         p[index0i]=fcomplex_sub(fcomplex_scale(vtilde[index0i],tempd),fcomplex_mul(mu,p[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. q=interaction_matrix^{T}omegatilde/xi-(gamma0*mu/xi)*q */
      /* Check for iterative.breakdown */
      if(xi<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","xi<iterative.breakdown"); 
      }
      tempdc=fcomplex_scale(mu,gamma0/xi);
      tempd=1.0F/xi;
      dftmatvec_S(omegatilde);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*omegatilde=M^{-1}interaction_matrix*omegatilde */
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
            q[array2]=fcomplex_sub(fcomplex_scale(zero_padded_vector[array1+index0i],tempd),fcomplex_mul(tempdc,q[array2]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. vtilde=interaction_matrix*p-(tau/gamma0)vtilde */
      tempdc=fcomplex_scale(tau,(1.0F/gamma0));
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
            vtilde[array2]=fcomplex_sub(zero_padded_vector[array1+index0i],fcomplex_mul(tempdc,vtilde[array2]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. omegatilde=q-(tau/xi)omegatilde */
      tempdc=fcomplex_scale(tau,(1.0F/xi));
      for(index0i=0;index0i<target.Nd3;index0i++){
         omegatilde[index0i]=fcomplex_sub(q[index0i],fcomplex_mul(tempdc,omegatilde[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. Check stopping criterion */
/* r=incident_E_field-interaction_matrixdipole_polarisation */
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
         break;
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. gamma0=||vtilde||
      xi=||omegatilde||
      rho=[omegatilde]^{T}vtilde
      epsilon=(interaction_matrix^{T}omegatilde)^{T}vtilde */
      gamma1=gamma0;
      xi0=xi;
      rho0=rho;
      gamma0=fcomplex_2vectornorm(vtilde);
      xi=fcomplex_2vectornorm(omegatilde);
      rho=zero;
      for(index0i=0;index0i<target.Nd3;index0i++){
         rho=fcomplex_add(rho,fcomplex_mul(omegatilde[index0i],vtilde[index0i]));
      }
      epsilon=zero;
      dftmatvec_S(omegatilde);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*omegatilde=M^{-1}interaction_matrix*omegatilde */
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
            epsilon=fcomplex_add(epsilon,fcomplex_mul(zero_padded_vector[array1+index0i],vtilde[array0+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. mu=(gamma1*xi0*rho)/(gamma0*tau*rho0) */
      /* Check for iterative.breakdown */
      if(gamma0<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","gamma0<iterative.breakdown"); 
      }
      if(fcomplex_abs(tau)<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","|tau|<iterative.breakdown"); 
      }
      if(fcomplex_abs(rho0)<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","|rho0|<iterative.breakdown"); 
      }
      mu=fcomplex_scale(fcomplex_div(rho,fcomplex_mul(tau,rho0)),(gamma1*xi0)/gamma0);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. tau=(epsilon/rho)-gamma0*mu */
      if(fcomplex_abs(rho)<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","|rho|<iterative.breakdown"); 
      }
      tau0=tau;
      tau=fcomplex_sub(fcomplex_div(epsilon,rho),fcomplex_scale(mu,gamma0));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. theta=(|tau0|^{2}(1-lambda))/(lambda*|tau0|^{2}+|gamma0|^{2}) */
      if(lambda<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","lambda<iterative.breakdown"); 
      }
      if(fcomplex_abs(tau)<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","|tau0|<iterative.breakdown");
      }
      denom=(lambda*fcomplex_abs2(tau0))+(gamma0*gamma0);
      theta=(fcomplex_abs2(tau0)*(1.0F-lambda))/denom;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. kappa=(-gamma1*CONJ(tau0)*kappa0)/(lambda*|tau0|^{2}+|gamma0|^{2}) */
      kappa0=kappa;
      tempdc=tau0;
      tempdc.dat[1]*=-1.0F;
      kappa=fcomplex_scale(fcomplex_scale(fcomplex_mul(tempdc,kappa0),-gamma1),(1.0F/denom));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. lambda=(lambda*|tau0|^{2})/(lambda*|tau|^{2}+|gamma0|^{2}) */
      lambda=(lambda*fcomplex_abs2(tau0))/denom;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. d=theta*d+kappa*p */
      for(index0i=0;index0i<target.Nd3;index0i++){
         d[index0i]=fcomplex_add(fcomplex_scale(d[index0i],theta),fcomplex_mul(kappa,p[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. s=theta*s+kappa*interaction_matrix*p */
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
            s[array2]=fcomplex_add(fcomplex_scale(s[array2],theta),fcomplex_mul(kappa,zero_padded_vector[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. dipole_polarisation=dipole_polarisation+d */
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=fcomplex_add(dipole_polarisation[index0i],d[index0i]);
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 20. r=r-s */
      for(index0i=0;index0i<target.Nd3;index0i++){
         r[index0i]=fcomplex_sub(r[index0i],s[index0i]);
      }
      iteration++;
   }
   if(timing.enabled){ /* Finalise timing */
      finish=walltime();
      fprintf(output.iter,"\n  Time for convergence: %.*gs\n\n",DBLP,finish-start);
      timing.iterative_solution+=finish-start;
   }
}
