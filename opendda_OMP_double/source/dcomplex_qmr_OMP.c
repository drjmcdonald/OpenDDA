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
   lambda,theta,denom,real,imag,v0,v1: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dcomplex_qmr_OMP.h"

void dcomplex_qmr_OMP(void){

   int vc,iteration=0;
   long int index0i,array0,array1,array2;
   dcomplex tempdc={{0.0,0.0}},kappa={{-1.0,0.0}},kappa0,rho,rho0;
   dcomplex epsilon,mu={{0.0,0.0}},tau,tau0,v0,v1;
   double ratio,norm0,norm1,gamma0,gamma1,xi,xi0,lambda=1.0;
   double theta=-1.0,denom,tempd,terminate,start,finish,real,imag;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=walltime();}

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 1. lambda=1,kappa=theta=-1 Set at variable declaration stage */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. r=incident_E_field-interaction_matrix(dipole_polarisation),omegatilde=vtilde=r */
   if(iterative.initial_guess>0){
      dftmatvec_OMP(dipole_polarisation);
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
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
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            r[index0i]=dcomplex_mul(r[index0i],point_jacobi[index0i]);
         }
      }
   }

   norm0=dcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0; 
   };

   /* omegatilde=vtilde=r */
   memcpy(omegatilde,r,(size_t)target.Nd3*sizeof(dcomplex));
   memcpy(vtilde,r,(size_t)target.Nd3*sizeof(dcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. p=q=d=s=0 Set to zero using calloc in memory allocation */ 

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. gamma0=||vtilde||
      xi=||omegatilde||
      rho=[omegatilde]^{T}vtilde
      epsilon=(interaction_matrix^{T}omegatilde)^{T}vtilde
      mu=0 set at variable declaration stage */
/*   gamma0=dcomplex_2vectornorm(vtilde);
   xi=dcomplex_2vectornorm(omegatilde); */
   gamma0=norm0;xi=norm0;
   real=0.0;imag=0.0;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
   for(index0i=0;index0i<target.Nd3;index0i++){
      v0=omegatilde[index0i];
      v1=vtilde[index0i];
      real+=(v0.dat[0]*v1.dat[0]-v0.dat[1]*v1.dat[1]);
      imag+=(v0.dat[0]*v1.dat[1]+v0.dat[1]*v1.dat[0]);
   }
   rho.dat[0]=real;rho.dat[1]=imag;
   dftmatvec_OMP(omegatilde);
   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*omegatilde=M^{-1}interaction_matrix*omegatilde */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=array1+index0i;
               zero_padded_vector[array2]=dcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
            }
         }
      }
   }
   real=0.0;imag=0.0;
   for(vc=0;vc<3;vc++){
      array0=vc*target.Nd;
      array1=vc*target.Nv;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd;index0i++){
         v0=zero_padded_vector[array1+index0i];
         v1=vtilde[array0+index0i];
         real+=(v0.dat[0]*v1.dat[0]-v0.dat[1]*v1.dat[1]);
         imag+=(v0.dat[0]*v1.dat[1]+v0.dat[1]*v1.dat[0]);
      }
   }
   epsilon.dat[0]=real;epsilon.dat[1]=imag;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. tau=epsilon/rho    */
   /* Check for iterative.breakdown */
   if(dcomplex_abs(rho)<iterative.breakdown){
      print_error("Iterative scheme error","QMR() iterative.breakdown","|rho|<iterative.breakdown"); 
   }
   tau=dcomplex_div(epsilon,rho);

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

   while(iteration<iterative.maximum){ /* Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. p=vtilde/gamma0-mu*p */
      /* Check for iterative.breakdown */
      if(gamma0<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","gamma0<iterative.breakdown"); 
      }
      tempd=1.0/gamma0;
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         p[index0i]=dcomplex_sub(dcomplex_scale(vtilde[index0i],tempd),dcomplex_mul(mu,p[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. q=interaction_matrix^{T}omegatilde/xi-(gamma0*mu/xi)*q */
      /* Check for iterative.breakdown */
      if(xi<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","xi<iterative.breakdown"); 
      }
      tempdc=dcomplex_scale(mu,gamma0/xi);
      tempd=1.0/xi;
      dftmatvec_OMP(omegatilde);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*omegatilde=M^{-1}interaction_matrix*omegatilde */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=dcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            q[array2]=dcomplex_sub(dcomplex_scale(zero_padded_vector[array1+index0i],tempd),dcomplex_mul(tempdc,q[array2]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. vtilde=interaction_matrix*p-(tau/gamma0)vtilde */
      tempdc=dcomplex_scale(tau,(1.0/gamma0));
      dftmatvec_OMP(p);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*p=M^{-1}interaction_matrix*p */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=dcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            vtilde[array2]=dcomplex_sub(zero_padded_vector[array1+index0i],dcomplex_mul(tempdc,vtilde[array2]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. omegatilde=q-(tau/xi)omegatilde */
      tempdc=dcomplex_scale(tau,(1.0/xi));
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         omegatilde[index0i]=dcomplex_sub(q[index0i],dcomplex_mul(tempdc,omegatilde[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. Check stopping criterion */
/* r=incident_E_field-interaction_matrixdipole_polarisation */
      dftmatvec_OMP(dipole_polarisation);
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            r[array2]=dcomplex_sub(incident_E_field[array2],zero_padded_vector[array1+index0i]);
         }
      }

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* r=M^{-1}r */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
#pragma omp parallel for schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd3;index0i++){
               r[index0i]=dcomplex_mul(r[index0i],point_jacobi[index0i]);
            }
         }
      }

      norm1=dcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      fprintf(output.iter,"  At iteration %d the convergence ratio is %.*g and the residual norm is %.*g\n",iteration,DBLP,ratio,DBLP,norm1);
      if(ratio<(terminate+DBL_EPSILON)){
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
      gamma0=dcomplex_2vectornorm(vtilde);
      xi=dcomplex_2vectornorm(omegatilde);
      real=0.0;imag=0.0;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         v0=omegatilde[index0i];
         v1=vtilde[index0i];
         real+=(v0.dat[0]*v1.dat[0]-v0.dat[1]*v1.dat[1]);
         imag+=(v0.dat[0]*v1.dat[1]+v0.dat[1]*v1.dat[0]);
      }
      rho.dat[0]=real;rho.dat[1]=imag;
      dftmatvec_OMP(omegatilde);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*omegatilde=M^{-1}interaction_matrix*omegatilde */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=dcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }
      real=0.0;imag=0.0;
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            v0=zero_padded_vector[array1+index0i];
            v1=vtilde[array0+index0i];
            real+=(v0.dat[0]*v1.dat[0]-v0.dat[1]*v1.dat[1]);
            imag+=(v0.dat[0]*v1.dat[1]+v0.dat[1]*v1.dat[0]);
         }
      }
      epsilon.dat[0]=real;epsilon.dat[1]=imag;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. mu=(gamma1*xi0*rho)/(gamma0*tau*rho0) */
      /* Check for iterative.breakdown */
      if(gamma0<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","gamma0<iterative.breakdown"); 
      }
      if(dcomplex_abs(tau)<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","|tau|<iterative.breakdown"); 
      }
      if(dcomplex_abs(rho0)<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","|rho0|<iterative.breakdown"); 
      }
      mu=dcomplex_scale(dcomplex_div(rho,dcomplex_mul(tau,rho0)),(gamma1*xi0)/gamma0);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. tau=(epsilon/rho)-gamma0*mu */
      if(dcomplex_abs(rho)<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","|rho|<iterative.breakdown"); 
      }
      tau0=tau;
      tau=dcomplex_sub(dcomplex_div(epsilon,rho),dcomplex_scale(mu,gamma0));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. theta=(|tau0|^{2}(1-lambda))/(lambda*|tau0|^{2}+|gamma0|^{2}) */
      if(lambda<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","lambda<iterative.breakdown"); 
      }
      if(dcomplex_abs(tau)<iterative.breakdown){
         print_error("Iterative scheme error","QMR() iterative.breakdown","|tau0|<iterative.breakdown");
      }
      denom=(lambda*dcomplex_abs2(tau0))+(gamma0*gamma0);
      theta=(dcomplex_abs2(tau0)*(1.0-lambda))/denom;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. kappa=(-gamma1*CONJ(tau0)*kappa0)/(lambda*|tau0|^{2}+|gamma0|^{2}) */
      kappa0=kappa;
      tempdc=tau0;
      tempdc.dat[1]*=-1.0;
      kappa=dcomplex_scale(dcomplex_scale(dcomplex_mul(tempdc,kappa0),-gamma1),(1.0/denom));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. lambda=(lambda*|tau0|^{2})/(lambda*|tau|^{2}+|gamma0|^{2}) */
      lambda=(lambda*dcomplex_abs2(tau0))/denom;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. d=theta*d+kappa*p */
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         d[index0i]=dcomplex_add(dcomplex_scale(d[index0i],theta),dcomplex_mul(kappa,p[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. s=theta*s+kappa*interaction_matrix*p */
      dftmatvec_OMP(p);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*p=M^{-1}interaction_matrix*p */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=dcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            s[array2]=dcomplex_add(dcomplex_scale(s[array2],theta),dcomplex_mul(kappa,zero_padded_vector[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. dipole_polarisation=dipole_polarisation+d */
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=dcomplex_add(dipole_polarisation[index0i],d[index0i]);
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 20. r=r-s */
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         r[index0i]=dcomplex_sub(r[index0i],s[index0i]);
      }
      iteration++;
   }
   if(timing.enabled){ /* Finalise timing */
      finish=walltime();
      ratio=finish-start;
      fprintf(output.iter,"\n  Time for convergence: %.*gs\n\n",DBLP,ratio);
      timing.iterative_solution+=ratio;
   }
}
