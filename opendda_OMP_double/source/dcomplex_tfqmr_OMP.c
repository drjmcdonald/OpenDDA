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
   tau0,theta0,c,kappa,real,imag,v0,v1: Local variables for iterative scheme
   error,ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dcomplex_tfqmr_OMP.h"

void dcomplex_tfqmr_OMP(void){

   int m,vc,iteration=0,exitflag=0;
   long int index0i,array0,array1,array2;
   dcomplex eta={{0.0,0.0}},eta0,rho,rhok,sigma,alpha,beta,temp,*ptr,v0,v1;
   double ratio,norm0,norm1,tau,theta=0.0,tau0,theta0,c,kappa,error;
   double terminate,start,finish,real,imag;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=walltime();}

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 1. r=incident_E_field-interaction_matrix(dipole_polarisation) */
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
   ratio=1.0;
   error=iterative.tolerance*norm0;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. omega=yb=r */
   memcpy(omega,r,(size_t)target.Nd3*sizeof(dcomplex));
   memcpy(yb,r,(size_t)target.Nd3*sizeof(dcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. v=interaction_matrix*yb */
   dftmatvec_OMP(yb);
   for(vc=0;vc<3;vc++){
      array0=vc*target.Nd;
      array1=vc*target.Nv;
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd;index0i++){
         v[array0+index0i]=zero_padded_vector[array1+index0i];
      }
   }

   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* v=M^{-1}v */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            v[index0i]=dcomplex_mul(v[index0i],point_jacobi[index0i]);
         }
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. d=0 Set to zero using calloc in memory allocation */ 

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. tau=||r|| */
   tau=dcomplex_2vectornorm(r);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. theta=0.0,eta=0.0+0.0i Set at variable declaration stage */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. rtilde=r */
   memcpy(rtilde,r,(size_t)target.Nd3*sizeof(dcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. rhok=rtilde^{T}r */
   real=0.0;imag=0.0;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
   for(index0i=0;index0i<target.Nd3;index0i++){
      v0=rtilde[index0i];
      v1=r[index0i];
      real+=(v0.dat[0]*v1.dat[0]-v0.dat[1]*v1.dat[1]);
      imag+=(v0.dat[0]*v1.dat[1]+v0.dat[1]*v1.dat[0]);
   }
   rhok.dat[0]=real;rhok.dat[1]=imag;

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
/* 9. sigma=rtilde^{T}v */
      real=0.0;imag=0.0;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         v0=rtilde[index0i];
         v1=v[index0i];
         real+=(v0.dat[0]*v1.dat[0]-v0.dat[1]*v1.dat[1]);
         imag+=(v0.dat[0]*v1.dat[1]+v0.dat[1]*v1.dat[0]);
      }
      sigma.dat[0]=real;sigma.dat[1]=imag;
      /* Check for iterative.breakdown */
      if(dcomplex_abs(sigma)<iterative.breakdown){
         print_error("Iterative scheme error","TFQMR() iterative.breakdown","|sigma|<iterative.breakdown");
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. alpha=rhok/sigma */
      alpha=dcomplex_div(rhok,sigma);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. ya=yb,yb=ya-alpha*v */
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         ya[index0i]=dcomplex_sub(yb[index0i],dcomplex_mul(alpha,v[index0i]));
      }
      /* ya=yb,yb=ya-alpha*v */
      ptr=ya;ya=yb;yb=ptr;

      for(m=2*iteration-1;m<=2*iteration;m++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. omega=omega-alpha*interaction_matrix*ya */
         dftmatvec_OMP(ya);
         if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*ya=M^{-1}interaction_matrix*ya */
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
               omega[array2]=dcomplex_sub(omega[array2],dcomplex_mul(alpha,zero_padded_vector[array1+index0i]));
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
         theta=dcomplex_2vectornorm(omega)/tau0;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. c=(1)/(sqrt(1+theta^{2})) */
         c=1.0/(sqrt(1.0+theta*theta));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. tau=tau0*theta*c */
         tau=tau0*theta*c;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. eta=c^{2}*alpha */
         eta0=eta;
         eta=dcomplex_scale(alpha,c*c);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. d=ya+(theta0^{2}*eta0/alpha)*d */
         temp=dcomplex_scale(dcomplex_div(eta0,alpha),theta0*theta0);
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            d[index0i]=dcomplex_add(ya[index0i],dcomplex_mul(temp,d[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. dipole_polarisation=dipole_polarisation+eta*d */
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            dipole_polarisation[index0i]=dcomplex_add(dipole_polarisation[index0i],dcomplex_mul(eta,d[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. kappa=tau*sqrt(m+1) */
         kappa=tau*sqrt((double)m+1.0);
         
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 20. Check stopping criterion */
         if(kappa<error){
/* r=incident_E_field-interaction_matrix*dipole_polarisation)*/
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
      real=0.0;imag=0.0;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         v0=rtilde[index0i];
         v1=omega[index0i];
         real+=(v0.dat[0]*v1.dat[0]-v0.dat[1]*v1.dat[1]);
         imag+=(v0.dat[0]*v1.dat[1]+v0.dat[1]*v1.dat[0]);
      }
      rhok.dat[0]=real;rhok.dat[1]=imag;
      /* Check for iterative.breakdown */
      if(dcomplex_abs(rho)<iterative.breakdown){
         print_error("Iterative scheme error","TFQMR() iterative.breakdown","|rho|<iterative.breakdown"); 
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 22. beta=rhok/rho   */
      beta=dcomplex_div(rhok,rho);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 23. yb=omega+beta*ya */
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         yb[index0i]=dcomplex_add(omega[index0i],dcomplex_mul(beta,ya[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 24. v=interaction_matrix*yb+beta*(interaction_matrix*ya+beta*v) */
      dftmatvec_OMP(ya);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*ya=M^{-1}interaction_matrix*ya */
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
            array2=array1+index0i;
            zero_padded_vector[array2]=dcomplex_add(zero_padded_vector[array2],dcomplex_mul(beta,v[array0+index0i]));
         }
      }
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            v[array0+index0i]=dcomplex_mul(beta,zero_padded_vector[array1+index0i]);
         }
      }
      dftmatvec_OMP(yb);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*yb=M^{-1}interaction_matrix*yb */
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
            v[array2]=dcomplex_add(v[array2],zero_padded_vector[array1+index0i]);
         }
      }
      iteration++;
   } /* End while loop */
   if(timing.enabled){ /* Finalise timing */
      finish=walltime();
      ratio=finish-start;
      fprintf(output.iter,"\n  Time for convergence: %.*gs\n\n",DBLP,ratio);
      timing.iterative_solution+=ratio;
   }
}
