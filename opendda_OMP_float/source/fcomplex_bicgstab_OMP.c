/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   Stabilised version of the BiConjugate-Gradients: BiCGSTAB

   Bi-CGSTAB : A fast and smoothly converging variant of Bi-CG
   for the solution of nonsymmetric linear systems, van der Vorst, H.A.,
   SIAM J. Sci. Stat. Comput., Vol. 13, No. 2, Pages 631-644, 1992.

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
   v[3],s[3]: Local vectors for iterative scheme
   rho,rhok,alpha,omega,
   beta,denom,numer,real,imag,v0,v1: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "fcomplex_bicgstab_OMP.h"

void fcomplex_bicgstab_OMP(void){

   int vc,iteration=0;
   long int index0i,array0,array1,array2;
   fcomplex rho,rhok,alpha,omega,beta,denom,numer,v0,v1;
   float ratio,norm0,norm1,terminate,real,imag;
   double start,finish;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=walltime();}

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 1. r=incident_E_field-interaction_matrix*dipole_polarisation */
   if(iterative.initial_guess>0){
      dftmatvec_OMP(dipole_polarisation);
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
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
#pragma omp parallel for schedule(dynamic,target.M)
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

   while((ratio>(terminate+FLT_EPSILON))&&(iteration<iterative.maximum)){ /* 17. Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. rhok=conj(rtilde^{T})r */
      real=0.0F;imag=0.0F;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         v0=rtilde[index0i];
         v1=r[index0i];
         real+=(v0.dat[0]*v1.dat[0]+v0.dat[1]*v1.dat[1]);
         imag+=(v0.dat[0]*v1.dat[1]-v0.dat[1]*v1.dat[0]);
      }
      rhok.dat[0]=real;rhok.dat[1]=imag;

      if(iteration==0){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. p=r */
         memcpy(p,r,(size_t)target.Nd3*sizeof(fcomplex));
      }
      else{
         /* Check for iterative.breakdown */
         if((fcomplex_abs(rho))<iterative.breakdown){
            print_error("Iterative scheme error","BiCGSTAB() iterative.breakdown","|rho|<iterative.breakdown");
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. beta=(rhok*alpha)/(rho*omega) */
         beta=fcomplex_div(fcomplex_mul(alpha,rhok),fcomplex_mul(omega,rho));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. p=r+beta(p-omega*v) */
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            p[index0i]=fcomplex_add(r[index0i],fcomplex_mul(beta,fcomplex_sub(p[index0i],fcomplex_mul(omega,v[index0i]))));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. v=interaction_matrix*p */
      dftmatvec_OMP(p);
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
               v[index0i]=fcomplex_mul(v[index0i],point_jacobi[index0i]);
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. alpha=rhok/(conj(rtilde^{T})v) */
      real=0.0F;imag=0.0F;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         v0=rtilde[index0i];
         v1=v[index0i];
         real+=(v0.dat[0]*v1.dat[0]+v0.dat[1]*v1.dat[1]);
         imag+=(v0.dat[0]*v1.dat[1]-v0.dat[1]*v1.dat[0]);
      }
      alpha.dat[0]=real;alpha.dat[1]=imag;
      /* Check for iterative.breakdown */
      if((fcomplex_abs(alpha))<iterative.breakdown){
         print_error("Iterative scheme error","BiCGSTAB() iterative.breakdown","|alpha|<iterative.breakdown"); 
      }
      alpha=fcomplex_div(rhok,alpha);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. s=r-alpha*v */
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         s[index0i]=fcomplex_sub(r[index0i],fcomplex_mul(alpha,v[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. Check for soft iterative.breakdown ||s||<iterative.breakdown*/
      norm1=fcomplex_2vectornorm(s);
      if(norm1<iterative.breakdown){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* dipole_polarisation=dipole_polarisation+alpha*p */
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            dipole_polarisation[index0i]=fcomplex_add(dipole_polarisation[index0i],fcomplex_mul(alpha,p[index0i]));
         }
         fprintf(stderr,"BiCGSTAB() soft iterative.breakdown (||s||=%.*g<iterative.breakdown=%.*g)...\n",FLTP,norm1,FLTP,iterative.breakdown);
         break;
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. t=interaction_matrix*s */
/* t is stored in zero_padded_vector */
      dftmatvec_OMP(s);

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* t=M^{-1}t */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=fcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. omega=(conj(t^{T})s)/(conj(t^{T})t) */
/* t is stored in zero_padded_vector */
      real=0.0F;imag=0.0F;
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            v0=s[array0+index0i];
            v1=zero_padded_vector[array1+index0i];
            real+=(v1.dat[0]*v0.dat[0]+v1.dat[1]*v0.dat[1]);
            imag+=(v1.dat[0]*v0.dat[1]-v1.dat[1]*v0.dat[0]);
         }
      }
      numer.dat[0]=real;numer.dat[1]=imag;
      real=0.0F;
      for(vc=0;vc<3;vc++){
         array1=vc*target.Nv;
#pragma omp parallel for reduction(+:real) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            real+=fcomplex_abs2(zero_padded_vector[array1+index0i]);
         }
      }
      denom.dat[0]=real;denom.dat[1]=0.0;
      /* Check for iterative.breakdown */
      if((fcomplex_abs(denom))<iterative.breakdown){
         print_error("Iterative scheme error","BiCGSTAB() iterative.breakdown","|t^{T}t|<iterative.breakdown");
      }
      omega=fcomplex_div(numer,denom);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. dipole_polarisation=dipole_polarisation+alpha*p+omega*s */
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=fcomplex_add(dipole_polarisation[index0i],fcomplex_add(fcomplex_mul(alpha,p[index0i]),fcomplex_mul(omega,s[index0i])));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. r=s-omega*t */
/* t is stored in zero_padded_vector */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            r[array2]=fcomplex_sub(s[array2],fcomplex_mul(omega,zero_padded_vector[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. Check stopping criterion */
      norm1=fcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      fprintf(output.iter,"  At iteration %d the convergence ratio is %.*g and the residual norm is %.*g\n",iteration,FLTP,ratio,FLTP,norm1);
      rho=rhok;
      iteration++;
   }
   if(timing.enabled){ /* Finalise timing */
      finish=walltime();
      fprintf(output.iter,"\n  Time for convergence: %.*gs\n\n",DBLP,finish-start);
      timing.iterative_solution+=finish-start;
   }
}
