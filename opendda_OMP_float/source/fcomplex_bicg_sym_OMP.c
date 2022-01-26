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
   alpha,beta,real,imag,v0,v1: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "fcomplex_bicg_sym_OMP.h"

void fcomplex_bicg_sym_OMP(void){

   int vc,iteration=0;
   long int index0i,array0,array1,array2;
   fcomplex rho,rhok,alpha,beta,delta,v0,v1;
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

   while(iteration<iterative.maximum){ /* 11. Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. rhok=r^{T}r */
      real=0.0F;imag=0.0F;
#pragma omp parallel for private(v0) reduction(+:real,imag) schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         v0=r[index0i];
         real+=(v0.dat[0]*v0.dat[0]-v0.dat[1]*v0.dat[1]);
         imag+=(2.0F*(v0.dat[0]*v0.dat[1]));
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
            print_error("Iterative scheme error","BiCG_sym() iterative.breakdown","|rho|<iterative.breakdown"); 
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. beta=rhok/rho */
         beta=fcomplex_div(rhok,rho);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. p=r+beta*p */
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            p[index0i]=fcomplex_add(r[index0i],fcomplex_mul(beta,p[index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. interaction_matrix*p */
/* interaction_matrix*p is stored in zero_padded_vector */
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
                  zero_padded_vector[array2]=fcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. alpha=p^{T}interaction_matrix*p */
/* interaction_matrix*p is stored in zero_padded_vector */
      real=0.0F;imag=0.0F;
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            v0=p[array0+index0i];
            v1=zero_padded_vector[array1+index0i];
            real+=(v0.dat[0]*v1.dat[0]-v0.dat[1]*v1.dat[1]);
            imag+=(v0.dat[0]*v1.dat[1]+v0.dat[1]*v1.dat[0]);
         }
      }
      alpha.dat[0]=real;alpha.dat[1]=imag;
      /* Check for iterative.breakdown */
      if((fcomplex_abs(alpha))<iterative.breakdown){
         print_error("Iterative scheme error","BiCG_sym() iterative.breakdown","|alpha|<iterative.breakdown"); 
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. delta=rhok/alpha */
      delta=fcomplex_div(rhok,alpha);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. dipole_polarisation=dipole_polarisation+delta*p */
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=fcomplex_add(dipole_polarisation[index0i],fcomplex_mul(delta,p[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. r=r-delta*interaction_matrix*p */
/* interaction_matrix*p is stored in zero_padded_vector */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            r[array2]=fcomplex_sub(r[array2],fcomplex_mul(delta,zero_padded_vector[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. Check stopping criterion */
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
