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
   beta,denom,numer: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "ldcomplex_bicgstab_S.h"

void ldcomplex_bicgstab_S(void){

   int vc,iteration=0;
   long int index0i,array0,array1,array2;
   ldcomplex rho,rhok,alpha,omega,beta,denom,numer;
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
   ratio=1.0L;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. rtilde=r */
   memcpy(rtilde,r,(size_t)target.Nd3*sizeof(ldcomplex));

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

   while((ratio>(terminate+LDBL_EPSILON))&&(iteration<iterative.maximum)){ /* 17. Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. rhok=conj(rtilde^{T})r */
      rhok=zero;
      for(index0i=0;index0i<target.Nd3;index0i++){
         rhok=ldcomplex_add(rhok,ldcomplex_mul(ldcomplex_conj(rtilde[index0i]),r[index0i]));
      }

      if(iteration==0){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. p=r */
         memcpy(p,r,(size_t)target.Nd3*sizeof(ldcomplex));
      }
      else{
         /* Check for iterative.breakdown */
         if((ldcomplex_abs(rho))<iterative.breakdown){
            print_error("Iterative scheme error","BiCGSTAB() iterative.breakdown","|rho|<iterative.breakdown");
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. beta=(rhok*alpha)/(rho*omega) */
         beta=ldcomplex_div(ldcomplex_mul(alpha,rhok),ldcomplex_mul(omega,rho));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. p=r+beta(p-omega*v) */
         for(index0i=0;index0i<target.Nd3;index0i++){
            p[index0i]=ldcomplex_add(r[index0i],ldcomplex_mul(beta,ldcomplex_sub(p[index0i],ldcomplex_mul(omega,v[index0i]))));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. v=interaction_matrix*p */
      dftmatvec_S(p);
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
               v[index0i]=ldcomplex_mul(v[index0i],point_jacobi[index0i]);
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. alpha=rhok/(conj(rtilde^{T})v) */
      alpha=zero;
      for(index0i=0;index0i<target.Nd3;index0i++){
         alpha=ldcomplex_add(alpha,ldcomplex_mul(ldcomplex_conj(rtilde[index0i]),v[index0i]));
      }
      /* Check for iterative.breakdown */
      if((ldcomplex_abs(alpha))<iterative.breakdown){
         print_error("Iterative scheme error","BiCGSTAB() iterative.breakdown","|alpha|<iterative.breakdown"); 
      }
      alpha=ldcomplex_div(rhok,alpha);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. s=r-alpha*v */
      for(index0i=0;index0i<target.Nd3;index0i++){
         s[index0i]=ldcomplex_sub(r[index0i],ldcomplex_mul(alpha,v[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. Check for soft iterative.breakdown ||s||<iterative.breakdown*/
      norm1=ldcomplex_2vectornorm(s);
      if(norm1<iterative.breakdown){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* dipole_polarisation=dipole_polarisation+alpha*p */
         for(index0i=0;index0i<target.Nd3;index0i++){
            dipole_polarisation[index0i]=ldcomplex_add(dipole_polarisation[index0i],ldcomplex_mul(alpha,p[index0i]));
         }
         fprintf(stderr,"BiCGSTAB() soft iterative.breakdown (||s||=%.*Lg<iterative.breakdown=%.*Lg)...\n",LDBLP,norm1,LDBLP,iterative.breakdown);
         break;
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. t=interaction_matrix*s */
/* t is stored in zero_padded_vector */
      dftmatvec_S(s);

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* t=M^{-1}t */
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
/* 14. omega=(conj(t^{T})s)/(conj(t^{T})t) */
/* t is stored in zero_padded_vector */
      numer=zero;
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            numer=ldcomplex_add(numer,ldcomplex_mul(ldcomplex_conj(zero_padded_vector[array1+index0i]),s[array0+index0i]));
         }
      }
      denom=zero;
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            denom.dat[0]+=ldcomplex_abs2(zero_padded_vector[array0+index0i]);
         }
      }
      /* Check for iterative.breakdown */
      if((ldcomplex_abs(denom))<iterative.breakdown){
         print_error("Iterative scheme error","BiCGSTAB() iterative.breakdown","|t^{T}t|<iterative.breakdown");
      }
      omega=ldcomplex_div(numer,denom);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. dipole_polarisation=dipole_polarisation+alpha*p+omega*s */
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=ldcomplex_add(dipole_polarisation[index0i],ldcomplex_add(ldcomplex_mul(alpha,p[index0i]),ldcomplex_mul(omega,s[index0i])));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. r=s-omega*t */
/* t is stored in zero_padded_vector */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            r[array2]=ldcomplex_sub(s[array2],ldcomplex_mul(omega,zero_padded_vector[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. Check stopping criterion */
      norm1=ldcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      fprintf(output.iter,"  At iteration %d the convergence ratio is %.*Lg and the residual norm is %.*Lg\n",iteration,LDBLP,ratio,LDBLP,norm1);
      rho=rhok;
      iteration++;
   }
   if(timing.enabled){ /* Finalise timing */
      finish=walltime();
      fprintf(output.iter,"\n  Time for convergence: %.*gs\n\n",DBLP,finish-start);
      timing.iterative_solution+=finish-start;
   }
}
