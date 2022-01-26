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
   beta,denom,numer,temp,
   result,start,finish,element: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   reduce: Used to facilitate the MPI_Allreduce function
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dcomplex_bicgstab_MPI.h"

void dcomplex_bicgstab_MPI(void){

   int vc,iteration=0;
   
   long int index0i,array0,array1,array2;
   dcomplex rho,rhok,alpha,omega,beta,denom,numer,reduce;
   double ratio,norm0,norm1,terminate,temp,result,start,finish;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=MPI_Wtime();}

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 1. r=incident_E_field-interaction_matrix*dipole_polarisation */
   if(iterative.initial_guess>0){
      dftmatvec_MPI(dipole_polarisation);
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            r[array2]=dcomplex_sub(incident_E_field[array2],zero_padded_vector[array1+index0i]);
         }
      }
   }
   else{ /* r=incident_E_field */
      memcpy(r,incident_E_field,(size_t)parallel.alloc_vector3*sizeof(dcomplex));
   }

   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* r=M^{-1}r */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            r[index0i]=dcomplex_mul(r[index0i],point_jacobi[index0i]);
         }
      }
   }

   norm0=dcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0; 
   };
   ratio=1.0;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. rtilde=r */
   memcpy(rtilde,r,(size_t)parallel.alloc_vector3*sizeof(dcomplex));

   if(parallel.myrank==0){ /* Restrict to master */
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
   }

   terminate=iterative.tolerance;

   while((ratio>(terminate+DBL_EPSILON))&&(iteration<iterative.maximum)){ /* 17. Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. rhok=conj(rtilde^{T})r */
      reduce=zero;
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(rtilde[index0i]),r[index0i]));
      }
      MPI_Allreduce(&reduce,&rhok,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

      if(iteration==0){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. p=r */
         memcpy(p,r,(size_t)parallel.alloc_vector3*sizeof(dcomplex));
      }
      else{
         /* Check for iterative.breakdown */
         if((dcomplex_abs(rho))<iterative.breakdown){
            print_error("Iterative scheme error","BiCGSTAB() iterative.breakdown","|rho|<iterative.breakdown");
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. beta=(rhok*alpha)/(rho*omega) */
         beta=dcomplex_div(dcomplex_mul(alpha,rhok),dcomplex_mul(omega,rho));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. p=r+beta(p-omega*v) */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            p[index0i]=dcomplex_add(r[index0i],dcomplex_mul(beta,dcomplex_sub(p[index0i],dcomplex_mul(omega,v[index0i]))));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. v=interaction_matrix*p */
      dftmatvec_MPI(p);
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            v[array0+index0i]=zero_padded_vector[array1+index0i];
         }
      }
      
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* v=M^{-1}v */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
               v[index0i]=dcomplex_mul(v[index0i],point_jacobi[index0i]);
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. alpha=rhok/(conj(rtilde^{T})v) */
      reduce=zero;
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(rtilde[index0i]),v[index0i]));
      }
      MPI_Allreduce(&reduce,&alpha,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
      /* Check for iterative.breakdown */
      if((dcomplex_abs(alpha))<iterative.breakdown){
         print_error("Iterative scheme error","BiCGSTAB() iterative.breakdown","|alpha|<iterative.breakdown"); 
      }
      alpha=dcomplex_div(rhok,alpha);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. s=r-alpha*v */
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         s[index0i]=dcomplex_sub(r[index0i],dcomplex_mul(alpha,v[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. Check for soft iterative.breakdown ||s||<iterative.breakdown*/
      norm1=dcomplex_2vectornorm(s);
      if(norm1<iterative.breakdown){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* dipole_polarisation=dipole_polarisation+alpha*p */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            dipole_polarisation[index0i]=dcomplex_add(dipole_polarisation[index0i],dcomplex_mul(alpha,p[index0i]));
         }
         fprintf(stderr,"BiCGSTAB() soft iterative.breakdown (||s||=%.*g<iterative.breakdown=%.*g)...\n",DBLP,norm1,DBLP,iterative.breakdown);
         break;
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. t=interaction_matrix*s */
/* t is stored in zero_padded_vector */
      dftmatvec_MPI(s);

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* t=M^{-1}t */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=vc*parallel.alloc_padded;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=dcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. omega=(conj(t^{T})s)/(conj(t^{T})t) */
/* t is stored in zero_padded_vector */
      reduce=zero;
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(zero_padded_vector[array1+index0i]),s[array0+index0i]));
         }
      }
      MPI_Allreduce(&reduce,&numer,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
      reduce=zero;denom.dat[1]=0.0;
      for(vc=0;vc<3;vc++){
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            reduce.dat[0]+=dcomplex_abs2(zero_padded_vector[array1+index0i]);
         }
      }
      MPI_Allreduce(&reduce.dat[0],&denom.dat[0],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      /* Check for iterative.breakdown */
      if((dcomplex_abs(denom))<iterative.breakdown){
         print_error("Iterative scheme error","BiCGSTAB() iterative.breakdown","|t^{T}t|<iterative.breakdown");
      }
      omega=dcomplex_div(numer,denom);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. dipole_polarisation=dipole_polarisation+alpha*p+omega*s */
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         dipole_polarisation[index0i]=dcomplex_add(dipole_polarisation[index0i],dcomplex_add(dcomplex_mul(alpha,p[index0i]),dcomplex_mul(omega,s[index0i])));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. r=s-omega*t */
/* t is stored in zero_padded_vector */
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            r[array2]=dcomplex_sub(s[array2],dcomplex_mul(omega,zero_padded_vector[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. Check stopping criterion */
      norm1=dcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      if(parallel.myrank==0){ /* Restrict to master */
         fprintf(output.iter,"  At iteration %d the convergence ratio is %.*g and the residual norm is %.*g\n",iteration,DBLP,ratio,DBLP,norm1);
      }
      rho=rhok;
      iteration++;
   }
   if(timing.enabled){ /* Finalise timing */
      finish=MPI_Wtime();
      temp=finish-start;
      result=0.0;
      MPI_Reduce(&temp,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      if(parallel.myrank==0){ /* Restrict to master */   
         fprintf(output.iter,"\n  Time for convergence [Averaged over the %d processors]: %.*gs\n\n",parallel.np,DBLP,(result/parallel.np));
      }
      timing.iterative_solution+=result;
   }
}
