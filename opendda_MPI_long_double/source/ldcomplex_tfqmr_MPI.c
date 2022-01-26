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
   alpha,beta,tempdc,tau,theta,
   tau0,theta0,c,kappa,temp,result,element: Local variables for iterative scheme
   error,ratio,norm0,norm1,terminate: Check stopping criterion variables
   reduce: Used to facilitate the MPI_Allreduce function
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "ldcomplex_tfqmr_MPI.h"

void ldcomplex_tfqmr_MPI(void){

   int m,vc,iteration=0,exitflag=0;
   long int index0i,array0,array1,array2;
   ldcomplex eta={{0.0L,0.0L}},eta0,rho,rhok,sigma,alpha,beta,tempdc,*ptr,reduce;
   long double ratio,norm0,norm1,tau,theta=0.0L,tau0,theta0,c,kappa,error;
   long double terminate;
   double temp,result,start,finish;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=MPI_Wtime();}

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 1. r=incident_E_field-interaction_matrix(dipole_polarisation) */
   if(iterative.initial_guess>0){
      dftmatvec_MPI(dipole_polarisation);
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            r[array2]=ldcomplex_sub(incident_E_field[array2],zero_padded_vector[array1+index0i]);
         }
      }
   }
   else{ /* r=incident_E_field */
      memcpy(r,incident_E_field,(size_t)parallel.alloc_vector3*sizeof(ldcomplex));
   }

   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* r=M^{-1}r */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            r[index0i]=ldcomplex_mul(r[index0i],point_jacobi[index0i]);
         }
      }
   }

   norm0=ldcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0L; 
   };
   ratio=1.0L;
   error=iterative.tolerance*norm0;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. omega=yb=r */
   memcpy(omega,r,(size_t)parallel.alloc_vector3*sizeof(ldcomplex));
   memcpy(yb,r,(size_t)parallel.alloc_vector3*sizeof(ldcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. v=interaction_matrix*yb */
   dftmatvec_MPI(yb);
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
            v[index0i]=ldcomplex_mul(v[index0i],point_jacobi[index0i]);
         }
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. d=0 Set to zero using calloc in memory allocation */ 

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. tau=||r|| */
   tau=ldcomplex_2vectornorm(r);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. theta=0.0,eta=0.0+0.0i Set at variable declaration stage */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. rtilde=r */
   memcpy(rtilde,r,(size_t)parallel.alloc_vector3*sizeof(ldcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. rhok=rtilde^{T}r */
   reduce=zero;
   for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
      reduce=ldcomplex_add(reduce,ldcomplex_mul(rtilde[index0i],r[index0i]));
   }
   MPI_Allreduce(&reduce,&rhok,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);

   if(parallel.myrank==0){ /* Restrict to master */
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
   }

   terminate=iterative.tolerance;

   while(iteration<iterative.maximum){ /* Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. sigma=rtilde^{T}v */
      reduce=zero;
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         reduce=ldcomplex_add(reduce,ldcomplex_mul(rtilde[index0i],v[index0i]));
      }
      MPI_Allreduce(&reduce,&sigma,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);
      /* Check for iterative.breakdown */
      if(ldcomplex_abs(sigma)<iterative.breakdown){
         print_error("Iterative scheme error","TFQMR() iterative.breakdown","|sigma|<iterative.breakdown");
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. alpha=rhok/sigma */
      alpha=ldcomplex_div(rhok,sigma);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. ya=yb,yb=ya-alpha*v */
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         ya[index0i]=ldcomplex_sub(yb[index0i],ldcomplex_mul(alpha,v[index0i]));
      }
      /* ya=yb,yb=ya-alpha*v */
      ptr=ya;ya=yb;yb=ptr;

      for(m=2*iteration-1;m<=2*iteration;m++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. omega=omega-alpha*interaction_matrix*ya */
         dftmatvec_MPI(ya);
         if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*ya=M^{-1}interaction_matrix*ya */
            if(iterative.precond==1){
               /* Point-Jacobi Preconditioning (Divide by the diagonal) */
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=vc*parallel.alloc_padded;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     array2=array1+index0i;
                     zero_padded_vector[array2]=ldcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
                  }
               }
            }
         }
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=vc*parallel.alloc_padded;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               omega[array2]=ldcomplex_sub(omega[array2],ldcomplex_mul(alpha,zero_padded_vector[array1+index0i]));
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
         theta=ldcomplex_2vectornorm(omega)/tau0;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. c=(1)/(sqrt(1+theta^{2})) */
         c=1.0L/(sqrtl(1.0L+theta*theta));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. tau=tau0*theta*c */
         tau=tau0*theta*c;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. eta=c^{2}*alpha */
         eta0=eta;
         eta=ldcomplex_scale(alpha,c*c);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. d=ya+(theta0^{2}*eta0/alpha)*d */
         tempdc=ldcomplex_scale(ldcomplex_div(eta0,alpha),theta0*theta0);
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            d[index0i]=ldcomplex_add(ya[index0i],ldcomplex_mul(tempdc,d[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. dipole_polarisation=dipole_polarisation+eta*d */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            dipole_polarisation[index0i]=ldcomplex_add(dipole_polarisation[index0i],ldcomplex_mul(eta,d[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. kappa=tau*sqrt(m+1) */
         kappa=tau*sqrtl((long double)m+1.0L);
         
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 20. Check stopping criterion */
         if(kappa<error){
/* r=incident_E_field-interaction_matrix*dipole_polarisation)*/
            dftmatvec_MPI(dipole_polarisation);
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=vc*parallel.alloc_padded;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  r[array2]=ldcomplex_sub(incident_E_field[array2],zero_padded_vector[array1+index0i]);
               }
            }

            if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* r=M^{-1}r */
               if(iterative.precond==1){
                  /* Point-Jacobi Preconditioning (Divide by the diagonal) */
                  for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
                     r[index0i]=ldcomplex_mul(r[index0i],point_jacobi[index0i]);
                  }
               }
            }

            norm1=ldcomplex_2vectornorm(r);
            ratio=norm1/norm0;
            if(parallel.myrank==0){ /* Restrict to master */
               fprintf(output.iter,"  At iteration %d the convergence ratio is %.*Lg and the residual norm is %.*Lg\n",iteration,LDBLP,ratio,LDBLP,norm1);
            }
            if(ratio<(terminate+LDBL_EPSILON)){
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
      reduce=zero;
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         reduce=ldcomplex_add(reduce,ldcomplex_mul(rtilde[index0i],omega[index0i]));
      }
      MPI_Allreduce(&reduce,&rhok,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);
      /* Check for iterative.breakdown */
      if(ldcomplex_abs(rho)<iterative.breakdown){
         print_error("Iterative scheme error","TFQMR() iterative.breakdown","|rho|<iterative.breakdown"); 
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 22. beta=rhok/rho   */
      beta=ldcomplex_div(rhok,rho);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 23. yb=omega+beta*ya */
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         yb[index0i]=ldcomplex_add(omega[index0i],ldcomplex_mul(beta,ya[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 24. v=interaction_matrix*yb+beta*(interaction_matrix*ya+beta*v) */
      dftmatvec_MPI(ya);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*ya=M^{-1}interaction_matrix*ya */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=vc*parallel.alloc_padded;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=ldcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array1+index0i;
            zero_padded_vector[array2]=ldcomplex_add(zero_padded_vector[array2],ldcomplex_mul(beta,v[array0+index0i]));
         }
      }
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            v[array0+index0i]=ldcomplex_mul(beta,zero_padded_vector[array1+index0i]);
         }
      }

      dftmatvec_MPI(yb);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*yb=M^{-1}interaction_matrix*yb */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=vc*parallel.alloc_padded;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=ldcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            v[array2]=ldcomplex_add(v[array2],zero_padded_vector[array1+index0i]);
         }
      }
      iteration++;
   } /* End while loop */
   if(timing.enabled){ /* Finalise timing */
      finish=MPI_Wtime();
      temp=finish-start;
      result=0.0;
      MPI_Reduce(&temp,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      if(parallel.myrank==0){ /* Restrict to master */   
         fprintf(output.iter,"\n  Time for convergence [Averaged over the %d processors]: %.*gs\n\n",parallel.np,DBLP,(result/parallel.dnp));
      }
      timing.iterative_solution+=result;
   }
}
