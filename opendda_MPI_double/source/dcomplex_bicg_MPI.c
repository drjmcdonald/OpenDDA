/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   BiConjugate-Gradients: BiCG

   Conjugate gradient methods for indefinite systems, Fletcher, R.,
   Proceedings of the Dundee Biennial Conference on Numerical Analysis,
   Ed. Watson, G.~A., Pages 73-89, 1975.

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
   ptilde[3],q[3]: Local vectors for iterative scheme
   rho,rhok,alpha,beta,temp,result,element: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   reduce: Used to facilitate the MPI_Allreduce function
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dcomplex_bicg_MPI.h"

void dcomplex_bicg_MPI(void){

   int vc,iteration=0;
   long int index0i,array0,array1,array2;
   dcomplex rho,rhok,alpha,beta,reduce;
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

   while(iteration<iterative.maximum){ /* 17. Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. rhok=r^{T}rtilde */
      reduce=zero;
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         reduce=dcomplex_add(reduce,dcomplex_mul(r[index0i],rtilde[index0i]));
      }
      MPI_Allreduce(&reduce,&rhok,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

      if(iteration==0){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. p=r */
         memcpy(p,r,(size_t)parallel.alloc_vector3*sizeof(dcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. ptilde=rtilde */
         memcpy(ptilde,rtilde,(size_t)parallel.alloc_vector3*sizeof(dcomplex));
      }
      else{
         /* Check for iterative.breakdown */
         if((dcomplex_abs(rho))<iterative.breakdown){
            print_error("Iterative scheme error","BiCG() iterative.breakdown","|rho|<iterative.breakdown"); 
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. beta=rhok/rho */
         beta=dcomplex_div(rhok,rho);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. p=r+beta*p */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            p[index0i]=dcomplex_add(r[index0i],dcomplex_mul(beta,p[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. ptilde=rtilde+beta*ptilde */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            ptilde[index0i]=dcomplex_add(rtilde[index0i],dcomplex_mul(beta,ptilde[index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. q=interaction_matrix*p */
      dftmatvec_MPI(p);
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            q[array0+index0i]=zero_padded_vector[array1+index0i];
         }
      }

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* q=M^{-1}q */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            q[index0i]=dcomplex_mul(q[index0i],point_jacobi[index0i]);
         }
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. qtilde=interaction_matrix^{T}ptilde */
/* qtilde is stored in zero_padded_vector */
      dftmatvec_MPI(ptilde);

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix^{T}ptilde=M^{-1}interaction_matrix^{T}ptilde */
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
/* 13. alpha=rhok/(ptilde^{T}q) */
      reduce=zero;
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         reduce=dcomplex_add(reduce,dcomplex_mul(ptilde[index0i],q[index0i]));
      }
      MPI_Allreduce(&reduce,&alpha,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
      /* Check for iterative.breakdown */
      if((dcomplex_abs(alpha))<iterative.breakdown){
         print_error("Iterative scheme error","BiCG() iterative.breakdown","|alpha|<iterative.breakdown"); 
      }
      alpha=dcomplex_div(rhok,alpha);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. dipole_polarisation=dipole_polarisation+alpha*p */
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         dipole_polarisation[index0i]=dcomplex_add(dipole_polarisation[index0i],dcomplex_mul(alpha,p[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. r=r-alpha*q */
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         r[index0i]=dcomplex_sub(r[index0i],dcomplex_mul(alpha,q[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. rtilde=rtilde-alpha*qtilde */
/* qtilde is stored in zero_padded_vector */
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            rtilde[array2]=dcomplex_sub(rtilde[array2],dcomplex_mul(alpha,zero_padded_vector[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. Check stopping criterion */
      norm1=dcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      if(parallel.myrank==0){ /* Restrict to master */
         fprintf(output.iter,"  At iteration %d the convergence ratio is %.*g and the residual norm is %.*g\n",iteration,DBLP,ratio,DBLP,norm1);
      }
      if(ratio<(terminate+DBL_EPSILON)){
         break;
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
