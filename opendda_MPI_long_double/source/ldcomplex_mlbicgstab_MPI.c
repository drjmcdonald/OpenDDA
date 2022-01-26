/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   A variant of the stabilised version of the
   BiConjugate-Gradients (BiCGSTAB) based on
   multiple Lanczos starting vectors: ML(n)BiCGSTAB_{new}

   ML(n)BiCGSTAB: reformulations and implementations,
   Yeung, M. and Boley, D., In press.
   Appendix. Algorithm 2.2' ML(n)BiCGSTAB with iterative.preconditioning.

   Uses SIMD oriented Fast Mersenne Twister RNG for initialisation
   of the q vectors:
   Mutsuo Saito, Makoto Matsumoto, Hiroshima University.

   Uses Modified Gram-Schmidt Orthogonalisation for the randomly
   generated q vectors:
   Matrix Computations, Golub, G.H. and Van Loan, C.F.,
   The Johns Hopkins University Press, Baltimore,
   3rd Ed., 1996, ISBN 0-8018-5414-8

   i,ii,s,t,vc,index0i,array*: Loop and array index control variables
   iteration: Iteration number
   iterative.precond: 0: None
                      1: Point-Jacobi Preconditioning
   r[3],q[iterative.vec][3],g[iterative.vec][3],
   gtilde[3],u[3],utilde[3],
   d[iterative.vec-1][3],omega[iterative.vec][3],
   zd[3],zg[3],zomega[3],c: Local vectors for iterative scheme
   alpha,rho,minusone,
   tempdc,tempd,tempe,
   tempf,beta,temp,result,element: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   reduce: Used to facilitate the MPI_Allreduce function
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "ldcomplex_mlbicgstab_MPI.h"

void ldcomplex_mlbicgstab_MPI(void){

   int i,ii,s,t,vc,iteration=0;
   long int index0i,array0,array1,array2;
   ldcomplex alpha,rho,tempdc,tempd,tempe,tempf,beta,reduce;
   long double ratio,norm0,norm1,terminate;
   double temp,result,start,finish;

   if(iterative.vec>50){
      print_error("Iterative scheme error","MLBiCGSTAB() error","Too many starting vectors");
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=MPI_Wtime();}

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 1. Fill q vectors with random numbers */
   /* Mersenne Twister RNG Initialisation is done by passing an array of seeds
      which are generated using `lrand48()' which is in turn seeded from the clock */
   srand48((long int)time(NULL)); 
   for(i=0;i<iterative.vec;i++){
      init[i]=(unsigned int)lrand48();
   }
   init_by_array(init,(unsigned int)iterative.vec);
   for(s=0;s<iterative.vec;s++){
      for(vc=0;vc<3;vc++){
         array0=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            q[array2].dat[0]=(long double)(genrand_real1()-genrand_real1());
            q[array2].dat[1]=(long double)(genrand_real1()-genrand_real1());            
         }
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Modified Gram-Schmidt Orthogonalisation */
   for(s=0;s<iterative.vec;s++){
      norm0=1.0L/ldcomplex_2vectornorm(&q[s*parallel.alloc_vector3]);
      for(vc=0;vc<3;vc++){ /* q[s]=q[s]/||q[s]|| */
         array0=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            q[array2]=ldcomplex_scale(q[array2],norm0);
         }
      }
      for(t=s+1;t<iterative.vec;t++){
         reduce=zero;
         for(vc=0;vc<3;vc++){ /* q[s]^{*}q[t] */
            array0=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
            array1=t*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[array0+index0i]),q[array1+index0i]));
            }
         }
         MPI_Allreduce(&reduce,&tempdc,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);
         for(vc=0;vc<3;vc++){ /* q[t]=q[t]-(q[s]^{*}q[t])q[s] */
            array0=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
            array1=t*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array1+index0i;
               q[array2]=ldcomplex_sub(q[array2],ldcomplex_mul(tempdc,q[array0+index0i]));
            }
         }
      }
   }

/*   printf("Check dot products are 0 for s!=t and 1 for s==t and check ||q[s]||=1\n");
   for(s=0;s<iterative.vec;s++){
      norm0=ldcomplex_2vectornorm(&q[s]);
      fprintf(stdout,"norm for q[%d] = %lf\n",s,norm0);
      for(t=0;t<iterative.vec;t++){
         reduce=zero;
         for(vc=0;vc<3;vc++){
            array0=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
            array1=t*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[array0+index0i]),q[array1+index0i]));
            }
         }
         MPI_Allreduce(&reduce,&tempdc,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);   
         if(parallel.myrank==0){printf("vc=%d\ts=%d\tii=%d\t%+.*Lg%+.*Lgi\n",vc,s,t,LDBLP,tempdc.dat[0],LDBLP,tempdc.dat[1]);}
      }
   } */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. r=incident_E_field-interaction_matrix*dipole_polarisation */
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

   norm0=ldcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0L; 
   };

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. g[iterative.vec-1]=r */
   memcpy(&g[(iterative.vec-1)*parallel.alloc_vector3],r,(size_t)parallel.alloc_vector3*sizeof(ldcomplex));

   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. gtilde=M^{-1}r */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            gtilde[index0i]=ldcomplex_mul(r[index0i],point_jacobi[index0i]);
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. omega[iterative.vec-1]=interaction_matrix*gtilde */
      dftmatvec_MPI(gtilde);
      for(vc=0;vc<3;vc++){
         array0=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            omega[array0+index0i]=zero_padded_vector[array1+index0i];
         }
      }
   }
   else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. omega[iterative.vec-1]=interaction_matrix*r */
      dftmatvec_MPI(r);
      for(vc=0;vc<3;vc++){
         array0=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            omega[array0+index0i]=zero_padded_vector[array1+index0i];
         }
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. c[iterative.vec-1]=q[0]^{H}omega[iterative.vec-1] */
   reduce=zero;
   for(vc=0;vc<3;vc++){
      array0=vc*parallel.alloc_vector;
      array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
      for(index0i=0;index0i<parallel.alloc_vector;index0i++){
         reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[array0+index0i]),omega[array1+index0i]));
      }
   }
   MPI_Allreduce(&reduce,&c[iterative.vec-1],1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* tempe=q[0]^{H}r */
   reduce=zero;
   for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
      reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[index0i]),r[index0i]));      
   }
   MPI_Allreduce(&reduce,&tempe,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);

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
/* 7. alpha=(q[0]^{H}r)/(c[iterative.vec-1]} */
      /* Check for iterative.breakdown */
      if(ldcomplex_abs(c[iterative.vec-1])<iterative.breakdown){
         print_error("Iterative scheme error","MLBiCGSTAB() iterative.breakdown","|c[iterative.vec-1]|<iterative.breakdown");
      }
      alpha=ldcomplex_div(tempe,c[iterative.vec-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. u=r-alpha*omega[iterative.vec-1] */
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            u[array2]=ldcomplex_sub(r[array2],ldcomplex_mul(alpha,omega[array1+index0i]));
         }
      }

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. utilde=M^{-1}u */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
               utilde[index0i]=ldcomplex_mul(u[index0i],point_jacobi[index0i]);
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. rho=(-[interaction_matrix*utilde]^{H}u]/(||interaction_matrix*utilde||^2)
      interaction_matrix*utilde */
         dftmatvec_MPI(utilde);
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=vc*parallel.alloc_padded;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               zd[array0+index0i]=zero_padded_vector[array1+index0i]; /* Using zd as a scratch space here */
            }
         }

         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(zd[index0i]),zd[index0i]));
         }
         MPI_Allreduce(&reduce,&tempdc,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);
         /* Check for iterative.breakdown */
         if(ldcomplex_abs(tempdc)<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB() iterative.breakdown","||interaction_matrix*utilde||<iterative.breakdown"); 
         }

/* rho=[interaction_matrix*utilde]^{H}u */
         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(zd[index0i]),u[index0i]));
         }
         MPI_Allreduce(&reduce,&rho,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);

/* rho=-[interaction_matrix*utilde]^{H}u */
         rho=ldcomplex_mul(minusone,rho);

/* rho=(-[interaction_matrix*utilde]^{H}u]/(||interaction_matrix*utilde||^2) */
         rho=ldcomplex_div(rho,tempdc);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. dipole_polarisation=dipole_polarisation-rho*utilde+alpha*gtilde */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            dipole_polarisation[index0i]=ldcomplex_sub(dipole_polarisation[index0i],ldcomplex_mul(rho,utilde[index0i]));
            dipole_polarisation[index0i]=ldcomplex_add(dipole_polarisation[index0i],ldcomplex_mul(alpha,gtilde[index0i]));            
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. r=rho*interaction_matrix*utilde+u */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            r[index0i]=ldcomplex_add(ldcomplex_mul(rho,zd[index0i]),u[index0i]);
         }
      }
      else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. rho=(-[interaction_matrix*u]^{H}u]/(||interaction_matrix*u||^2)
      interaction_matrix*u */
         dftmatvec_MPI(u);
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=vc*parallel.alloc_padded;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               zd[array0+index0i]=zero_padded_vector[array1+index0i]; /* Using zd as a scratch space here */
            }
         }

         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(zd[index0i]),zd[index0i]));
         }
         MPI_Allreduce(&reduce,&tempdc,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);
         /* Check for iterative.breakdown */
         if(ldcomplex_abs(tempdc)<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB() iterative.breakdown","||A*u||<iterative.breakdown");
         }

/* rho=[interaction_matrix*u]^{H}u */
         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(zd[index0i]),u[index0i]));
         }
         MPI_Allreduce(&reduce,&rho,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);

/* rho=-[interaction_matrix*u]^{H}u */
         rho=ldcomplex_mul(minusone,rho);

/* rho=(-[interaction_matrix*u]^{H}u]/(||interaction_matrix*u||^2) */
         rho=ldcomplex_div(rho,tempdc);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. dipole_polarisation=dipole_polarisation-rho*u+alpha*g[iterative.vec-1] */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               dipole_polarisation[array2]=ldcomplex_sub(dipole_polarisation[array2],ldcomplex_mul(rho,u[array2]));
               dipole_polarisation[array2]=ldcomplex_add(dipole_polarisation[array2],ldcomplex_mul(alpha,g[array1+index0i]));               
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. r=rho*interaction_matrix*u+u */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            r[index0i]=ldcomplex_add(ldcomplex_mul(rho,zd[index0i]),u[index0i]);
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Check stopping criterion */
      norm1=ldcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      if(parallel.myrank==0){ /* Restrict to master */
         fprintf(output.iter,"  At iteration %d the convergence ratio is %.*Lg and the residual norm is %.*Lg\n",iteration,LDBLP,ratio,LDBLP,norm1);
      }
      if(ratio<(terminate+LDBL_EPSILON)){
         break;
      }

      tempd=ldcomplex_mul(rho,c[iterative.vec-1]);
      /* Check for iterative.breakdown */
      if(ldcomplex_abs(tempd)<iterative.breakdown){
         print_error("Iterative scheme error","MLBiCGSTAB() iterative.breakdown","|rho*c[iterative.vec-1]|<iterative.breakdown");
      }

      for(ii=1;ii<iterative.vec;ii++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. zd=u,zg=0,zomega=0 */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            zg[index0i]=zero;
         }
         memcpy(zd,u,(size_t)parallel.alloc_vector3*sizeof(ldcomplex));
         memcpy(zomega,zg,(size_t)parallel.alloc_vector3*sizeof(ldcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* tempf=q[ii]^{H}u */
         reduce=zero;
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=ii*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[array1+index0i]),u[array0+index0i]));
            }
         }
         MPI_Allreduce(&reduce,&tempf,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);

         if(iteration>0){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* beta=(-q[ii]^{H}zd)(c[ii-1]) */
            /* Check for iterative.breakdown */
            if(ldcomplex_abs(c[ii-1])<iterative.breakdown){
               print_error("Iterative scheme error","MLBiCGSTAB() iterative.breakdown","|c[ii-1]|<iterative.breakdown");
            }
            beta=ldcomplex_div(ldcomplex_mul(minusone,tempf),c[ii-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* zd=zd+beta*d[ii-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  zd[array2]=ldcomplex_add(zd[array2],ldcomplex_mul(beta,d[array1+index0i]));
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* zg=zg+beta*g[ii-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  zg[array2]=ldcomplex_add(zg[array2],ldcomplex_mul(beta,g[array1+index0i]));
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* zomega=zomega+beta*omega[ii-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  zomega[array2]=ldcomplex_add(zomega[array2],ldcomplex_mul(beta,omega[array1+index0i]));
               }
            }

            for(s=ii+1;s<iterative.vec;s++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. beta=(-q[s]^{H}zd)/(c[s-1]) */
               reduce=zero;
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[array1+index0i]),zd[array0+index0i]));
                  }
               }
               MPI_Allreduce(&reduce,&beta,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);   
               /* Check for iterative.breakdown */
               if(ldcomplex_abs(c[s-1])<iterative.breakdown){
                  print_error("Iterative scheme error","MLBiCGSTAB() iterative.breakdown","|c[s-1]|<iterative.breakdown");
               }
               beta=ldcomplex_div(ldcomplex_mul(minusone,beta),c[s-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. zd=zd+beta*d[s-1] */
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     array2=array0+index0i;
                     zd[array2]=ldcomplex_add(zd[array2],ldcomplex_mul(beta,d[array1+index0i]));
                  }
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. zg=zg+beta*g[s-1] */
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     array2=array0+index0i;
                     zg[array2]=ldcomplex_add(zg[array2],ldcomplex_mul(beta,g[array1+index0i]));
                  }
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. zomega=zomega+beta*omega[s-1] */
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     array2=array0+index0i;
                     zomega[array2]=ldcomplex_add(zomega[array2],ldcomplex_mul(beta,omega[array1+index0i]));
                  }
               }
            } /* End s for loop */
         } /* End if iteration>0 */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. zomega=r+rho*zomega */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            zomega[index0i]=ldcomplex_add(r[index0i],ldcomplex_mul(rho,zomega[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. beta=(-q[0]^{H}zomega)(rho*c[iterative.vec-1]) */
         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[index0i]),zomega[index0i]));
         }
         MPI_Allreduce(&reduce,&beta,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);
         beta=ldcomplex_div(ldcomplex_mul(minusone,beta),tempd);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 20. zomega=zomega+rho*beta*omega[iterative.vec-1] */
         tempdc=ldcomplex_mul(rho,beta);
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               zomega[array2]=ldcomplex_add(zomega[array2],ldcomplex_mul(tempdc,omega[array1+index0i]));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 21. zg=zg+zomega+beta*g[iterative.vec-1] */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               zg[array2]=ldcomplex_add(zg[array2],ldcomplex_add(zomega[array2],ldcomplex_mul(beta,g[array1+index0i])));
            }
         }

         for(s=1;s<ii;s++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 22. beta=(-q[s]^{H}zomega)/(c[s-1]) */
            reduce=zero;
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[array1+index0i]),zomega[array0+index0i]));
               }
            }
            MPI_Allreduce(&reduce,&beta,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);
            /* Check for iterative.breakdown */
            if(ldcomplex_abs(c[s-1])<iterative.breakdown){
               print_error("Iterative scheme error","MLBiCGSTAB() iterative.breakdown","|c[s-1]|<iterative.breakdown");
            }
            beta=ldcomplex_div(ldcomplex_mul(minusone,beta),c[s-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 23. zg=zg+beta*g[s-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  zg[array2]=ldcomplex_add(zg[array2],ldcomplex_mul(beta,g[array1+index0i]));
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 24. zomega=zomega+beta*d[s-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  zomega[array2]=ldcomplex_add(zomega[array2],ldcomplex_mul(beta,d[array1+index0i]));
               }
            }
         } /* End s for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 25. d[ii-1]=zomega-u */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               d[array1+index0i]=ldcomplex_sub(zomega[array2],u[array2]);
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 26. c[ii-1]=q[ii]^{H}d[ii-1] */
         reduce=zero;
         for(vc=0;vc<3;vc++){
            array0=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            array1=ii*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[array1+index0i]),d[array0+index0i]));
            }
         }
         MPI_Allreduce(&reduce,&c[ii-1],1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 27. g[ii-1]=zg */
         memcpy(&g[(ii-1)*parallel.alloc_vector3],zg,(size_t)parallel.alloc_vector3*sizeof(ldcomplex));

         if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 28. gtilde=M^{-1}zg */
            if(iterative.precond==1){
               /* Point-Jacobi Preconditioning (Divide by the diagonal) */
               for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
                  gtilde[index0i]=ldcomplex_mul(zg[index0i],point_jacobi[index0i]);
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 29. omega[ii-1]=interaction_matrix*gtilde */
            dftmatvec_MPI(gtilde);
            for(vc=0;vc<3;vc++){
               array0=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               array1=vc*parallel.alloc_padded;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  omega[array0+index0i]=zero_padded_vector[array1+index0i];
               }
            }
         }
         else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 29. omega[ii-1]=interaction_matrix*zg */
            dftmatvec_MPI(zg);
            for(vc=0;vc<3;vc++){
               array0=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               array1=vc*parallel.alloc_padded;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  omega[array0+index0i]=zero_padded_vector[array1+index0i];
               }
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 30. alpha=(q[ii]^{H}u)/(c[ii-1]) */
         /* Check for iterative.breakdown */
         if(ldcomplex_abs(c[ii-1])<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB() iterative.breakdown","|c[ii-1]|<iterative.breakdown");
         }
         alpha=ldcomplex_div(tempf,c[ii-1]);

         if(ii<iterative.vec-1){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 31. u=u-alpha*d[ii-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  u[array2]=ldcomplex_sub(u[array2],ldcomplex_mul(alpha,d[array1+index0i]));
               }
            }
         } /* End if ii<iterative.vec-1 */

         if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 32. dipole_polarisation=dipole_polarisation+rho*alpha*gtilde */
            tempdc=ldcomplex_mul(rho,alpha);
            for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
               dipole_polarisation[index0i]=ldcomplex_add(dipole_polarisation[index0i],ldcomplex_mul(tempdc,gtilde[index0i]));
            }
         }
         else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 32. dipole_polarisation=dipole_polarisation+rho*alpha*zg */
            tempdc=ldcomplex_mul(rho,alpha);
            for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
               dipole_polarisation[index0i]=ldcomplex_add(dipole_polarisation[index0i],ldcomplex_mul(tempdc,zg[index0i]));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 33. r=r-rho*alpha*omega[ii-1] */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               r[array2]=ldcomplex_sub(r[array2],ldcomplex_mul(tempdc,omega[array1+index0i]));
            }
         }
      } /* End ii for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* tempe=q[0]^{H}r */
      reduce=zero;
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[index0i]),r[index0i]));
      }
      MPI_Allreduce(&reduce,&tempe,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);
/* 34. beta=(-q[0]^{H}r)/(rho*c[iterative.vec-1]) */
      beta=ldcomplex_div(ldcomplex_mul(minusone,tempe),tempd);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 35. zomega=r+rho*beta*omega[iterative.vec-1] */
      tempdc=ldcomplex_mul(rho,beta);
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            zomega[array2]=ldcomplex_add(r[array2],ldcomplex_mul(tempdc,omega[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 36. zg=zomega+beta*g[iterative.vec-1] */
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            zg[array2]=ldcomplex_add(zomega[array2],ldcomplex_mul(beta,g[array1+index0i]));
         }
      }

      for(s=1;s<iterative.vec;s++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 37. beta=(-q[s]^{H}zomega)/(c[s-1]) */
         reduce=zero;
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[array1+index0i]),zomega[array0+index0i]));
            }
         }
         MPI_Allreduce(&reduce,&beta,1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);
         /* Check for iterative.breakdown */
         if(ldcomplex_abs(c[s-1])<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB() iterative.breakdown","|c[s-1]|<iterative.breakdown");
         }
         beta=ldcomplex_div(ldcomplex_mul(minusone,beta),c[s-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 38. zg=zg+beta*g[s-1] */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               zg[array2]=ldcomplex_add(zg[array2],ldcomplex_mul(beta,g[array1+index0i]));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 39. zomega=zomega+beta*d[s-1] */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               zomega[array2]=ldcomplex_add(zomega[array2],ldcomplex_mul(beta,d[array1+index0i]));
            }
         }
      } /* End s for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 40. g[iterative.vec-1]=zg */
      memcpy(&g[(iterative.vec-1)*parallel.alloc_vector3],zg,(size_t)parallel.alloc_vector3*sizeof(ldcomplex));

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 41. gtilde=M^{-1}zg */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
               gtilde[index0i]=ldcomplex_mul(zg[index0i],point_jacobi[index0i]);
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 42. omega[iterative.vec-1]=interaction_matrix*gtilde */
         dftmatvec_MPI(gtilde);
         for(vc=0;vc<3;vc++){
            array0=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            array1=vc*parallel.alloc_padded;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               omega[array0+index0i]=zero_padded_vector[array1+index0i];
            }
         }
      }
      else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 42. omega[iterative.vec-1]=interaction_matrix*zg */
         dftmatvec_MPI(zg);
         for(vc=0;vc<3;vc++){
            array0=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            array1=vc*parallel.alloc_padded;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               omega[array0+index0i]=zero_padded_vector[array1+index0i];
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 43. c[iterative.vec-1]=q[0]^{H}omega[iterative.vec-1] */
      reduce=zero;
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            reduce=ldcomplex_add(reduce,ldcomplex_mul(ldcomplex_conj(q[array0+index0i]),omega[array1+index0i]));
         }
      }
      MPI_Allreduce(&reduce,&c[iterative.vec-1],1,ldcomplex_type,add_ldcomplex,MPI_COMM_WORLD);
      iteration++;
   } /* while(iteration<iterative.maximum) */
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
