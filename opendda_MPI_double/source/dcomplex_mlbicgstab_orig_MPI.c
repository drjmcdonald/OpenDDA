/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   A variant of the stabilised version of the
   BiConjugate-Gradients (BiCGSTAB) based on
   multiple Lanczos starting vectors: ML(n)BiCGSTAB_{orig}

   ML(k)BiCGSTAB: A BiCGSTAB variant based on multiple Lanczos
   starting vectors, Yeung, M. and Chan, T.F., SIAM J. Sci. Comput.,
   Vol. 21, No. 4, Pages 1263–1290, 1999, Algorithm 3.

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
   r[3],g[iterative.vec][3],gtilde[3],
   d[iterative.vec-1][3],omega[iterative.vec][3],
   q[iterative.vec][3],u[3],utilde[3],
   zd[3],zg[3],zomega[3],c: Local vectors for iterative scheme
   alpha,rho,beta,minusone,
   tempdc,temp,result,element: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   reduce: Used to facilitate the MPI_Allreduce function
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dcomplex_mlbicgstab_orig_MPI.h"

void dcomplex_mlbicgstab_orig_MPI(void){

   int i,ii,s,t,vc,iteration=0;
   long int index0i,array0,array1,array2;
   dcomplex alpha,rho,beta,tempdc,reduce;
   double ratio,norm0,norm1,temp,result,terminate,start,finish;

   if(iterative.vec>50){
      print_error("Iterative scheme error","MLBiCGSTAB_orig() error","Too many starting vectors");
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
            q[array2].dat[0]=genrand_real1()-genrand_real1();
            q[array2].dat[1]=genrand_real1()-genrand_real1();            
         }
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Modified Gram-Schmidt Orthogonalisation */
   for(s=0;s<iterative.vec;s++){
      norm0=1.0/dcomplex_2vectornorm(&q[s*parallel.alloc_vector3]);
      for(vc=0;vc<3;vc++){ /* q[s]=q[s]/||q[s]|| */
         array0=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            q[array2]=dcomplex_scale(q[array2],norm0);
         }
      }
      for(t=s+1;t<iterative.vec;t++){
         reduce=zero;
         for(vc=0;vc<3;vc++){ /* q[s]^{*}q[t] */
            array0=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
            array1=t*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array0+index0i]),q[array1+index0i]));
            }
         }
         MPI_Allreduce(&reduce,&tempdc,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
         for(vc=0;vc<3;vc++){ /* q[t]=q[t]-(q[s]^{*}q[t])q[s] */
            array0=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
            array1=t*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array1+index0i;
               q[array2]=dcomplex_sub(q[array2],dcomplex_mul(tempdc,q[array0+index0i]));
            }
         }
      }
   }

/*   printf("Check dot products are 0 for s!=t and 1 for s==t and check ||q[s]||=1\n");
   for(s=0;s<iterative.vec;s++){
      norm0=dcomplex_2vectornorm(&q[s]);
      fprintf(stdout,"norm for q[%d] = %lf\n",s,norm0);
      for(t=0;t<iterative.vec;t++){
         reduce=zero;
         for(vc=0;vc<3;vc++){
            array0=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
            array1=t*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array0+index0i]),q[array1+index0i]));
            }
         }
         MPI_Allreduce(&reduce,&tempdc,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);   
         if(parallel.myrank==0){printf("vc=%d\ts=%d\tii=%d\t%+.*g%+.*gi\n",vc,s,t,DBLP,tempdc.dat[0],DBLP,tempdc.dat[1]);}
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
            r[array2]=dcomplex_sub(incident_E_field[array2],zero_padded_vector[array1+index0i]);
         }
      }
   }
   else{ /* r=incident_E_field */
      memcpy(r,incident_E_field,(size_t)parallel.alloc_vector3*sizeof(dcomplex));
   }

   norm0=dcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0; 
   };

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. g[iterative.vec-1]=r */
   memcpy(&g[(iterative.vec-1)*parallel.alloc_vector3],r,(size_t)parallel.alloc_vector3*sizeof(dcomplex));
         
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

   while(iteration<iterative.maximum){ /* Check stopping criterion */
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. gtilde=M^{-1}g[iterative.vec-1] */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  gtilde[array2]=dcomplex_mul(g[array1+index0i],point_jacobi[array2]);
               }
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
/* 5. omega[iterative.vec-1]=interaction_matrix*g[iterative.vec-1] */
         dftmatvec_MPI(&g[(iterative.vec-1)*parallel.alloc_vector3]);
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
            reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array0+index0i]),omega[array1+index0i]));
         }
      }
      MPI_Allreduce(&reduce,&c[iterative.vec-1],1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. alpha=(q[0]^{H}r)/(c[iterative.vec-1]) */
      reduce=zero;
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[index0i]),r[index0i]));
      }
      MPI_Allreduce(&reduce,&alpha,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
      /* Check for iterative.breakdown */
      if(dcomplex_abs(c[iterative.vec-1])<iterative.breakdown){
         print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","|c[iterative.vec-1]|<iterative.breakdown"); 
      }
      alpha=dcomplex_div(alpha,c[iterative.vec-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. u=r-alpha*omega[iterative.vec-1] */
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            u[array2]=dcomplex_sub(r[array2],dcomplex_mul(alpha,omega[array1+index0i]));
         }
      }

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. utilde=M^{-1}u */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
               utilde[index0i]=dcomplex_mul(u[index0i],point_jacobi[index0i]);
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. rho=(-[interaction_matrix*utilde]^{H}u)/(||interaction_matrix*utilde||^2)
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
            reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(zd[index0i]),zd[index0i]));
         }
         MPI_Allreduce(&reduce,&tempdc,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
         /* Check for iterative.breakdown */
         if(dcomplex_abs(tempdc)<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","||interaction_matrix*utilde||<iterative.breakdown");
         }

/* rho=[interaction_matrix*utilde]^{H}u */
         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(zd[index0i]),u[index0i]));
         }
         MPI_Allreduce(&reduce,&rho,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

/* rho=-[interaction_matrix*utilde]^{H}u */
         rho=dcomplex_mul(minusone,rho);

/* rho=(-[interaction_matrix*utilde]^{H}u)/(||interaction_matrix*utilde||^2) */
         rho=dcomplex_div(rho,tempdc);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. dipole_polarisation=dipole_polarisation-rho*utilde+alpha*gtilde */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            dipole_polarisation[index0i]=dcomplex_sub(dipole_polarisation[index0i],dcomplex_mul(rho,utilde[index0i]));
            dipole_polarisation[index0i]=dcomplex_add(dipole_polarisation[index0i],dcomplex_mul(alpha,gtilde[index0i]));   
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. r=rho*interaction_matrix*utilde+u */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            r[index0i]=dcomplex_add(dcomplex_mul(rho,zd[index0i]),u[index0i]);
         }
      }
      else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. rho=(-[interaction_matrix*u]^{H}u)/(||interaction_matrix*u||^2)
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
            reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(zd[index0i]),zd[index0i]));
         }
         MPI_Allreduce(&reduce,&tempdc,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
         /* Check for iterative.breakdown */
         if(dcomplex_abs(tempdc)<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","||interaction_matrix*u||<iterative.breakdown"); 
         }

/* rho=[interaction_matrix*u]^{H}u */
         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(zd[index0i]),u[index0i]));
         }
         MPI_Allreduce(&reduce,&rho,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

/* rho=-[interaction_matrix*u]^{H}u */
         rho=dcomplex_mul(minusone,rho);

/* rho=(-[interaction_matrix*u]^{H}u)/(||interaction_matrix*u||^2) */
         rho=dcomplex_div(rho,tempdc);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. dipole_polarisation=dipole_polarisation-rho*u+alpha*g[iterative.vec-1] */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               dipole_polarisation[array2]=dcomplex_sub(dipole_polarisation[array2],dcomplex_mul(rho,u[array2]));
               dipole_polarisation[array2]=dcomplex_add(dipole_polarisation[array2],dcomplex_mul(alpha,g[array1+index0i]));               
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. r=rho*interaction_matrix*u+u */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            r[index0i]=dcomplex_add(dcomplex_mul(rho,zd[index0i]),u[index0i]);
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Check stopping criterion */
      norm1=dcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      if(parallel.myrank==0){ /* Restrict to master */
         fprintf(output.iter,"  At iteration %d the convergence ratio is %.*g and the residual norm is %.*g\n",iteration,DBLP,ratio,DBLP,norm1);
      }
      if(ratio<(terminate+DBL_EPSILON)){
         break;
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Start loop over ii */
      for(ii=1;ii<=iterative.vec;ii++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13.  zd=u,zg=r,zomega=0 */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            zomega[index0i]=zero;
         }
         memcpy(zd,u,(size_t)parallel.alloc_vector3*sizeof(dcomplex));
         memcpy(zg,r,(size_t)parallel.alloc_vector3*sizeof(dcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Start loop over s and check if iteration>0 */
         if(iteration>0){
            for(s=ii;s<iterative.vec;s++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. beta=(-q[s]^{H}zd)/(c[s-1]) */
               reduce=zero;
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array1+index0i]),zd[array0+index0i]));
                  }
               }
               MPI_Allreduce(&reduce,&beta,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
               /* Check for iterative.breakdown */
               if(dcomplex_abs(c[s-1])<iterative.breakdown){
                  print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","|c[s-1]|<iterative.breakdown");
               }
               beta=dcomplex_div(dcomplex_mul(minusone,beta),c[s-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. zd=zd+beta*d[s-1] */
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     array2=array0+index0i;
                     zd[array2]=dcomplex_add(zd[array2],dcomplex_mul(beta,d[array1+index0i]));
                  }
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. zg=zg+beta*g[s-1] */
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     array2=array0+index0i;
                     zg[array2]=dcomplex_add(zg[array2],dcomplex_mul(beta,g[array1+index0i]));
                  }
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. zomega=zomega+beta*omega[s-1] */
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     array2=array0+index0i;
                     zomega[array2]=dcomplex_add(zomega[array2],dcomplex_mul(beta,omega[array1+index0i]));
                  }
               }
            } /* End s for loop */
         } /* End if iteration>0 */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. beta=(-q[0]^{H}(r+rho*zomega))/(rho*c[iterative.vec-1]) */
         tempdc=dcomplex_mul(rho,c[iterative.vec-1]);
         /* Check for iterative.breakdown */
         if(dcomplex_abs(tempdc)<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","|rho*c[iterative.vec-1]|<iterative.breakdown"); 
         }
         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[index0i]),dcomplex_add(r[index0i],dcomplex_mul(rho,zomega[index0i]))));
         }
         MPI_Allreduce(&reduce,&beta,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
         beta=dcomplex_div(dcomplex_mul(minusone,beta),tempdc);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. zg=zg+beta*g[iterative.vec-1] */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               zg[array2]=dcomplex_add(zg[array2],dcomplex_mul(beta,g[array1+index0i]));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 20. zomega=rho(zomega+beta*omega[iterative.vec-1]) */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               zomega[array2]=dcomplex_mul(rho,dcomplex_add(zomega[array2],dcomplex_mul(beta,omega[array1+index0i])));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 21. zd=r+zomega */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            zd[index0i]=dcomplex_add(r[index0i],zomega[index0i]);
         }

         for(s=1;s<ii;s++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 22. beta=(-q[s]^{H}zd)/(c[s-1]) */
            reduce=zero;
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array1+index0i]),zd[array0+index0i]));
               }
            }
            MPI_Allreduce(&reduce,&beta,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);   
            /* Check for iterative.breakdown */
            if(dcomplex_abs(c[s-1])<iterative.breakdown){
               print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","|c[s-1]|<iterative.breakdown");
            }
            beta=dcomplex_div(dcomplex_mul(minusone,beta),c[s-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 23. zd=zd+beta*d[s-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  zd[array2]=dcomplex_add(zd[array2],dcomplex_mul(beta,d[array1+index0i]));
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 24. zg=zg+beta*g[s-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(s-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  zg[array2]=dcomplex_add(zg[array2],dcomplex_mul(beta,g[array1+index0i]));
               }
            }
         } /* End s for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 25. g[ii-1]=zg+zomega */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               g[array1+index0i]=dcomplex_add(zg[array2],zomega[array2]);
            }
         }

         if(ii<iterative.vec){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 26. d[ii-1]=zd-u */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  d[array1+index0i]=dcomplex_sub(zd[array2],u[array2]);   
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 27. c[ii-1]=q[ii]^{H}d[ii-1] */
            reduce=zero;
            for(vc=0;vc<3;vc++){
               array0=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               array1=ii*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array1+index0i]),d[array0+index0i]));   
               }
            }
            MPI_Allreduce(&reduce,&c[ii-1],1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 28. alpha=(q[ii]^{H}u)/(c[ii-1]) */
            reduce=zero;
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=ii*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array1+index0i]),u[array0+index0i]));
               }
            }
            MPI_Allreduce(&reduce,&alpha,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
            /* Check for iterative.breakdown */
            if(dcomplex_abs(c[ii-1])<iterative.breakdown){
               print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","|c[ii-1]|<iterative.breakdown");
            }
            alpha=dcomplex_div(alpha,c[ii-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 29. u=u-alpha*d[ii-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  u[array2]=dcomplex_sub(u[array2],dcomplex_mul(alpha,d[array1+index0i]));
               }
            }

            if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 30. gtilde=M^{-1}g[ii-1] */
               if(iterative.precond==1){
                  /* Point-Jacobi Preconditioning (Divide by the diagonal) */
                  for(vc=0;vc<3;vc++){
                     array0=vc*parallel.alloc_vector;
                     array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
                     for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                        array2=array0+index0i;
                        gtilde[array2]=dcomplex_mul(g[array1+index0i],point_jacobi[array2]);
                     }
                  }
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 31. dipole_polarisation=dipole_polarisation+rho*alpha*gtilde */
               tempdc=dcomplex_mul(rho,alpha);
               for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
                  dipole_polarisation[index0i]=dcomplex_add(dipole_polarisation[index0i],dcomplex_mul(tempdc,gtilde[index0i]));
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 32. omega[ii-1]=interaction_matrix*gtilde */
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
/* 31. dipole_polarisation=dipole_polarisation+rho*alpha*g[ii-1] */
               tempdc=dcomplex_mul(rho,alpha);
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     array2=array0+index0i;
                     dipole_polarisation[array2]=dcomplex_add(dipole_polarisation[array2],dcomplex_mul(tempdc,g[array1+index0i]));
                  }
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 32. omega[ii-1]=interaction_matrix*g[ii-1] */
               dftmatvec_MPI(&g[(ii-1)*parallel.alloc_vector3]);
               for(vc=0;vc<3;vc++){
                  array0=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  array1=vc*parallel.alloc_padded;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     omega[array0+index0i]=zero_padded_vector[array1+index0i];
                  }
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 33. r=r-rho*alpha*omega[ii-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  r[array2]=dcomplex_sub(r[array2],dcomplex_mul(tempdc,omega[array1+index0i]));
               }
            }
         } /* End if ii<iterative.vec */
      } /* End ii for loop */
      iteration++;
   } /* while(iteration<iterative.maximum) */
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
