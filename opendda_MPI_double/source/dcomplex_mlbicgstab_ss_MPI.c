/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   A variant of the stabilised version of the
   BiConjugate-Gradients (BiCGSTAB) based on
   multiple Lanczos starting vectors: ML(n)BiCGSTAB_{new_ss}

   ML(n)BiCGSTAB: reformulations and implementations,
   Yeung, M. and Boley, D., In press.
   Appendix. Algorithm 2.3' ML(n)BiCGSTAB with
   iterative.preconditioning (storage-saving).

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
   omega[iterative.vec][3],zg[3],
   zomega[3],c: Local vectors for iterative scheme
   alpha,rho,minusone,
   tempdc,tempe,beta,temp,result,element: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   reduce: Used to facilitate the MPI_Allreduce function
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dcomplex_mlbicgstab_ss_MPI.h"

void dcomplex_mlbicgstab_ss_MPI(void){

   int i,ii,s,t,vc,iteration=0;
   long int index0i,array0,array1,array2;
   dcomplex alpha,rho,tempdc,tempe,beta,reduce;
   double ratio,norm0,norm1,terminate,temp,result,start,finish;

   if(iterative.vec>50){
      print_error("Iterative scheme error","MLBiCGSTAB_ss() error","Too many starting vectors");
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

   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. g[0]=M^{-1}r */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            g[index0i]=dcomplex_mul(r[index0i],point_jacobi[index0i]);
         }
      }
   }
   else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. g[0]=r */
      memcpy(g,r,(size_t)parallel.alloc_vector3*sizeof(dcomplex));
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. omega[0]=interaction_matrix*g[0] */
   dftmatvec_MPI(g);
   for(vc=0;vc<3;vc++){
      array0=vc*parallel.alloc_vector;
      array1=vc*parallel.alloc_padded;
      for(index0i=0;index0i<parallel.alloc_vector;index0i++){
         omega[array0+index0i]=zero_padded_vector[array1+index0i];
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. c[0]=q[0]^{H}omega[0] */
   reduce=zero;
   for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
      reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[index0i]),omega[index0i]));
   }
   MPI_Allreduce(&reduce,&c[0],1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* tempe=q[0]^{H}r */
   reduce=zero;
   for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
      reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[index0i]),r[index0i]));      
   }
   MPI_Allreduce(&reduce,&tempe,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

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
      for(ii=1;ii<iterative.vec;ii++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. alpha=(q[ii-1]^{H}r)/(c[ii-1]) */
         /* Check for iterative.breakdown */
         if(dcomplex_abs(c[ii-1])<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB_ss() iterative.breakdown","|c[ii-1]|<iterative.breakdown"); 
         }
         alpha=dcomplex_div(tempe,c[ii-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. dipole_polarisation=dipole_polarisation+alpha*g[ii-1] */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               dipole_polarisation[array2]=dcomplex_add(dipole_polarisation[array2],dcomplex_mul(alpha,g[array1+index0i]));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. r=r-alpha*omega[ii-1] */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=(ii-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               r[array2]=dcomplex_sub(r[array2],dcomplex_mul(alpha,omega[array1+index0i]));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. zomega=r,zg=0 */
         memcpy(zomega,r,(size_t)parallel.alloc_vector3*sizeof(dcomplex));
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            zg[index0i]=zero;
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* tempe=q[ii]^{H}r */
         reduce=zero;
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=ii*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array1+index0i]),r[array0+index0i]));
            }
         }
         MPI_Allreduce(&reduce,&tempe,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

         if(iteration>0){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* beta=(-q[ii]^{H}zomega)/(c[ii]) */
            /* Check for iterative.breakdown */
            if(dcomplex_abs(c[ii])<iterative.breakdown){
               print_error("Iterative scheme error","MLBiCGSTAB_ss() iterative.breakdown","|c[ii]|<iterative.breakdown"); 
            }
            beta=dcomplex_div(dcomplex_mul(minusone,tempe),c[ii]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* zomega=zomega+rho*beta*omega[ii] */
            tempdc=dcomplex_mul(rho,beta);
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=ii*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  zomega[array2]=dcomplex_add(zomega[array2],dcomplex_mul(tempdc,omega[array1+index0i]));
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* zg=zg+beta*g[ii] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=ii*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  zg[array2]=dcomplex_add(zg[array2],dcomplex_mul(beta,g[array1+index0i]));
               }
            }

            for(s=ii+1;s<iterative.vec;s++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. beta=(-q[s]^{H}zomega)/(c[s]) */
               reduce=zero;
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     array2=array0+index0i;
                     reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array1+index0i]),zomega[array0+index0i]));
                  }
               }
               MPI_Allreduce(&reduce,&beta,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
               /* Check for iterative.breakdown */
               if(dcomplex_abs(c[s])<iterative.breakdown){
                  print_error("Iterative scheme error","MLBiCGSTAB_ss() iterative.breakdown","|c[s]|<iterative.breakdown");
               }
               beta=dcomplex_div(dcomplex_mul(minusone,beta),c[s]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. zomega=zomega+rho*beta*omega[s] */
               tempdc=dcomplex_mul(rho,beta);
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     array2=array0+index0i;
                     zomega[array2]=dcomplex_add(zomega[array2],dcomplex_mul(tempdc,omega[array1+index0i]));
                  }
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. zg=zg+beta*g[s] */
               for(vc=0;vc<3;vc++){
                  array0=vc*parallel.alloc_vector;
                  array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
                  for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                     array2=array0+index0i;
                     zg[array2]=dcomplex_add(zg[array2],dcomplex_mul(beta,g[array1+index0i]));
                  }
               }
            } /* End s for loop */
         } /* End if iteration>0 */

         if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. zg=zg+M^{-1}zomega */
            if(iterative.precond==1){
               /* Point-Jacobi Preconditioning (Divide by the diagonal) */
               for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
                  zg[index0i]=dcomplex_add(zg[index0i],dcomplex_mul(zomega[index0i],point_jacobi[index0i]));
               }
            }
         }
         else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. zg=zg+zomega */
            for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
               zg[index0i]=dcomplex_add(zg[index0i],zomega[index0i]);
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. zomega=interaction_matrix*zg */
         dftmatvec_MPI(zg);
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=vc*parallel.alloc_padded;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               zomega[array0+index0i]=zero_padded_vector[array1+index0i];
            }
         }

         for(s=0;s<ii;s++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. beta=(-q[s]^{H}zomega)/(c[s]) */
            reduce=zero;
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array1+index0i]),zomega[array0+index0i]));
               }
            }
            MPI_Allreduce(&reduce,&beta,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
            /* Check for iterative.breakdown */
            if(dcomplex_abs(c[s])<iterative.breakdown){
               print_error("Iterative scheme error","MLBiCGSTAB_ss() iterative.breakdown","|c[s]|<iterative.breakdown");
            }
            beta=dcomplex_div(dcomplex_mul(minusone,beta),c[s]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. zomega=zomega+beta*omega[s] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  zomega[array2]=dcomplex_add(zomega[array2],dcomplex_mul(beta,omega[array1+index0i]));
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. zg=zg+beta*g[s] */
            for(vc=0;vc<3;vc++){
               array0=vc*parallel.alloc_vector;
               array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
               for(index0i=0;index0i<parallel.alloc_vector;index0i++){
                  array2=array0+index0i;
                  zg[array2]=dcomplex_add(zg[array2],dcomplex_mul(beta,g[array1+index0i]));
               }
            }
         } /* End s for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. g[ii]=zg */
         memcpy(&g[ii*parallel.alloc_vector3],zg,(size_t)parallel.alloc_vector3*sizeof(dcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. omega[ii]=zomega */
         memcpy(&omega[ii*parallel.alloc_vector3],zomega,(size_t)parallel.alloc_vector3*sizeof(dcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 20. c[ii]=q[ii]^{H}zomega */
         reduce=zero;
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=ii*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array1+index0i]),zomega[array0+index0i]));
            }
         }
         MPI_Allreduce(&reduce,&c[ii],1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
      } /* End ii for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 21. alpha=(q[iterative.vec-1]^{H}r)/(c[iterative.vec-1]) */
      /* Check for iterative.breakdown */
      if(dcomplex_abs(c[iterative.vec-1])<iterative.breakdown){
         print_error("Iterative scheme error","MLBiCGSTAB_ss() iterative.breakdown","|c[iterative.vec-1]|<iterative.breakdown");
      }
      alpha=dcomplex_div(tempe,c[iterative.vec-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 22. dipole_polarisation=dipole_polarisation+alpha*g[iterative.vec-1] */
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            dipole_polarisation[array2]=dcomplex_add(dipole_polarisation[array2],dcomplex_mul(alpha,g[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 23. r=r-alpha*omega[iterative.vec-1] */
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=(iterative.vec-1)*parallel.alloc_vector3+vc*parallel.alloc_vector;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            array2=array0+index0i;
            r[array2]=dcomplex_sub(r[array2],dcomplex_mul(alpha,omega[array1+index0i]));
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

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 24. zg=M^{-1}r */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
               zg[index0i]=dcomplex_mul(r[index0i],point_jacobi[index0i]);
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 25. rho=(-[interaction_matrix*M^{-1}r]^{H}r)/(||interaction_matrix*M^{-1}r||^2)=
           (-[interaction_matrix*zg]^{H}r)/(||interaction_matrix*zg||^2)
         interaction_matrix*M^{-1}r */
         dftmatvec_MPI(zg);
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=vc*parallel.alloc_padded;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               zomega[array0+index0i]=zero_padded_vector[array1+index0i]; /* Using zomega as a scratch space here */
            }
         }

         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(zomega[index0i]),zomega[index0i]));
         }
         MPI_Allreduce(&reduce,&tempdc,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
         /* Check for iterative.breakdown */
         if(dcomplex_abs(tempdc)<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB_ss() iterative.breakdown","||interaction_matrix*M^{-1}r||<iterative.breakdown");
         }

/* rho=[interaction_matrix*M^{-1}r]^{H}r */
         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(zomega[index0i]),r[index0i]));   
         }
         MPI_Allreduce(&reduce,&rho,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

/* rho=-[interaction_matrix*M^{-1}r]^{H}r */
         rho=dcomplex_mul(minusone,rho);

/* rho=(-[interaction_matrix*M^{-1}r]^{H}r)/(||interaction_matrix*M^{-1}r||^2) */
         rho=dcomplex_div(rho,tempdc);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 26. dipole_polarisation=dipole_polarisation-rho*zg */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            dipole_polarisation[index0i]=dcomplex_sub(dipole_polarisation[index0i],dcomplex_mul(rho,zg[index0i]));   
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 27. r=rho*interaction_matrix*M^{-1}r+r=rho*zomega+r= */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            r[index0i]=dcomplex_add(dcomplex_mul(rho,zomega[index0i]),r[index0i]);
         }
      }
      else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 25. rho=(-[interaction_matrix*r]^{H}r)/(||interaction_matrix*r||^2)
         interaction_matrix*r */
         dftmatvec_MPI(r);
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=vc*parallel.alloc_padded;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               zomega[array0+index0i]=zero_padded_vector[array1+index0i]; /* Using zomega as a scratch space here */
            }
         }

         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(zomega[index0i]),zomega[index0i]));
         }
         MPI_Allreduce(&reduce,&tempdc,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
         /* Check for iterative.breakdown */
         if(dcomplex_abs(tempdc)<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB_ss() iterative.breakdown","||interaction_matrix*r||<iterative.breakdown");
         }

/* rho=[interaction_matrix*r]^{H}r */
         reduce=zero;
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(zomega[index0i]),r[index0i]));      
         }
         MPI_Allreduce(&reduce,&rho,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

/* rho=-[interaction_matrix*r]^{H}r */
         rho=dcomplex_mul(minusone,rho);

/* rho=(-[interaction_matrix*r]^{H}r)/(||interaction_matrix*r||^2) */
         rho=dcomplex_div(rho,tempdc);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 26. dipole_polarisation=dipole_polarisation-rho*r */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            dipole_polarisation[index0i]=dcomplex_sub(dipole_polarisation[index0i],dcomplex_mul(rho,r[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 27. r=rho*interaction_matrix*r+r=rho*zomega+r */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            r[index0i]=dcomplex_add(dcomplex_mul(rho,zomega[index0i]),r[index0i]);
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 28. zomega=r,zg=0 */
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         zg[index0i]=zero;
      }
      memcpy(zomega,r,(size_t)parallel.alloc_vector3*sizeof(dcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* tempe=q[0]^{H}r */
      reduce=zero;
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[index0i]),r[index0i]));
      }
      MPI_Allreduce(&reduce,&tempe,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

      for(i=0;i<iterative.vec;i++){
         c[i]=dcomplex_mul(c[i],rho);
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* beta=(-q[0]^{H}zomega)/(rho*c[0]) */
      /* Check for iterative.breakdown */
      if(dcomplex_abs(c[0])<iterative.breakdown){
         print_error("Iterative scheme error","MLBiCGSTAB_ss() iterative.breakdown","|rho*c[0]|<iterative.breakdown"); 
      }
      beta=dcomplex_div(dcomplex_mul(minusone,tempe),c[0]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* zomega=zomega+rho*beta*omega[0] */
      tempdc=dcomplex_mul(rho,beta);
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         zomega[index0i]=dcomplex_add(zomega[index0i],dcomplex_mul(tempdc,omega[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* zg=zg+beta*g[0] */
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         zg[index0i]=dcomplex_add(zg[index0i],dcomplex_mul(beta,g[index0i]));
      }

      for(s=1;s<iterative.vec;s++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 29. beta=(-q[s]^{H}zomega)/(rho*c[s]) */
         reduce=zero;
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[array1+index0i]),zomega[array0+index0i]));
            }
         }
         MPI_Allreduce(&reduce,&beta,1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);
         /* Check for iterative.breakdown */
         if(dcomplex_abs(c[s])<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB_ss() iterative.breakdown","|rho*c[s]|<iterative.breakdown");
         }
         beta=dcomplex_div(dcomplex_mul(minusone,beta),c[s]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 30. zomega=zomega+rho*beta*omega[s] */
         tempdc=dcomplex_mul(rho,beta);
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               zomega[array2]=dcomplex_add(zomega[array2],dcomplex_mul(tempdc,omega[array1+index0i]));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 31. zg=zg+beta*g[s] */
         for(vc=0;vc<3;vc++){
            array0=vc*parallel.alloc_vector;
            array1=s*parallel.alloc_vector3+vc*parallel.alloc_vector;
            for(index0i=0;index0i<parallel.alloc_vector;index0i++){
               array2=array0+index0i;
               zg[array2]=dcomplex_add(zg[array2],dcomplex_mul(beta,g[array1+index0i]));
            }
         }
      } /* End s for loop */

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 32. zg=zg+M^{-1}zomega */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
               zg[index0i]=dcomplex_add(zg[index0i],dcomplex_mul(zomega[index0i],point_jacobi[index0i]));
            }
         }
      }
      else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 32. zg=zg+zomega */
         for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
            zg[index0i]=dcomplex_add(zg[index0i],zomega[index0i]);
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 33. zomega=interaction_matrix*zg */
      dftmatvec_MPI(zg);
      for(vc=0;vc<3;vc++){
         array0=vc*parallel.alloc_vector;
         array1=vc*parallel.alloc_padded;
         for(index0i=0;index0i<parallel.alloc_vector;index0i++){
            zomega[array0+index0i]=zero_padded_vector[array1+index0i];
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 34. c[0]=q[0]^{H}zomega */
      reduce=zero;
      for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
         reduce=dcomplex_add(reduce,dcomplex_mul(dcomplex_conj(q[index0i]),zomega[index0i]));
      }
      MPI_Allreduce(&reduce,&c[0],1,dcomplex_type,add_dcomplex,MPI_COMM_WORLD);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 35. g[0]=zg,omega[0]=zomega */
      memcpy(g,zg,(size_t)parallel.alloc_vector3*sizeof(dcomplex));
      memcpy(omega,zomega,(size_t)parallel.alloc_vector3*sizeof(dcomplex));
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
