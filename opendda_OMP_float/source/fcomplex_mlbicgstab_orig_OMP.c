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

   i,ii,s,t,vc,index0*,array*: Loop and array index control variables
   iteration: Iteration number
   iterative.precond: 0: None
                      1: Point-Jacobi Preconditioning
   r[3],g[iterative.vec][3],gtilde[3],
   d[iterative.vec-1][3],omega[iterative.vec][3],
   q[iterative.vec][3],u[3],utilde[3],
   zd[3],zg[3],zomega[3],c: Local vectors for iterative scheme
   alpha,rho,beta,minusone,
   tempdc,real,imag,v0,v1,v2: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "fcomplex_mlbicgstab_orig_OMP.h"

void fcomplex_mlbicgstab_orig_OMP(void){

   int i,ii,s,t,vc,iteration=0;
   long int index0i,array0,array1,array2;
   fcomplex alpha,rho,beta,tempdc,v0,v1,v2;
   float ratio,norm0,norm1,terminate,real,imag;
   double start,finish;

   if(iterative.vec>50){
      print_error("Iterative scheme error","MLBiCGSTAB_orig() error","Too many starting vectors");
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=walltime();}

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
         array0=s*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=index0i+array0;
            q[array2].dat[0]=(float)(genrand_real1()-genrand_real1());
            q[array2].dat[1]=(float)(genrand_real1()-genrand_real1());   
         }
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Modified Gram-Schmidt Orthogonalisation */
   for(s=0;s<iterative.vec;s++){
      norm0=1.0F/fcomplex_2vectornorm(&q[s*target.Nd3]);
      for(vc=0;vc<3;vc++){ /* q[s]=q[s]/||q[s]|| */
         array0=s*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=index0i+array0;
            q[array2]=fcomplex_scale(q[array2],norm0);
         }
      }
      for(t=s+1;t<iterative.vec;t++){
         real=0.0F;imag=0.0F;
         for(vc=0;vc<3;vc++){ /* q[s]^{*}q[t] */
            array0=s*target.Nd3+vc*target.Nd;
            array1=t*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic)
            for(index0i=0;index0i<target.Nd;index0i++){
               v0=q[array0+index0i];
               v1=q[array1+index0i];
               real+=(v0.dat[0]*v1.dat[0]+v0.dat[1]*v1.dat[1]);
               imag+=(v0.dat[0]*v1.dat[1]-v0.dat[1]*v1.dat[0]);
            }
         }
         tempdc.dat[0]=real;tempdc.dat[1]=imag;
         for(vc=0;vc<3;vc++){ /* q[ii]=q[ii]-(q[s]^{*}q[t])q[s] */
            array0=s*target.Nd3+vc*target.Nd;
            array1=t*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=index0i+array1;
               q[array2]=fcomplex_sub(q[array2],fcomplex_mul(tempdc,q[array0+index0i]));
            }
         }
      }
   }

/*   printf("Check dot products are 0 for s!=t and 1 for s==t and check ||q[s]||=1\n");
   for(s=0;s<iterative.vec;s++){
      norm0=fcomplex_2vectornorm(&q[s]);
      fprintf(stdout,"norm for q[%d] = %lf\n",s,norm0);
      for(t=0;t<iterative.vec;t++){
         tempdc=zero;
         for(vc=0;vc<3;vc++){
            array0=s*target.Nd3+vc*target.Nd;
            array1=t*target.Nd3+vc*target.Nd;
            for(index0i=0;index0i<target.Nd;index0i++){
               tempdc=fcomplex_add(tempdc,fcomplex_mul(fcomplex_conj(q[array0+index0i]),q[array1+index0i]));
            }
         }
         printf("vc=%d\ts=%d\tii=%d\t%+.*g%+.*gi\n",vc,s,t,FLTP,tempdc.dat[0],FLTP,tempdc.dat[1]);
      }
   } */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. r=incident_E_field-interaction_matrix*dipole_polarisation */
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

   norm0=fcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0F; 
   };

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. g[iterative.vec-1]=r */
   memcpy(&g[(iterative.vec-1)*target.Nd3],r,(size_t)target.Nd3*sizeof(fcomplex));

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

   while(iteration<iterative.maximum){ /* Check stopping criterion */
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. gtilde=M^{-1}g[iterative.vec-1] */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=(iterative.vec-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=index0i+array0;
                  gtilde[array2]=fcomplex_mul(g[array1+index0i],point_jacobi[array2]);
               }
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. omega[iterative.vec-1]=interaction_matrix*gtilde */
         dftmatvec_OMP(gtilde);
         for(vc=0;vc<3;vc++){
            array0=(iterative.vec-1)*target.Nd3+vc*target.Nd;
            array1=vc*target.Nv;
#pragma omp parallel for schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               omega[array0+index0i]=zero_padded_vector[array1+index0i];
            }
         }
      }
      else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. omega[iterative.vec-1]=interaction_matrix*g[iterative.vec-1] */
         dftmatvec_OMP(&g[(iterative.vec-1)*target.Nd3]);
         for(vc=0;vc<3;vc++){
            array0=(iterative.vec-1)*target.Nd3+vc*target.Nd;
            array1=vc*target.Nv;
#pragma omp parallel for schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               omega[array0+index0i]=zero_padded_vector[array1+index0i];
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. c[iterative.vec-1]=q[0]^{H}omega[iterative.vec-1] */
      real=0.0F;imag=0.0F;
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=(iterative.vec-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            v0=q[array0+index0i];
            v1=omega[array1+index0i];
            real+=(v0.dat[0]*v1.dat[0]+v0.dat[1]*v1.dat[1]);
            imag+=(v0.dat[0]*v1.dat[1]-v0.dat[1]*v1.dat[0]);
         }
      }
      c[iterative.vec-1].dat[0]=real;c[iterative.vec-1].dat[1]=imag;


/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. alpha=(q[0]^{H}r)/(c[iterative.vec-1]) */
      real=0.0F;imag=0.0F;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         v0=q[index0i];
         v1=r[index0i];
         real+=(v0.dat[0]*v1.dat[0]+v0.dat[1]*v1.dat[1]);
         imag+=(v0.dat[0]*v1.dat[1]-v0.dat[1]*v1.dat[0]);
      }
      alpha.dat[0]=real;alpha.dat[1]=imag;
      /* Check for iterative.breakdown */
      if(fcomplex_abs(c[iterative.vec-1])<iterative.breakdown){
         print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","|c[iterative.vec-1]|<iterative.breakdown"); 
      }
      alpha=fcomplex_div(alpha,c[iterative.vec-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. u=r-alpha*omega[iterative.vec-1] */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=(iterative.vec-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=index0i+array0;
            u[array2]=fcomplex_sub(r[array2],fcomplex_mul(alpha,omega[array1+index0i]));
         }
      }

      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. utilde=M^{-1}u */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
#pragma omp parallel for schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd3;index0i++){
               utilde[index0i]=fcomplex_mul(u[index0i],point_jacobi[index0i]);
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. rho=(-[interaction_matrix*utilde]^{H}u)/(||interaction_matrix*utilde||^2)
      interaction_matrix*utilde */
         dftmatvec_OMP(utilde);
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=vc*target.Nv;
#pragma omp parallel for schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               zd[array0+index0i]=zero_padded_vector[array1+index0i]; /* Using zd as a scratch space here */
            }
         }

         real=0.0F;
#pragma omp parallel for private(v0) reduction(+:real) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            v0=zd[index0i];
            real+=(v0.dat[0]*v0.dat[0]+v0.dat[1]*v0.dat[1]);
         }
         tempdc.dat[0]=real;tempdc.dat[1]=0.0;
         /* Check for iterative.breakdown */
         if(fcomplex_abs(tempdc)<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","||interaction_matrix*utilde||<iterative.breakdown");
         }

/* rho=[interaction_matrix*utilde]^{H}u */
         real=0.0F;imag=0.0F;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            v0=zd[index0i];
            v1=u[index0i];
            real+=(v0.dat[0]*v1.dat[0]+v0.dat[1]*v1.dat[1]);
            imag+=(v0.dat[0]*v1.dat[1]-v0.dat[1]*v1.dat[0]);
         }
         rho.dat[0]=real;rho.dat[1]=imag;

/* rho=-[interaction_matrix*utilde]^{H}u */
         rho=fcomplex_mul(minusone,rho);

/* rho=(-[interaction_matrix*utilde]^{H}u)/(||interaction_matrix*utilde||^2) */
         rho=fcomplex_div(rho,tempdc);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. dipole_polarisation=dipole_polarisation-rho*utilde+alpha*gtilde */
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            dipole_polarisation[index0i]=fcomplex_sub(dipole_polarisation[index0i],fcomplex_mul(rho,utilde[index0i]));
            dipole_polarisation[index0i]=fcomplex_add(dipole_polarisation[index0i],fcomplex_mul(alpha,gtilde[index0i]));
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. r=rho*interaction_matrix*utilde+u */
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            r[index0i]=fcomplex_add(fcomplex_mul(rho,zd[index0i]),u[index0i]);
         }
      }
      else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. rho=(-[interaction_matrix*u]^{H}u)/(||interaction_matrix*u||^2)
      interaction_matrix*u */
         dftmatvec_OMP(u);
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=vc*target.Nv;
#pragma omp parallel for schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               zd[array0+index0i]=zero_padded_vector[array1+index0i]; /* Using zd as a scratch space here */
            }
         }

         real=0.0F;
#pragma omp parallel for private(v0) reduction(+:real) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            v0=zd[index0i];
            real+=(v0.dat[0]*v0.dat[0]+v0.dat[1]*v0.dat[1]);
         }
         tempdc.dat[0]=real;tempdc.dat[1]=0.0;
         /* Check for iterative.breakdown */
         if(fcomplex_abs(tempdc)<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","||interaction_matrix*u||<iterative.breakdown"); 
         }

/* rho=[interaction_matrix*u]^{H}u */
         real=0.0F;imag=0.0F;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            v0=zd[index0i];
            v1=u[index0i];
            real+=(v0.dat[0]*u[index0i].dat[0]+v0.dat[1]*v1.dat[1]);
            imag+=(v0.dat[0]*u[index0i].dat[1]-v0.dat[1]*v1.dat[0]);
         }
         rho.dat[0]=real;rho.dat[1]=imag;

/* rho=-[interaction_matrix*u]^{H}u */
         rho=fcomplex_mul(minusone,rho);

/* rho=(-[interaction_matrix*u]^{H}u)/(||interaction_matrix*u||^2) */
         rho=fcomplex_div(rho,tempdc);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. dipole_polarisation=dipole_polarisation-rho*u+alpha*g[iterative.vec-1] */
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=(iterative.vec-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=index0i+array0;
               dipole_polarisation[array2]=fcomplex_sub(dipole_polarisation[array2],fcomplex_mul(rho,u[array2]));
               dipole_polarisation[array2]=fcomplex_add(dipole_polarisation[array2],fcomplex_mul(alpha,g[array1+index0i]));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. r=rho*interaction_matrix*u+u */
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            r[index0i]=fcomplex_add(fcomplex_mul(rho,zd[index0i]),u[index0i]);
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Check stopping criterion */
      norm1=fcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      fprintf(output.iter,"  At iteration %d the convergence ratio is %.*g and the residual norm is %.*g\n",iteration,FLTP,ratio,FLTP,norm1);
      if(ratio<(terminate+FLT_EPSILON)){
         break;
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Start loop over ii */
      for(ii=1;ii<=iterative.vec;ii++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13.  zd=u,zg=r,zomega=0 */
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            zomega[index0i]=zero;
         }
         memcpy(zd,u,(size_t)target.Nd3*sizeof(fcomplex));
         memcpy(zg,r,(size_t)target.Nd3*sizeof(fcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Start loop over s and check if iteration>0 */
         if(iteration>0){
            for(s=ii;s<iterative.vec;s++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. beta=(-q[s]^{H}zd)/(c[s-1]) */
               real=0.0F;imag=0.0F;
               for(vc=0;vc<3;vc++){
                  array0=vc*target.Nd;
                  array1=s*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
                  for(index0i=0;index0i<target.Nd;index0i++){
                     v0=zd[array0+index0i];
                     v1=q[array1+index0i];
                     real+=(v1.dat[0]*v0.dat[0]+v1.dat[1]*v0.dat[1]);
                     imag+=(v1.dat[0]*v0.dat[1]-v1.dat[1]*v0.dat[0]);
                  }
               }
               beta.dat[0]=real;beta.dat[1]=imag;
               /* Check for iterative.breakdown */
               if(fcomplex_abs(c[s-1])<iterative.breakdown){
                  print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","|c[s-1]|<iterative.breakdown");
               }
               beta=fcomplex_div(fcomplex_mul(minusone,beta),c[s-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. zd=zd+beta*d[s-1] */
               for(vc=0;vc<3;vc++){
                  array0=vc*target.Nd;
                  array1=(s-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
                  for(index0i=0;index0i<target.Nd;index0i++){
                     array2=index0i+array0;
                     zd[array2]=fcomplex_add(zd[array2],fcomplex_mul(beta,d[array1+index0i]));
                  }
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. zg=zg+beta*g[s-1] */
               for(vc=0;vc<3;vc++){
                  array0=vc*target.Nd;
                  array1=(s-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
                  for(index0i=0;index0i<target.Nd;index0i++){
                     array2=index0i+array0;
                     zg[array2]=fcomplex_add(zg[array2],fcomplex_mul(beta,g[array1+index0i]));
                  }
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. zomega=zomega+beta*omega[s-1] */
               for(vc=0;vc<3;vc++){
                  array0=vc*target.Nd;
                  array1=(s-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
                  for(index0i=0;index0i<target.Nd;index0i++){
                     array2=index0i+array0;
                     zomega[array2]=fcomplex_add(zomega[array2],fcomplex_mul(beta,omega[array1+index0i]));
                  }
               }
            } /* End s for loop */
         } /* End if iteration>0 */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. beta=(-q[0]^{H}(r+rho*zomega))/(rho*c[iterative.vec-1]) */
         tempdc=fcomplex_mul(rho,c[iterative.vec-1]);
         /* Check for iterative.breakdown */
         if(fcomplex_abs(tempdc)<iterative.breakdown){
            print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","|rho*c[iterative.vec-1]|<iterative.breakdown"); 
         }
         real=0.0F;imag=0.0F;
#pragma omp parallel for private(v0,v1,v2) reduction(+:real,imag) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            v0=q[index0i];
            v1=r[index0i];
            v2=zomega[index0i];
            real+=(v0.dat[0]*(v1.dat[0]+rho.dat[0]*v2.dat[0]-rho.dat[1]*v2.dat[1])+v0.dat[1]*(v1.dat[1]+rho.dat[0]*v2.dat[1]-rho.dat[1]*v2.dat[0]));
            imag+=(v0.dat[0]*(v1.dat[1]+rho.dat[0]*v2.dat[1]-rho.dat[1]*v2.dat[0])-v0.dat[1]*(v1.dat[0]+rho.dat[0]*v2.dat[0]-rho.dat[1]*v2.dat[1]));
         }
         beta.dat[0]=real;beta.dat[1]=imag;
         beta=fcomplex_div(fcomplex_mul(minusone,beta),tempdc);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. zg=zg+beta*g[iterative.vec-1] */
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=(iterative.vec-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=index0i+array0;
               zg[array2]=fcomplex_add(zg[array2],fcomplex_mul(beta,g[array1+index0i]));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 20. zomega=rho(zomega+beta*omega[iterative.vec-1]) */
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=(iterative.vec-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=index0i+array0;
               zomega[array2]=fcomplex_mul(rho,fcomplex_add(zomega[array2],fcomplex_mul(beta,omega[array1+index0i])));
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 21. zd=r+zomega */
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            zd[index0i]=fcomplex_add(r[index0i],zomega[index0i]);
         }

         for(s=1;s<ii;s++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 22. beta=(-q[s]^{H}zd)/(c[s-1]) */
            real=0.0F;imag=0.0F;
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=s*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  v0=zd[array0+index0i];
                  v1=q[array1+index0i];
                  real+=(v1.dat[0]*v0.dat[0]+v1.dat[1]*v0.dat[1]);
                  imag+=(v1.dat[0]*v0.dat[1]-v1.dat[1]*v0.dat[0]);
               }
            }
            beta.dat[0]=real;beta.dat[1]=imag;
            /* Check for iterative.breakdown */
            if(fcomplex_abs(c[s-1])<iterative.breakdown){
               print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","|c[s-1]|<iterative.breakdown");
            }
            beta=fcomplex_div(fcomplex_mul(minusone,beta),c[s-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 23. zd=zd+beta*d[s-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=(s-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=index0i+array0;
                  zd[array2]=fcomplex_add(zd[array2],fcomplex_mul(beta,d[array1+index0i]));
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 24. zg=zg+beta*g[s-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=(s-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=index0i+array0;
                  zg[array2]=fcomplex_add(zg[array2],fcomplex_mul(beta,g[array1+index0i]));
               }
            }
         } /* End s for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 25. g[ii-1]=zg+zomega */
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=(ii-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=index0i+array0;
               g[array1+index0i]=fcomplex_add(zg[array2],zomega[array2]);
            }
         }

         if(ii<iterative.vec){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 26. d[ii-1]=zd-u */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=(ii-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=index0i+array0;
                  d[array1+index0i]=fcomplex_sub(zd[array2],u[array2]);
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 27. c[ii-1]=q[ii]^{H}d[ii-1] */
            real=0.0F;imag=0.0F;
            for(vc=0;vc<3;vc++){
               array0=(ii-1)*target.Nd3+vc*target.Nd;
               array1=ii*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  v0=d[array0+index0i];
                  v1=q[array1+index0i];
                  real+=(v1.dat[0]*v0.dat[0]+v1.dat[1]*v0.dat[1]);
                  imag+=(v1.dat[0]*v0.dat[1]-v1.dat[1]*v0.dat[0]);
               }
            }
            c[ii-1].dat[0]=real;c[ii-1].dat[1]=imag;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 28. alpha=(q[ii]^{H}u)/(c[ii-1]) */
            real=0.0F;imag=0.0F;
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=ii*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  v0=u[array0+index0i];
                  v1=q[array1+index0i];
                  real+=(v1.dat[0]*v0.dat[0]+v1.dat[1]*v0.dat[1]);
                  imag+=(v1.dat[0]*v0.dat[1]-v1.dat[1]*v0.dat[0]);
               }
            }
            alpha.dat[0]=real;alpha.dat[1]=imag;
            /* Check for iterative.breakdown */
            if(fcomplex_abs(c[ii-1])<iterative.breakdown){
               print_error("Iterative scheme error","MLBiCGSTAB_orig() iterative.breakdown","|c[ii-1]|<iterative.breakdown");
            }
            alpha=fcomplex_div(alpha,c[ii-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 29. u=u-alpha*d[ii-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=(ii-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=index0i+array0;
                  u[array2]=fcomplex_sub(u[array2],fcomplex_mul(alpha,d[array1+index0i]));
               }
            }

            if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 30. gtilde=M^{-1}g[ii-1] */
               if(iterative.precond==1){
                  /* Point-Jacobi Preconditioning (Divide by the diagonal) */
                  for(vc=0;vc<3;vc++){
                     array0=vc*target.Nd;
                     array1=(ii-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
                     for(index0i=0;index0i<target.Nd;index0i++){
                        array2=index0i+array0;
                        gtilde[array2]=fcomplex_mul(g[array1+index0i],point_jacobi[array2]);
                     }
                  }
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 31. dipole_polarisation=dipole_polarisation+rho*alpha*gtilde */
               tempdc=fcomplex_mul(rho,alpha);
#pragma omp parallel for schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd3;index0i++){
                  dipole_polarisation[index0i]=fcomplex_add(dipole_polarisation[index0i],fcomplex_mul(tempdc,gtilde[index0i]));
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 32. omega[ii-1]=interaction_matrix*gtilde */
               dftmatvec_OMP(gtilde);
               for(vc=0;vc<3;vc++){
                  array0=(ii-1)*target.Nd3+vc*target.Nd;
                  array1=vc*target.Nv;
#pragma omp parallel for schedule(dynamic,target.M)
                  for(index0i=0;index0i<target.Nd;index0i++){
                     omega[array0+index0i]=zero_padded_vector[array1+index0i];
                  }
               }
            }
            else{ /* Preconditioning disabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 31. dipole_polarisation=dipole_polarisation+rho*alpha*g[ii-1] */
               tempdc=fcomplex_mul(rho,alpha);
               for(vc=0;vc<3;vc++){
                  array0=vc*target.Nd;
                  array1=(ii-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
                  for(index0i=0;index0i<target.Nd;index0i++){
                     array2=index0i+array0;
                     dipole_polarisation[array2]=fcomplex_add(dipole_polarisation[array2],fcomplex_mul(tempdc,g[array1+index0i]));
                  }
               }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 32. omega[ii-1]=interaction_matrix*g[ii-1] */
               dftmatvec_OMP(&g[(ii-1)*target.Nd3]);
               for(vc=0;vc<3;vc++){
                  array0=(ii-1)*target.Nd3+vc*target.Nd;
                  array1=vc*target.Nv;
#pragma omp parallel for schedule(dynamic,target.M)
                  for(index0i=0;index0i<target.Nd;index0i++){
                     omega[array0+index0i]=zero_padded_vector[array1+index0i];
                  }
               }
            }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 33. r=r-rho*alpha*omega[ii-1] */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=(ii-1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=index0i+array0;
                  r[array2]=fcomplex_sub(r[array2],fcomplex_mul(tempdc,omega[array1+index0i]));
               }
            }
         } /* End if ii<iterative.vec */
      } /* End ii for loop */
      iteration++;
   } /* while(iteration<iterative.maximum) */
   if(timing.enabled){ /* Finalise timing */
      finish=walltime();
      fprintf(output.iter,"\n  Time for convergence: %.*gs\n\n",DBLP,finish-start);
      timing.iterative_solution+=finish-start;
   }
}
