/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   Restarted, stabilised version of the BiConjugate-Gradients: RBiCGSTAB

   BiCGSTAB(L) for linear equations involving unsymmetric matrices
   with complex spectrum, Sleijpen, G.L.G. and Fokkema, D.R.,
   Elect. Trans. Numer. Anal., Vol. 1, Pages 11-3, 1993.

   Adapted from:
   PIM 2.2: The parallel iterative methods package for systems of linear
   equations. User's Guide (Fortran 77 version)}, da Cunha, R.D. and
   Hopkins, T., 1997.

   i,ii,jj,vc,index0i,array*: Loop and array index control variables
   iteration: Iteration number
   iterative.precond: 0: None
                      1: Point-Jacobi Preconditioning
   r[iterative.vec+1][3],rtilde[3],
   u[iterative.vec+1][3]: Local vectors for iterative scheme
   minusone,rhozero,rho,
   alpha,omega,beta,
   sigma[iterative.vec],gam[iterative.vec],
   gampp[iterative.vec],
   gamp[iterative.vec],
   tau[iterative.vec][iterative.vec],
   real,imag,v0,v1: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "ldcomplex_rbicgstab_OMP.h"

void ldcomplex_rbicgstab_OMP(void){

   int i,ii,jj,vc,iteration=0;
   long int index0i,array0,array1,array2;
   ldcomplex rhozero={{1.0L,0.0L}},rho,alpha={{0.0L,0.0L}},omega={{1.0L,0.0L}},beta,xi,v0,v1;
   long double ratio,norm0,norm1,terminate,real,imag;
   double start,finish;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=walltime();}

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 1. r=incident_E_field-interaction_matrix(dipole_polarisation) */
   if(iterative.initial_guess>0){
      dftmatvec_OMP(dipole_polarisation);
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
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
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            r[index0i]=ldcomplex_mul(r[index0i],point_jacobi[index0i]);
         }
      }
   }

   norm0=ldcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0L; 
   };

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. rtilde=r[0] */
   memcpy(rtilde,r,(size_t)target.Nd3*sizeof(ldcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. u[0]=0 Set to zero using calloc in memory allocation */ 

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. rhozero=1,alpha=0,omega=1 Set at variable declaration stage */

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

   while(iteration<iterative.maximum){ /* Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. rhozero=-omega*rhozero */
      rhozero=ldcomplex_mul(ldcomplex_mul(minusone,omega),rhozero);

      for(jj=0;jj<iterative.vec;jj++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. rho=r[jj]^{T}rtilde */
         real=0.0L;imag=0.0L;
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=jj*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               v0=rtilde[array0+index0i];
               v1=r[array1+index0i];
               real+=(v1.dat[0]*v0.dat[0]-v1.dat[1]*v0.dat[1]);
               imag+=(v1.dat[0]*v0.dat[1]+v1.dat[1]*v0.dat[0]);
            }
         }
         rho.dat[0]=real;rho.dat[1]=imag;
         /* Check for iterative.breakdown */
         if(ldcomplex_abs(rhozero)<iterative.breakdown){
            print_error("Iterative scheme error","RBiCGSTAB() iterative.breakdown","|rhozero|<iterative.breakdown");
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. beta=(alpha*rho)/(rhozero) */
         beta=ldcomplex_div(ldcomplex_mul(alpha,rho),rhozero);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. rhozero=rho */
         rhozero=rho;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. u[ii]=r[ii]-beta*u[ii] */
         for(ii=0;ii<=jj;ii++){
            for(vc=0;vc<3;vc++){
               array0=ii*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array0+index0i;
                  u[array2]=ldcomplex_sub(r[array2],ldcomplex_mul(beta,u[array2]));
               }
            }
         } /* End ii for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. u[jj+1]=interaction_matrix*u[jj] */
         dftmatvec_OMP(&u[jj*target.Nd3]);
         for(vc=0;vc<3;vc++){
            array0=(jj+1)*target.Nd3+vc*target.Nd;
            array1=vc*target.Nv;
#pragma omp parallel for schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               u[array0+index0i]=zero_padded_vector[array1+index0i];
            }
         }

         if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* u[jj+1]=M^{-1}u[jj+1] */
            if(iterative.precond==1){
               /* Point-Jacobi Preconditioning (Divide by the diagonal) */
               for(vc=0;vc<3;vc++){
                  array0=vc*target.Nd;
                  array1=(jj+1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
                  for(index0i=0;index0i<target.Nd;index0i++){
                     array2=index0i+array1;
                     u[array2]=ldcomplex_mul(u[array2],point_jacobi[array0+index0i]);
                  }
               }
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. xi=u[jj+1]^{T}rtilde */
         real=0.0L;imag=0.0L;
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=(jj+1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               v0=rtilde[array0+index0i];
               v1=u[array1+index0i];
               real+=(v1.dat[0]*v0.dat[0]-v1.dat[1]*v0.dat[1]);
               imag+=(v1.dat[0]*v0.dat[1]+v1.dat[1]*v0.dat[0]);
            }
         }
         xi.dat[0]=real;xi.dat[1]=imag;
         /* Check for iterative.breakdown */
         if(ldcomplex_abs(xi)<iterative.breakdown){
            print_error("Iterative scheme error","RBiCGSTAB() iterative.breakdown","|xi|<iterative.breakdown");
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. alpha=rhozero/xi */
         alpha=ldcomplex_div(rhozero,xi);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. r[ii]=r[ii]-alpha*u[ii+1] */
         for(ii=0;ii<=jj;ii++){
            for(vc=0;vc<3;vc++){
               array0=(ii+1)*target.Nd3+vc*target.Nd;
               array1=ii*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array1+index0i;
                  r[array2]=ldcomplex_sub(r[array2],ldcomplex_mul(alpha,u[array0+index0i]));
               }
            }
         } /* End ii for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. r[jj+1]=interaction_matrix*r[jj] */
         dftmatvec_OMP(&r[jj*target.Nd3]);
         for(vc=0;vc<3;vc++){
            array0=(jj+1)*target.Nd3+vc*target.Nd;
            array1=vc*target.Nv;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               r[array0+index0i]=zero_padded_vector[array1+index0i];
            }
         }

         if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* r[jj+1]=M^{-1}r[jj+1] */
            if(iterative.precond==1){
               /* Point-Jacobi Preconditioning (Divide by the diagonal) */
               for(vc=0;vc<3;vc++){
                  array0=vc*target.Nd;
                  array1=(jj+1)*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
                  for(index0i=0;index0i<target.Nd;index0i++){
                     array2=index0i+array1;
                     r[array2]=ldcomplex_mul(r[array2],point_jacobi[array0+index0i]);
                  }
               }
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. dipole_polarisation=dipole_polarisation+alpha*u[0] */
#pragma omp parallel for schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd3;index0i++){
            dipole_polarisation[index0i]=ldcomplex_add(dipole_polarisation[index0i],ldcomplex_mul(alpha,u[index0i]));
         }
      } /* End jj for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. Check stopping criterion */
      norm1=ldcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      fprintf(output.iter,"  At iteration %d the convergence ratio is %.*Lg and the residual norm is %.*Lg\n",iteration,LDBLP,ratio,LDBLP,norm1);
      if(ratio<(terminate+LDBL_EPSILON)){
         break;
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. sigma[0]=r[1]^{T}r[1], gamp[0]=(r[0]^{T}r[1])/(sigma[0]) */
      real=0.0L;imag=0.0L;
      for(vc=0;vc<3;vc++){
         array0=target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0) reduction(+:real,imag) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            v0=r[index0i+array0];
            real+=(v0.dat[0]*v0.dat[0]-v0.dat[1]*v0.dat[1]);
            imag+=(2.0L*(v0.dat[0]*v0.dat[1]));
         }
      }
      sigma[0].dat[0]=real;sigma[0].dat[1]=imag;
      /* Check for iterative.breakdown */
      if(ldcomplex_abs(sigma[0])<iterative.breakdown){
         print_error("Iterative scheme error","RBiCGSTAB() iterative.breakdown","|sigma[0]|<iterative.breakdown"); 
      }
      real=0.0L;imag=0.0L;
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            v0=r[array0+index0i];
            v1=r[array1+index0i];
            real+=(v0.dat[0]*v1.dat[0]-v0.dat[1]*v1.dat[1]);
            imag+=(v0.dat[0]*v1.dat[1]+v0.dat[1]*v1.dat[0]);
         
         }
      }
      gamp[0].dat[0]=real;gamp[0].dat[1]=imag;
      gamp[0]=ldcomplex_div(gamp[0],sigma[0]);

      for(jj=2;jj<=iterative.vec;jj++){
         for(ii=1;ii<jj;ii++){

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. tau[ii-1][jj-1]=(r[jj]^{T}r[ii])/(sigma[ii-1]) */
            /* Check for iterative.breakdown */
            if(ldcomplex_abs(sigma[ii-1])<iterative.breakdown){
               print_error("Iterative scheme error","RBiCGSTAB() iterative.breakdown","|sigma[ii-1]|<iterative.breakdown");
            }
            real=0.0L;imag=0.0L;
            for(vc=0;vc<3;vc++){
               array0=jj*target.Nd3+vc*target.Nd;
               array1=ii*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  v0=r[array0+index0i];
                  v1=r[array1+index0i];
                  real+=(v0.dat[0]*v1.dat[0]-v0.dat[1]*v1.dat[1]);
                  imag+=(v0.dat[0]*v1.dat[1]+v0.dat[1]*v1.dat[0]);
               }
            }
            tau[ii-1][jj-1].dat[0]=real;tau[ii-1][jj-1].dat[1]=imag;
            tau[ii-1][jj-1]=ldcomplex_div(tau[ii-1][jj-1],sigma[ii-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* r[jj]=r[jj]-tau[ii-1][jj-1]*r[ii] */
            for(vc=0;vc<3;vc++){
               array0=jj*target.Nd3+vc*target.Nd;
               array1=ii*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array0+index0i;
                  r[array2]=ldcomplex_sub(r[array2],ldcomplex_mul(tau[ii-1][jj-1],r[array1+index0i]));
               }
            }
         } /* End ii for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. sigma[jj-1]=(r[jj])^{T}r[jj], gamp[jj-1]=(r[0]^{T}r[jj])/(sigma[jj-1]) */
         real=0.0L;imag=0.0L;
         for(vc=0;vc<3;vc++){
            array0=jj*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0) reduction(+:real,imag) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               v0=r[index0i+array0];
               real+=(v0.dat[0]*v0.dat[0]-v0.dat[1]*v0.dat[1]);
               imag+=(2.0L*(v0.dat[0]*v0.dat[1]));
            }
         }
         sigma[jj-1].dat[0]=real;sigma[jj-1].dat[1]=imag;
         /* Check for iterative.breakdown */
         if(ldcomplex_abs(sigma[jj-1])<iterative.breakdown){
            print_error("Iterative scheme error","RBiCGSTAB() iterative.breakdown","|sigma[jj-1]|<iterative.breakdown"); 
         }
         real=0.0L;imag=0.0L;
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=jj*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(v0,v1) reduction(+:real,imag) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               v0=r[array0+index0i];
               v1=r[array1+index0i];
               real+=(v0.dat[0]*v1.dat[0]-v0.dat[1]*v1.dat[1]);
               imag+=(v0.dat[0]*v1.dat[1]+v0.dat[1]*v1.dat[0]);
            }
         }
         gamp[jj-1].dat[0]=real;gamp[jj-1].dat[1]=imag;
         gamp[jj-1]=ldcomplex_div(gamp[jj-1],sigma[jj-1]);
      } /* End jj for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 20. gam[iterative.vec-1]=omega=gamp[iterative.vec-1] */
      gam[iterative.vec-1]=gamp[iterative.vec-1];
      omega=gam[iterative.vec-1];

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* gam[jj]=gamp[jj]-sum_{i=jj+1}^{iterative.vec}tau[jj][i]*gam[i] jj=iterative.vec-2...0 */
      for(jj=iterative.vec-2;jj>=0;jj--){
         gam[jj]=zero;
         for(i=jj+1;i<iterative.vec;i++){
            gam[jj]=ldcomplex_add(gam[jj],ldcomplex_mul(tau[jj][i],gam[i]));
         }
         gam[jj]=ldcomplex_sub(gamp[jj],gam[jj]);
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 21. gampp[jj]=gam[jj+1]+sum_{i=jj+1}^{iterative.vec-1}tau[jj][i]*gam[i+1] jj=0...iterative.vec-2 */
      for(jj=0;jj<iterative.vec-1;jj++){
         gampp[jj]=zero;
         for(i=jj+1;i<iterative.vec-1;i++){
            gampp[jj]=ldcomplex_add(gampp[jj],ldcomplex_mul(tau[jj][i],gam[i+1]));
         }
         gampp[jj]=ldcomplex_add(gam[jj+1],gampp[jj]);
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 22. dipole_polarisation=dipole_polarisation+gam[0]r[0] */
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=ldcomplex_add(dipole_polarisation[index0i],ldcomplex_mul(gam[0],r[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 23. r[0]=r[0]-gamp[iterative.vec-1]r[iterative.vec] */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=iterative.vec*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            r[array2]=ldcomplex_sub(r[array2],ldcomplex_mul(gamp[iterative.vec-1],r[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 24. u[0]=u[0]-gam[iterative.vec-1]u[iterative.vec] */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=iterative.vec*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            u[array2]=ldcomplex_sub(u[array2],ldcomplex_mul(gam[iterative.vec-1],u[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 25. Start loop over jj and u[0]=u[0]-gam[jj-1]*u[jj] 
   26.                       dipole_polarisation[0]=dipole_polarisation[0]+gampp[jj-1]*r[jj]
   27.                       r[0]=r[0]-gamp[jj-1]*r[jj] */
      for(jj=1;jj<iterative.vec;jj++){
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=jj*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=array0+index0i;
               u[array2]=ldcomplex_sub(u[array2],ldcomplex_mul(gam[jj-1],u[array1+index0i]));
            }
         }
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=jj*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=array0+index0i;
               dipole_polarisation[array2]=ldcomplex_add(dipole_polarisation[array2],ldcomplex_mul(gampp[jj-1],r[array1+index0i]));
            }
         }
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=jj*target.Nd3+vc*target.Nd;
#pragma omp parallel for private(array2) schedule(dynamic,target.M)
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=array0+index0i;
               r[array2]=ldcomplex_sub(r[array2],ldcomplex_mul(gamp[jj-1],r[array1+index0i]));
            }
         }
      }
      iteration++;
   } /* End while loop */
   if(timing.enabled){ /* Finalise timing */
      finish=walltime();
      fprintf(output.iter,"\n  Time for convergence: %.*gs\n\n",DBLP,finish-start);
      timing.iterative_solution+=finish-start;
   }
}
