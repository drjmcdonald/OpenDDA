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
   tau[iterative.vec][iterative.vec]: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "fcomplex_rbicgstab_S.h"

void fcomplex_rbicgstab_S(void){

   int i,ii,jj,vc,iteration=0;
   long int index0i,array0,array1,array2;
   fcomplex rhozero={{1.0F,0.0F}},rho,alpha={{0.0F,0.0F}},omega={{1.0F,0.0F}},beta,xi;
   float ratio,norm0,norm1,terminate;
   double start,finish;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Initialise timing */
   if(timing.enabled){start=walltime();}

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 1. r=incident_E_field-interaction_matrix(dipole_polarisation) */
   if(iterative.initial_guess>0){
      dftmatvec_S(dipole_polarisation);
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
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
         for(index0i=0;index0i<target.Nd3;index0i++){
            r[index0i]=fcomplex_mul(r[index0i],point_jacobi[index0i]);
         }
      }
   }

   norm0=fcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0F; 
   };

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. rtilde=r[0] */
   memcpy(rtilde,r,(size_t)target.Nd3*sizeof(fcomplex));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. u[0]=0 Set to zero using calloc in memory allocation */ 

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. rhozero=1,alpha=0,omega=1 Set at variable declaration stage */

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
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. rhozero=-omega*rhozero */
      rhozero=fcomplex_mul(fcomplex_mul(minusone,omega),rhozero);

      for(jj=0;jj<iterative.vec;jj++){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. rho=r[jj]^{T}rtilde */
         rho=zero;
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=jj*target.Nd3+vc*target.Nd;
            for(index0i=0;index0i<target.Nd;index0i++){
               rho=fcomplex_add(rho,fcomplex_mul(r[array1+index0i],rtilde[array0+index0i]));
            }
         }
         /* Check for iterative.breakdown */
         if(fcomplex_abs(rhozero)<iterative.breakdown){
            print_error("Iterative scheme error","RBiCGSTAB() iterative.breakdown","|rhozero|<iterative.breakdown");
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. beta=(alpha*rho)/(rhozero) */
         beta=fcomplex_div(fcomplex_mul(alpha,rho),rhozero);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. rhozero=rho */
         rhozero=rho;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. u[ii]=r[ii]-beta*u[ii] */
         for(ii=0;ii<=jj;ii++){
            for(vc=0;vc<3;vc++){
               array0=ii*target.Nd3+vc*target.Nd;
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array0+index0i;
                  u[array2]=fcomplex_sub(r[array2],fcomplex_mul(beta,u[array2]));
               }
            }
         } /* End ii for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. u[jj+1]=interaction_matrix*u[jj] */
         dftmatvec_S(&u[jj*target.Nd3]);
         for(vc=0;vc<3;vc++){
            array0=(jj+1)*target.Nd3+vc*target.Nd;
            array1=vc*target.Nv;
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
                  for(index0i=0;index0i<target.Nd;index0i++){
                     array2=index0i+array1;
                     u[array2]=fcomplex_mul(u[array2],point_jacobi[array0+index0i]);
                  }
               }
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. xi=u[jj+1]^{T}rtilde */
         xi=zero;
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=(jj+1)*target.Nd3+vc*target.Nd;
            for(index0i=0;index0i<target.Nd;index0i++){
               xi=fcomplex_add(xi,fcomplex_mul(u[array1+index0i],rtilde[array0+index0i]));
            }
         }
         /* Check for iterative.breakdown */
         if(fcomplex_abs(xi)<iterative.breakdown){
            print_error("Iterative scheme error","RBiCGSTAB() iterative.breakdown","|xi|<iterative.breakdown");
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. alpha=rhozero/xi */
         alpha=fcomplex_div(rhozero,xi);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. r[ii]=r[ii]-alpha*u[ii+1] */
         for(ii=0;ii<=jj;ii++){
            for(vc=0;vc<3;vc++){
               array0=(ii+1)*target.Nd3+vc*target.Nd;
               array1=ii*target.Nd3+vc*target.Nd;
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array1+index0i;
                  r[array2]=fcomplex_sub(r[array2],fcomplex_mul(alpha,u[array0+index0i]));
               }
            }
         } /* End ii for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. r[jj+1]=interaction_matrix*r[jj] */
         dftmatvec_S(&r[jj*target.Nd3]);
         for(vc=0;vc<3;vc++){
            array0=(jj+1)*target.Nd3+vc*target.Nd;
            array1=vc*target.Nv;
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
                  for(index0i=0;index0i<target.Nd;index0i++){
                     array2=index0i+array1;
                     r[array2]=fcomplex_mul(r[array2],point_jacobi[array0+index0i]);
                  }
               }
            }
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. dipole_polarisation=dipole_polarisation+alpha*u[0] */
         for(index0i=0;index0i<target.Nd3;index0i++){
            dipole_polarisation[index0i]=fcomplex_add(dipole_polarisation[index0i],fcomplex_mul(alpha,u[index0i]));
         }
      } /* End jj for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. Check stopping criterion */
      norm1=fcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      fprintf(output.iter,"  At iteration %d the convergence ratio is %.*g and the residual norm is %.*g\n",iteration,FLTP,ratio,FLTP,norm1);
      if(ratio<(terminate+FLT_EPSILON)){
         break;
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. sigma[0]=r[1]^{T}r[1], gamp[0]=(r[0]^{T}r[1])/(sigma[0]) */
      sigma[0]=zero;
      for(vc=0;vc<3;vc++){
         array0=target.Nd3+vc*target.Nd;
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=index0i+array0;
            sigma[0]=fcomplex_add(sigma[0],fcomplex_mul(r[array2],r[array2]));
         }
      }
      /* Check for iterative.breakdown */
      if(fcomplex_abs(sigma[0])<iterative.breakdown){
         print_error("Iterative scheme error","RBiCGSTAB() iterative.breakdown","|sigma[0]|<iterative.breakdown"); 
      }
      gamp[0]=zero;
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=target.Nd3+vc*target.Nd;
         for(index0i=0;index0i<target.Nd;index0i++){
            gamp[0]=fcomplex_add(gamp[0],fcomplex_mul(r[array0+index0i],r[array1+index0i]));
         }
      }
      gamp[0]=fcomplex_div(gamp[0],sigma[0]);

      for(jj=2;jj<=iterative.vec;jj++){
         for(ii=1;ii<jj;ii++){

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. tau[ii-1][jj-1]=(r[jj]^{T}r[ii])/(sigma[ii-1]) */
            /* Check for iterative.breakdown */
            if(fcomplex_abs(sigma[ii-1])<iterative.breakdown){
               print_error("Iterative scheme error","RBiCGSTAB() iterative.breakdown","|sigma[ii-1]|<iterative.breakdown");
            }
            tau[ii-1][jj-1]=zero;
            for(vc=0;vc<3;vc++){
               array0=jj*target.Nd3+vc*target.Nd;
               array1=ii*target.Nd3+vc*target.Nd;
               for(index0i=0;index0i<target.Nd;index0i++){
                  tau[ii-1][jj-1]=fcomplex_add(tau[ii-1][jj-1],fcomplex_mul(r[array0+index0i],r[array1+index0i]));
               }
            }
            tau[ii-1][jj-1]=fcomplex_div(tau[ii-1][jj-1],sigma[ii-1]);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* r[jj]=r[jj]-tau[ii-1][jj-1]*r[ii] */
            for(vc=0;vc<3;vc++){
               array0=jj*target.Nd3+vc*target.Nd;
               array1=ii*target.Nd3+vc*target.Nd;
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array0+index0i;
                  r[array2]=fcomplex_sub(r[array2],fcomplex_mul(tau[ii-1][jj-1],r[array1+index0i]));
               }
            }
         } /* End ii for loop */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. sigma[jj-1]=(r[jj])^{T}r[jj], gamp[jj-1]=(r[0]^{T}r[jj])/(sigma[jj-1]) */
         sigma[jj-1]=zero;
         for(vc=0;vc<3;vc++){
            array0=jj*target.Nd3+vc*target.Nd;
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=index0i+array0;
               sigma[jj-1]=fcomplex_add(sigma[jj-1],fcomplex_mul(r[array2],r[array2]));
            }
         }
         /* Check for iterative.breakdown */
         if(fcomplex_abs(sigma[jj-1])<iterative.breakdown){
            print_error("Iterative scheme error","RBiCGSTAB() iterative.breakdown","|sigma[jj-1]|<iterative.breakdown"); 
         }
         gamp[jj-1]=zero;
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=jj*target.Nd3+vc*target.Nd;
            for(index0i=0;index0i<target.Nd;index0i++){
               gamp[jj-1]=fcomplex_add(gamp[jj-1],fcomplex_mul(r[array0+index0i],r[array1+index0i]));
               
            }
         }
         gamp[jj-1]=fcomplex_div(gamp[jj-1],sigma[jj-1]);
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
            gam[jj]=fcomplex_add(gam[jj],fcomplex_mul(tau[jj][i],gam[i]));
         }
         gam[jj]=fcomplex_sub(gamp[jj],gam[jj]);
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 21. gampp[jj]=gam[jj+1]+sum_{i=jj+1}^{iterative.vec-1}tau[jj][i]*gam[i+1] jj=0...iterative.vec-2 */
      for(jj=0;jj<iterative.vec-1;jj++){
         gampp[jj]=zero;
         for(i=jj+1;i<iterative.vec-1;i++){
            gampp[jj]=fcomplex_add(gampp[jj],fcomplex_mul(tau[jj][i],gam[i+1]));
         }
         gampp[jj]=fcomplex_add(gam[jj+1],gampp[jj]);
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 22. dipole_polarisation=dipole_polarisation+gam[0]r[0] */
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=fcomplex_add(dipole_polarisation[index0i],fcomplex_mul(gam[0],r[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 23. r[0]=r[0]-gamp[iterative.vec-1]r[iterative.vec] */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=iterative.vec*target.Nd3+vc*target.Nd;
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            r[array2]=fcomplex_sub(r[array2],fcomplex_mul(gamp[iterative.vec-1],r[array1+index0i]));
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 24. u[0]=u[0]-gam[iterative.vec-1]u[iterative.vec] */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=iterative.vec*target.Nd3+vc*target.Nd;
         for(index0i=0;index0i<target.Nd;index0i++){
            array2=array0+index0i;
            u[array2]=fcomplex_sub(u[array2],fcomplex_mul(gam[iterative.vec-1],u[array1+index0i]));
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
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=array0+index0i;
               u[array2]=fcomplex_sub(u[array2],fcomplex_mul(gam[jj-1],u[array1+index0i]));
            }
         }
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=jj*target.Nd3+vc*target.Nd;
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=array0+index0i;
               dipole_polarisation[array2]=fcomplex_add(dipole_polarisation[array2],fcomplex_mul(gampp[jj-1],r[array1+index0i]));
            }
         }
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=jj*target.Nd3+vc*target.Nd;
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=array0+index0i;
               r[array2]=fcomplex_sub(r[array2],fcomplex_mul(gamp[jj-1],r[array1+index0i]));
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
