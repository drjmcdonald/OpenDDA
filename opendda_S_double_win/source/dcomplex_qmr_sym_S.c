/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   Quasi-minimal residual for symmetric systems: QMR_sym

   Conjugate gradient-type methods for linear systems with
   complex symmetric coefficient matrices, Freund, R. W.,
   SIAM Journal on Scientific and Statistical Computing,
   Special issue on iterative methods in numerical linear
   algebra, Vol. 13, No. 1, Pages 425-448, 1992.

   vc,index0i,array*: Loop and array index control variables
   iteration: Iteration number
   iterative.precond: 0: None
                      1: Point-Jacobi Preconditioning
   r[3],v[3],vtilde[3],
   p[3],pold[3]: Local vectors for iterative scheme
   alpha,beta,tau,tautilde,
   s,sold,tempdc0,tempdc1,
   theta,eta,zeta,zetatilde,
   one,minusone,zetamod,
   omega,omegaold,c,cold: Local variables for iterative scheme
   ratio,norm0,norm1,terminate: Check stopping criterion variables
   start,finish: Timing variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dcomplex_qmr_sym_S.h"

void dcomplex_qmr_sym_S(void){

   int vc,iteration=0;
   long int index0i,array0,array1,array2;
   dcomplex alpha,beta,tau,tautilde,s={{0.0,0.0}},sold={{0.0,0.0}},tempdc0,tempdc1;
   dcomplex theta,eta,zeta,zetatilde,*ptr;
   double ratio,norm0,norm1,zetamod,omega,omegaold=0.0;
   double c=1.0,cold=1.0,terminate,start,finish;

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
            r[array2]=dcomplex_sub(incident_E_field[array2],zero_padded_vector[array1+index0i]);
         }
      }
   }
   else{ /* r=incident_E_field */
      memcpy(r,incident_E_field,(size_t)target.Nd3*sizeof(dcomplex));
   }

   if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* r=M^{-1}r */
      if(iterative.precond==1){
         /* Point-Jacobi Preconditioning (Divide by the diagonal) */
         for(index0i=0;index0i<target.Nd3;index0i++){
            r[index0i]=dcomplex_mul(r[index0i],point_jacobi[index0i]);
         }
      }
   }

   norm0=dcomplex_2vectornorm(r);
   if(norm0<iterative.breakdown){
      norm0=1.0; 
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 2. v=p=pold=0 set to zero using calloc in memory allocation */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 3. beta=sqrt[vtilde^{T}vtilde]=sqrt[r^{T}r],vtilde=r */
   beta=zero;
   for(index0i=0;index0i<target.Nd3;index0i++){
      beta=dcomplex_add(beta,dcomplex_mul(r[index0i],r[index0i]));
   }
   beta=dcomplex_sqrt(beta);
   /* Check for iterative.breakdown */
   if(dcomplex_abs(beta)<iterative.breakdown){
      print_error("Iterative scheme error","QMR_sym() iterative.breakdown","|beta|<iterative.breakdown");
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 4. tautilde=omega*beta,omega=||vtilde||/|beta|=||r||/|beta| */
/* omega_{0}=||v_{0}||=0.0 set at declaration stage */
   omega=norm0/dcomplex_abs(beta);
   tautilde=dcomplex_scale(beta,omega);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 5. c_{0}=c_{-1}=c=cold=1.0,s_{0}=s_{-1}=s=sold=0.0+0.0i set at declaration stage */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 6. v=vtilde/beta=r/beta */
   tempdc0=dcomplex_div(one,beta); /* (tempdc0=1/beta) */
   for(index0i=0;index0i<target.Nd3;index0i++){
      v[index0i]=dcomplex_mul(r[index0i],tempdc0);
   }

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

   terminate=iterative.tolerance;

   while(iteration<iterative.maximum){ /* 21. Check stopping criterion */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 7. alpha=v^{T}interaction_matrix*v */
/* interaction_matrix*v is stored in zero_padded_vector */
      dftmatvec_S(v);
      if(iterative.precond!=0){ /* Preconditioning enabled */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix*v=M^{-1}interaction_matrix*v */
         if(iterative.precond==1){
            /* Point-Jacobi Preconditioning (Divide by the diagonal) */
            for(vc=0;vc<3;vc++){
               array0=vc*target.Nd;
               array1=vc*target.Nv;
               for(index0i=0;index0i<target.Nd;index0i++){
                  array2=array1+index0i;
                  zero_padded_vector[array2]=dcomplex_mul(zero_padded_vector[array2],point_jacobi[array0+index0i]);
               }
            }
         }
      }
      alpha=zero;
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
         array1=vc*target.Nv;
         for(index0i=0;index0i<target.Nd;index0i++){
            alpha=dcomplex_add(alpha,dcomplex_mul(v[array0+index0i],zero_padded_vector[array1+index0i]));
         }
      }

      if(iteration==0){ /* v_{0} initialised to zero */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. vtilde=interaction_matrix*v-alpha*v */
/* interaction_matrix*v is stored in zero_padded_vector */
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=vc*target.Nv;
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=array0+index0i;
               vtilde[array2]=dcomplex_sub(zero_padded_vector[array1+index0i],dcomplex_mul(alpha,v[array2]));
            }
         }
      }
      else{
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 8. vtilde=interaction_matrix*v-alpha*v-beta*vold (vtilde=vold) */
/* interaction_matrix*v is stored in zero_padded_vector */
         for(vc=0;vc<3;vc++){
            array0=vc*target.Nd;
            array1=vc*target.Nv;
            for(index0i=0;index0i<target.Nd;index0i++){
               array2=array0+index0i;
               vtilde[array2]=dcomplex_sub(zero_padded_vector[array1+index0i],dcomplex_add(dcomplex_mul(alpha,
                              v[array2]),dcomplex_mul(vtilde[array2],beta)));
            }
         }
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 9. theta=conj(sold)*omegaold*beta */
      tempdc0=dcomplex_scale(beta,omegaold); /* tempdc0=omegaold*beta */
      theta=dcomplex_mul(dcomplex_conj(sold),tempdc0);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 10. eta=c*cold*omegaold*beta+conj(s)*omega*alpha */
      tempdc1=dcomplex_scale(alpha,omega); /* tempdc1=omega*alpha */
      eta=dcomplex_add(dcomplex_scale(tempdc0,c*cold),dcomplex_mul(dcomplex_conj(s),tempdc1));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 11. zetatilde=c*omega*alpha-s*cold*omegaold*beta */
      zetatilde=dcomplex_sub(dcomplex_scale(tempdc1,c),dcomplex_mul(s,dcomplex_scale(tempdc0,cold)));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 12. beta=sqrt[vtilde^{T}vtilde],omega=||vtilde||/|beta| */
      beta=zero;
      for(index0i=0;index0i<target.Nd3;index0i++){
         beta=dcomplex_add(beta,dcomplex_mul(vtilde[index0i],vtilde[index0i]));
      }
      beta=dcomplex_sqrt(beta);
      /* Check for iterative.breakdown */
      if(dcomplex_abs(beta)<iterative.breakdown){
         print_error("Iterative scheme error","QMR_sym() iterative.breakdown","|beta|<iterative.breakdown");
      }
      omegaold=omega;
      omega=dcomplex_2vectornorm(vtilde)/dcomplex_abs(beta);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 13. zetamod=sqrt(|zetatilde|^2+(omega^2)|beta|^2) */
      zetamod=sqrt(dcomplex_abs2(zetatilde)+omega*omega*dcomplex_abs2(beta));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 14. zeta=zetamod*zetatilde/|zetatilde| if(zetatilde!=0.0)
            zetamod                       if(zetatilde==0.0) */
      if(dcomplex_abs(zetatilde)<iterative.breakdown){
         zeta.dat[0]=zetamod;
         zeta.dat[1]=0.0;
      }
      else{
         /* Check for iterative.breakdown */
         if(dcomplex_abs(zetatilde)<iterative.breakdown){
            print_error("Iterative scheme error","QMR_sym() iterative.breakdown","|zetatilde|<iterative.breakdown");
         }
         zeta=dcomplex_scale(zetatilde,(zetamod/dcomplex_abs(zetatilde)));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 15. c=zetatilde/zeta */
      cold=c;
      /* Check for iterative.breakdown */
      if(dcomplex_abs(zeta)<iterative.breakdown){
         print_error("Iterative scheme error","QMR_sym() iterative.breakdown","|zeta|<iterative.breakdown");
      }
      tempdc0=dcomplex_div(one,zeta); /* tempdc0=1/zeta */
      c=dcomplex_abs(dcomplex_mul(zetatilde,tempdc0));

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 16. s=omega*beta/zeta (tempdc0=1/zeta) */
      sold=s;
      s=dcomplex_mul(dcomplex_scale(beta,omega),tempdc0);

      if(iteration==0){ /* p_{0} and p_{-1} initialised to zero */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. p=v/zeta (tempdc0=1/zeta) */
         for(index0i=0;index0i<target.Nd3;index0i++){
            p[index0i]=dcomplex_mul(v[index0i],tempdc0);
         }
      }
      else if(iteration==1){ /* p_{0} initialised to zero */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. p=(v-eta*p)/zeta (tempdc0=1/zeta) */
         for(index0i=0;index0i<target.Nd3;index0i++){
            /* pold=(v-eta*p)/zeta */
            pold[index0i]=dcomplex_mul(dcomplex_sub(v[index0i],dcomplex_mul(eta,p[index0i])),tempdc0);
         }
         /* pold=p,p=(v-eta*p)/zeta */
         ptr=pold;pold=p;p=ptr;
      }
      else{
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 17. p=(v-eta*p-theta*pold)/zeta (tempdc0=1/zeta) */
         for(index0i=0;index0i<target.Nd3;index0i++){
            /* pold=(v-eta*p-theta*pold)/zeta */
            pold[index0i]=dcomplex_mul(dcomplex_sub(v[index0i],
                           dcomplex_add(dcomplex_mul(eta,p[index0i]),dcomplex_mul(theta,pold[index0i]))),tempdc0);
         }
         /* pold=p,p=(v-eta*p-theta*pold)/zeta */
         ptr=pold;pold=p;p=ptr;
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 18. tau=c*tautilde,tautilde=-s*tautilde */
      tau=dcomplex_scale(tautilde,c);
      tautilde=dcomplex_mul(dcomplex_mul(minusone,s),tautilde);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 19. dipole_polarisation=dipole_polarisation+tau*p */
      for(index0i=0;index0i<target.Nd3;index0i++){
         dipole_polarisation[index0i]=dcomplex_add(dipole_polarisation[index0i],dcomplex_mul(tau,p[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 20. r=(|s|^2)r+(c*tautilde/omega)*v */
/* v=vtilde/beta */
      /* vtilde=vtilde/beta */
      tempdc0=dcomplex_div(one,beta); /* (tempdc0=1/beta) */ 
      for(index0i=0;index0i<target.Nd3;index0i++){
         vtilde[index0i]=dcomplex_mul(vtilde[index0i],tempdc0);
      }
      /* vtilde=vold,v=vtilde/beta */
      ptr=vtilde;vtilde=v;v=ptr;
      tempdc0=dcomplex_scale(tautilde,c/omega); /* tempdc0=(c*tautilde)/omega */
      ratio=dcomplex_abs2(s); /* ratio used as temp variable for |s| */
      for(index0i=0;index0i<target.Nd3;index0i++){
         r[index0i]=dcomplex_add(dcomplex_scale(r[index0i],ratio),dcomplex_mul(tempdc0,v[index0i]));
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* 21. Check stopping criterion */
      norm1=dcomplex_2vectornorm(r);
      ratio=norm1/norm0;
      fprintf(output.iter,"  At iteration %d the convergence ratio is %.*g and the residual norm is %.*g\n",iteration,DBLP,ratio,DBLP,norm1);
      if(ratio<(terminate+DBL_EPSILON)){
         break;
      }
      iteration++;
   }
   if(timing.enabled){ /* Finalise timing */
      finish=walltime();
      ratio=finish-start;
      fprintf(output.iter,"\n  Time for convergence: %.*gs\n\n",DBLP,ratio);
      timing.iterative_solution+=ratio;
   }
}
