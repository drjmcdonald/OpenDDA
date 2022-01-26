/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Allocate memory

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "memory_allocate_OMP.h"

void memory_allocate_OMP(void){

   int i;

   dcomplex_malloc(&incident_E_field,(size_t)target.Nd3); /* incident_E_field[3][Nd] */
   dcomplex_fftw_malloc(&zero_padded_vector,(size_t)(target.Nv*3)); /* zero_padded_vector[3][P][Jpp][K] */
   dcomplex_malloc(&interaction_matrix,(size_t)(target.Na*6)); /* interaction_matrix[6][Pa][Ja][Ka] */
   dcomplex_malloc(&dipole_polarisation,(size_t)target.Nd3); /* dipole_polarisation[3][Nd] */

   if(iterative.precond==1){ /* Point-Jacobi Preconditioning */
      dcomplex_malloc(&point_jacobi,(size_t)(target.Nd3)); /* point_jacobi[3][Nd] */
   }

   /* To conserve memory array elements that correspond to vacant lattice sites
      and are thus zero are not stored. The populated array stores the starting
      index for the non-zero elements for each of the P xy-planes */
   int_malloc(&target.populated,(size_t)target.P);

   /* 1D scratch arrays for the DFT of the interaction_matrix [1 for each thread] */
   dcomplex_fftw_malloc_ptr(&xdft,(size_t)openmp_threads);
   dcomplex_fftw_malloc_ptr(&ydft,(size_t)openmp_threads);
   dcomplex_fftw_malloc_ptr(&zdft,(size_t)openmp_threads);
   for(i=0;i<openmp_threads;i++){
      dcomplex_fftw_malloc(&xdft[i],(size_t)target.Kpp); /* xdft[openmp_threads][Kpp] */
      dcomplex_fftw_malloc(&ydft[i],(size_t)target.Jpp); /* ydft[openmp_threads][Jpp] */
      dcomplex_fftw_malloc(&zdft[i],(size_t)target.Ppp); /* zdft[openmp_threads][Ppp] */
   }

   /* Scratch arrays for the DFT, FD multiplication & iDFT of the 3 zero-padded vector components [1 for each thread] */
   dcomplex_fftw_malloc_ptr(&xzscratch,(size_t)openmp_threads); /* xzscratch[3*Qpp] */
   for(i=0;i<openmp_threads;i++){
      dcomplex_fftw_malloc(&xzscratch[i],(size_t)target.Qpp*3); /* xzscratch[openmp_threads][3*Qpp] */
   }

   /* The diagonal of interaction_matrix */
   dcomplex_malloc(&interaction_matrix_diagonal,(size_t)(target.Nd3)); /* interaction_matrix_diagonal[3][Nd] */

   /* Iterative solution vectors */
   if(strcmp(iterative.scheme,"bicg")==0){
      dcomplex_malloc(&r,(size_t)target.Nd3);
      dcomplex_malloc(&rtilde,(size_t)target.Nd3);
      dcomplex_malloc(&p,(size_t)target.Nd3);
      dcomplex_malloc(&ptilde,(size_t)target.Nd3);
      dcomplex_malloc(&q,(size_t)target.Nd3);
   }
   else if(strcmp(iterative.scheme,"bicg_sym")==0){
      dcomplex_malloc(&r,(size_t)target.Nd3);
      dcomplex_malloc(&p,(size_t)target.Nd3);
   }
   else if(strcmp(iterative.scheme,"bicgstab")==0){
      dcomplex_malloc(&r,(size_t)target.Nd3);
      dcomplex_malloc(&rtilde,(size_t)target.Nd3);
      dcomplex_malloc(&p,(size_t)target.Nd3);
      dcomplex_malloc(&s,(size_t)target.Nd3);
      dcomplex_malloc(&v,(size_t)target.Nd3);
   }
   else if(strcmp(iterative.scheme,"cg")==0){
      dcomplex_malloc(&r,(size_t)target.Nd3);
      dcomplex_malloc(&p,(size_t)target.Nd3);
   }
   else if(strcmp(iterative.scheme,"cgs")==0){
      dcomplex_malloc(&r,(size_t)target.Nd3);
      dcomplex_malloc(&rtilde,(size_t)target.Nd3);
      dcomplex_malloc(&p,(size_t)target.Nd3);
      dcomplex_malloc(&u,(size_t)target.Nd3);
      dcomplex_malloc(&q,(size_t)target.Nd3);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab_orig")==0){
      uint_malloc(&init,(size_t)iterative.vec); /* From SIMD oriented Fast Mersenne Twister initialisation */
      dcomplex_malloc(&g,(size_t)(iterative.vec*target.Nd3));
      dcomplex_malloc(&omega,(size_t)(iterative.vec*target.Nd3));
      dcomplex_malloc(&q,(size_t)(iterative.vec*target.Nd3));
      if(iterative.vec>1){
         dcomplex_malloc(&d,(size_t)((iterative.vec-1)*target.Nd3));
      }
      dcomplex_malloc(&r,(size_t)target.Nd3);
      dcomplex_malloc(&u,(size_t)target.Nd3);
      dcomplex_malloc(&zd,(size_t)target.Nd3);
      dcomplex_malloc(&zg,(size_t)target.Nd3);
      dcomplex_malloc(&zomega,(size_t)target.Nd3);
      if(iterative.precond!=0){ /* Preconditioning enabled */
         dcomplex_malloc(&gtilde,(size_t)target.Nd3);
         dcomplex_malloc(&utilde,(size_t)target.Nd3);
      }
      dcomplex_malloc(&c,(size_t)iterative.vec);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab")==0){
      uint_malloc(&init,(size_t)iterative.vec); /* From SIMD oriented Fast Mersenne Twister initialisation */
      dcomplex_malloc(&g,(size_t)(iterative.vec*target.Nd3));
      dcomplex_malloc(&omega,(size_t)(iterative.vec*target.Nd3));
      dcomplex_malloc(&q,(size_t)(iterative.vec*target.Nd3));
      if(iterative.vec>1){
         dcomplex_malloc(&d,(size_t)((iterative.vec-1)*target.Nd3));
      }
      dcomplex_malloc(&r,(size_t)target.Nd3);
      dcomplex_malloc(&u,(size_t)target.Nd3);
      dcomplex_malloc(&zd,(size_t)target.Nd3);
      dcomplex_malloc(&zg,(size_t)target.Nd3);
      dcomplex_malloc(&zomega,(size_t)target.Nd3);
      if(iterative.precond!=0){ /* Preconditioning enabled */
         dcomplex_malloc(&gtilde,(size_t)target.Nd3);
         dcomplex_malloc(&utilde,(size_t)target.Nd3);
      }
      dcomplex_malloc(&c,(size_t)iterative.vec);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab_ss")==0){
      uint_malloc(&init,(size_t)iterative.vec); /* From SIMD oriented Fast Mersenne Twister initialisation */
      dcomplex_malloc(&g,(size_t)(iterative.vec*target.Nd3));
      dcomplex_malloc(&omega,(size_t)(iterative.vec*target.Nd3));
      dcomplex_malloc(&q,(size_t)(iterative.vec*target.Nd3));
      dcomplex_malloc(&r,(size_t)target.Nd3);
      dcomplex_malloc(&zg,(size_t)target.Nd3);
      dcomplex_malloc(&zomega,(size_t)target.Nd3);
      dcomplex_malloc(&c,(size_t)iterative.vec);
   }
   else if(strcmp(iterative.scheme,"qmr")==0){
      dcomplex_malloc(&r,(size_t)target.Nd3);
      dcomplex_malloc(&p,(size_t)target.Nd3);
      dcomplex_malloc(&q,(size_t)target.Nd3);
      dcomplex_malloc(&d,(size_t)target.Nd3);
      dcomplex_malloc(&s,(size_t)target.Nd3);
      dcomplex_malloc(&omegatilde,(size_t)target.Nd3);
      dcomplex_malloc(&vtilde,(size_t)target.Nd3);
   }
   else if(strcmp(iterative.scheme,"qmr_sym")==0){
      dcomplex_malloc(&r,(size_t)target.Nd3);
      dcomplex_malloc(&p,(size_t)target.Nd3);
      dcomplex_malloc(&pold,(size_t)target.Nd3);
      dcomplex_malloc(&v,(size_t)target.Nd3);
      dcomplex_malloc(&vtilde,(size_t)target.Nd3);
   }
   else if(strcmp(iterative.scheme,"rbicgstab")==0){
      dcomplex_malloc(&u,(size_t)((iterative.vec+1)*target.Nd3));
      dcomplex_malloc(&r,(size_t)((iterative.vec+1)*target.Nd3));
      dcomplex_malloc(&rtilde,(size_t)target.Nd3);
      dcomplex_malloc(&gam,(size_t)iterative.vec);
      dcomplex_malloc(&gamp,(size_t)iterative.vec);
      dcomplex_malloc(&gampp,(size_t)iterative.vec);
      dcomplex_malloc(&sigma,(size_t)iterative.vec);
      dcomplex_malloc_ptr(&tau,(size_t)iterative.vec);
      for(i=0;i<iterative.vec;i++){
         dcomplex_malloc(&tau[i],(size_t)iterative.vec);
      }
   }
   else if(strcmp(iterative.scheme,"tfqmr")==0){
      dcomplex_malloc(&r,(size_t)target.Nd3);
      dcomplex_malloc(&rtilde,(size_t)target.Nd3);
      dcomplex_malloc(&d,(size_t)target.Nd3);
      dcomplex_malloc(&v,(size_t)target.Nd3);
      dcomplex_malloc(&ya,(size_t)target.Nd3);
      dcomplex_malloc(&yb,(size_t)target.Nd3);
      dcomplex_malloc(&omega,(size_t)target.Nd3);
   }
}
