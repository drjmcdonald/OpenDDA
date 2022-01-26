/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Allocate memory

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "memory_allocate_MPI.h"

void memory_allocate_MPI(void){

   int i;

   ldcomplex_malloc(&incident_E_field,(size_t)parallel.alloc_vector3); /* incident_E_field[3][Nd] */
   ldcomplex_fftwl_malloc(&zero_padded_vector,(size_t)(parallel.alloc_padded*3)); /* zero_padded_vector[3][P][Jpp][K] */
   ldcomplex_malloc(&interaction_matrix,(size_t)(parallel.alloc_tensor*6)); /* interaction_matrix[6][Pa][Jpp][Ka] */
   ldcomplex_malloc(&dipole_polarisation,(size_t)parallel.alloc_vector3); /* dipole_polarisation[3][Nd] */

   if(iterative.precond==1){ /* Point-Jacobi Preconditioning */
      ldcomplex_malloc(&point_jacobi,(size_t)parallel.alloc_vector3); /* point_jacobi[3][Nd] */
   }

   /* To conserve memory array elements that correspond to vacant lattice sites
      and are thus zero are not stored. The populated array stores the starting
      index for the non-zero elements for each of the parallel.planesY_vector
      local xy-planes on each node */
   int_malloc(&target.populated,(size_t)parallel.planesY_vector);

   /* 1D scratch arrays for the DFT of the interaction_matrix */
   ldcomplex_fftwl_malloc(&xdft,(size_t)target.Kpp); /* xdft[Kpp] */
   ldcomplex_fftwl_malloc(&zdft,(size_t)target.Ppp); /* zdft[Ppp] */

   /* Scratch array for the DFT, FD multiplication & iDFT of the 3 zero-padded vector components */
   ldcomplex_fftwl_malloc(&xzscratch,(size_t)target.Qpp*3); /* xzscratch[3*Qpp] */

   /* The diagonal of interaction_matrix */
   ldcomplex_malloc(&interaction_matrix_diagonal,(size_t)parallel.alloc_vector3); /* interaction_matrix_diagonal[3][Nd] */

   if(parallel.myrank!=0){ /* Restrict to slaves Note: Master node already allocated this in read_parameters() */
      /* Allocate memory to store the current complex refractive index values for materials 0 to (number_of_materials-1) */
      ldcomplex_malloc(&refractive_index.value,(size_t)refractive_index.number_of_materials);
   }

   ldcomplex_malloc(&parallel.local_tensor,(size_t)target.Ka); /* Temporary scratch array for local transpose */
   ldcomplex_malloc(&parallel.local_padded,(size_t)target.K); /* Temporary scratch array for local transpose */
   parallel.tf_ctrl_tensor=(int)((double)(parallel.alloc_tensor-1)/(32.0*(double)target.Ka)+1.0);
   int_malloc(&parallel.transposed_flag_tensor,(size_t)parallel.tf_ctrl_tensor); /* Bit-fielded array for local transpose */
   parallel.tf_ctrl_padded=(int)((double)(parallel.alloc_padded-1)/(32.0*(double)target.K)+1.0);
   int_malloc(&parallel.transposed_flag_padded,(size_t)parallel.tf_ctrl_padded); /* Bit-fielded array for local transpose */

   /* Iterative solution vectors */
   if(strcmp(iterative.scheme,"bicg")==0){
      ldcomplex_malloc(&r,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&rtilde,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&p,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&ptilde,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&q,(size_t)parallel.alloc_vector3);
   }
   else if(strcmp(iterative.scheme,"bicg_sym")==0){
      ldcomplex_malloc(&r,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&p,(size_t)parallel.alloc_vector3);
   }
   else if(strcmp(iterative.scheme,"bicgstab")==0){
      ldcomplex_malloc(&r,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&rtilde,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&p,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&s,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&v,(size_t)parallel.alloc_vector3);
   }
   else if(strcmp(iterative.scheme,"cg")==0){
      ldcomplex_malloc(&r,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&p,(size_t)parallel.alloc_vector3);
   }
   else if(strcmp(iterative.scheme,"cgs")==0){
      ldcomplex_malloc(&r,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&rtilde,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&p,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&u,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&q,(size_t)parallel.alloc_vector3);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab_orig")==0){
      uint_malloc(&init,(size_t)iterative.vec); /* From SIMD oriented Fast Mersenne Twister initialisation */
      ldcomplex_malloc(&g,(size_t)(iterative.vec*parallel.alloc_vector3));
      ldcomplex_malloc(&omega,(size_t)(iterative.vec*parallel.alloc_vector3));
      ldcomplex_malloc(&q,(size_t)(iterative.vec*parallel.alloc_vector3));
      if(iterative.vec>1){
         ldcomplex_malloc(&d,(size_t)((iterative.vec-1)*parallel.alloc_vector3));
      }
      ldcomplex_malloc(&r,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&u,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&zd,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&zg,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&zomega,(size_t)parallel.alloc_vector3);
      if(iterative.precond!=0){ /* Preconditioning enabled */
         ldcomplex_malloc(&gtilde,(size_t)parallel.alloc_vector3);
         ldcomplex_malloc(&utilde,(size_t)parallel.alloc_vector3);
      }
      ldcomplex_malloc(&c,(size_t)iterative.vec);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab")==0){
      uint_malloc(&init,(size_t)iterative.vec); /* From SIMD oriented Fast Mersenne Twister initialisation */
      ldcomplex_malloc(&g,(size_t)(iterative.vec*parallel.alloc_vector3));
      ldcomplex_malloc(&omega,(size_t)(iterative.vec*parallel.alloc_vector3));
      ldcomplex_malloc(&q,(size_t)(iterative.vec*parallel.alloc_vector3));
      if(iterative.vec>1){
         ldcomplex_malloc(&d,(size_t)((iterative.vec-1)*parallel.alloc_vector3));
      }
      ldcomplex_malloc(&r,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&u,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&zd,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&zg,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&zomega,(size_t)parallel.alloc_vector3);
      if(iterative.precond!=0){ /* Preconditioning enabled */
         ldcomplex_malloc(&gtilde,(size_t)parallel.alloc_vector3);
         ldcomplex_malloc(&utilde,(size_t)parallel.alloc_vector3);
      }
      ldcomplex_malloc(&c,(size_t)iterative.vec);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab_ss")==0){
      uint_malloc(&init,(size_t)iterative.vec); /* From SIMD oriented Fast Mersenne Twister initialisation */
      ldcomplex_malloc(&g,(size_t)(iterative.vec*parallel.alloc_vector3));
      ldcomplex_malloc(&omega,(size_t)(iterative.vec*parallel.alloc_vector3));
      ldcomplex_malloc(&q,(size_t)(iterative.vec*parallel.alloc_vector3));
      ldcomplex_malloc(&r,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&zg,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&zomega,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&c,(size_t)iterative.vec);
   }
   else if(strcmp(iterative.scheme,"qmr")==0){
      ldcomplex_malloc(&r,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&p,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&q,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&d,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&s,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&omegatilde,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&vtilde,(size_t)parallel.alloc_vector3);
   }
   else if(strcmp(iterative.scheme,"qmr_sym")==0){
      ldcomplex_malloc(&r,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&p,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&pold,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&v,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&vtilde,(size_t)parallel.alloc_vector3);
   }
   else if(strcmp(iterative.scheme,"rbicgstab")==0){
      ldcomplex_malloc(&u,(size_t)((iterative.vec+1)*parallel.alloc_vector3));
      ldcomplex_malloc(&r,(size_t)((iterative.vec+1)*parallel.alloc_vector3));
      ldcomplex_malloc(&rtilde,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&gam,(size_t)iterative.vec);
      ldcomplex_malloc(&gamp,(size_t)iterative.vec);
      ldcomplex_malloc(&gampp,(size_t)iterative.vec);
      ldcomplex_malloc(&sigma,(size_t)iterative.vec);
      ldcomplex_malloc_ptr(&tau,(size_t)iterative.vec);
      for(i=0;i<iterative.vec;i++){
         ldcomplex_malloc(&tau[i],(size_t)iterative.vec);
      }
   }
   else if(strcmp(iterative.scheme,"tfqmr")==0){
      ldcomplex_malloc(&r,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&rtilde,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&d,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&v,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&ya,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&yb,(size_t)parallel.alloc_vector3);
      ldcomplex_malloc(&omega,(size_t)parallel.alloc_vector3);
   }
}
