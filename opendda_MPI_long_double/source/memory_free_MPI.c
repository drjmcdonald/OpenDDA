/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Free allocated memory

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "memory_free_MPI.h"

void memory_free_MPI(void){

   int i;

   uint_free(&target.occupied); /* Occupied/Unoccupied bitfielded flag array */
   if(refractive_index.number_of_materials>1){ /* Composition array is only required
      if more than one material is used for target construction */
      for(i=0;i<3;i++){
         uint_free(&target.composition[i]); /* Composition bitfielded array */
      }
   }
   ldcomplex_free(&dipole_polarisation); /* Dipole polarisations */
   ldcomplex_free(&incident_E_field); /* Incident electric field */
   ldcomplex_fftwl_free(&zero_padded_vector); /* Zero-padded vector for DFT matrix-vector multiplication */
   ldcomplex_free(&interaction_matrix); /* 1st column of the interaction matrix with zero-diagonal */
   ldcomplex_free(&interaction_matrix_diagonal); /* Diagonal of the interaction matrix */

   /* To conserve memory array elements that correspond to vacant lattice sites
      and are thus zero are not stored. The populated array stores the starting
      index for the non-zero elements for each of the parallel.planesY_vector
      local xy-planes on each node */
   int_free(&target.populated);

   /* 1D scratch arrays for the DFT of the interaction_matrix */
   ldcomplex_fftwl_free(&xdft);  /* xdft[Kpp] */
   ldcomplex_fftwl_free(&zdft);  /* zdft[Ppp] */

   /* Scratch array for the DFT, FD multiplication & iDFT of the 3 zero-padded vector components */
   ldcomplex_fftwl_free(&xzscratch); /* xzscratch[3*Qpp] */

   if(iterative.precond==1){ /* Point-Jacobi Preconditioning */
      ldcomplex_free(&point_jacobi);
   }

   if(parallel.myrank==0){ /* Restrict to master */
      /* If data list was read in from file then memory was allocated to store it */   
      if(wavelength.flag==1){ /* List was read from file */
         long_double_free(&wavelength.list); /* Wavelengths */
      }
      if(radius.flag==1){ /* List was read from file */
         long_double_free(&radius.list); /* Effective radii */
      }
      if(euler_phi.flag==1){ /* List was read from file */
         long_double_free(&euler_phi.list); /* Euler phi */
      }
      if(euler_theta.flag==1){ /* List was read from file */
         long_double_free(&euler_theta.list); /* Euler theta */
      }
      if(euler_psi.flag==1){ /* List was read from file */
         long_double_free(&euler_psi.list); /* Euler psi */
      }
      if(phi.flag==1){ /* List was read from file */
         long_double_free(&phi.list); /* Scattering angle phi */
      }
      if(theta.flag==1){ /* List was read from file */
         long_double_free(&theta.list); /* Scattering angle theta */
      }
   }

   if(parallel.myrank==0){ /* Restrict to master */
      /* For each constituent material a refractive index versus wavelength table
      can be read in for subsequent interpolation. Scan through the materials and
      if data was read in from file, then memory was allocated, so deallocate it now */
      for(i=0;i<refractive_index.number_of_materials;i++){
         if(refractive_index.flag[i]==1){ /* Data was read from file for interpolation */
            long_double_free(&refractive_index.wavelengths[i]);
            ldcomplex_free(&refractive_index.refractive_indices[i]);
         }
      }
      long_double_free_ptr(&refractive_index.wavelengths); /* Free pointers */
      ldcomplex_free_ptr(&refractive_index.refractive_indices); /* Free pointers */

      int_free(&refractive_index.flag); /* Refractive index flag: 0 - control.input file 1 - read from file */
      int_free(&refractive_index.list_length); /* Length of the refractive index versus wavelength table */
      file_free_ptr(&output.refindex_vs_wave); /* Output file pointers for the various materials */
   }
   ldcomplex_free(&refractive_index.value); /* Current refractive index values */

   ldcomplex_free(&parallel.local_tensor); /* Temporary scratch array for local transpose */
   ldcomplex_free(&parallel.local_padded); /* Temporary scratch array for local transpose */
   int_free(&parallel.transposed_flag_tensor); /* Bit-fielded array for local transpose */
   int_free(&parallel.transposed_flag_padded); /* Bit-fielded array for local transpose */

   if(strcmp(iterative.scheme,"bicg")==0){
      ldcomplex_free(&r);
      ldcomplex_free(&rtilde);
      ldcomplex_free(&p);
      ldcomplex_free(&ptilde);
      ldcomplex_free(&q);
   }
   else if(strcmp(iterative.scheme,"bicg_sym")==0){
      ldcomplex_free(&r);
      ldcomplex_free(&p);
   }
   else if(strcmp(iterative.scheme,"bicgstab")==0){
      ldcomplex_free(&r);
      ldcomplex_free(&rtilde);
      ldcomplex_free(&p);
      ldcomplex_free(&s);
      ldcomplex_free(&v);
   }
   else if(strcmp(iterative.scheme,"cg")==0){
      ldcomplex_free(&r);
      ldcomplex_free(&p);
   }
   else if(strcmp(iterative.scheme,"cgs")==0){
      ldcomplex_free(&r);
      ldcomplex_free(&rtilde);
      ldcomplex_free(&p);
      ldcomplex_free(&u);
      ldcomplex_free(&q);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab_orig")==0){
      uint_free(&init); /* From SIMD oriented Fast Mersenne Twister initialisation */
      ldcomplex_free(&g);
      ldcomplex_free(&omega);
      ldcomplex_free(&q);
      if(iterative.vec>1){
         ldcomplex_free(&d);
      }
      ldcomplex_free(&r);
      ldcomplex_free(&u);
      ldcomplex_free(&zd);
      ldcomplex_free(&zg);
      ldcomplex_free(&zomega);
      if(iterative.precond!=0){ /* Preconditioning enabled */
         ldcomplex_free(&gtilde);
         ldcomplex_free(&utilde);
      }
      ldcomplex_free(&c);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab")==0){
      uint_free(&init); /* From SIMD oriented Fast Mersenne Twister initialisation */
      ldcomplex_free(&g);
      ldcomplex_free(&omega);
      ldcomplex_free(&q);
      if(iterative.vec>1){
         ldcomplex_free(&d);
      }
      ldcomplex_free(&r);
      ldcomplex_free(&u);
      ldcomplex_free(&zd);
      ldcomplex_free(&zg);
      ldcomplex_free(&zomega);
      if(iterative.precond!=0){ /* Preconditioning enabled */
         ldcomplex_free(&gtilde);
         ldcomplex_free(&utilde);
      }
      ldcomplex_free(&c);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab_ss")==0){
      uint_free(&init); /* From SIMD oriented Fast Mersenne Twister initialisation */
      ldcomplex_free(&g);
      ldcomplex_free(&omega);
      ldcomplex_free(&q);
      ldcomplex_free(&r);
      ldcomplex_free(&zg);
      ldcomplex_free(&zomega);
      ldcomplex_free(&c);
   }
   else if(strcmp(iterative.scheme,"qmr")==0){
      ldcomplex_free(&r);
      ldcomplex_free(&d);
      ldcomplex_free(&p);
      ldcomplex_free(&q);
      ldcomplex_free(&s);
      ldcomplex_free(&vtilde);
      ldcomplex_free(&omegatilde);
   }
   else if(strcmp(iterative.scheme,"qmr_sym")==0){
      ldcomplex_free(&r);
      ldcomplex_free(&p);
      ldcomplex_free(&pold);
      ldcomplex_free(&v);
      ldcomplex_free(&vtilde);
   }
   else if(strcmp(iterative.scheme,"rbicgstab")==0){
      ldcomplex_free(&u);
      ldcomplex_free(&r);
      ldcomplex_free(&rtilde);
      for(i=0;i<iterative.vec;i++){
         ldcomplex_free(&tau[i]);
      }
      ldcomplex_free_ptr(&tau);
      ldcomplex_free(&gam);
      ldcomplex_free(&gamp);
      ldcomplex_free(&gampp);
      ldcomplex_free(&sigma);
   }
   else if(strcmp(iterative.scheme,"tfqmr")==0){
      ldcomplex_free(&r);
      ldcomplex_free(&rtilde);
      ldcomplex_free(&d);
      ldcomplex_free(&omega);
      ldcomplex_free(&v);
      ldcomplex_free(&ya);
      ldcomplex_free(&yb);
   }
}
