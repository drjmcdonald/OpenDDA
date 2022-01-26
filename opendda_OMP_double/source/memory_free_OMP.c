/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Free allocated memory

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "memory_free_OMP.h"

void memory_free_OMP(void){

   int i;

   uint_free(&target.occupied); /* Occupied/Unoccupied bitfielded flag array */
   if(refractive_index.number_of_materials>1){ /* Composition array is only required
      if more than one material is used for target construction */
      for(i=0;i<3;i++){
         uint_free(&target.composition[i]); /* Composition bitfielded array */
      }
   }
   dcomplex_free(&dipole_polarisation); /* Dipole polarisations */
   dcomplex_free(&incident_E_field); /* Incident electric field */
   dcomplex_fftw_free(&zero_padded_vector); /* Zero-padded vector for DFT matrix-vector multiplication */
   dcomplex_free(&interaction_matrix); /* 1st column of the interaction matrix with zero-diagonal */
   dcomplex_free(&interaction_matrix_diagonal); /* Diagonal of the interaction matrix */

   /* To conserve memory array elements that correspond to vacant lattice sites
      and are thus zero are not stored. The populated array stores the starting
      index for the non-zero elements for each of the P xy-planes */
   int_free(&target.populated);

   /* 1D scratch arrays for the DFT of the interaction_matrix */
   for(i=0;i<openmp_threads;i++){
      dcomplex_fftw_free(&xdft[i]);
      dcomplex_fftw_free(&ydft[i]);
      dcomplex_fftw_free(&zdft[i]);
   }
   dcomplex_fftw_free_ptr(&xdft);  /* xdft[openmp_threads][Kpp] */
   dcomplex_fftw_free_ptr(&ydft);  /* ydft[openmp_threads][Jpp] */
   dcomplex_fftw_free_ptr(&zdft);  /* zdft[openmp_threads][Ppp] */

   /* Scratch arrays for the DFT, FD multiplication & iDFT of the 3 zero-padded vector components [1 for each thread] */
   for(i=0;i<openmp_threads;i++){
      dcomplex_fftw_free(&xzscratch[i]);
   }
   dcomplex_fftw_free_ptr(&xzscratch);  /* xzscratch[openmp_threads][3*Qpp] */

   if(iterative.precond==1){ /* Point-Jacobi Preconditioning */
      dcomplex_free(&point_jacobi);
   }

   /* If data list was read in from file then memory was allocated to store it */
   if(wavelength.flag==1){ /* List was read from file */
      double_free(&wavelength.list); /* Wavelengths */
   }
   if(radius.flag==1){ /* List was read from file */
      double_free(&radius.list); /* Effective radii */
   }
   if(euler_phi.flag==1){ /* List was read from file */
      double_free(&euler_phi.list); /* Euler phi */
   }
   if(euler_theta.flag==1){ /* List was read from file */
      double_free(&euler_theta.list); /* Euler theta */
   }
   if(euler_psi.flag==1){ /* List was read from file */
      double_free(&euler_psi.list); /* Euler psi */
   }
   if(phi.flag==1){ /* List was read from file */
      double_free(&phi.list); /* Scattering angle phi */
   }
   if(theta.flag==1){ /* List was read from file */
      double_free(&theta.list); /* Scattering angle theta */
   }

   /* For each constituent material a refractive index versus wavelength table
   can be read in for subsequent interpolation. Scan through the materials and
   if data was read in from file, then memory was allocated, so deallocate it now */
   for(i=0;i<refractive_index.number_of_materials;i++){
      if(refractive_index.flag[i]==1){ /* Data was read from file for interpolation */
         double_free(&refractive_index.wavelengths[i]);
         dcomplex_free(&refractive_index.refractive_indices[i]);
      }
   }
   double_free_ptr(&refractive_index.wavelengths); /* Free pointers */
   dcomplex_free_ptr(&refractive_index.refractive_indices); /* Free pointers */

   dcomplex_free(&refractive_index.value); /* Current refractive index values */
   int_free(&refractive_index.flag); /* Refractive index flag: 0 - control.input file 1 - read from file */
   int_free(&refractive_index.list_length); /* Length of the refractive index versus wavelength table */
   file_free_ptr(&output.refindex_vs_wave); /* Output file pointers for the various materials */

   if(strcmp(iterative.scheme,"bicg")==0){
      dcomplex_free(&r);
      dcomplex_free(&rtilde);
      dcomplex_free(&p);
      dcomplex_free(&ptilde);
      dcomplex_free(&q);
   }
   else if(strcmp(iterative.scheme,"bicg_sym")==0){
      dcomplex_free(&r);
      dcomplex_free(&p);
   }
   else if(strcmp(iterative.scheme,"bicgstab")==0){
      dcomplex_free(&r);
      dcomplex_free(&rtilde);
      dcomplex_free(&p);
      dcomplex_free(&s);
      dcomplex_free(&v);
   }
   else if(strcmp(iterative.scheme,"cg")==0){
      dcomplex_free(&r);
      dcomplex_free(&p);
   }
   else if(strcmp(iterative.scheme,"cgs")==0){
      dcomplex_free(&r);
      dcomplex_free(&rtilde);
      dcomplex_free(&p);
      dcomplex_free(&u);
      dcomplex_free(&q);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab_orig")==0){
      uint_free(&init); /* From SIMD oriented Fast Mersenne Twister initialisation */
      dcomplex_free(&g);
      dcomplex_free(&omega);
      dcomplex_free(&q);
      if(iterative.vec>1){
         dcomplex_free(&d);
      }
      dcomplex_free(&r);
      dcomplex_free(&u);
      dcomplex_free(&zd);
      dcomplex_free(&zg);
      dcomplex_free(&zomega);
      if(iterative.precond!=0){ /* Preconditioning enabled */
         dcomplex_free(&gtilde);
         dcomplex_free(&utilde);
      }
      dcomplex_free(&c);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab")==0){
      uint_free(&init); /* From SIMD oriented Fast Mersenne Twister initialisation */
      dcomplex_free(&g);
      dcomplex_free(&omega);
      dcomplex_free(&q);
      if(iterative.vec>1){
         dcomplex_free(&d);
      }
      dcomplex_free(&r);
      dcomplex_free(&u);
      dcomplex_free(&zd);
      dcomplex_free(&zg);
      dcomplex_free(&zomega);
      if(iterative.precond!=0){ /* Preconditioning enabled */
         dcomplex_free(&gtilde);
         dcomplex_free(&utilde);
      }
      dcomplex_free(&c);
   }
   else if(strcmp(iterative.scheme,"mlbicgstab_ss")==0){
      uint_free(&init); /* From SIMD oriented Fast Mersenne Twister initialisation */
      dcomplex_free(&g);
      dcomplex_free(&omega);
      dcomplex_free(&q);
      dcomplex_free(&r);
      dcomplex_free(&zg);
      dcomplex_free(&zomega);
      dcomplex_free(&c);
   }
   else if(strcmp(iterative.scheme,"qmr")==0){
      dcomplex_free(&r);
      dcomplex_free(&d);
      dcomplex_free(&p);
      dcomplex_free(&q);
      dcomplex_free(&s);
      dcomplex_free(&vtilde);
      dcomplex_free(&omegatilde);
   }
   else if(strcmp(iterative.scheme,"qmr_sym")==0){
      dcomplex_free(&r);
      dcomplex_free(&p);
      dcomplex_free(&pold);
      dcomplex_free(&v);
      dcomplex_free(&vtilde);
   }
   else if(strcmp(iterative.scheme,"rbicgstab")==0){
      dcomplex_free(&u);
      dcomplex_free(&r);
      dcomplex_free(&rtilde);
      for(i=0;i<iterative.vec;i++){
         dcomplex_free(&tau[i]);
      }
      dcomplex_free_ptr(&tau);
      dcomplex_free(&gam);
      dcomplex_free(&gamp);
      dcomplex_free(&gampp);
      dcomplex_free(&sigma);
   }
   else if(strcmp(iterative.scheme,"tfqmr")==0){
      dcomplex_free(&r);
      dcomplex_free(&rtilde);
      dcomplex_free(&d);
      dcomplex_free(&omega);
      dcomplex_free(&v);
      dcomplex_free(&ya);
      dcomplex_free(&yb);
   }
}
