/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   External definitions header file

   dft_interaction_matrix_flag: Flag for the DFT of the 6 independent tensor components
   polarisations: The number of incident polarisations
                  1 - Just specified polarisation is used
                  2 - Calculations are also performed for the orthonormal state
   polarisation_state: The current polarisation state
                  0 - Main specified incident polarisation
                  1 - Orthonormal polarisation state
   DBLP: Number of significant digits for the program output for type double
   FLTP: Number of significant digits for the program output for type float
   openmp_threads: The number of threads to use for parallel regions
   threads_available: Maximum number of threads available by the architecture
   interaction_matrix[6][Ka*Ja*Pa],
   dipole_polarisation[3][Nd],
   incident_E_field[3][Nd]: Stores the 6 independent tensor components of the
                   interaction matrix, the unknown polarisation vector and the
                   incident electric field components for the linear system Ax=b
                   i.e., interaction_matrix*dipole_polarisation=incident_E_field
   zero_padded_vector[3][Kp*Jp*Pp]: Stores the zero-padded vector components for
                  the multiplication and acts as a scratch space for
                  the FD multiplications
   interaction_matrix_diagonal: Stores the polarisability tensors i.e. self-induction
                  terms for the Nd occupied lattice sites
   zero: zero=0.0+0.0i
   one: one=1.0+0.0i
   onei: onei=0.0+1.0i
   two: two=2.0+0.0i
   minusone: minusone=-1.0+0.0i
   xdft,ydft,zdft: Since the 3D DFT of the 6 independent tensor components
                   of the interaction matrix are performed as ensembles
                   of 1D transforms and due to the storage saving techniques
                   employed, these 3 arrays are needed as scratch spaces to
                   perform the DFTs
   xzscratch: Scratch array for the xdft, zdft, xidft & zidft of the three zero-padded
              vector components and the simultaneous Fourier domain multiplication
   plan_vcxdft,
   plan_vcydft,
   plan_vczdft: FFTW plans for the DFT of the 3 zero-padded vector components
                in dftmatvec_OMP, zero_padded_vector[3][K*Jpp*P]
   plan_vcxdfti,
   plan_vcydfti,
   plan_vczdfti: FFTW plans for the iDFT of the 3 zero-padded vector components
                 in dftmatvec_OMP, zero_padded_vector[3][K*Jpp*P]
   plan_tcxdft,
   plan_tcydft,
   plan_tczdft: FFTW plans for the DFT of the 6 tensor components in dftmatvec_OMP
                interaction_matrix[6][Ka*Ja*Pa]
   fpi: PI
   degree: the degree of the polynomials for the interpolating function,
           the default is cubic i.e. degree=3 to preserve the accuracy of a
           3rd-degree polynomial, however, higher degree polynomials, which
           can reduce undulations, can be used at the expense of the 3rd-order
           accuracy. Should be used prudently and sparingly.
   size_parameter: size_parameter~wavenumber*effective_radius
                                 =[2*PI*effective_radius]/wavelength
                   where the effective radius is the radius of a sphere of equal volume
   wavenumber: (2*PI)/wavelength
   wavenumber_previous: The 1st column of the interaction matrix only needs to be
                        recalculated if the wavenumber changes so this stores the
                        previous value for comparison
   wavelength,
   radius,
   euler_phi,
   euler_theta,
   euler_psi,
   phi,
   theta: Control structures for the quantities
          Defined in attributes_type.h

   iterative: Control structure for the iterative scheme
              Defined in attributes_type.h

   refractive_index: Control structure for the complex refractive index
                     Defined in attributes_type.h

   target: Control structure for the target description
           Defined in attributes_type.h

   output: Control structure for data output
           Defined in attributes_type.h

   incident_LF,
   incident_TF: Control structure for incident polarisation in the lab
                frame and target frame, respectively.
                Defined in attributes_type.h

   scattering_LF,
   scattering_TF: Control structure for scattering vectors in the lab
                  frame and target frame, respectively.
                  Defined in attributes_type.h

   cross_section: Control structure for the cross sections
                  Defined in attributes_type.h

   efficiency: Control structure for effeciencies
               Defined in attributes_type.h

   timing: Control structure for timing
           Defined in attributes_type.h

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#ifndef _OPENDDA_DEFS_
#define _OPENDDA_DEFS_
#include <fftw3.h>
#include "attributes_type.h"
#include "fcomplex_type.h"
#define fpi 3.14159265358979323846264338327950288419716939937510F
extern int wl,er,f0,t0,p0,degree,dft_interaction_matrix_flag;
extern int polarisation_state,polarisations,DBLP,FLTP,openmp_threads,threads_available;
extern float size_parameter,wavenumber,wavenumber_previous;
extern fcomplex *dipole_polarisation,*incident_E_field,*zero_padded_vector;
extern fcomplex *interaction_matrix,*interaction_matrix_diagonal,*point_jacobi;
extern fcomplex **xdft,**ydft,**zdft,**xzscratch,two,one,onei,zero,minusone;
extern fftwf_plan plan_vcxdft,plan_vcydft,plan_vczdft,plan_vcxdfti,plan_vcydfti,plan_vczdfti;
extern fftwf_plan plan_tcxdft,plan_tcydft,plan_tczdft;
extern attributes_general wavelength,radius,euler_phi,euler_theta,euler_psi,phi,theta;
extern attributes_refractive_index refractive_index;
extern attributes_iterative iterative;
extern attributes_target target;
extern attributes_output output;
extern attributes_incident_polarisation incident_LF,incident_TF;
extern attributes_scattering_vectors scattering_LF,scattering_TF;
extern attributes_cross_sections cross_section;
extern attributes_efficiency efficiency;
extern attributes_timing timing;
/* Iterative vectors */
extern fcomplex *c,*d,*g,*gam,*gamp,*gampp,*gtilde,*omega,*omegatilde,*p,*pold,*ptilde,*q;
extern fcomplex *r,*rtilde,*s,*sigma,**tau,*u,*utilde,*v,*vtilde,*ya,*yb,*zd,*zg,*zomega;
/* For SIMD oriented Fast Mersenne Twister RNG initialisation */
extern unsigned int *init;
#endif
