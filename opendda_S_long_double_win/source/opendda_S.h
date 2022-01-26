/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for opendda and global variables
   Note: see definitions.h for an explanation of the variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <float.h>
#include <math.h>
#include "attributes_type.h"
#include "bitfielding.h"
#include "build_target.h"
#include "cross_sections_S.h"
#include "ldcomplex_type.h"
#include "dda_interaction_matrix_S.h"
#include "definitions.h"
#include "dft_plan.h"
#include "dipole_polarisabilities_S.h"
#include "efficiencies.h"
#include "find_dft_size.h"
#include "incident_electric_field_S.h"
#include "initialise_output.h"
#include "interpolation.h"
#include "iterative_solution.h"
#include "memory_allocate_S.h"
#include "memory_free_S.h"
#include "read_parameters.h"
#include "rotation.h"
#include "set_domain.h"
#include "set_initial_guess.h"
#include "timing_info.h"
#include "walltime.h"

/* Function prototypes */
void initialise(int flag,int *list_length,long double *limits);

/* Global variables */
int wl,er,f0,t0,p0,degree,dft_interaction_matrix_flag=0;
int polarisation_state,polarisations,DBLP,LDBLP;
long double size_parameter,wavenumber,wavenumber_previous=-1.0L;
ldcomplex *dipole_polarisation,*incident_E_field,*zero_padded_vector;
ldcomplex *interaction_matrix,*interaction_matrix_diagonal,*point_jacobi;
ldcomplex *xdft,*ydft,*zdft,*xzscratch,two={{2.0L,0.0L}},one={{1.0L,0.0L}};
ldcomplex onei={{0.0L,1.0L}},zero={{0.0L,0.0L}},minusone={{-1.0L,0.0L}};
fftwl_plan plan_vcxdft,plan_vcydft,plan_vczdft,plan_vcxdfti,plan_vcydfti,plan_vczdfti;
fftwl_plan plan_tcxdft,plan_tcydft,plan_tczdft;
attributes_general wavelength={0},radius={0},euler_phi={0},euler_theta={0},euler_psi={0},phi,theta={0};
attributes_refractive_index refractive_index={0};
attributes_iterative iterative={0};
attributes_target target={0};
attributes_output output={0};
attributes_incident_polarisation incident_LF={0},incident_TF={0};
attributes_scattering_vectors scattering_LF={0},scattering_TF={0};
attributes_cross_sections cross_section={0};
attributes_efficiency efficiency={0};
attributes_timing timing={0};
/* Iterative vectors */
ldcomplex *c,*d,*g,*gam,*gamp,*gampp,*gtilde,*omega,*omegatilde,*p,*pold,*ptilde,*q;
ldcomplex *r,*rtilde,*s,*sigma,**tau,*u,*utilde,*v,*vtilde,*ya,*yb,*zd,*zg,*zomega;
/* For SIMD oriented Fast Mersenne Twister RNG initialisation */
unsigned int *init;
