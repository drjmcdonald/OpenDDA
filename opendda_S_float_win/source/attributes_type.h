/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Type definitions for control of

   attributes_general:
      wavelength: Incident wavelength
      radius: Effective radius of the target
      euler_phi: Target rotation Euler angle phi
      euler_theta: Target rotation Euler angle theta
      euler_psi: Target rotation Euler angle psi
      phi: Scattering direction angle phi
      theta: Scattering direction angle theta

   attributes_refractive_index:
      Complex refractive index control

   attributes_iterative:
      Iterative scheme control

   attributes_target:
      Target description

   attributes_output:
      Output file control

   attributes_incident_polarisation:
      Incident polarisation control

   attributes_scattering_vectors:
      Scattering vector control

   attributes_cross_sections:
      Cross section control

   attributes_effeciency:
      Efficiency control

   Copyright (C) 2006 James Mc Donald
   Computational Astrophysics Laboratory
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include "fcomplex_type.h"

/* General attribute control type definition */
#ifndef _OPENDDA_ATTRIBUTES_
#define _OPENDDA_ATTRIBUTES_
typedef struct {
/* Type definition for a structure to hold the control information for

   wavelength: Incident wavelength
   radius: Effective radius of the target
   euler_phi: Target rotation Euler angle phi
   euler_theta: Target rotation Euler angle theta
   euler_psi: Target rotation Euler angle psi
   phi: Scattering direction angle phi
   theta: Scattering direction angle theta

   flag: 0 if data is governed by minimum,increment,maximum from file control.input
         1 if data is read in from a file
   list_length: number in data entries
   value: variable to store current value
   limits: minimum,increment,maximum from file control.input. Only meaningful if flag=0
   *list: storage array for data read from file. Only meaningful if flag=1, otherwise *list='NULL' */
   int flag;
   int list_length;
   float value;
   float limits[3];
   float *list;
} attributes_general;
#endif

/* Complex refractive index control type definition */
#ifndef _OPENDDA_REFRACTIVE_INDEX_
#define _OPENDDA_REFRACTIVE_INDEX_
typedef struct {
/* flag: for each material 'j',j=0..(number_of_dielectric_materials-1), the flag
         indicates whether the value from the control.input file is being used
         (flag[j]=0) OR data for refractive index interpolation has been read from
         the input file INPUT_refractive_index_vs_wavelength_'j'.input (flag[j]=1)
   list_length: if the data for refractive index interpolation has been read from
         the input file INPUT_refractive_index_vs_wavelength_'j'.input (flag[j]=1)
         then list_length[j] = length of the wavelength and complex refractive index
         arrays for this material. However, if the value from the control.input file
         is being used (flag[j]=0) then list_length[j]=0   
   number_of_materials: number of dielectric materials for the target [MAX=8]
   *value: value[j] stores the refractive index[j] for the current wavelength
           If (flag[j]=1) then this value has been interpolated from the data for
           refractive index interpolation that was read from the input file
           INPUT_refractive_index_vs_wavelength_'j'.input
   **wavelengths,
   **refractive_indices: if (flag[j]=1) then the arrays wavelengths[j] and refractive_indices[j]
                        are of length list_length[j] and store the wavelength and complex
                        refractive index interpolation data for material 'j'. However, if
                        (flag[j]=0) then the value for the complex refractive index is the
                        value from the control.input file and no interpolation is required.
                        For this case, list_length[j]=0 and the pointers for wavelengths[j]
                        and refractive_indices[j] are 'NULL' */
   int *flag;
   int *list_length;
   int number_of_materials;
   fcomplex *value;
   float **wavelengths;
   fcomplex **refractive_indices;
} attributes_refractive_index;
#endif

/* Iterative control type definition */
#ifndef _OPENDDA_ITERATIVE_
#define _OPENDDA_ITERATIVE_
typedef struct {
/* Type definition for a structure to hold the control information for
   the iterative solution scheme

   maximum: Maximum number of iterations for the iterative solvers
   precond: Preconditioning type 0: None 1: Point-Jacobi
   initial_guess: Initial guess for the iterative solvers
                  0: set x=1
                  1: set x=b
                  2: set x=(1/alpha) i.e. (1/polarisability)
   vec: Number of starting vectors which is only relevant for
        'mlbicgstab_orig','mlbicgstab','mlbicgstab_ss' & 'rbicgstab'
   tolerance: Solution tolerance for the iterative solvers
   breakdown: Breakdown tolerance for the iterative solvers
   scheme[15]: Solution scheme
      bicg            BiConjugate-Gradients
      bicg_sym        BiConjugate-Gradients for symmetric systems
      bicgstab        Stabilised version of the BiConjugate-Gradients
      cg              Conjugate-Gradients
      cgs             Conjugate-Gradients Squared
      mlbicgstab_orig A variant of BiCGSTAB based on multiple Lanczos starting vectors
      (Author's original algorithm)
      mlbicgstab      A variant of BiCGSTAB based on multiple Lanczos starting vectors
                      (Author's reformulation)
      mlbicgstab_ss   A variant of BiCGSTAB based on multiple Lanczos starting vectors
                      (Author's space saving algorithm)
      qmr             Quasi-minimal residual with coupled two-term recurrences
      qmr_sym         Quasi-minimal residual for symmetric systems
      rbicgstab       Restarted, stabilised version of the BiConjugate-Gradients      
      tfqmr           Transpose-free quasi-minimal residual */
   int maximum;
   int precond;
   int initial_guess;
   int vec;
   float tolerance;
   float breakdown;
   char scheme[15];
} attributes_iterative;
#endif

/* Target description type definition */
#ifndef _OPENDDA_TARGET_
#define _OPENDDA_TARGET_
typedef struct {
/* Type definition for a structure to hold the target description

   K: number of sites in the x direction i.e. number of sites in line
   J: number of sites in the y direction i.e. number of lines in an xy-plane
   P: number of sites in the z direction i.e. number of xy-planes planes in the array
   M: M=K*J number of sites in an xy-plane
   N: N=K*J*P number of sites in the array
   Ka,kf: Ka=K+kf where `kf' is the zero-fill in the x-direction for to permit the use
          of fast DFT algorithms
   Ja,jf: Ja=J+jf where `jf' is the zero-fill in the y-direction for to permit the use
          of fast DFT algorithms
   Pa,pf: Pa=P+pf where `pf' is the zero-fill in the z-direction for to permit the use
          of fast DFT algorithms
   Ma: Ma=Ka*Ja
   Na: Ka*Ja*Pa

   Note: The Toeplitz-to-circulant conversion required to permit the use of DFTs for
         the multiplication of the interaction matrix and a given vector extends the
         6 independent tensor components arrays from [Ka*Ja*Pa] to [Kpp*Jpp*Ppp]. It
         also requires the 3 vector components to be zero-padded from [K*J*P] to
         [Kpp*Jpp*Ppp]

   Kp: [K'] Kp=2K-1
   Jp: [J'] Jp=2J-1
   Pp: [P'] Pp=2P-1
   Kpp: [K''] Kpp=Kp+kf
   Jpp: [J''] Jpp=Jp+jf
   Ppp: [P''] Ppp=Pp+pf
   Mpp: [M''] Mpp=Kpp*Jpp
   Mv: Mv=K*Jpp
   Qpp: [Q''] Qpp=Kpp*Ppp
   Npp: [N''] Npp=Kpp*Jpp*Ppp
   Nv: Nv=Kpp*Jpp*P
   *occupied: A bitfielded array that contains 0 or 1 flags to indicate whether lattice
              sites are occupied or not
              0: unoccupied
              1: occupied
   construction_parameters[6]: target construction parameters
   Nd: number of lattice sites actually occupied by dipoles
   composition[3]: Bit-fielded arrays which specifies the x, y and z material properties
                for the Nd dipoles that constitute the actual target. There are
                (refractive_index.number_of_materials) materials to choose from
   mask: Bitfielding mask for the composition array i.e. for 5 disparate materials in
         the target, 0..4=000..101 in binary, the extraction mask is 7 = 111 in binary 
   mask_size: Size of the bitfielding mask for composition array, which corresponds to
              number of bits required to store the integer range of allowable materials
              i.e. if there are 5 disparate materials in the target then 3 bits are
              required 0..4=000..101 in binary
   masks_per_integer: Number of data entries that can stored within a single array
   N3: N3=N*3
   Nd3: Nd3=Nd*3
   dipoles_per_wavelength: number of dipoles per wavelength
   dipole_spacing: inter-dipole spacing
   shape[50]: target description string
   from_file: target data read from file flag
                  0: Values from the control.input file were used
                  1: Data was read from file
   output_xyz_vtk: VTK target xyz coordinate data output flag.
                   0: No output
                   1: target xyz coordinates are output to a vtk file
                      for subsequent 3D visualisation */
   int K; /* number of sites in the x direction */
   int J; /* number of sites in the y direction */
   int P; /* number of sites in the z direction */
   int M; /* M=K*J number of sites in an xy-plane */
   long int N; /* N=K*J*P number of sites in the array */
   int Ka; /* Ka=K+kf */
   int Ja; /* Ja=J+jf */
   int Pa; /* Pa=P+pf */
   int Ma; /* Ma=Ka*Ja */
   long int Na; /* Ka*Ja*Pa */
   int Kp; /* Kp=2K-1 */
   int Jp; /* Jp=2J-1 */
   int Pp; /* Pp=2P-1 */
   int Kpp; /* Kpp=Kp+kf */
   int Jpp; /* Jpp=Jp+jf */
   int Ppp; /* Ppp=Pp+pf */
   int Mpp; /* Mpp=Kpp*Jpp */
   int Mv; /* Mv=K*Jpp */
   int Qpp; /* Qpp=Kpp*Ppp */
   long int Npp; /* Npp=Kpp*Jpp*Ppp */
   long int Nv; /* Nv=Kpp*Jpp*P */
   int kf; /* zero-fill in the x-direction */
   int jf; /* zero-fill in the y-direction */
   int pf; /* zero-fill in the z-direction */
   unsigned int *occupied; /* bitfielded array of 0/1 flags for whether a site is occupied by a dipole or not */
   int construction_parameters[6]; /* target construction parameters */
   long int Nd; /* number of lattice sites actually occupied by dipoles */
   unsigned int *composition[3]; /* bitfielded x, y and z material property for the Nd dipoles */
   int mask; /* Bitfielding mask for composition array */
   int mask_size; /* Size of the bitfielding mask for composition array */
   int masks_per_integer; /* Number of data entries that can stored within a single array */
   long int N3; /* N3=N*3 */
   long int Nd3; /* Nd3=Nd*3 */
   float dipoles_per_wavelength; /* number of dipoles per wavelength */
   float dipole_spacing; /* inter-dipole spacing */
   char shape[50]; /* target description string */
   int from_file; /* target data read from file flag */
   int output_xyz_vtk; /* VTK target xyz coordinate data output flag */
} attributes_target;
#endif

/* Output file type definition */
#ifndef _OPENDDA_OUTPUT_FILES_
#define _OPENDDA_OUTPUT_FILES_
typedef struct {
/* Type definition for a structure to control data output */
   FILE *eff_radius; /* Effective radius */
   FILE *incident_wave; /* Incident wave */
   FILE *iter; /* Iterative scheme */
   FILE **refindex_vs_wave; /* Complex refractive index versus wavelength */
   FILE *scatt_angle; /* Scattering angles */
   FILE *target_orien; /* Target orientations */
   FILE *target_descrip; /* Target description */
   FILE *time; /* Timing for the iterative kernel and DFT engine */
} attributes_output;
#endif

/* Incident polarisation type definition */
#ifndef _OPENDDA_INCIDENT_POLARISATION_
#define _OPENDDA_INCIDENT_POLARISATION_
typedef struct {
/* Type definition for a structure to control incident polarisation

   n_inc: incident direction
   polarisation[3]: Incident polarisation state
   orthonormal[3]: Orthonormal polarisation state */
   float n_inc[3]; /* Incident direction */
   fcomplex polarisation[3]; /* Incident polarisation state */
   fcomplex orthonormal[3]; /* Orthonormal polarisation state */
} attributes_incident_polarisation;
#endif

/* Scattering vector type definition */
#ifndef _OPENDDA_SCATTERING_VECTORS_
#define _OPENDDA_SCATTERING_VECTORS_
typedef struct {
/* Type definition for a structure to control scattering vectors

   n_sca: current scattering direction
   phi_sca: current scattering vector perpendicular to
            the scattering plane
   theta_sca: current scattering vector parallel to the
              scattering plane */
   float n_sca[3]; /* Scattering direction */
   float phi_sca[3]; /* Scattering vector perpendicular to the scattering plane */
   float theta_sca[3]; /* Scattering vector parallel to the scattering plane */
} attributes_scattering_vectors;
#endif

/* Cross section type definition */
#ifndef _OPENDDA_CROSS_SECTIONS_
#define _OPENDDA_CROSS_SECTIONS_
typedef struct {
/* Type definition for a structure to control cross sections

   extinction: extinction cross section
   absorption: absorption cross section
   scattering: scattering cross section
   geometrical: geometrical cross section */
   float extinction; /* Extinction cross section */
   float absorption; /* Absorption cross section */
   float scattering; /* Scattering cross section */
   float geometrical; /* Geometrical cross section */
} attributes_cross_sections;
#endif

/* Efficiency type definition */
#ifndef _OPENDDA_EFFICIENCY_
#define _OPENDDA_EFFICIENCY_
typedef struct {
/* Type definition for a structure to control efficiencies

   extinction: extinction efficiency
   absorption: absorption efficiency
   scattering: scattering efficiency */
   float extinction; /* Extinction efficiency */
   float absorption; /* Absorption efficiency */
   float scattering; /* Scattering efficiency */
} attributes_efficiency;
#endif

/* Timing type definition */
#ifndef _OPENDDA_TIMING_
#define _OPENDDA_TIMING_
typedef struct {
/* Type definition for a structure to control timing
   
   overall[2]: Overall runtime of the code, start and finish
   enabled: Whether timing is enabled or not 1: enabled 0: disabled
   iterative_kernel_count: Number of times the iterative kernel was executed
   dft_engine_count: Number of times the DFT engine was executed
   dft_tensor_components_count: Number of times the DFT was applied to the 6
                                independent tensor components. Note: the interaction
                                matrix only changes if the wave vector changes which
                                means that although this component is part of the
                                DFT engine it may not be executed every time the
                                engine is                    
   zero_padded_vector: Time to initialise the zero padded vector
   dft_tensor_components: Time for the forward DFT of the 6 independent tensor
                          components of the interaction matrix
   zdft_zidft_and_fd_mult: Time for the DFT of the 3 vector components, the
                           subsequent Fourier domain multplication and the iDFT
                           of the resultant
   norm_diag_compress: Time for the normalisation of the resultant of the iDFT + 
                       Time to include the diagonal contributions +
                       Time to compress data to occupied sites only
   iterative_solution: Time for the iterative solution */
   double overall[2];
   int enabled;
   int iterative_kernel_count;
   int dft_engine_count;
   int dft_tensor_components_count;
   double zero_padded_vector;
   double dft_tensor_components;
   double zdft_zidft_and_fd_mult;
   double norm_diag_compress;
   double iterative_solution;
} attributes_timing;
#endif
