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
   Q: Q=K*P number of sites in an xz-plane
   N: N=K*J*P number of sites in the array
   Ka,kf: Ka=K+kf where `kf' is the zero-fill in the x-direction for to permit the use
          of fast DFT algorithms
   Ja,jf: Ja=J+jf where `jf' is the zero-fill in the y-direction for to permit the use
          of fast DFT algorithms
   Pa,pf: Pa=P+pf where `pf' is the zero-fill in the z-direction for to permit the use
          of fast DFT algorithms
   Ma: Ma=Ka*Ja
   Qa: Qa=Ka*Pa
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
   Qpp: [M''] Qpp=Kpp*Ppp
   Npp: [N''] Npp=Kpp*Jpp*Ppp
   KaJpp: KaJpp=Ka*Jpp
   PKpp: PKpp=P*Kpp
   occ_ctrl: Number of integers required to store `occupied' bit-fielded array
   occupied: A bitfielded array that contains 0 or 1 flags to indicate whether lattice
              sites are occupied or not
              0: unoccupied
              1: occupied
   populated: To conserve memory array elements that correspond to vacant lattice sites
              and are thus zero are not stored. The populated array stores the starting
              index for the non-zero elements for each of the P xy-planes
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
   int Q; /* Q=K*P number of sites in an xz-plane */
   long int N; /* N=K*J*P number of sites in the array */
   int Ka; /* Ka=K+kf */
   int Ja; /* Ja=J+jf */
   int Pa; /* Pa=P+pf */
   int Ma; /* Ma=Ka*Ja */
   int Qa; /* Qa=Ka*Pa */
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
   int KaJpp; /* KaJpp=Ka*Jpp */
   int PKpp; /* PKpp=P*Kpp */
   int kf; /* zero-fill in the x-direction */
   int jf; /* zero-fill in the y-direction */
   int pf; /* zero-fill in the z-direction */
   int occ_ctrl; /* Number of integers required to store `occupied' bit-fielded array */
   unsigned int *occupied; /* bitfielded array of 0/1 flags for whether a site is occupied by a dipole or not */
   int *populated; /* Stores how many lattice sites are occupied in each of the P xy-planes */
   int construction_parameters[6]; /* target construction parameters */
   long int Nd; /* number of lattice sites actually occupied by dipoles */
   int mat_ctrl; /* Number of integers required to store one component of the `composition' bit-fielded array */
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
   transpose_tensor_forward: Time for the forward-distributed transpose of the 6
                             independent tensor components
   transpose_vector_forward: Time for the forward-distributed transpose of the 3
                             zero-padded vector components
   transpose_vector_reverse: Time for the reverse-distributed transpose of the 3
                             zero-padded vector components
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
   double transpose_tensor_forward;
   double transpose_vector_forward;
   double transpose_vector_reverse;
   double zdft_zidft_and_fd_mult;
   double norm_diag_compress;
   double iterative_solution;
} attributes_timing;
#endif

/* Parallel type definition */
#ifndef _OPENDDA_PARALLEL_
#define _OPENDDA_PARALLEL_
typedef struct {
/* Type definition for a structure to control the parallel
   aspects of the code execution

   alloc_tensor: The 6 independent tensor components of the interaction
                 matrix are stored in parallel so `alloc_tensor' is the
                 memory allocated on each processor to store its portion of
                 the interaction matrix before OR after the distributed
                 xy-to-xz plane transpose 
   alloc_padded: The 3 components of the zero-padded vector components
                 are stored in parallel so `alloc_alloc_padded' is the
                 memory allocated on each processor to store its
                 portion of the zero-padded verctor before OR after the
                 distributed xy-to-xz plane transpose
   alloc_vector: The 3 components of all the basic vector components
                 are stored in parallel so `alloc_vector' is the
                 memory allocated on each processor to store its
                 portion of the vector. N.B. Only stores information for
                 occupied sites
   alloc_vector3: alloc_vector*3
   np,
   dnp: `np' is the number of processors and `dnp' is it cast to a double
   myrank: MPI processor rank 0(master) to np-1

   BEFORE XY-TO-XZ DISTRIBUTED TRANSPOSE
   Znp_tensor: Base number of xy-planes per proc before the transpose.
               For the interaction matrix (gets transposed)
   Znp_padded: Base number of xy-planes per proc before the transpose.
               For the zero padded vector (gets transposed)
   Znp_vector: Base number of xy-planes per proc before the transpose.
               For all other vectors (does not get transposed)
   AFTER XY-TO-XZ DISTRIBUTED TRANSPOSE
   Ynp: Base number of xz-planes per proc after the transpose.
        For the interaction matrix & the zero padded vector
   plane: Since the system is `plane' decomposed across the
          available processors, plane k, k=0..P, is no longer
          plane k but is plane ((int)((double)k/dnp)) on
          processor k mod np.
   planesY_padded: Number of xy-planes locally on the processor before the transpose
                   i.e., Znp_padded+1 if myrank < P%np, Znp_padded if myrank >= P%np
   limitY_padded: Number of lines in the x-direction locally before the transpose
                  i.e., planesY_padded*Jpp
   planesY_tensor: Number of xy-planes locally on the processor before the transpose
                   i.e., Znp_tensor+1 if myrank < Pa%np, Znp_tensor if myrank >= Pa%np
   limitY_tensor: Number of lines in the x-direction locally before the
                  transpose i.e., planesY_tensor*Jpp, N.B. For this MPI version the
                  6 independent tensor components of the interaction matrix are stored
                  in interaction_matrix[6][Ka*Jpp*Pa] rather than the
                  interaction_matrix[6][Ka*Ja*Pa] stored for the serial and OpenMP cases.
                  This is to facilitate the distributed FD tensor-vector multiplication.
   planesZ_padded: Number of xz-planes locally on the processor after the
                   transpose i.e., Ynp+1 if myrank < Jpp%np, Ynp if myrank >= Jpp%np
   limitZ_padded: Number of lines in the x-direction locally after the transpose i.e.,
                  planesZ_padded*P for zero_padded_vector
   planesZ_tensor: Number of xz-planes locally on the processor after the transpose i.e.,
                   Ynp+1 if myrank < Jpp%np, Ynp if myrank >= Jpp%np
   limitZ_tensor: Number of lines in the x-direction locally before the transpose i.e.,
                  planesY_tensor*Pa, N.B. For this MPI version the 6 independent tensor
                  components of the interaction matrix are stored in
                  interaction_matrix[6][Ka*Jpp*Pa] rather than the
                  interaction_matrix[6][Ka*Ja*Pa] stored for the serial and OpenMP cases.
                  This is to facilitate the distributed FD tensor-vector multiplication.
   planesY_vector: Number of xy-planes locally on the processor for all non-transposed
                   vectors i.e., Znp_vector+1 if myrank < P%np, Znp_vector if myrank >= P%np
   tensor_transpose: Transpose algorithm to use for the 6 independent tensor components, 0..5
   padded_transpose: Transpose algorithm to use for the zero-padded vector components 0..5
   local_tensor: Temporary scratch array for local transposes of the 6 independent tensor components
                 of the interaction matrix
   tf_ctrl_tensor: Required size of the `transposed_flag_tensor' array
   transposed_flag_tensor: Bit-fielded array to specify if a line has been locally transposed or not
                           for the 6 independent tensor components of the interaction matrix
   local_padded: Temporary scratch array for local transposes of the zero-padded vector array
   tf_ctrl_padded: Required size of the `transposed_flag_padded' array
   transposed_flag_padded: Bit-fielded array to specify if a line has been locally transposed or not
                           for the zero-padded vector array */
   int alloc_tensor;
   int alloc_padded;
   int alloc_vector;
   int alloc_vector3;
   int np;
   double dnp;
   int myrank;
   int Ynp;
   int Znp_tensor;
   int Znp_padded;
   int Znp_vector;
   int plane;
   int limitY_tensor;
   int limitZ_tensor;
   int planesY_tensor;
   int planesZ_tensor;
   int limitY_padded;
   int limitZ_padded;
   int planesY_padded;
   int planesZ_padded;
   int planesY_vector;
   int tensor_transpose;
   int padded_transpose;
   fcomplex *local_tensor;
   fcomplex *local_padded;
   int tf_ctrl_tensor;
   int tf_ctrl_padded;
   int *transposed_flag_tensor;
   int *transposed_flag_padded;
} attributes_parallel;
#endif
