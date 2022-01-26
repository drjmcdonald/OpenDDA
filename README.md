# OpenDDA
OpenDDA is a highly optimised computational framework, written in the C language, for the Discrete Dipole Approximation, a numerical method for calculating the optical properties associated with a target of arbitrary geometry that is widely used in atmospheric, astrophysical and industrial simulations.

Core optimisations include the bit-fielding of integer data and iterative methods that complement a new Discrete Fourier Transform (DFT) kernel, which efficiently calculates the matrix-vector products required by these iterative solution schemes. The new kernel performs the requisite 3D DFTs as ensembles of 1D transforms, and by doing so, is able to reduce the number of constituent 1D transforms by 60% and the memory by over 80%. The optimisations also facilitate the use of parallel techniques to further enhance the performance. Complete OpenMP-based shared-memory and MPI-based distributed-memory implementations have been created to take full advantage of the various architectures.

OpenDDA is free to use for everyone. The only request pertaining to its use is that any subsequent work, use, modification and / or augmentation of the available software, should cite the associated publication: OpenDDA: A Novel High-Performance Computational Framework for the Discrete Dipole Approximation, Mc Donald, J., Golden, A., Jennings, S. G., International Journal of High Performance Computing Applications (IJHPCA), Volume 23, No. 1, 42-61, 2009.

The are two associated documents that are useful when working with the OpenDDA framework:

1. The associated publication:
OpenDDA: A Novel High-Performance Computational Framework for the Discrete Dipole Approximation, Mc Donald, J., Golden, A., Jennings, S. G., International Journal of High Performance Computing Applications (IJHPCA), Volume 23, No. 1, 42-61, 2009.
http://hpc.sagepub.com/cgi/content/abstract/23/1/42

2. The associated Ph.D. Thesis:
OpenDDA: A Novel High-Performance Computational Framework for the Discrete Dipole Approximation, Mc Donald, J., Ph.D. Thesis, School of Physics, National University of Ireland, Galway, Ireland, September 2007.

There are 9 different LINUX versions of the OpenDDA framework

3 different architecture specific versions (each with 3 precision specific variants, i.e., float, double and long double):

OpenDDA_S is the basic serial implementation
OpenDDA_OMP is the OpenMP-based shared-memory implementation
OpenDDA_MPI is the MPI-based distributed-memory implementation.

The WINDOWS versions (Visual Studio 2008 project) exist in SERIAL implmentations only (float, double and long double). N.B. Note that these WINDOWS versions are completely untested, and have been supplied, including fftw-3.3.2-dll32, for convenience as a starting point.

For the WINDOWS versions, the executable is NOT in the source directory. It is created in the Debug or Release directory as appropriate. The ".input" files MUST be modified in this Debug or Release directory as appropriate.
