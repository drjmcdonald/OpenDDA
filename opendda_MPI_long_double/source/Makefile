CC = mpicc
CFLAGS = -O3 -ipo -long_double -no-gcc -Wall -wd981 -wd869 -wd810 -mp
INCFFTW = -I/opt/fftw-3.1.2-intel-long-double/include
LIBFFTW = -L/opt/fftw-3.1.2-intel-long-double/lib
INCMPI = -I/opt/mpich2-1.0.6p1-intel/include
LIBMPI = -L/opt/mpich2-1.0.6p1-intel/lib
LDFLAGS = -lfftw3l -lmpich -lm

.c.o:
	$(CC) $(CFLAGS) -c $(INCFFTW) $(INCMPI) $<

SRCS0 = opendda_MPI.c \
			bitfielding.c \
			broadcast_parameters.c \
			build_target.c \
			cross_sections_MPI.c \
			dda_interaction_matrix_MPI.c \
			degree_radian_conversion.c \
			dft_plan.c \
			dftmatvec_MPI.c \
			dipole_polarisabilities_MPI.c \
			efficiencies.c \
			file_alloc.c \
			find_dft_size.c \
			incident_electric_field_MPI.c \
			initialise_output.c \
			int_alloc.c \
			interpolation.c \
			iterative_solution.c \
			ldcomplex_alloc.c \
			ldcomplex_bicg_MPI.c \
			ldcomplex_bicg_sym_MPI.c \
			ldcomplex_bicgstab_MPI.c \
			ldcomplex_cg_MPI.c \
			ldcomplex_cgs_MPI.c \
			ldcomplex_math_MPI.c \
			ldcomplex_mlbicgstab_MPI.c \
			ldcomplex_mlbicgstab_orig_MPI.c \
			ldcomplex_mlbicgstab_ss_MPI.c \
			ldcomplex_qmr_MPI.c \
			ldcomplex_qmr_sym_MPI.c \
			ldcomplex_rbicgstab_MPI.c \
			ldcomplex_tfqmr_MPI.c \
			local_transpose_tensor_component.c \
			local_transpose_vector_component.c \
			long_double_alloc.c \
			memory_allocation_sizes.c \
			mpi_environment.c \
			mpi_request_alloc.c \
			memory_allocate_MPI.c \
			memory_free_MPI.c \
			preconditioner_MPI.c \
			print_details.c \
			read_parameters.c \
			reset_string.c \
			rotation.c \
			set_domain.c \
			set_initial_guess.c \
			sfmt_rng.c \
			target_cuboid.c \
			target_ellipsoid.c \
			target_from_file.c \
			target_xyz_vtk.c \
			timing_info.c \
			transpose_tensor_components.c \
			transpose_tensor_components_6in1.c \
			transpose_vector_components.c \
			transpose_vector_components_3in1.c \
			uint_alloc.c \
			ulong_alloc.c

SRCS1 = guesstimate_MPI.c \
			find_dft_size.c \
			memory_allocation_sizes.c \
			mpi_environment.c \
			print_details.c \
			reset_string.c \
			uint_alloc.c

OBJS0=$(SRCS0:.c=.o)
OBJS1=$(SRCS1:.c=.o)

opendda_MPI: $(OBJS0) Makefile
	$(CC) $(CFLAGS) $(OBJS0) $(INCFFTW) $(LIBFFTW) $(INCMPI) $(LIBMPI) $(LDFLAGS) -o $@
	if [ -e control.input ]; then cp control.input control.input.copy; fi

guesstimate_MPI: $(OBJS1) Makefile
	$(CC) $(CFLAGS) $(OBJS1) -o $@ 

all: opendda_MPI guesstimate_MPI

infoclean:; rm -f *.info *.vtk

.PHONY: clean
clean:; rm -f *.o opendda_MPI guesstimate_MPI *.info *.vtk
