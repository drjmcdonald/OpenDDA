CC = icc
CFLAGS = -O3 -ipo -static -long_double -no-gcc -openmp -openmp-report0 -Wall -wd981 -wd810 -mp
INCFFTW = -I/opt/fftw-3.1.2-intel-long-double/include
LIBFFTW = -L/opt/fftw-3.1.2-intel-long-double/lib
LDFLAGS = -lfftw3l -lm

.c.o:
	$(CC) $(CFLAGS) -c $(INCFFTW) $<

SRCS0 = opendda_OMP.c \
			bitfielding.c \
			build_target.c \
			cross_sections_OMP.c \
			dda_interaction_matrix_OMP.c \
			degree_radian_conversion.c \
			dft_plan.c \
			dftmatvec_OMP.c \
			dipole_polarisabilities_OMP.c \
			efficiencies.c \
			file_alloc.c \
			find_dft_size.c \
			incident_electric_field_OMP.c \
			initialise_output.c \
			int_alloc.c \
			interpolation.c \
			iterative_solution.c \
			ldcomplex_alloc.c \
			ldcomplex_bicg_OMP.c \
			ldcomplex_bicg_sym_OMP.c \
			ldcomplex_bicgstab_OMP.c \
			ldcomplex_cg_OMP.c \
			ldcomplex_cgs_OMP.c \
			ldcomplex_math_OMP.c \
			ldcomplex_mlbicgstab_OMP.c \
			ldcomplex_mlbicgstab_orig_OMP.c \
			ldcomplex_mlbicgstab_ss_OMP.c \
			ldcomplex_qmr_OMP.c \
			ldcomplex_qmr_sym_OMP.c \
			ldcomplex_rbicgstab_OMP.c \
			ldcomplex_tfqmr_OMP.c \
			long_double_alloc.c \
			memory_allocate_OMP.c \
			memory_free_OMP.c \
			preconditioner_OMP.c \
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
			uint_alloc.c \
			ulong_alloc.c \
			walltime.c

SRCS1 = guesstimate_OMP.c \
			find_dft_size.c \
			print_details.c \
			reset_string.c

OBJS0=$(SRCS0:.c=.o)
OBJS1=$(SRCS1:.c=.o)

opendda_OMP: $(OBJS0) Makefile
	$(CC) $(CFLAGS) $(OBJS0) $(INCFFTW) $(LIBFFTW) $(LDFLAGS) -o $@
	if [ -e control.input ]; then cp control.input control.input.copy; fi

guesstimate_OMP: $(OBJS1) Makefile
	$(CC) $(CFLAGS) $(OBJS1) -o $@ 

all: opendda_OMP guesstimate_OMP

infoclean:; rm -f *.info *.vtk

.PHONY: clean
clean:; rm -f *.o opendda_OMP guesstimate_OMP *.info *.vtk
