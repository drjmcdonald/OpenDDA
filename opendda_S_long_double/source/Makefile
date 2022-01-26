CC = icc
CFLAGS = -O3 -ipo -static -long_double -no-gcc -Wall -wd981 -wd810 -mp
INCFFTW = -I/opt/fftw-3.1.2-intel-long-double/include
LIBFFTW = -L/opt/fftw-3.1.2-intel-long-double/lib
LDFLAGS = -lfftw3l -lm

.c.o:
	$(CC) $(CFLAGS) -c $(INCFFTW) $<

SRCS0 = opendda_S.c \
			bitfielding.c \
			build_target.c \
			cross_sections_S.c \
			dda_interaction_matrix_S.c \
			degree_radian_conversion.c \
			dft_plan.c \
			dftmatvec_S.c \
			dipole_polarisabilities_S.c \
			efficiencies.c \
			file_alloc.c \
			find_dft_size.c \
			incident_electric_field_S.c \
			initialise_output.c \
			int_alloc.c \
			interpolation.c \
			iterative_solution.c \
			ldcomplex_alloc.c \
			ldcomplex_bicg_S.c \
			ldcomplex_bicg_sym_S.c \
			ldcomplex_bicgstab_S.c \
			ldcomplex_cg_S.c \
			ldcomplex_cgs_S.c \
			ldcomplex_math_S.c \
			ldcomplex_mlbicgstab_S.c \
			ldcomplex_mlbicgstab_orig_S.c \
			ldcomplex_mlbicgstab_ss_S.c \
			ldcomplex_qmr_S.c \
			ldcomplex_qmr_sym_S.c \
			ldcomplex_rbicgstab_S.c \
			ldcomplex_tfqmr_S.c \
			long_double_alloc.c \
			memory_allocate_S.c \
			memory_free_S.c \
			preconditioner_S.c \
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

SRCS1 = guesstimate_S.c \
			find_dft_size.c \
			print_details.c \
			reset_string.c

OBJS0=$(SRCS0:.c=.o)
OBJS1=$(SRCS1:.c=.o)

opendda_S: $(OBJS0) Makefile
	$(CC) $(CFLAGS) $(OBJS0) $(INCFFTW) $(LIBFFTW) $(LDFLAGS) -o $@
	if [ -e control.input ]; then cp control.input control.input.copy; fi

guesstimate_S: $(OBJS1) Makefile
	$(CC) $(CFLAGS) $(OBJS1) -o $@ 

all: opendda_S guesstimate_S

infoclean:; rm -f *.info *.vtk

.PHONY: clean
clean:; rm -f *.o opendda_S guesstimate_S *.info *.vtk
