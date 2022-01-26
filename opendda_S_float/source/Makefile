CC = icc
CFLAGS = -O3 -ipo -static -no-gcc -Wall -wd981 -wd810 -mp
INCFFTW = -I/opt/fftw-3.1.2-intel-float/include
LIBFFTW = -L/opt/fftw-3.1.2-intel-float/lib
LDFLAGS = -lfftw3f -lm

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
			fcomplex_alloc.c \
			fcomplex_bicg_S.c \
			fcomplex_bicg_sym_S.c \
			fcomplex_bicgstab_S.c \
			fcomplex_cg_S.c \
			fcomplex_cgs_S.c \
			fcomplex_math_S.c \
			fcomplex_mlbicgstab_S.c \
			fcomplex_mlbicgstab_orig_S.c \
			fcomplex_mlbicgstab_ss_S.c \
			fcomplex_qmr_S.c \
			fcomplex_qmr_sym_S.c \
			fcomplex_rbicgstab_S.c \
			fcomplex_tfqmr_S.c \
			file_alloc.c \
			find_dft_size.c \
			float_alloc.c \
			incident_electric_field_S.c \
			initialise_output.c \
			int_alloc.c \
			interpolation.c \
			iterative_solution.c \
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
