/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for dftmatvec_MPI

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <string.h>
#include <fftw3.h>
#include <mpi.h>
#include "ldcomplex_type.h"
#include "ldcomplex_math_MPI.h"
#include "definitions.h"
#include "transpose_tensor_components.h"
#include "transpose_vector_components.h"
#include "transpose_tensor_components_6in1.h"
#include "transpose_vector_components_3in1.h"

void dftmatvec_MPI(ldcomplex *vector);
