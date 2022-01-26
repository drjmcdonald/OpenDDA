/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for iterative_solution

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "definitions.h"
#include "dcomplex_bicg_MPI.h"
#include "dcomplex_bicg_sym_MPI.h"
#include "dcomplex_bicgstab_MPI.h"
#include "dcomplex_cg_MPI.h"
#include "dcomplex_cgs_MPI.h"
#include "dcomplex_mlbicgstab_MPI.h"
#include "dcomplex_mlbicgstab_orig_MPI.h"
#include "dcomplex_mlbicgstab_ss_MPI.h"
#include "dcomplex_qmr_MPI.h"
#include "dcomplex_qmr_sym_MPI.h"
#include "dcomplex_rbicgstab_MPI.h"
#include "dcomplex_tfqmr_MPI.h"
#include "preconditioner_MPI.h"

void iterative_solution(void);
