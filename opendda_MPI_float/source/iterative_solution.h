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
#include "fcomplex_bicg_MPI.h"
#include "fcomplex_bicg_sym_MPI.h"
#include "fcomplex_bicgstab_MPI.h"
#include "fcomplex_cg_MPI.h"
#include "fcomplex_cgs_MPI.h"
#include "fcomplex_mlbicgstab_MPI.h"
#include "fcomplex_mlbicgstab_orig_MPI.h"
#include "fcomplex_mlbicgstab_ss_MPI.h"
#include "fcomplex_qmr_MPI.h"
#include "fcomplex_qmr_sym_MPI.h"
#include "fcomplex_rbicgstab_MPI.h"
#include "fcomplex_tfqmr_MPI.h"
#include "preconditioner_MPI.h"

void iterative_solution(void);
