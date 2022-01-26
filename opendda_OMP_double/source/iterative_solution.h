/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for iterative_solution

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <string.h>
#include "definitions.h"
#include "walltime.h"
#include "dcomplex_bicg_OMP.h"
#include "dcomplex_bicg_sym_OMP.h"
#include "dcomplex_bicgstab_OMP.h"
#include "dcomplex_cg_OMP.h"
#include "dcomplex_cgs_OMP.h"
#include "dcomplex_mlbicgstab_OMP.h"
#include "dcomplex_mlbicgstab_orig_OMP.h"
#include "dcomplex_mlbicgstab_ss_OMP.h"
#include "dcomplex_qmr_OMP.h"
#include "dcomplex_qmr_sym_OMP.h"
#include "dcomplex_rbicgstab_OMP.h"
#include "dcomplex_tfqmr_OMP.h"
#include "preconditioner_OMP.h"

void iterative_solution(void);
