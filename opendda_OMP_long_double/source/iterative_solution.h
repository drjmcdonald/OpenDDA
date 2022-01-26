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
#include "ldcomplex_bicg_OMP.h"
#include "ldcomplex_bicg_sym_OMP.h"
#include "ldcomplex_bicgstab_OMP.h"
#include "ldcomplex_cg_OMP.h"
#include "ldcomplex_cgs_OMP.h"
#include "ldcomplex_mlbicgstab_OMP.h"
#include "ldcomplex_mlbicgstab_orig_OMP.h"
#include "ldcomplex_mlbicgstab_ss_OMP.h"
#include "ldcomplex_qmr_OMP.h"
#include "ldcomplex_qmr_sym_OMP.h"
#include "ldcomplex_rbicgstab_OMP.h"
#include "ldcomplex_tfqmr_OMP.h"
#include "preconditioner_OMP.h"

void iterative_solution(void);
