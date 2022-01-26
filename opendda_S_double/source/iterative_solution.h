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
#include "dcomplex_bicg_S.h"
#include "dcomplex_bicg_sym_S.h"
#include "dcomplex_bicgstab_S.h"
#include "dcomplex_cg_S.h"
#include "dcomplex_cgs_S.h"
#include "dcomplex_mlbicgstab_S.h"
#include "dcomplex_mlbicgstab_orig_S.h"
#include "dcomplex_mlbicgstab_ss_S.h"
#include "dcomplex_qmr_S.h"
#include "dcomplex_qmr_sym_S.h"
#include "dcomplex_rbicgstab_S.h"
#include "dcomplex_tfqmr_S.h"
#include "preconditioner_S.h"

void iterative_solution(void);
