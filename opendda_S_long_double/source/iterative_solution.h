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
#include "ldcomplex_bicg_S.h"
#include "ldcomplex_bicg_sym_S.h"
#include "ldcomplex_bicgstab_S.h"
#include "ldcomplex_cg_S.h"
#include "ldcomplex_cgs_S.h"
#include "ldcomplex_mlbicgstab_S.h"
#include "ldcomplex_mlbicgstab_orig_S.h"
#include "ldcomplex_mlbicgstab_ss_S.h"
#include "ldcomplex_qmr_S.h"
#include "ldcomplex_qmr_sym_S.h"
#include "ldcomplex_rbicgstab_S.h"
#include "ldcomplex_tfqmr_S.h"
#include "preconditioner_S.h"

void iterative_solution(void);
