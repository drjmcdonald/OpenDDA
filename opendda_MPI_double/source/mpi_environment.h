/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for mpi_environment

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <mpi.h>
#include "dcomplex_type.h"
#include "definitions.h"

void mpi_initialise(int *ac,char ***av);
void mpi_finalise(void);
int dcomplex_sum(dcomplex *in,dcomplex *inout,int *len,MPI_Datatype *dtptr);
