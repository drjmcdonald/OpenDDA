/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Complex type definition

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#ifndef _FCOMPLEX_
#define _FCOMPLEX_
#include <mpi.h>
typedef struct{
   float dat[2];
} fcomplex;
MPI_Datatype fcomplex_type; /* MPI data type for fcomplex */
#endif
