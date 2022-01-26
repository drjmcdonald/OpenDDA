/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Complex type definition

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#ifndef _DCOMPLEX_
#define _DCOMPLEX_
#include <mpi.h>
typedef struct{
   double dat[2];
} dcomplex;
MPI_Datatype dcomplex_type; /* MPI data type for dcomplex */
#endif
