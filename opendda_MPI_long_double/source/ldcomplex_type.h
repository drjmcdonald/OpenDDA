/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Complex type definition

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#ifndef _LDCOMPLEX_
#define _LDCOMPLEX_
#include <mpi.h>
typedef struct{
   long double dat[2];
} ldcomplex;
MPI_Datatype ldcomplex_type; /* MPI data type for ldcomplex */
#endif
