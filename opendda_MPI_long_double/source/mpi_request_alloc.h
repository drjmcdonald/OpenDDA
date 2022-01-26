/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for mpi_request_alloc

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "print_details.h"

void mpi_request_malloc(MPI_Request **p2p2p,size_t order);
void mpi_request_calloc(MPI_Request **p2p2p,size_t order);
void mpi_request_free(MPI_Request **p2p2p);
