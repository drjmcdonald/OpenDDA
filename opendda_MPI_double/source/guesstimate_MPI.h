/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for guesstimate_MPI

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <ctype.h>
#include "attributes_type.h"
#include "dcomplex_type.h"
#include "find_dft_size.h"
#include "memory_allocation_sizes.h"
#include "mpi_environment.h"
#include "print_details.h"
#include "reset_string.h"
#include "uint_alloc.h"
