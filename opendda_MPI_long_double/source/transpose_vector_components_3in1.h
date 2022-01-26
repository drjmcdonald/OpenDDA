/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for transpose_vector_component_3in1

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "ldcomplex_type.h"
#include "ldcomplex_alloc.h"
#include "definitions.h"
#include "local_transpose_vector_component.h"
#include "mpi_request_alloc.h"

void transpose_vector_component_3in1_0(int before,int after,int limitbefore,int limitafter);
void transpose_vector_component_3in1_1(int before,int after,int planesbefore,int planesafter,int recv);
void transpose_vector_component_3in1_2(int before,int after,int planesbefore,int planesafter,int recv);
