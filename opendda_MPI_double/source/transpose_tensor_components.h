/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for transpose_tensor_component

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "dcomplex_type.h"
#include "dcomplex_alloc.h"
#include "definitions.h"
#include "local_transpose_tensor_component.h"
#include "mpi_request_alloc.h"

void transpose_tensor_component_0(int before,int after,int limitbefore,int limitafter,int comp);
void transpose_tensor_component_1(int before,int after,int planesbefore,int planesafter,int recv,int comp);
void transpose_tensor_component_2(int before,int after,int planesbefore,int planesafter,int recv,int comp);
