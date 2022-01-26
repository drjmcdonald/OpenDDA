/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for broadcast_parameters

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "dcomplex_alloc.h"
#include "definitions.h"
#include "int_alloc.h"
#include "uint_alloc.h"

void broadcast_control_parameters(void);
void broadcast_length_parameters(void);
void broadcast_physical_parameters(void);
void broadcast_target_from_file(void);
