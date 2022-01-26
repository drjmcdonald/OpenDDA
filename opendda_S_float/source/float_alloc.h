/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for float_alloc

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include "print_details.h"

void float_malloc(float **p2p,size_t order);
void float_calloc(float **p2p,size_t order);
void float_malloc_ptr(float ***p2p2p,size_t order);
void float_free(float **p2p);
void float_free_ptr(float ***p2p2p);
