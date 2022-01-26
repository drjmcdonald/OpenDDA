/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for int_alloc

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include "print_details.h"

void int_malloc(int **p2p,size_t order);
void int_calloc(int **p2p,size_t order);
void int_malloc_ptr(int ***p2p2p,size_t order);
void int_free(int **p2p);
void int_free_ptr(int ***p2p2p);
