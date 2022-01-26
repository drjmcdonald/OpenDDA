/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for uint_alloc

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include "print_details.h"

void uint_malloc(unsigned int **p2p,size_t order);
void uint_calloc(unsigned int **p2p,size_t order);
void uint_malloc_ptr(unsigned int ***p2p2p,size_t order);
void uint_free(unsigned int **p2p);
void uint_free_ptr(unsigned int ***p2p2p);
