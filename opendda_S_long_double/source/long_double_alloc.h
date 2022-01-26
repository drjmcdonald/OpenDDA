/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for long_double_alloc

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include "print_details.h"

void long_double_malloc(long double **p2p,size_t order);
void long_double_calloc(long double **p2p,size_t order);
void long_double_malloc_ptr(long double ***p2p2p,size_t order);
void long_double_free(long double **p2p);
void long_double_free_ptr(long double ***p2p2p);
