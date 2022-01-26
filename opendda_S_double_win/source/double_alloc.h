/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for double_alloc

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include "print_details.h"

void double_malloc(double **p2p,size_t order);
void double_calloc(double **p2p,size_t order);
void double_malloc_ptr(double ***p2p2p,size_t order);
void double_free(double **p2p);
void double_free_ptr(double ***p2p2p);
