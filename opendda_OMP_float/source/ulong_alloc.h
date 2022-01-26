/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for ulong_alloc

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include "print_details.h"

void ulong_malloc(unsigned long **p2p,size_t order);
void ulong_calloc(unsigned long **p2p,size_t order);
void ulong_free(unsigned long **p2p);
