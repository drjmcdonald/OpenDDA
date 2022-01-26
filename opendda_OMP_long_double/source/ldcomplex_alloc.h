/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for ldcomplex_alloc

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include "ldcomplex_type.h"
#include "print_details.h"

void ldcomplex_malloc(ldcomplex **p2p,size_t order);
void ldcomplex_fftwl_malloc(ldcomplex **p2p,size_t order);
void ldcomplex_calloc(ldcomplex **p2p,size_t order);
void ldcomplex_malloc_ptr(ldcomplex ***p2p2p,size_t order);
void ldcomplex_fftwl_malloc_ptr(ldcomplex ***p2p2p,size_t order);
void ldcomplex_free(ldcomplex **p2p);
void ldcomplex_fftwl_free(ldcomplex **p2p);
void ldcomplex_free_ptr(ldcomplex ***p2p2p);
void ldcomplex_fftwl_free_ptr(ldcomplex ***p2p2p);
