/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for dcomplex_alloc

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include "dcomplex_type.h"
#include "print_details.h"

void dcomplex_malloc(dcomplex **p2p,size_t order);
void dcomplex_fftw_malloc(dcomplex **p2p,size_t order);
void dcomplex_calloc(dcomplex **p2p,size_t order);
void dcomplex_malloc_ptr(dcomplex ***p2p2p,size_t order);
void dcomplex_fftw_malloc_ptr(dcomplex ***p2p2p,size_t order);
void dcomplex_free(dcomplex **p2p);
void dcomplex_fftw_free(dcomplex **p2p);
void dcomplex_free_ptr(dcomplex ***p2p2p);
void dcomplex_fftw_free_ptr(dcomplex ***p2p2p);
