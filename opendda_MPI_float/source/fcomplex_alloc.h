/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for fcomplex_alloc

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include "fcomplex_type.h"
#include "print_details.h"

void fcomplex_malloc(fcomplex **p2p,size_t order);
void fcomplex_fftwf_malloc(fcomplex **p2p,size_t order);
void fcomplex_calloc(fcomplex **p2p,size_t order);
void fcomplex_malloc_ptr(fcomplex ***p2p2p,size_t order);
void fcomplex_fftwf_malloc_ptr(fcomplex ***p2p2p,size_t order);
void fcomplex_free(fcomplex **p2p);
void fcomplex_fftwf_free(fcomplex **p2p);
void fcomplex_free_ptr(fcomplex ***p2p2p);
void fcomplex_fftwf_free_ptr(fcomplex ***p2p2p);
