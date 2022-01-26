/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   dcomplex memory allocation functions

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dcomplex_alloc.h"

/* Note: **p2p is a pointer the original pointer i.e.

   p2p==&a   (Address where the original pointer is stored)
   *p2p==a   (Original pointer)
   **p2p==*a (Element that `*p2p==a' points to)

   Remember `int **p2p' is just declaring `p2p'
   as a pointer to a pointer to an element */

void dcomplex_malloc(dcomplex **p2p,size_t order){

   *p2p=(dcomplex*)malloc(sizeof(dcomplex)*order);
   if(*p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void dcomplex_fftw_malloc(dcomplex **p2p,size_t order){

   *p2p=(dcomplex*)fftw_malloc(sizeof(dcomplex)*order);
   if(*p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void dcomplex_calloc(dcomplex **p2p,size_t order){

   *p2p=(dcomplex*)calloc(order,sizeof(dcomplex));
   if(*p2p==NULL){
      print_error("Memory allocation error","calloc() failure while trying to allocate space",NULL);
   }
}

void dcomplex_malloc_ptr(dcomplex ***p2p2p,size_t order){

   *p2p2p=(dcomplex**)malloc(sizeof(dcomplex*)*order);
   if(*p2p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void dcomplex_fftw_malloc_ptr(dcomplex ***p2p2p,size_t order){

   *p2p2p=(dcomplex**)fftw_malloc(sizeof(dcomplex*)*order);
   if(*p2p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void dcomplex_free(dcomplex **p2p){

   if(*p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p);
      *p2p=NULL;
   }
}

void dcomplex_fftw_free(dcomplex **p2p){

   if(*p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      fftw_free(*p2p);
      *p2p=NULL;
   }
}

void dcomplex_free_ptr(dcomplex ***p2p2p){

   if(*p2p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p2p);
      *p2p2p=NULL;
   }
}

void dcomplex_fftw_free_ptr(dcomplex ***p2p2p){

   if(*p2p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      fftw_free(*p2p2p);
      *p2p2p=NULL;
   }
}
