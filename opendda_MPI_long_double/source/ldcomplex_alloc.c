/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   ldcomplex memory allocation functions

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "ldcomplex_alloc.h"

/* Note: **p2p is a pointer the original pointer i.e.

   p2p==&a   (Address where the original pointer is stored)
   *p2p==a   (Original pointer)
   **p2p==*a (Element that `*p2p==a' points to)

   Remember `int **p2p' is just declaring `p2p'
   as a pointer to a pointer to an element */

void ldcomplex_malloc(ldcomplex **p2p,size_t order){

   *p2p=(ldcomplex*)malloc(sizeof(ldcomplex)*order);
   if(*p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void ldcomplex_fftwl_malloc(ldcomplex **p2p,size_t order){

   *p2p=(ldcomplex*)fftwl_malloc(sizeof(ldcomplex)*order);
   if(*p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void ldcomplex_calloc(ldcomplex **p2p,size_t order){

   *p2p=(ldcomplex*)calloc(order,sizeof(ldcomplex));
   if(*p2p==NULL){
      print_error("Memory allocation error","calloc() failure while trying to allocate space",NULL);
   }
}

void ldcomplex_malloc_ptr(ldcomplex ***p2p2p,size_t order){

   *p2p2p=(ldcomplex**)malloc(sizeof(ldcomplex*)*order);
   if(*p2p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void ldcomplex_fftwl_malloc_ptr(ldcomplex ***p2p2p,size_t order){

   *p2p2p=(ldcomplex**)fftwl_malloc(sizeof(ldcomplex*)*order);
   if(*p2p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void ldcomplex_free(ldcomplex **p2p){

   if(*p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p);
      *p2p=NULL;
   }
}

void ldcomplex_fftwl_free(ldcomplex **p2p){

   if(*p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      fftwl_free(*p2p);
      *p2p=NULL;
   }
}

void ldcomplex_free_ptr(ldcomplex ***p2p2p){

   if(*p2p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p2p);
      *p2p2p=NULL;
   }
}

void ldcomplex_fftwl_free_ptr(ldcomplex ***p2p2p){

   if(*p2p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      fftwl_free(*p2p2p);
      *p2p2p=NULL;
   }
}
