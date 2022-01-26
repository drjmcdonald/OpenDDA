/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   long double memory allocation functions

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "long_double_alloc.h"

/* Note: **p2p is a pointer the original pointer i.e.

   p2p==&a     (Address where the original pointer is stored)
   *p2p==a     (Original pointer)
   **p2p==*a   (Element that `*p2p==a' points to)

   Remember `int **p2p' is just declaring `p2p'
   as a pointer to a pointer to an element */

void long_double_malloc(long double **p2p,size_t order){

   *p2p=(long double*)malloc(sizeof(long double)*order);
   if(*p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void long_double_calloc(long double **p2p,size_t order){

   *p2p=(long double*)calloc(order,sizeof(long double));
   if(*p2p==NULL){
      print_error("Memory allocation error","calloc() failure while trying to allocate space",NULL);
   }
}

void long_double_malloc_ptr(long double ***p2p2p,size_t order){

   *p2p2p=(long double**)malloc(sizeof(long double*)*order);
   if(*p2p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void long_double_free(long double **p2p){

   if(*p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p);
      *p2p=NULL;
   }
}

void long_double_free_ptr(long double ***p2p2p){

   if(*p2p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p2p);
      *p2p2p=NULL;
   }
}
