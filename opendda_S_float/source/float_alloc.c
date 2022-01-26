/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   float memory allocation functions

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "float_alloc.h"

/* Note: **p2p is a pointer the original pointer i.e.

   p2p==&a     (Address where the original pointer is stored)
   *p2p==a     (Original pointer)
   **p2p==*a   (Element that `*p2p==a' points to)

   Remember `int **p2p' is just declaring `p2p'
   as a pointer to a pointer to an element */

void float_malloc(float **p2p,size_t order){

   *p2p=(float*)malloc(sizeof(float)*order);
   if(*p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void float_calloc(float **p2p,size_t order){

   *p2p=(float*)calloc(order,sizeof(float));
   if(*p2p==NULL){
      print_error("Memory allocation error","calloc() failure while trying to allocate space",NULL);
   }
}

void float_malloc_ptr(float ***p2p2p,size_t order){

   *p2p2p=(float**)malloc(sizeof(float*)*order);
   if(*p2p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void float_free(float **p2p){

   if(*p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p);
      *p2p=NULL;
   }
}

void float_free_ptr(float ***p2p2p){

   if(*p2p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p2p);
      *p2p2p=NULL;
   }
}
