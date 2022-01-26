/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   unsigned int memory allocation functions

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "uint_alloc.h"

void uint_malloc(unsigned int **p2p,size_t order){

   *p2p=(unsigned int*)malloc(sizeof(unsigned int)*order);
   if(*p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void uint_calloc(unsigned int **p2p,size_t order){

   *p2p=(unsigned int*)calloc(order,sizeof(unsigned int));
   if(*p2p==NULL){
      print_error("Memory allocation error","calloc() failure while trying to allocate space",NULL);
   }
}

void uint_malloc_ptr(unsigned int ***p2p2p,size_t order){

   *p2p2p=(unsigned int**)malloc(sizeof(unsigned int*)*order);
   if(*p2p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void uint_free(unsigned int **p2p){

   if(*p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p);
      *p2p=NULL;
   }
}

void uint_free_ptr(unsigned int ***p2p2p){

   if(*p2p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p2p);
      *p2p2p=NULL;
   }
}
