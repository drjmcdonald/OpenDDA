/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   int memory allocation functions

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "int_alloc.h"

void int_malloc(int **p2p,size_t order){

   *p2p=(int*)malloc(sizeof(int)*order);
   if(*p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void int_calloc(int **p2p,size_t order){

   *p2p=(int*)calloc(order,sizeof(int));
   if(*p2p==NULL){
      print_error("Memory allocation error","calloc() failure while trying to allocate space",NULL);
   }
}

void int_malloc_ptr(int ***p2p2p,size_t order){

   *p2p2p=(int**)malloc(sizeof(int*)*order);
   if(*p2p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void int_free(int **p2p){

   if(*p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p);
      *p2p=NULL;
   }
}

void int_free_ptr(int ***p2p2p){

   if(*p2p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p2p);
      *p2p2p=NULL;
   }
}
