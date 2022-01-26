/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   FILE memory allocation functions

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "file_alloc.h"

void file_malloc_ptr(FILE ***p2p2p,size_t order){

   *p2p2p=(FILE**)malloc(sizeof(FILE*)*order);
   if(*p2p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void file_free_ptr(FILE ***p2p2p){

   if(*p2p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p2p);
      *p2p2p=NULL;
   }
}
