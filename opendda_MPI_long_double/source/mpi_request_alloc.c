/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj
   
   MPI_Request memory allocation functions

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "mpi_request_alloc.h"

void mpi_request_malloc(MPI_Request **p2p2p,size_t order){

   *p2p2p=(MPI_Request*)malloc(sizeof(MPI_Request)*order);
   if(*p2p2p==NULL){
      print_error("Memory allocation error","malloc() failure while trying to allocate space",NULL);
   }
}

void mpi_request_calloc(MPI_Request **p2p2p,size_t order){

   *p2p2p=(MPI_Request*)calloc(order,sizeof(MPI_Request));
   if(*p2p2p==NULL){
      print_error("Memory allocation error","calloc() failure while trying to allocate space",NULL);
   }
}

void mpi_request_free(MPI_Request **p2p2p){

   if(*p2p2p==NULL){
      print_error("Memory deallocation error","Trying to free unallocated memory",NULL);
   }
   else{
      free(*p2p2p);
      *p2p2p=NULL;
   }
}
