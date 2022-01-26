/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Initialise and finalise MPI

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "mpi_environment.h"

void mpi_initialise(int *ac,char ***av){

   MPI_Init(ac,av);
   MPI_Comm_rank(MPI_COMM_WORLD,&parallel.myrank);
   MPI_Comm_size(MPI_COMM_WORLD,&parallel.np);
   parallel.dnp=(double)parallel.np; /* Double cast of the number of processors */

/* Create an MPI data type for ldcomplex */
   MPI_Type_contiguous(2,MPI_LONG_DOUBLE,&ldcomplex_type);
   MPI_Type_commit(&ldcomplex_type); /* Commit the datatype */

/* Create a user-defined MPI function to add complex numbers */
   MPI_Op_create((MPI_User_function*)ldcomplex_sum,1,&add_ldcomplex);
}

void mpi_finalise(void){

/* Free the MPI data type for ldcomplex */
   MPI_Type_free(&ldcomplex_type);

/* Free the user-defined MPI function to add complex numbers */
   MPI_Op_free(&add_ldcomplex);

   MPI_Finalize();
}

int ldcomplex_sum(ldcomplex *in,ldcomplex *inout,int *len,MPI_Datatype *dtptr){
/* Function to add ldcomplex for MPI_Allreduce */

   int i;
   ldcomplex temp;

   for(i=0;i<(*len);i++){
      temp.dat[0]=(*in).dat[0]+(*inout).dat[0];
      temp.dat[1]=(*in).dat[1]+(*inout).dat[1];
      *inout=temp;
      in++;inout++;
   }

   return 0;
}
