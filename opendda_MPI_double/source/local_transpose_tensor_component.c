/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Performs the local transpose of the tensor components of the
   interaction matrix required by the distributed algorithms 1 & 2

   comp: Tensor component to transpose
   local_tensor: Temporary scratch array for local transpose

   before: Jpp
   planesbefore: planesY_tensor
   i,j,k,j0,k0,
   initialsource,
   initialdest,
   source,dest,
   index0,index1: Loop and array index control variables
   transposed_flag_tensor: Bitfielded flag array to specify if a line
                    has been transposed or not 
   alloc_tensor: The 6 independent tensor components of the interaction
            matrix are stored in parallel so `alloc_tensor' is the
            memory allocated on each processor to store its portion of
            the interaction matrix before OR after the distributed
            xy-to-xz plane transpose

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "local_transpose_tensor_component.h"

void local_transpose_tensor_component(int before,int planesbefore,int comp){

   int i,j,k,j0,k0,initialsource,initialdest,array0;
   int source=-1,dest,index0,index1;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Transpose local data */
   array0=comp*parallel.alloc_tensor;
   if(before==planesbefore){ /* Square matrix transpose */
      for(k=1;k<planesbefore;k++){
         for(j=0;j<k;j++){
            source=k*before+j;
            dest=j*planesbefore+k;
            /* parallel.local_tensor=interaction_matrix[source];
               interaction_matrix[source]=interaction_matrix[dest];
               interaction_matrix[dest]=parallel.local_tensor; */
            index0=source*target.Ka+array0;
            index1=dest*target.Ka+array0;
            for(i=0;i<target.Ka;i++){
               parallel.local_tensor[i]=interaction_matrix[index0+i];
            }
            for(i=0;i<target.Ka;i++){
               interaction_matrix[index0+i]=interaction_matrix[index1+i];
            }
            for(i=0;i<target.Ka;i++){
               interaction_matrix[index1+i]=parallel.local_tensor[i];
            }
         }
      }
   }
   else{  /* Non-square matrix transpose */
      for(k=0;k<parallel.tf_ctrl_tensor;k++){
         parallel.transposed_flag_tensor[k]=0;
      }
      for(k=0;k<planesbefore;k++){
         for(j=0;j<before;j++){
            initialsource=k*before+j;
            initialdest=j*planesbefore+k;
            if(initialsource==initialdest){
               continue;
            }
            if((parallel.transposed_flag_tensor[(int)((double)initialsource/32.0)]&(1<<(initialsource%32)))==0){ /* If not already transposed */
               /* Move 'initialsource' to BUFFER */
               /* parallel.local_tensor=interaction_matrix[initialsource]; */
               index0=initialsource*target.Ka+array0;
               for(i=0;i<target.Ka;i++){
                  parallel.local_tensor[i]=interaction_matrix[index0+i];
               }
               dest=initialsource;
               while(source!=initialdest){
                  /* Treat source as a destination(j*planesbefore+k) and work out its source(k*before+j) */
                  j0=(int)((double)dest/(double)planesbefore);
                  k0=dest-(j0*planesbefore);
                  source=k0*before+j0;
                  /* Move 'source' to 'dest' */
                  /* interaction_matrix[dest]=interaction_matrix[source]; */
                  index0=source*target.Ka+array0;
                  index1=dest*target.Ka+array0;
                  for(i=0;i<target.Ka;i++){
                     interaction_matrix[index1+i]=interaction_matrix[index0+i];
                  }               
                  parallel.transposed_flag_tensor[(int)((double)source/32.0)]|=(1<<(source%32)); /* Mark 'source' as transposed */
                  dest=source;
               }
               /* Move BUFFER to 'dest' */
               /*   interaction_matrix[dest]=parallel.local_tensor; */
               index0=dest*target.Ka+array0;
               for(i=0;i<target.Ka;i++){
                  interaction_matrix[index0+i]=parallel.local_tensor[i];
               }
            }
         }
      }
   }
}
