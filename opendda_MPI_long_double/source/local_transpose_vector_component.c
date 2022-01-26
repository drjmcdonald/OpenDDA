/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Performs the local transpose of the zero-padded vector components
   required by the distributed algorithms 1 & 2

   comp: Vector component to transpose
   local_padded: Temporary scratch array for local transpose

   before: Forward transpose - Jp
           Reverse transpose - Pp
   planesbefore: Forward transpose - planesY_padded
                 Reverse transpose - planesZ_padded
   i,j,k,j0,k0,
   initialsource,
   initialdest,
   source,dest,
   index0,index1: Loop and array index control variables
   transposed_flag_padded: Bitfielded flag array to specify if a line
                    has been transposed or not 
   alloc_padded: The 3 components of the zero-padded vector components
                 are stored in parallel so `alloc_alloc_padded' is the
                 memory allocated on each processor to store its
                 portion of the zero-padded verctor before OR after the
                 distributed xy-to-xz plane transpose

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "local_transpose_vector_component.h"

void local_transpose_vector_component(int before,int planesbefore,int comp){

   int i,j,k,j0,k0,initialsource,initialdest,array0;
   int source=-1,dest,index0,index1;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Transpose local data */
   array0=comp*parallel.alloc_padded;
   if(before==planesbefore){ /* Square matrix transpose */
      for(k=1;k<planesbefore;k++){
         for(j=0;j<k;j++){
            source=k*before+j;
            dest=j*planesbefore+k;
/*            parallel.local_padded=zero_padded_vector[source];
            zero_padded_vector[source]=zero_padded_vector[dest];
            zero_padded_vector[dest]=parallel.local_padded; */
            index0=source*target.K+array0;
            index1=dest*target.K+array0;
            for(i=0;i<target.K;i++){
               parallel.local_padded[i]=zero_padded_vector[index0+i];
            }
            for(i=0;i<target.K;i++){
               zero_padded_vector[index0+i]=zero_padded_vector[index1+i];
            }
            for(i=0;i<target.K;i++){
               zero_padded_vector[index1+i]=parallel.local_padded[i];
            }
         }
      }
   }
   else{  /* Non-square matrix transpose */
      for(k=0;k<parallel.tf_ctrl_padded;k++){
         parallel.transposed_flag_padded[k]=0;
      }
      for(k=0;k<planesbefore;k++){
         for(j=0;j<before;j++){
            initialsource=k*before+j;
            initialdest=j*planesbefore+k;
            if(initialsource==initialdest){
               continue;
            }
            if((parallel.transposed_flag_padded[(int)((double)initialsource/32.0)]&(1<<(initialsource%32)))==0){ /* If not already transposed */
               /* Move 'initialsource' to BUFFER */
/*               parallel.local_padded=zero_padded_vector[initialsource]; */
               index0=initialsource*target.K+array0;
               for(i=0;i<target.K;i++){
                  parallel.local_padded[i]=zero_padded_vector[index0+i];
               }
               dest=initialsource;
               while(source!=initialdest){
                  /* Treat source as a destination(j*planesbefore+k) and work out its source(k*before+j) */
                  j0=(int)((double)dest/(double)planesbefore);
                  k0=dest-(j0*planesbefore);
                  source=k0*before+j0;
                  /* Move 'source' to 'dest' */
/*                  zero_padded_vector[dest]=zero_padded_vector[source]; */
                  index0=source*target.K+array0;
                  index1=dest*target.K+array0;
                  for(i=0;i<target.K;i++){
                     zero_padded_vector[index1+i]=zero_padded_vector[index0+i];
                  }               
                  parallel.transposed_flag_padded[(int)((double)source/32.0)]|=(1<<(source%32)); /* Mark 'source' as transposed */
                  dest=source;
               }
               /* Move BUFFER to 'dest' */
/*               zero_padded_vector[dest]=parallel.local_padded; */
               index0=dest*target.K+array0;
               for(i=0;i<target.K;i++){
                  zero_padded_vector[index0+i]=parallel.local_padded[i];
               }
            }
         }
      }
   }
}
