/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Calculates the memory allocation required for each node for the 6 independent
   tensor components of the interaction matrix, the zero-padded vector components
   and for all other vectors

   K: Number of sites in each line
   J: Number of lines in a plane
   P: Number of planes in array
   Kp: Kp=2K-1
   Jp: Jp=2J-1
   Pp: Pp=2P-1
   np,dnp: `np' is the number of processors and `dnp' is it cast to a double

   alloc_tensor: The 6 independent tensor components of the interaction
                 matrix are stored in parallel so `alloc_tensor' is the
                 memory allocated on each processor to store its portion of
                 the interaction matrix before OR after the distributed
                 xy-to-xz plane transpose 
   alloc_padded: The 3 components of the zero-padded vector components
                 are stored in parallel so `alloc_alloc_padded' is the
                 memory allocated on each processor to store its
                 portion of the zero-padded verctor before OR after the
                 distributed xy-to-xz plane transpose
   alloc_vector: The 3 components of all the basic vector components
                 are stored in parallel so `alloc_vector' is the
                 memory allocated on each processor to store its
                 portion of the vector before OR after the distributed
                 xy-to-xz plane transpose
   alloc_vector3: alloc_vector*3
   myrank: MPI processor rank 0(master) to np-1

   BEFORE XY-TO-XZ DISTRIBUTED TRANSPOSE
   Znp_tensor: Base number of xy-planes per proc before the transpose.
               For the interaction matrix (gets transposed)
   Znp_padded: Base number of xy-planes per proc before the transpose.
               For the zero padded vector (gets transposed)
   Znp_vector: Base number of xy-planes per proc before the transpose.
               For all other vectors (does not get transposed)
   AFTER XY-TO-XZ DISTRIBUTED TRANSPOSE
   Ynp: Base number of xz-planes per proc after the transpose.
        For the interaction matrix & the zero padded vector
   plane: Since the system is `plane' decomposed across the
          available processors, plane k, k=0..P, is no longer
          plane k but is plane ((int)((double)k/dnp)) on
          processor k mod np.
   planesY_padded: Number of xy-planes locally on the processor before the transpose
                   i.e., Znp_padded+1 if myrank < P%np, Znp_padded if myrank >= P%np
   limitY_padded: Number of lines in the x-direction locally before the transpose
                  i.e., planesY_padded*Jpp
   planesY_tensor: Number of xy-planes locally on the processor before the transpose
                   i.e., Znp_tensor+1 if myrank < Pa%np, Znp_tensor if myrank >= Pa%np
   limitY_tensor: Number of lines in the x-direction locally before the
                  transpose i.e., planesY_tensor*Jpp, N.B. For this MPI version the
                  6 independent tensor components of the interaction matrix are stored
                  in interaction_matrix[6][Ka*Jpp*Pa] rather than the
                  interaction_matrix[6][Ka*Ja*Pa] stored for the serial and OpenMP cases.
                  This is to facilitate the distributed FD tensor-vector multiplication.
   planesZ_padded: Number of xz-planes locally on the processor after the
                   transpose i.e., Ynp+1 if myrank < Jpp%np, Ynp if myrank >= Jpp%np
   limitZ_padded: Number of lines in the x-direction locally after the transpose i.e.,
                  planesZ_padded*P for zero_padded_vector
   planesZ_tensor: Number of xz-planes locally on the processor after the transpose i.e.,
                   Ynp+1 if myrank < Jpp%np, Ynp if myrank >= Jpp%np
   limitZ_tensor: Number of lines in the x-direction locally before the transpose i.e.,
                  planesY_tensor*Pa, N.B. For this MPI version the 6 independent tensor
                  components of the interaction matrix are stored in
                  interaction_matrix[6][Ka*Jpp*Pa] rather than the
                  interaction_matrix[6][Ka*Ja*Pa] stored for the serial and OpenMP cases.
                  This is to facilitate the distributed FD tensor-vector multiplication.

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "memory_allocation_sizes.h"

void memory_allocation_sizes(void){

   parallel.Ynp=(int)((double)target.Jpp/parallel.dnp);
   parallel.Znp_padded=(int)((double)target.P/parallel.dnp);
   parallel.Znp_tensor=(int)((double)target.Pa/parallel.dnp);
   parallel.Znp_vector=(int)((double)target.P/parallel.dnp);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* interaction_matrix:
   Before transpose if myrank < Pa%np=Pa%np then processor has Znp_tensor+1 Ka*Jpp planes
                              >= Pa%np=Pa%np then processor has Znp_tensor Ka*Jpp planes
   After transpose if myrank < Jpp%np then processor has Ynp+1 Qa=Ka*Pa planes
                             >= Jpp%np then processor has Ynp  Qa=Ka*Pa planes */
   if(parallel.myrank<target.Pa%parallel.np){ /* Before transpose: myrank < Pa%np */
      parallel.planesY_tensor=parallel.Znp_tensor+1;
      parallel.limitY_tensor=parallel.planesY_tensor*target.Jpp;
   }
   else{ /* Before transpose: myrank >= Pa%np */
      parallel.planesY_tensor=parallel.Znp_tensor;
      parallel.limitY_tensor=parallel.planesY_tensor*target.Jpp;
   }
   if(parallel.myrank<target.Jpp%parallel.np){ /* After transpose: myrank < Jpp%np */
      parallel.planesZ_tensor=parallel.Ynp+1;
      parallel.limitZ_tensor=parallel.planesZ_tensor*target.Pa;
   }
   else{ /* After transpose: myrank >= Jpp%np */
      parallel.planesZ_tensor=parallel.Ynp;
      parallel.limitZ_tensor=parallel.planesZ_tensor*target.Pa;
   }
   parallel.alloc_tensor=parallel.limitZ_tensor>parallel.limitY_tensor?
      parallel.limitZ_tensor*target.Ka:parallel.limitY_tensor*target.Ka;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* zero_padded_vector:
   Before transpose if myrank < P%np then processor has Znp_padded+1 K*Jpp=Mv planes
                              >= P%np then processor has Znp_padded K*Jpp=Mv planes
   After transpose if myrank < Jpp%np then processor has Ynp+1 K*P=M planes
                             >= Jpp%np then processor has Ynp K*P=M planes */
   if(parallel.myrank<target.P%parallel.np){ /* Before transpose: myrank < P%np */
      parallel.planesY_padded=parallel.Znp_padded+1;
      parallel.limitY_padded=parallel.planesY_padded*target.Jpp;
   }
   else{ /* Before transpose: myrank >= P%np */
      parallel.planesY_padded=parallel.Znp_padded;
      parallel.limitY_padded=parallel.Znp_padded*target.Jpp;
   }
   if(parallel.myrank<target.Jpp%parallel.np){ /* After transpose: myrank < Jpp%np*/
      parallel.planesZ_padded=parallel.Ynp+1;
      parallel.limitZ_padded=parallel.planesZ_padded*target.P;
   }
   else{ /* After transpose: myrank >= Jpp%np */
      parallel.planesZ_padded=parallel.Ynp;
      parallel.limitZ_padded=parallel.planesZ_padded*target.P;
   }
   parallel.alloc_padded=parallel.limitZ_padded>parallel.limitY_padded?
      parallel.limitZ_padded*target.K:parallel.limitY_padded*target.K;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* All other vectors (not transposed)
   If myrank < P%np then processor has Znp_vector+1 xy-planes planes
             >= P%np then processor has Znp_vector xy-planes planes */
   if(parallel.myrank<target.P%parallel.np){ /* myrank < P%np */
      parallel.planesY_vector=parallel.Znp_vector+1;
   }
   else{ /* myrank >= P%np */
      parallel.planesY_vector=parallel.Znp_vector;
   }
}
