/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Transpose the 3 zero-padded vector components from xy-planes
   to xz-planes OR vice versa (Distributed transpose)

   Transposes the 3 zero-padded vector components separately

   transpose_vector_component_0: Transposes zero_padded_vector[3][K*Jpp*P]
                                 line-by-line (lines in the x-direction)
   transpose_vector_component_1: Transposes zero_padded_vector[3][K*Jpp*P]
                                 plane-by-plane (includes a local transpose
                                 and only works if message size is below the
                                 MPI EAGER limit)
   transpose_vector_component_2: Transposes zero_padded_vector[3][K*Jpp*P]
                                 plane-by-plane (includes a local transpose) 

   Note: Due to the EAGER limit, transposes

         transpose_vector_component_0()
         transpose_vector_component_1()

   are only compliant with MPICH1 and NOT MPICH2 and OPENMPI in their default
   configurations. Transpose

         transpose_vector_component_2()

   is fully compliant and are the recommended choice

   before: Forward transpose - Jpp
           Reverse transpose - P
   after: Forward transpose - P
          Reverse transpose - Jpp
   limitbefore: Forward transpose - limitY_padded
                Reverse transpose - limitZ_padded
   limitafter: Forward transpose - limitZ_padded
               Reverse transpose - limitY_padded
   recv: Forward transpose - Znp_padded
         Reverse transpose - Ynp
   comp: Vector component to be transposed.

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "transpose_vector_components.h"

void transpose_vector_component_0(int before,int after,int limitbefore,int limitafter,int comp){

   int i,j,k,j1,k1,cr,limit;
   MPI_Request *rqrecv,rqsend;

/* Allocate memory for the MPI_Request for Irecv.
   The max number of these is the max number of receives that a
   processor will perform which is of course equal to the number
   of valid memory locations on that proc after transposition
   which is `limitafter' */
   mpi_request_malloc(&rqrecv,(size_t)limitafter);

   cr=0;
   limit=limitbefore>limitafter?limitbefore:limitafter;
   for(i=0;i<limit;i++){
      if(i==0&parallel.myrank==0){continue;} /* Skip proc 0 element 0 as it does not move */
      if(i<limitbefore){
/* Treat (proc,element)=(myrank,i) as a source and work out the destination */
         j1=i%before; /* j1=[Send_element]%before */
         k1=(int)(((double)(i-j1)/(double)before))*parallel.np+parallel.myrank; /* k1=(([Send_element]-j1)/before)*np+[Sending_processor_rank] */
         MPI_Isend(&zero_padded_vector[comp*parallel.alloc_padded+(i*target.K)],target.K,fcomplex_type,j1%parallel.np,k1*before+j1,MPI_COMM_WORLD,&rqsend); /* Non-blocking send */
         MPI_Wait(&rqsend,MPI_STATUS_IGNORE);
      }
      if(i<limitafter){
/* Treat (proc,element)=(myrank,i) as a destination and work out the source */
         k=i%after; /* k=[Receive_element]%after */
         j=(int)(((double)(i-k)/(double)after))*parallel.np+parallel.myrank; /* j=(([Receive_element]-k)/after)*np+[Receiving_processor_rank] */
         MPI_Irecv(&zero_padded_vector[comp*parallel.alloc_padded+(i*target.K)],target.K,fcomplex_type,k%parallel.np,k*before+j,MPI_COMM_WORLD,&rqrecv[cr]); /* Non-blocking receive */
         cr++;
      }
   }
   MPI_Waitall(cr,rqrecv,MPI_STATUS_IGNORE); /* Wait for the non-blocking receives to complete */
   
 /* Free allocated memory */
   mpi_request_free(&rqrecv);
}

void transpose_vector_component_1(int before,int after,int planesbefore,int planesafter,int recv,int comp){

   int i,j,limit,cs=0,cr=0,index0,index1,index2;
   MPI_Request *rqrecv,*rqsend;
   MPI_Datatype vector_recv0,vector_recv1;
   int count_recv0=recv+1,count_recv1=recv,blocklen_recv=target.K,stride_recv=parallel.np*target.K;

/* Create and commit the MPI_Type_vector for receiving data */
/* Forward: If receiving from proc with `myrank' < Pp%np, each plane has `Znp+1' x-dirn lines
   Reverse: If receiving from proc with `myrank' < Jp%np, each plane has `Ynp+1' x-dirn lines */
   MPI_Type_vector(count_recv0,blocklen_recv,stride_recv,fcomplex_type,&vector_recv0);
/* Forward: If receiving from proc with `myrank' >= Pp%np, each plane has `Znp' x-dirn lines
   Reverse: If receiving from proc with `myrank' >= Jp%np, each plane has `Ynp' x-dirn lines */
   MPI_Type_vector(count_recv1,blocklen_recv,stride_recv,fcomplex_type,&vector_recv1);
   MPI_Type_commit(&vector_recv0); /* Commit the vector */
   MPI_Type_commit(&vector_recv1); /* Commit the vector */

/* Allocate memory for MPI_Irecv requests */
   mpi_request_malloc(&rqrecv,(size_t)(planesafter*parallel.np)); /* There will be `planesafter' planes on the local
   proc after the transpose and there is one message sent by each proc (including the local proc)
   per plane so there will be the local proc will have to make `planesafter*np' receives i.e.,
   each of the received `planesafter' planes on each local proc is made up of a single
   contribution from all the `np' processor. This means that each proc will perform
   `np*planesafter' receives */
   mpi_request_malloc(&rqsend,(size_t)before); /* Max is the number of planes after the local transpose
   which is equal the the number of planes on the local proc after the transpose which is in turn
   equal to the length of the planes before the transpose i.e., after the local transpose, each
   processor has `before' planes. Thus during the communication phase each processor will perform
   `before' sends. */

/* Performs a local transpose on its portion of the zero-padded vector component specified by `comp'. */
   local_transpose_vector_component(before,planesbefore,comp);

/* Send the locally transposed data to the appropriate processors */
   limit=after%parallel.np; /* Used to control how many elements to reveive depending on the sending processor's rank */
   index0=0;
   index1=planesbefore*target.K; /* The number of elements in each plane after local transpose and before communication */
   index2=parallel.myrank*before; /* After the local transpose BUT before the communication, each local proc
   has `before' planes of size `index' so `index2' is used in MPI tag generation */
   i=0;
   while((i<before)||(index0<planesafter)){
      cs=0;
      while((i*planesbefore)<((index0+1)*after)&&(i<before)){
 /* Send the planes that are impeding the in-place receive.
         Note that this algorithm only works when the message size i.e. index1*sizeof(fcomplex) < MPI_EAGER_LIMIT.
         The reason for this is that when the message is below this limit, MPI sends the entire message immediately,
         allowing for the in-place receive. When the message size is bigger than this limit, MPI uses the rendezvous protocal
         and the transpose requires expficit buffering to implement in-place functionality */
         MPI_Isend(&zero_padded_vector[comp*parallel.alloc_padded+(i*index1)],index1,fcomplex_type,i%parallel.np,index2+i,MPI_COMM_WORLD,&rqsend[cs]);
         cs++;i++;
      }
      MPI_Waitall(cs,rqsend,MPI_STATUS_IGNORE); /* Wait for the sends to finish so that the receives can be made */
      if(index0<planesafter){ /* Perform the receives */
         for(j=0;j<parallel.np;j++){
            if(j<limit){ /* Forward: If receiving from proc with `myrank' < Pp%np, each plane has `Znp+1' x-dirn lines
                           Reverse: If receiving from proc with `myrank' < Jp%np, each plane has `Ynp+1' x-dirn lines */
               MPI_Irecv(&zero_padded_vector[comp*parallel.alloc_padded+((index0*after+j)*target.K)],1,vector_recv0,j,j*before+index0*parallel.np+parallel.myrank,MPI_COMM_WORLD,&rqrecv[cr]);
               cr++; /* Increment the number of non-blocking receives that are currently active */
            }
            else{ /* Forward: If receiving from proc with `myrank' >= Pp%np, each plane has `Znp' x-dirn lines
                     Reverse: If receiving from proc with `myrank' >= Jp%np, each plane has `Ynp' x-dirn lines */
               MPI_Irecv(&zero_padded_vector[comp*parallel.alloc_padded+((index0*after+j)*target.K)],1,vector_recv1,j,j*before+index0*parallel.np+parallel.myrank,MPI_COMM_WORLD,&rqrecv[cr]);
               cr++; /* Increment the number of non-blocking receives that are currently active */
            }
         }
      }
      index0++; /* Go to the next receive plane */
   }
   MPI_Waitall(cr,rqrecv,MPI_STATUS_IGNORE); /* Wait for all receives to finish */

/* Free allocated memory */
   mpi_request_free(&rqsend);
   mpi_request_free(&rqrecv);

   MPI_Type_free(&vector_recv0);
   MPI_Type_free(&vector_recv1);
}

void transpose_vector_component_2(int before,int after,int planesbefore,int planesafter,int recv,int comp){

   int i,i0,j,k,limit,cr=0,cs,index0,index1,bufsize,buffind;
   fcomplex *buffer;
   MPI_Request *rqrecv,*rqsend;
   MPI_Datatype vector_recv0,vector_recv1;
   int count_recv0=recv+1,count_recv1=recv,blocklen_recv=target.K,stride_recv=parallel.np*target.K;

/* Create and commit the MPI_Type_vector for receiving data */
/* Forward: If receiving from proc with `myrank' < Pp%np, each plane has `Znp+1' x-dirn lines
   Reverse: If receiving from proc with `myrank' < Jp%np, each plane has `Ynp+1' x-dirn lines */
   MPI_Type_vector(count_recv0,blocklen_recv,stride_recv,fcomplex_type,&vector_recv0);
/* Forward: If receiving from proc with `myrank' >= Pp%np, each plane has `Znp' x-dirn lines
   Reverse: If receiving from proc with `myrank' >= Jp%np, each plane has `Ynp' x-dirn lines */
   MPI_Type_vector(count_recv1,blocklen_recv,stride_recv,fcomplex_type,&vector_recv1);
   MPI_Type_commit(&vector_recv0); /* Commit the vector */
   MPI_Type_commit(&vector_recv1); /* Commit the vector */

/* Allocate memory for MPI_Irecv requests */
   mpi_request_malloc(&rqrecv,(size_t)(parallel.np*planesafter)); /* There will be `planesafter' planes
   on the local proc after the transpose and there is one message sent by each proc
   (including the local proc) per plane so there will be the local proc will have to
   make `planesafter*np' receives i.e., each of the received `planesafter' planes on
   each local proc is made up of a single contribution from all the `np' processor.
   This means that each proc will perform `np*planesafter' receives */
   mpi_request_malloc(&rqsend,(size_t)parallel.np); /* For each communication round, each processor must send one
   plane of length `before' to all processors (including itself) that has a relevant plane for
   `index0=0..planesafter-1', i.e. upto `np' sends per communication round */

/* Performs a local transpose on its portion of the zero-padded vector component specified by `comp'. */
   local_transpose_vector_component(before,planesbefore,comp);

/* Find the buffer memory required
   For each receiving plane index0=0 to (planesafter-1), there are `np' send
   and receives. In order to facilitate the in-place receives, the way must
   be cleared for the receive by buffering the data that is occupying the
   data elements for the receive. Since the number of buffered planes may be
   more than the `np' sent planes in each round, the buffering requirement maybe
   slightly higher than `np' (a slight accumulation in the number of planes buffered)

   index0=0 to (planesafter-1): controls the receiving plane. After the local transpose,
   each processor has a block `before' high and `planesbefore' wide. */
   index0=0;i=0;i0=0;buffind=0;bufsize=0; /* i controls the buffering, i0 controls the
   sending and index0 controls the receiving plane */
   while(i<before||i0<before||index0<planesafter){
      while((i<((index0+1)*parallel.np)||(i*planesbefore)<((index0+1)*after))&&i<before){
         /* i<before:
               Needed because each processor has `before' planes, thus `before' is the
               maximum buffering that could ever be required
            (i*planesbefore)<((index0+1)*after):
               This is what calculates how many planes need to be buffered. After the local
               transpose, the planes on each local proc are of length `planesbefore'. If this
               is smaller than `after' then more than one plane of length `planesbefore' will
               need to be buffered in order to allow the receive to be carried out in-place
            i<((index0+1)*np):
               To facilitate the in-place transpose, the planes impeding the receives are
               buffered and thus the sending is performed out of the buffering memory. Since
               `np' planes are sent in each communication round, `index0'=0 to (planesafter-1),
               we need to make sure that the `np' planes have been buffered and are ready for sending.
               Thus the code ensures that even if it is not required to buffer the `np' planes to
               perform the receive that the `np' planes are buffered for sending */
         buffind++;
         i++;
      }
      for(j=0;j<parallel.np;j++){ /* During each communication round, the algorithm and thus each proc performs
               `np' sends, which frees up buffer space, thus the `number of planes that need to be
               buffered' can be decremented */
         if(i0<before){
            buffind--;
            i0++;
         }
      }
      index0++;
      bufsize=buffind>bufsize?buffind:bufsize;
   }
/* Buffering requirement is the bufsize to store the accumulated planes and
   `np' planes for the dynamic send and receive contribution */
   bufsize+=parallel.np;
   fcomplex_malloc(&buffer,(size_t)(bufsize*planesbefore*target.K)); /* Allocate the buffer */

/* Send the locally transposed data to the appropriate processors */
   limit=after%parallel.np; /* Used to control how many elements to reveive depending on the sending processor's rank */
   index0=0;
   index1=planesbefore*target.K; /* The number of elements in each plane after local transpose and before communication */
   i=0;i0=0; /* i controls the buffering, i0 controls the sending and index0 controls the receiving plane */
   while(i<before||i0<before||index0<planesafter){
      while((i<((index0+1)*parallel.np)||(i*planesbefore)<((index0+1)*after))&&i<before){
         for(k=0;k<index1;k++){ /* Buffer the plane */
            buffer[(i%bufsize)*index1+k]=zero_padded_vector[comp*parallel.alloc_padded+(i*index1+k)];
         }
         i++; /* Increment the plane number: Buffering */
      }
      if(index0<planesafter){ /* Perform the `np' receives */
         for(j=0;j<parallel.np;j++){
            if(j<limit){ /* Forward: If receiving from proc with `myrank' < Pp%np, each plane has `Znp+1' x-dirn lines
                           Reverse: If receiving from proc with `myrank' < Jp%np, each plane has `Ynp+1' x-dirn lines */
               MPI_Irecv(&zero_padded_vector[comp*parallel.alloc_padded+((index0*after+j)*target.K)],1,vector_recv0,j,j*before+index0*parallel.np+parallel.myrank,MPI_COMM_WORLD,&rqrecv[cr]);
               cr++; /* Increment the number of non-blocking receives that are currently active */
            }
            else{ /* Forward: If receiving from proc with `myrank' >= Pp%np, each plane has `Znp' x-dirn lines
                     Reverse: If receiving from proc with `myrank' >= Jp%np, each plane has `Ynp' x-dirn lines */
               MPI_Irecv(&zero_padded_vector[comp*parallel.alloc_padded+((index0*after+j)*target.K)],1,vector_recv1,j,j*before+index0*parallel.np+parallel.myrank,MPI_COMM_WORLD,&rqrecv[cr]);
               cr++; /* Increment the number of non-blocking receives that are currently active */
            }
         }
      }
		cs=0; /* Set the number of active sends to zero */
      for(j=0;j<parallel.np;j++){ /* Perform the sends */
         if(i0<before){
            MPI_Isend(&buffer[(i0%bufsize)*index1],index1,fcomplex_type,i0%parallel.np,parallel.myrank*before+i0,MPI_COMM_WORLD,&rqsend[cs]);
            i0++; /* Increment the plane number: Sending */
				cs++; /* Increment the number of non-blocking sends that are currently active. */
         }
      }
      MPI_Waitall(cs,rqsend,MPI_STATUS_IGNORE); /* Wait for cs sends to complete, cannot reuse buffering location with modulo fn until send has completed. Note: due to the if statement, if(i0<before), this may be < np */
      index0++; /* Go to the next receive plane */
   }
   MPI_Waitall(cr,rqrecv,MPI_STATUS_IGNORE); /* Wait for receives to complete */

/* Free allocated memory */
   mpi_request_free(&rqrecv);
   mpi_request_free(&rqsend);
   fcomplex_free(&buffer);

/* Free the vector datatypes */
   MPI_Type_free(&vector_recv0);
   MPI_Type_free(&vector_recv1);
}
