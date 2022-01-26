/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Output timing information

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "timing_info.h"

void timing_info(void){

   double temp,result;
   char string[100];

/* Overall runtime */
   if(parallel.myrank==0){ /* Restrict to master */
      reset_string(string);
      sprintf(string,"Overall runtime information [Averaged over the %d processors]",parallel.np);
      print_section_title(output.time,string);
   }
   timing.overall[1]=MPI_Wtime();
   temp=timing.overall[1]-timing.overall[0];
   MPI_Reduce(&temp,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      result/=parallel.dnp;   
      fprintf(output.time,"\n  Overall runtime: %.*gs [%d]\n",DBLP,result,1);
   }

/* Iterative kernel runtime */
   if(parallel.myrank==0){ /* Restrict to master */
      reset_string(string);
      sprintf(string,"Iterative kernel runtime information [Averaged over the %d processors]",parallel.np);
      print_section_title(output.time,string);
   }
   MPI_Reduce(&timing.iterative_solution,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      result/=parallel.dnp;
      fprintf(output.time,"\n  Iterative kernel average runtime: %.*gs [%d]\n",DBLP,
         (result/(double)timing.iterative_kernel_count),timing.iterative_kernel_count);
   }

/* DFT engine runtime */
   if(parallel.myrank==0){ /* Restrict to master */
      reset_string(string);
      sprintf(string,"DFT engine runtime information [Averaged over the %d processors]",parallel.np);
      print_section_title(output.time,string);
   }
   temp=timing.zero_padded_vector+timing.dft_tensor_components+timing.transpose_tensor_forward+
         timing.transpose_vector_forward+timing.zdft_zidft_and_fd_mult+timing.transpose_vector_reverse+
         timing.norm_diag_compress;

   MPI_Reduce(&temp,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      result/=parallel.dnp;
      fprintf(output.time,"\n  DFT engine average runtime: %.*gs [%d]\n",DBLP,
            (result/(double)timing.dft_engine_count),timing.dft_engine_count);
   }
   if(parallel.myrank==0){ /* Restrict to master */
      fprintf(output.time,"\n  ------ Breakdown [Averaged over the %d processors] ------\n\n",parallel.np);
   }
   MPI_Reduce(&timing.zero_padded_vector,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      result/=parallel.dnp;
      fprintf(output.time,"  Initialisation of the zero-padded vector    %.*gs [%d]\n\n",DBLP,
         (result/(double)timing.dft_engine_count),timing.dft_engine_count);
   }
   MPI_Reduce(&timing.dft_tensor_components,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      result/=parallel.dnp;
      fprintf(output.time,"  DFT of the 6 independent tensor components\n  [excluding communication]                   %.*gs [%d]\n",DBLP,
         (result/(double)timing.dft_tensor_components_count),timing.dft_tensor_components_count);
   }
   MPI_Reduce(&timing.transpose_tensor_forward,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      result/=parallel.dnp;
      fprintf(output.time,"  Distributed forward-transpose of the six\n  independent tensor components               %.*gs [%d]\n\n",DBLP,
         (result/(double)timing.dft_tensor_components_count),timing.dft_tensor_components_count);
   }
   MPI_Reduce(&timing.zdft_zidft_and_fd_mult,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      result/=parallel.dnp;
      fprintf(output.time,"  DFT of the 3 vector components &\n");
      fprintf(output.time,"  Fourier domain tensor-vector multiplication &\n");
      fprintf(output.time,"  iDFT of the 3 resulting vector components\n  [excluding communication]                   %.*gs [%d]\n",DBLP,
      (result/(double)timing.dft_engine_count),timing.dft_engine_count);
   }
   MPI_Reduce(&timing.transpose_vector_forward,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      result/=parallel.dnp;
      fprintf(output.time,"  Distributed forward-transpose of the three\n  zero-padded vector components               %.*gs [%d]\n",DBLP,
         (result/(double)timing.dft_engine_count),timing.dft_engine_count);
   }
   MPI_Reduce(&timing.transpose_vector_reverse,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      result/=parallel.dnp;
      fprintf(output.time,"  Distributed reverse-transpose of the three\n  zero-padded vector components               %.*gs [%d]\n",DBLP,
         (result/(double)timing.dft_engine_count),timing.dft_engine_count);
   }
   MPI_Reduce(&timing.norm_diag_compress,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   if(parallel.myrank==0){ /* Restrict to master */
      result/=parallel.dnp;
      fprintf(output.time,"  Normalisation of the resultant transform &\n");
      fprintf(output.time,"  addition of the diagonal contribution &\n");
      fprintf(output.time,"  data compression to occupied sites only     %.*gs [%d]\n",DBLP,
         (result/(double)timing.dft_engine_count),timing.dft_engine_count);
   }
   if(parallel.myrank==0){ /* Restrict to master */
      fclose(output.time);
   }
}
