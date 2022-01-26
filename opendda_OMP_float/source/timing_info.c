/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Output timing information

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "timing_info.h"

void timing_info(void){

/* Overall runtime */
   print_section_title(output.time,"Overall runtime information");
   timing.overall[1]=walltime();
   fprintf(output.time,"\n  Overall runtime: %.*gs [%d]\n",DBLP,timing.overall[1]-timing.overall[0],1);

/* Iterative kernel runtime */
   print_section_title(output.time,"Iterative kernel runtime information");
   fprintf(output.time,"\n  Iterative kernel average runtime: %.*gs [%d]\n",DBLP,
      (timing.iterative_solution/(double)timing.iterative_kernel_count),timing.iterative_kernel_count);

/* DFT engine runtime */
   print_section_title(output.time,"DFT engine runtime information");
   fprintf(output.time,"\n  DFT engine average runtime: %.*gs [%d]\n",DBLP,
      ((timing.zero_padded_vector+timing.dft_tensor_components+timing.zdft_zidft_and_fd_mult+
      timing.norm_and_diag+timing.compress)/(double)timing.dft_engine_count),timing.dft_engine_count);

   fprintf(output.time,"\n  ------ Breakdown ------\n");
   fprintf(output.time,"  Initialisation of the zero-padded vector    %.*gs [%d]\n",DBLP,
      (timing.zero_padded_vector/(double)timing.dft_engine_count),timing.dft_engine_count);
   fprintf(output.time,"  DFT of the 6 independent tensor components  %.*gs [%d]\n",DBLP,
      (timing.dft_tensor_components/(double)timing.dft_tensor_components_count),timing.dft_tensor_components_count);
   fprintf(output.time,"  DFT of the 3 vector components &\n");
   fprintf(output.time,"  Fourier domain tensor-vector multiplication &\n");
   fprintf(output.time,"  iDFT of the 3 resulting vector components   %.*gs [%d]\n",DBLP,
      (timing.zdft_zidft_and_fd_mult/(double)timing.dft_engine_count),timing.dft_engine_count);
   fprintf(output.time,"  Normalisation of the resultant transform &\n");
   fprintf(output.time,"  addition of the diagonal contribution       %.*gs [%d]\n",DBLP,
      (timing.norm_and_diag/(double)timing.dft_engine_count),timing.dft_engine_count);
   fprintf(output.time,"  Compress data to occupied sites only     %.*gs [%d]\n",DBLP,
      (timing.compress/(double)timing.dft_engine_count),timing.dft_engine_count);

   fclose(output.time);
}
