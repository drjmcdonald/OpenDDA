/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   FFTW plan creation/destruction

   Depending on the flag passed:

   If flag==1
   Performs the x, y & z FFTW plan creation for the DFT and iDFT of the
   3 zero-padded vector components

   If flag==0
   Destroys all created FFTW plans

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dft_plan.h"

void dft_plan(int flag){

   if(flag){
/* Plan creation for the DFT of the vector components */
      plan_vcydft=fftwl_plan_many_dft(1,&target.Jpp,target.K,(fftwl_complex *)(&zero_padded_vector[0]),NULL,target.K,1,(fftwl_complex *)(&zero_padded_vector[0]),NULL,target.K,1,-1,FFTW_MEASURE);

/* Plan creation for the iDFT of the vector components */
      plan_vcydfti=fftwl_plan_many_dft(1,&target.Jpp,target.K,(fftwl_complex *)(&zero_padded_vector[0]),NULL,target.K,1,(fftwl_complex *)(&zero_padded_vector[0]),NULL,target.K,1,+1,FFTW_MEASURE);
   }
   else{
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Destory FFTW plans */
      fftwl_destroy_plan(plan_tcxdft);
      fftwl_destroy_plan(plan_tcydft);
      fftwl_destroy_plan(plan_tczdft);
      fftwl_destroy_plan(plan_vcxdft);
      fftwl_destroy_plan(plan_vcydft);
      fftwl_destroy_plan(plan_vczdft);
      fftwl_destroy_plan(plan_vcxdfti);
      fftwl_destroy_plan(plan_vcydfti);
      fftwl_destroy_plan(plan_vczdfti);
      fftwl_cleanup();
   }
}
