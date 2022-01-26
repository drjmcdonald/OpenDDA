/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   FFTW plan creation/destruction

   Depending on the flag passed:

   If flag==1
   Performs the x, y & z FFTW plan creation for the DFT and iDFT of the
   3 zero-padded vector components and the DFT of the 6 independent
   tensor components

   If flag==0
   Destroys all created FFTW plans

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dft_plan.h"

void dft_plan(int flag){

   if(flag){
/* Plan creation for the DFT of the tensor components */
      plan_tcxdft=fftw_plan_dft(1,&target.Kpp,(fftw_complex *)(&xdft[0]),(fftw_complex *)(&xdft[0]),-1,FFTW_MEASURE);
      plan_tcydft=fftw_plan_dft(1,&target.Jpp,(fftw_complex *)(&ydft[0]),(fftw_complex *)(&ydft[0]),-1,FFTW_MEASURE);
      plan_tczdft=fftw_plan_dft(1,&target.Ppp,(fftw_complex *)(&zdft[0]),(fftw_complex *)(&zdft[0]),-1,FFTW_MEASURE);

/* Plan creation for the DFT of the vector components */
      plan_vcxdft=fftw_plan_many_dft(1,&target.Kpp,target.P,(fftw_complex *)(&xzscratch[0]),NULL,1,target.Kpp,(fftw_complex *)(&xzscratch[0]),NULL,1,target.Kpp,-1,FFTW_MEASURE);
      plan_vcydft=fftw_plan_many_dft(1,&target.Jpp,target.K,(fftw_complex *)(&zero_padded_vector[0]),NULL,target.K,1,(fftw_complex *)(&zero_padded_vector[0]),NULL,target.K,1,-1,FFTW_MEASURE);
      plan_vczdft=fftw_plan_many_dft(1,&target.Ppp,target.Kpp,(fftw_complex *)(&xzscratch[0]),NULL,target.Kpp,1,(fftw_complex *)(&xzscratch[0]),NULL,target.Kpp,1,-1,FFTW_MEASURE);

/* Plan creation for the iDFT of the vector components */
      plan_vcxdfti=fftw_plan_many_dft(1,&target.Kpp,target.P,(fftw_complex *)(&xzscratch[0]),NULL,1,target.Kpp,(fftw_complex *)(&xzscratch[0]),NULL,1,target.Kpp,+1,FFTW_MEASURE);
      plan_vcydfti=fftw_plan_many_dft(1,&target.Jpp,target.K,(fftw_complex *)(&zero_padded_vector[0]),NULL,target.K,1,(fftw_complex *)(&zero_padded_vector[0]),NULL,target.K,1,+1,FFTW_MEASURE);
      plan_vczdfti=fftw_plan_many_dft(1,&target.Ppp,target.Kpp,(fftw_complex *)(&xzscratch[0]),NULL,target.Kpp,1,(fftw_complex *)(&xzscratch[0]),NULL,target.Kpp,1,+1,FFTW_MEASURE);
   }
   else{
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Destory FFTW plans */
      fftw_destroy_plan(plan_tcxdft);
      fftw_destroy_plan(plan_tcydft);
      fftw_destroy_plan(plan_tczdft);
      fftw_destroy_plan(plan_vcxdft);
      fftw_destroy_plan(plan_vcydft);
      fftw_destroy_plan(plan_vczdft);
      fftw_destroy_plan(plan_vcxdfti);
      fftw_destroy_plan(plan_vcydfti);
      fftw_destroy_plan(plan_vczdfti);
      fftw_cleanup();
   }
}
