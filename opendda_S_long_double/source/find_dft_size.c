/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   FFTW is best at handling sizes of the form 2^{a}3^{b}5^{c}7^{d}11^{e}13^{f},
   where (e + f) is either 0 or 1, and the other exponents are arbitrary.
   Finds the smallest size >= `input' that fits the criteria.

   Find the sizes `kf', `jf' and `pf', the zero-fill sizes
   that extend the system to Kpp, Jpp And Ppp to permit the use
   of fast DFT algorithms

   K: no. of sites in x direction, no. of sites in each line
   J: no. of sites in y direction, no. of lines in a plane
   P: no. of sites in z direction, no. of planes in array
   
   Kpp: Kpp=Kp+kf
   Jpp: Jpp=Jp+jf
   Ppp:  Ppp=Pp+pf

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "find_dft_size.h"

void find_dft_size(int X,int *Xp,int *Xpp,int *xf){

   int size,tmp;

   *Xp=2*X-1; /* General Toeplitz to circulant conversion is X to 2X-1 */
   size=*Xp;*xf=0;
   while(size!=1){
      size=(*Xp)+(*xf);
      tmp=size;
      while((size&1)==0){ /* While even, factor 2 */
         size=size>>1; /* Bit-shift divide */
      }
      if(size==1){break;} /* Break if factored */
      while(size%3==0){ /* Factor 3 */
         size/=3;
      }
      if(size==1){break;} /* Break if factored */
      while(size%5==0){ /* Factor 5 */
         size/=5;
      }
      if(size==1){break;} /* Break if factored */
      while(size%7==0){ /* Factor 7 */
         size/=7;
      }
      if(size==1){break;} /* Break if factored */
      if(size%11==0){ /* Factor 11 (Only 1 of 11 or 13 allowed ) */
         size/=11;
      }
      else if(size%13==0){ /* Factor 13 (Only 1 of 11 or 13 allowed ) */
         size/=13;
      }
      if(size==1){break;} /* Break if factored */
      else{
         (*xf)++; /* If not factored increment zero-fill and try again */
      }
   }
   *Xpp=tmp; /* Set Xpp to the size */
}
