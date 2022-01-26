/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Functions to perform fcomplex maths

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "fcomplex_math_S.h"

float fcomplex_abs(fcomplex in){
/* Absolute magnitude
   in=a+ib
   |in|=sqrt(a^2+b^2) */

   float out;

   out=sqrtf(in.dat[0]*in.dat[0]+in.dat[1]*in.dat[1]);

   return out;
}

float fcomplex_abs2(fcomplex in){
/* in=a+ib
   |in|^2=a^2+b^2 */

   float out;

   out=in.dat[0]*in.dat[0]+in.dat[1]*in.dat[1];

   return out;
}

fcomplex fcomplex_add(fcomplex in0,fcomplex in1){
/* Add complex numbers
   in0=a+ib
   n1=c+id
   out=(a+c)+i(b+d) */

   fcomplex out;
                          
   out.dat[0]=in0.dat[0]+in1.dat[0];
   out.dat[1]=in0.dat[1]+in1.dat[1];

   return out;
}

fcomplex fcomplex_sub(fcomplex in0,fcomplex in1){
/* Subtract complex numbers
   in0=a+ib
   in1=c+id
   out=(a-c)+i(b-d) */

   fcomplex out;
                          
   out.dat[0]=in0.dat[0]-in1.dat[0];
   out.dat[1]=in0.dat[1]-in1.dat[1];

   return out;
}

fcomplex fcomplex_mul(fcomplex in0,fcomplex in1){
/* Multiply complex numbers
   in0=a+ib
   in1=c+id
   out=(ac-bd)+i(ad-bc) */

   fcomplex out;
                          
   out.dat[0]=in0.dat[0]*in1.dat[0]-in0.dat[1]*in1.dat[1];
   out.dat[1]=in0.dat[0]*in1.dat[1]+in0.dat[1]*in1.dat[0];

   return out;
}

fcomplex fcomplex_div(fcomplex in0,fcomplex in1){
/* Divide complex numbers
   in0=a+ib
   in1=c+id
   out=[(a+ib)/(c+id)]*[(c-id)/(c-id)]=[(ac+bd)+i(bc-ad)]/[c^2+d^2] */

   fcomplex out;
   float scale;

   scale=1.0F/fcomplex_abs2(in1);
   out.dat[0]=(in0.dat[0]*in1.dat[0]+in0.dat[1]*in1.dat[1])*scale;
   out.dat[1]=(in0.dat[1]*in1.dat[0]-in0.dat[0]*in1.dat[1])*scale;

   return out;
}

fcomplex fcomplex_scale(fcomplex in,float scale){
/* Scale by a float
   in=(a*scale)+i(b*scale) */
                 
   in.dat[0]*=scale;
   in.dat[1]*=scale;

   return in;
}

fcomplex fcomplex_no_r_scale(fcomplex in,float scale){
/* Scale by a float
   in=a+i(b*scale) */
                          
   in.dat[1]*=scale;

   return in;
}

float fcomplex_2vectornorm(fcomplex *in){
/* 2-vector-norm
   out=sqrt[sum_{i=0}^{n-1}(|in[i]|^2)] */

   long int index0i;
   float out=0.0F;

   for(index0i=0;index0i<target.Nd3;index0i++){
      out+=fcomplex_abs2(in[index0i]);
   }
   out=sqrtf(out);

   return out;
}

fcomplex fcomplex_exp(fcomplex in){
/* Complex exponential
   in=a+ib
   out=exp(in)=exp(a)*[cos(b)+isin(b)]*/

   float a=expf(in.dat[0]);
   float b=in.dat[1];
   fcomplex out;

   out.dat[0]=a*cosf(b);
   out.dat[1]=a*sinf(b);

   return out;
}

fcomplex fcomplex_no_r_exp(fcomplex in){
/* Complex Exponential where real part of 'in' is zero
   in=a+ib=0.0+ib
   out=exp(in)=exp(a)*[cos(b)+isin(b)]
              =exp(0.0)*[cos(b)+isin(b)]
              =1.0*[cos(b)+isin(b)]
              =cos(b)+isin(b) */                             

   float b=in.dat[1];
   fcomplex out;

   out.dat[0]=cosf(b);
   out.dat[1]=sinf(b);

   return out;
}

fcomplex fcomplex_conj(fcomplex in){
/* Complex conjugate of input
   in=a+ib
   out=a-ib */

   in.dat[1]*=-1.0F;

   return in;
}

fcomplex fcomplex_sqrt(fcomplex in){
/* sqrt[a+ib]=p+iq=
   p=[1/sqrt(2)]sqrt[sqrt(a^2+b^2)+a]
   q=[(sign of b)/sqrt(2)]sqrt[sqrt(a^2+b^2)-a]

   sqrt[(sqrt{a^2+b^2}+a)/(2)] +/- i*sqrt[(sqrt{a^2+b^2}-a)/(2)]
   where the sign is the same as `b' */

   fcomplex out;
   float sign,temp,fa,fb;

   fa=fabsf(in.dat[0]);
   fb=fabsf(in.dat[1]);
   if((fa<FLT_EPSILON)&&(fb<FLT_EPSILON)){ /* Check for zeroes */
      out.dat[0]=0.0F;
      out.dat[1]=0.0F;
   }
   else if(fb<FLT_EPSILON){ /* Imaginary part ~zero */
      out.dat[0]=sqrtf(in.dat[0]);
      out.dat[1]=0.0F;
   }
   else{
      sign=in.dat[1]/fb; /* sign of b */
      temp=fcomplex_abs(in); /* sqrt{a^2+b^2} */
      /* sqrt[(sqrt{a^2+b^2}+a)/(2)] */
      out.dat[0]=sqrtf((temp+in.dat[0])/2.0F);
      /* (sign of b)*sqrt[(sqrt{a^2+b^2}-a)/(2)] */
      out.dat[1]=sign*sqrtf((temp-in.dat[0])/2.0F);
   }

   return out;
}
