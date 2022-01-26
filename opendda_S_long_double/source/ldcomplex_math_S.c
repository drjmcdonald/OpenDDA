/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Functions to perform ldcomplex maths

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "ldcomplex_math_S.h"

long double ldcomplex_abs(ldcomplex in){
/* Absolute magnitude
   in=a+ib
   |in|=sqrt(a^2+b^2) */

   long double out;

   out=sqrtl(in.dat[0]*in.dat[0]+in.dat[1]*in.dat[1]);

   return out;
}

long double ldcomplex_abs2(ldcomplex in){
/* in=a+ib
   |in|^2=a^2+b^2 */

   long double out;

   out=in.dat[0]*in.dat[0]+in.dat[1]*in.dat[1];

   return out;
}

ldcomplex ldcomplex_add(ldcomplex in0,ldcomplex in1){
/* Add complex numbers
   in0=a+ib
   n1=c+id
   out=(a+c)+i(b+d) */

   ldcomplex out;
                          
   out.dat[0]=in0.dat[0]+in1.dat[0];
   out.dat[1]=in0.dat[1]+in1.dat[1];

   return out;
}

ldcomplex ldcomplex_sub(ldcomplex in0,ldcomplex in1){
/* Subtract complex numbers
   in0=a+ib
   in1=c+id
   out=(a-c)+i(b-d) */

   ldcomplex out;
                          
   out.dat[0]=in0.dat[0]-in1.dat[0];
   out.dat[1]=in0.dat[1]-in1.dat[1];

   return out;
}

ldcomplex ldcomplex_mul(ldcomplex in0,ldcomplex in1){
/* Multiply complex numbers
   in0=a+ib
   in1=c+id
   out=(ac-bd)+i(ad-bc) */

   ldcomplex out;
                          
   out.dat[0]=in0.dat[0]*in1.dat[0]-in0.dat[1]*in1.dat[1];
   out.dat[1]=in0.dat[0]*in1.dat[1]+in0.dat[1]*in1.dat[0];

   return out;
}

ldcomplex ldcomplex_div(ldcomplex in0,ldcomplex in1){
/* Divide complex numbers
   in0=a+ib
   in1=c+id
   out=[(a+ib)/(c+id)]*[(c-id)/(c-id)]=[(ac+bd)+i(bc-ad)]/[c^2+d^2] */

   ldcomplex out;
   long double scale;

   scale=1.0L/ldcomplex_abs2(in1);
   out.dat[0]=(in0.dat[0]*in1.dat[0]+in0.dat[1]*in1.dat[1])*scale;
   out.dat[1]=(in0.dat[1]*in1.dat[0]-in0.dat[0]*in1.dat[1])*scale;

   return out;
}

ldcomplex ldcomplex_scale(ldcomplex in,long double scale){
/* Scale by a long double
   in=(a*scale)+i(b*scale) */
                 
   in.dat[0]*=scale;
   in.dat[1]*=scale;

   return in;
}

ldcomplex ldcomplex_no_r_scale(ldcomplex in,long double scale){
/* Scale by a long double
   in=a+i(b*scale) */
                          
   in.dat[1]*=scale;

   return in;
}

long double ldcomplex_2vectornorm(ldcomplex *in){
/* 2-vector-norm
   out=sqrt[sum_{i=0}^{n-1}(|in[i]|^2)] */

   long int index0i;
   long double out=0.0L;

   for(index0i=0;index0i<target.Nd3;index0i++){
      out+=ldcomplex_abs2(in[index0i]);
   }
   out=sqrtl(out);

   return out;
}

ldcomplex ldcomplex_exp(ldcomplex in){
/* Complex exponential
   in=a+ib
   out=exp(in)=exp(a)*[cos(b)+isin(b)]*/

   long double a=expl(in.dat[0]);
   long double b=in.dat[1];
   ldcomplex out;

   out.dat[0]=a*cosl(b);
   out.dat[1]=a*sinl(b);

   return out;
}

ldcomplex ldcomplex_no_r_exp(ldcomplex in){
/* Complex Exponential where real part of 'in' is zero
   in=a+ib=0.0+ib
   out=exp(in)=exp(a)*[cos(b)+isin(b)]
              =exp(0.0)*[cos(b)+isin(b)]
              =1.0*[cos(b)+isin(b)]
              =cos(b)+isin(b) */                             

   long double b=in.dat[1];
   ldcomplex out;

   out.dat[0]=cosl(b);
   out.dat[1]=sinl(b);

   return out;
}

ldcomplex ldcomplex_conj(ldcomplex in){
/* Complex conjugate of input
   in=a+ib
   out=a-ib */

   in.dat[1]*=-1.0L;

   return in;
}

ldcomplex ldcomplex_sqrt(ldcomplex in){
/* sqrt[a+ib]=p+iq=
   p=[1/sqrt(2)]sqrt[sqrt(a^2+b^2)+a]
   q=[(sign of b)/sqrt(2)]sqrt[sqrt(a^2+b^2)-a]

   sqrt[(sqrt{a^2+b^2}+a)/(2)] +/- i*sqrt[(sqrt{a^2+b^2}-a)/(2)]
   where the sign is the same as `b' */

   ldcomplex out;
   long double sign,temp,fa,fb;

   fa=fabsl(in.dat[0]);
   fb=fabsl(in.dat[1]);
   if((fa<LDBL_EPSILON)&&(fb<LDBL_EPSILON)){ /* Check for zeroes */
      out.dat[0]=0.0L;
      out.dat[1]=0.0L;
   }
   else if(fb<LDBL_EPSILON){ /* Imaginary part ~zero */
      out.dat[0]=sqrtl(in.dat[0]);
      out.dat[1]=0.0L;
   }
   else{
      sign=in.dat[1]/fb; /* sign of b */
      temp=ldcomplex_abs(in); /* sqrt{a^2+b^2} */
      /* sqrt[(sqrt{a^2+b^2}+a)/(2)] */
      out.dat[0]=sqrtl((temp+in.dat[0])/2.0L);
      /* (sign of b)*sqrt[(sqrt{a^2+b^2}-a)/(2)] */
      out.dat[1]=sign*sqrtl((temp-in.dat[0])/2.0L);
   }

   return out;
}
