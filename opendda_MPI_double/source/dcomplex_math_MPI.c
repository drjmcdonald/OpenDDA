/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Functions to perform dcomplex maths

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dcomplex_math_MPI.h"

double dcomplex_abs(dcomplex in){
/* Absolute magnitude
   in=a+ib
   |in|=sqrt(a^2+b^2) */

   double out;

   out=sqrt(in.dat[0]*in.dat[0]+in.dat[1]*in.dat[1]);

   return out;
}

double dcomplex_abs2(dcomplex in){
/* in=a+ib
   |in|^2=a^2+b^2 */

   double out;

   out=in.dat[0]*in.dat[0]+in.dat[1]*in.dat[1];

   return out;
}

dcomplex dcomplex_add(dcomplex in0,dcomplex in1){
/* Add complex numbers
   in0=a+ib
   n1=c+id
   out=(a+c)+i(b+d) */

   dcomplex out;
                          
   out.dat[0]=in0.dat[0]+in1.dat[0];
   out.dat[1]=in0.dat[1]+in1.dat[1];

   return out;
}

dcomplex dcomplex_sub(dcomplex in0,dcomplex in1){
/* Subtract complex numbers
   in0=a+ib
   in1=c+id
   out=(a-c)+i(b-d) */

   dcomplex out;
                          
   out.dat[0]=in0.dat[0]-in1.dat[0];
   out.dat[1]=in0.dat[1]-in1.dat[1];

   return out;
}

dcomplex dcomplex_mul(dcomplex in0,dcomplex in1){
/* Multiply complex numbers
   in0=a+ib
   in1=c+id
   out=(ac-bd)+i(ad-bc) */

   dcomplex out;
                          
   out.dat[0]=in0.dat[0]*in1.dat[0]-in0.dat[1]*in1.dat[1];
   out.dat[1]=in0.dat[0]*in1.dat[1]+in0.dat[1]*in1.dat[0];

   return out;
}

dcomplex dcomplex_div(dcomplex in0,dcomplex in1){
/* Divide complex numbers
   in0=a+ib
   in1=c+id
   out=[(a+ib)/(c+id)]*[(c-id)/(c-id)]=[(ac+bd)+i(bc-ad)]/[c^2+d^2] */

   dcomplex out;
   double scale;

   scale=1.0/dcomplex_abs2(in1);
   out.dat[0]=(in0.dat[0]*in1.dat[0]+in0.dat[1]*in1.dat[1])*scale;
   out.dat[1]=(in0.dat[1]*in1.dat[0]-in0.dat[0]*in1.dat[1])*scale;

   return out;
}

dcomplex dcomplex_scale(dcomplex in,double scale){
/* Scale by a double
   in=(a*scale)+i(b*scale) */
                 
   in.dat[0]*=scale;
   in.dat[1]*=scale;

   return in;
}

dcomplex dcomplex_no_r_scale(dcomplex in,double scale){
/* Scale by a double
   in=a+i(b*scale) */
                          
   in.dat[1]*=scale;

   return in;
}

double dcomplex_2vectornorm(dcomplex *in){
/* 2-vector-norm
   out=sqrt[sum_{i=0}^{n-1}(|in[i]|^2)] */

   long int index0i;
   double out=0.0,result=0.0;

   for(index0i=0;index0i<parallel.alloc_vector3;index0i++){
      out+=dcomplex_abs2(in[index0i]);
   }
   MPI_Allreduce(&out,&result,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   result=sqrt(result);

   return result;
}

dcomplex dcomplex_exp(dcomplex in){
/* Complex exponential
   in=a+ib
   out=exp(in)=exp(a)*[cos(b)+isin(b)]*/

   double a=exp(in.dat[0]);
   double b=in.dat[1];
   dcomplex out;

   out.dat[0]=a*cos(b);
   out.dat[1]=a*sin(b);

   return out;
}

dcomplex dcomplex_no_r_exp(dcomplex in){
/* Complex Exponential where real part of 'in' is zero
   in=a+ib=0.0+ib
   out=exp(in)=exp(a)*[cos(b)+isin(b)]
              =exp(0.0)*[cos(b)+isin(b)]
              =1.0*[cos(b)+isin(b)]
              =cos(b)+isin(b) */                             

   double b=in.dat[1];
   dcomplex out;

   out.dat[0]=cos(b);
   out.dat[1]=sin(b);

   return out;
}

dcomplex dcomplex_conj(dcomplex in){
/* Complex conjugate of input
   in=a+ib
   out=a-ib */

   in.dat[1]*=-1.0;

   return in;
}

dcomplex dcomplex_sqrt(dcomplex in){
/* sqrt[a+ib]=p+iq=
   p=[1/sqrt(2)]sqrt[sqrt(a^2+b^2)+a]
   q=[(sign of b)/sqrt(2)]sqrt[sqrt(a^2+b^2)-a]

   sqrt[(sqrt{a^2+b^2}+a)/(2)] +/- i*sqrt[(sqrt{a^2+b^2}-a)/(2)]
   where the sign is the same as `b' */

   dcomplex out;
   double sign,temp,fa,fb;

   fa=fabs(in.dat[0]);
   fb=fabs(in.dat[1]);
   if((fa<DBL_EPSILON)&&(fb<DBL_EPSILON)){ /* Check for zeroes */
      out.dat[0]=0.0;
      out.dat[1]=0.0;
   }
   else if(fb<DBL_EPSILON){ /* Imaginary part ~zero */
      out.dat[0]=sqrt(in.dat[0]);
      out.dat[1]=0.0;
   }
   else{
      sign=in.dat[1]/fb; /* sign of b */
      temp=dcomplex_abs(in); /* sqrt{a^2+b^2} */
      /* sqrt[(sqrt{a^2+b^2}+a)/(2)] */
      out.dat[0]=sqrt((temp+in.dat[0])/2.0);
      /* (sign of b)*sqrt[(sqrt{a^2+b^2}-a)/(2)] */
      out.dat[1]=sign*sqrt((temp-in.dat[0])/2.0);
   }

   return out;
}
