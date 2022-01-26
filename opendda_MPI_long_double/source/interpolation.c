/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Interpolation and extrapolation

   If the number of data > 4

      [Point <= Range minimum OR Point >= Range maximum]
      Linear extrapolation: [DANGEROUS]
      Create a tangent line at the end of the known data and extend it beyond
      that limit y=y0+m(x-x0) where m is the slope i.e. first derivative

      [Range minimum < Point < Range maximum]
      **** Univariate interpolation using the improved Akima Method ****

   If the number of data = 2

      [All ranges]
      Linear interpolation and extrapolation: [DANGEROUS]
      Know (x0,y0) and (x1,y1) then at given x:
      y=y0(1-a)+a*y1 where a=(x-x0)/(x1-x0)


   If the number of data = 3

      [Point <= Range minimum OR Point >= Range maximum]
         Linear extrapolation: [DANGEROUS]
         Create a tangent line at the end of the known data and extend it beyond
         that limit y=y0+m(x-x0) where m is the slope i.e. first derivative
      [Range minimum < Point < Range maximum]
         Quadratic interpolation

   If the number of data = 4

      [Point <= Range minimum OR Point >= Range maximum]
         Linear extrapolation: [DANGEROUS]
         Create a tangent line at the end of the known data and extend it beyond
         that limit y=y0+m(x-x0) where m is the slope i.e. first derivative
      [Range minimum < Point < Range maximum]
         Cubic interpolation

   The function can be called to perform interpolation in one call at all points

      interpolation(degree,datasize,x,y,numinterp,xi,yi,component);

   OR

   The function can be called to perform interpolation `1 at a time' in repeated calls

      for(i=0;i<numinterp;i++){
         interpolation(degree,datasize,x,y,1,&xi[i],&yi[i],component);
      }

   where

      degree: degree of the polynomials for the interpolating function
      x: x data
      datasize: is the number of data points
      y: y data i.e. y=f(x)
      numinterp: is the number of points for which interpolation is desired
      xi: x points for which interpolation is desired
      yi: interpolated y values
      component: complex numbers 0: interpolate real part
                                 1: interpolate imaginary part

   Akima, H., A method of univariate interpolation that has the accuracy
   of a third-degree polynomial, Association for Computing Machinery,
   Transactions on Mathematical Software
   ACM TOMS, Vol. 17, No. 3, Sept. 1991, pp. 341-366.
   Note that the equation numbers given refer to this paper
   Adapted from http://www.netlib.org/toms/697 

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "interpolation.h"

void interpolation(int degree,int size,long double *xdata,ldcomplex *ydata,int points,long double *xinterp,ldcomplex *yinterp,int part){
/* degree   degree of the polynomials for the interpolating function
   size     number of data points (must be >= 2)
   xdata    x data
   ydata    y data i.e. y=f(x)
   points   the number of points for which interpolation is desired
   xinterp  x points for which interpolation is desired
   yinterp  interpolated y values
   part     complex numbers 0: interpolate real part
                            1: interpolate imaginary part */

   int i,pe,endpoint,left,right,mid,interval=-1,previous_interval=-1,id0,id1,id2,id3;
   long double first_derivative,primary_estimate,volatility_factor,distance_factor,weighting_factor;
   long double sum_primary_estimate_finite,sum_primary_estimate_infinite,final_estimate,final_estimate0,slope;
   long double sum_weighting_factor_finite,sum_weighting_factor_infinite,denom,d0,d1,d2,temp;
   long double x10,x13,x20,x21,x30,x32,b0,b1,v0,v1,v2,v3,sx,sy,sxy,sxs,sx2,a0,a1,a2,a3,u,v,a,b,c,d,e,f;
   char string[100];

/* Description of variables
   i: Loop control variable which loops over the desired interpolation points 0 to points-1

   pe: Loop over the 4 primary estimates i.e. there are 4 sets of 4 data points that include the data point P_{i}:
       P_{i-3} to P_{i}   given by [Equation 12] y'_{imm}=F(i,i-3,i-2,i-1)
       P_{i-2} to P_{i+1} given by [Equation 12] y'_{im}=F(i,i-2,i-1,i+1)
       P_{i-1} to P_{i+2} given by [Equation 12] y'_{ip}=F(i,i-1,i+1,i+2)
       P_{i} to P_{i+3}   given by [Equation 12] y'_{ipp}=F(i,i+1,i+2,i+3)
       The 4 primary estimates are calculated as the first derivative of a 3rd-degree polynomial fitted to each set
       Note that for F(i,j,k,l) the first index must be id0=i, which is the point in question, however, the
       remaining indices can be in arbitrary order i.e. id1,id2 & id3 can be in any order

       pe=1  i-3,i-2,i-1,i
       pe=2  i-2,i-1,i,i+1
       pe=3  i-1,i,i+1,i+2
       pe=4  i,i+1,i+2,i+3

   id0,id1,id2,id3: Store the indices for the 4 primary estimates pe=1..4

   left,right,mid,
   interval,
   previous_interval: Used in the binary search algorithm to locate and store the xdata interval that includes
                      the desired interpolation point

   first_derivative: Stores the first derivative at point P_{i} of a 3rd-degree polynomial fitted to a set of 4 data
                     points P_{i}, P_{j}, P_{k} & P_{l}, which is given by:
                     [Equation 8] first_derivative=F(i,j,k,l)=F(0,1,2,3)=[f1 + f2 + f3]/[f4]

                     f1=(y1-y0)*(x2-x0)^{2}*(x3-x0)^{2}*(x3-x2)
                     f2=(y2-y0)*(x3-x0)^{2}*(x1-x0)^{2}*(x1-x3)
                     f3=(y3-y0)*(x1-x0)^{2}*(x2-x0)^{2}*(x2-x1)
                     f4=(x1-x0)*(x2-x0)*(x3-x0)*(x2-x1)*(x3-x2)*(x3-x1)

   primary_estimate: The primary estimate for the first derivative of the 4 3rd-degree polynomials fitted to
                     the 4 sets of data points that include P_{i}
                     P_{i-3} to P_{i}   given by [Equation 12] y'_{imm}=F(i,i-3,i-2,i-1)
                     P_{i-2} to P_{i+1} given by [Equation 12] y'_{im}=F(i,i-2,i-1,i+1)
                     P_{i-1} to P_{i+2} given by [Equation 12] y'_{ip}=F(i,i-1,i+1,i+2)
                     P_{i} to P_{i+3}   given by [Equation 12] y'_{ipp}=F(i,i+1,i+2,i+3)

   volatility_factor,
   v0,v1,v2,v3: The volatility factor [Equation 9] V=sum[y-(b0+b1x)]^{2}
                v0=y0-(b0+b1*x0)
                v1=y1-(b0+b1*x1)
                v2=y2-(b0+b1*x2)
                v3=y3-(b0+b1*x3)
                volatility_factor=v0*v0+v1*v1+v2*v2+v3*v3

   b0,b1: The coefficients of the linear least-squares fit to the 4 data points
          [Equation 10] b0=(sum[x^{2}]sum[y]-sum[x]sum[y])/(4*sum[x^{2}]-sum[x]^{2})
          [Equation 10] b1=(4*sum[xy]-sum[x]sum[y])/(4*sum[x^{2}]-sum[x]^{2}

          x0=xdata[id0]  y0=ydata[id0]
          x1=xdata[id1]  y0=ydata[id0]
          x2=xdata[id2]  y2=ydata[id2]
          x3=xdata[id3]  y3=ydata[id3]

          b0=[sxs*sy-sx*sxy]/4*sxs-sx2
          b1=[4*sxy-sx*sy]/4*sxs-sx2

   sx,sy,sxy,
   sxs,sx2: Sums for the calculation of b0 & b1
            sx=sum[x]=(x0+x1+x2+x3)
            sy=sum[y]=(y0+y1+y2+y3)
            sxy=sum[xy]=(x0y0+x1y1+x2y2+x3y3)
            sxs=sum[x^{2}]=(x0x0+x1x1+x2x2+x3x3)
            sx2=sum[x]^{2}=sx*sx

   distance_factor,
   d0,d1,d2: The distance factor [Equation 11] D(0,1,2,3)=(x1-x0)^{2}+(x2-x0)^{2}+(x3-x0)^{2}
             d0=x1-x0
             d1=x2-x0
             d2=x3-x0
             distance_factor=d0*d0+d1*d1+d2*d2

   weighting_factor: Each of the 4 primary estimates is weighted using the reciprocal of the product
                     of the volatility factor and the distance factor

                     pe=1 P_{i-3} to P_{i}   y'_{imm}=F(id0,id1,id2,id3)=F(i,i-3,i-2,i-1)
                          corresponding weight [Equation 13] w_{imm}=1/[V(i,i-3,i-2,i-1)D(i,i-3,i-2,i-1)]
                     pe=2 P_{i-2} to P_{i+1} y'_{im}=F(id0,id2,id3,id1)=F(i,i-2,i-1,i+1)
                          corresponding weight [Equation 13] w_{im}=1/[V(i,i-2,i-1,i+1)D(i,i-2,i-1,i+1)]
                     pe=3 P_{i-1} to P_{i+2} y'_{ip}=F(id0,id3,id1,id2)=F(i,i-1,i+1,i+2)
                          corresponding weight [Equation 13] w_{ip}=1/[V(i,i-1,i+1,i+2)D(i,i-1,i+1,i+2)]
                     pe=4 P_{i} to P_{i+3}   y'_{ipp}=F(id0,id1,id2,id3)=F(i,i+1,i+2,i+3)
                          corresponding weight [Equation 13] w_{ipp}=1/[V(i,i+1,i+2,i+3)D(i,i+1,i+2,i+3)]

   final_estimate,
   final_estimate0: The final estimate of the first derivative of the interpolating function:

                    [Equation 14] y'_{i}=[y'_{imm}w_{imm}+
                                          y'_{im}w_{im}+
                                          y'_{ip}w_{ip}+
                                          y'_{ipp}w_{ipp}]/
                                          [w_{imm}+w_{im}+w_{ip}+w_{ipp}]
                    final_estimate0 stores the value for endpoint=0
                    final_estimate stores the value for endpoint=1

   sum_primary_estimate_finite,
   sum_primary_estimate_infinite,
   sum_weighting_factor_finite,
   sum_weighting_factor_infinite:
                     pe=1 sum_primary_estimate=y'_{imm}w_{imm}
                          sum_weighting_factor=w_{imm}
                     pe=2 sum_primary_estimate=y'_{imm}w_{imm}+y'_{im}w_{im}
                          sum_weighting_factor=w_{imm}+w_{im}
                     pe=3 sum_primary_estimate=y'_{imm}w_{imm}+y'_{im}w_{im}+y'_{ip}w_{ip}
                          sum_weighting_factor=w_{imm}+w_{im}+w_{ip}
                     pe=4 sum_primary_estimate=y'_{imm}w_{imm}+y'_{im}w_{im}+y'_{ip}w_{ip}+y'_{ipp}w_{ipp}
                          sum_weighting_factor=w_{imm}+w_{im}+w_{ip}+w_{ipp}
                     Note that the difference between finite and infinite is whether or not the
                     volatility factor is essentially zero
                     weighting_factor=1/[V*D] if V~0 then the weight is infinite else it is finite

                     For finite weight
                        weighting_factor=1/[V*D]
                        sum_primary_estimate_finite+=primary_estimate*weighting_factor
                        sum_weighting_factor_finite+=weighting_factor

                     For infinite weight [V~0]

                        sum_primary_estimate_infinite+=primary_estimate
                        sum_weighting_factor_infinite+=1

   slope
   a0,a1,
   a2,a3,
   u,v: Used in the calculation of the coefficients of the 3rd OR higher degree polynomial
        i.e. The 3rd-degree polynomial for the y value in the interval x_{i} and x_{i+1} is

             [Equation 3] y=a0+a1(x-x_{i})+a2(x-x_{i})^{2}+a3(x-x_{i})^{3}

             [Equation 4] a0=y_{i}
             [Equation 4] a1=y'_{i}
             [Equation 4] a2=-[2(y'_{i}-m_{i})+(y'_{i+1}-m_{i})]/[x_{i+1}-x_{i}]
             [Equation 4] a3=[(y'_{i}-m_{i})+(y'_{i+1}-m_{i})]/[(x_{i+1}-x_{i})^{2}]

             where m_{i} is the slope of the line segment connecting P_{i} and P_{i+1}

             [Equation 5] m_{i}=(y_{i+1}-y_{i})/(x_{i+1}-x_{i})

             The nth-degree polynomial is represented by

             [Equation 19] v(u)=A0[u^{n}-u]+A1[(1-u)^{n}-(1-u)]
             [Equation 20] A0=[v'_{0}+(n-1)v'_{1}]/[n(n-2)]
             [Equation 20] A1=-[(n-1)v'_{0}+v'_{1}]/[n(n-2)]

             u=0,v=0,v'_{0}=(y'_{i}-m_{i})/(x_{i+1}-x_{i}) @ P_{i}
             u=1,v=0,v'_{1}=(y'_{i+1}-m_{i})/(x_{i+1}-x_{i}) @ P_{i}

             where m_{i} is the slope of the line segment connecting P_{i} and P_{i+1}

             [Equation 5] m_{i}=(y_{i+1}-y_{i})/(x_{i+1}-x_{i})

   a,b,c,d,e,f,
   x10,x13,
   x20,x21,
   x30,x32: If the number of data points > 4 then the Improved Akima method is used
            If the number of data points = 2 then linear interpolation and extrapolation
            If the number of data points = 3 then quadratic interpolation and linear extrapolation
            If the number of data points = 4 then cubic interpolation and linear extrapolation
            a,b,c,d,e,f store the coefficients for the interpolation

            Number of data points = 2
               Linear interpolation and extrapolation
               Know (x0,y0) and (x1,y1) then at given x:
               y=y0(1-a)+a*y1
               a=(x-x0)/(x1-x0)

            Number of data points = 3
               If Point <= Range minimum
                  Linear extrapolation
                  Create a tangent line at the end of the known data and extend it beyond that limit
                  y=y0+m(x-x0) where m is the slope i.e. first derivative

                  first_derivative=F(0,1,2)=[f1 + f2]/[f3 + f4]

                  f1=(y1-y0)*(x0-x2)^{2}=(ydata[1]-ydata[0])*x20*x20
                  f2=(y0-y2)*(x1-x0)^{2}=(ydata[0]-ydata[2])*x10*x10
                  f3=(x0-x2)^{2}(x1-x0)=x20*x20*x10
                  f4=(x1-x0)^{2}(x0-x2)=x10*x10*x20
                  x10=xdata[1]-xdata[0]
                  x20=xdata[0]-xdata[2]

               If Range minimum < Point < Range maximum
                  Quadratic Interpolation
                  Know (x0,y0), (x1,y1) and (x2,y2) then at given x:
                  y=y0*c(1-a)+y1*a(1-b)+y2*b(1-c)
                  a=(x-x0)/(x1-x0)
                  b=(x-x1)/(x2-x1)
                  c=(x-x2)/(x0-x2)

               If Point >= Range maximum
                  Linear extrapolation
                  Create a tangent line at the end of the known data and extend it beyond that limit
                  y=y0+m(x-x0) where m is the slope i.e. first derivative

                  first_derivative=F(size-1,size-2,size-3)=F(0,1,2)=[f1 + f2]/[f3 + f4]

                  f1=(y1-y0)*(x0-x2)^{2}=(ydata[size-2]-ydata[size-1])*x20*x20
                  f2=(y0-y2)*(x1-x0)^{2}=(ydata[size-1]-ydata[2])*x10*x10
                  f3=(x0-x2)^{2}(x1-x0)=x20*x20*x10
                  f4=(x1-x0)^{2}(x0-x2)=x10*x10*x20
                  x10=xdata[size-2]-xdata[size-1];
                  x20=xdata[size-1]-xdata[size-3];

            Number of data points = 4
               If Point <= Range minimum
                  Linear extrapolation
                  Create a tangent line at the end of the known data and extend it beyond that limit
                  y=y0+m(x-x0) where m is the slope i.e. first derivative

                  first_derivative=F(0,1,2,3)=[f1 + f2 + f3]/[f4]

                  f1=(y1-y0)*(x2-x0)^{2}*(x3-x0)^{2}*(x3-x2)=(ydata[1]-ydata[0])*x20*x20*x30*x30*x32
                  f2=(y2-y0)*(x3-x0)^{2}*(x1-x0)^{2}*(x1-x3)=(ydata[2]-ydata[0])*x30*x30*x10*x10*x13
                  f3=(y3-y0)*(x1-x0)^{2}*(x2-x0)^{2}*(x2-x1)=(ydata[3]-ydata[0])*x10*x10*x20*x20*x21
                  f4=(x1-x0)*(x2-x0)*(x3-x0)*(x2-x1)*(x3-x2)*(x3-x1)=x10*x20*x30*x21*x32*-x13
                  x10=xdata[1]-xdata[0];x13=xdata[1]-xdata[3];
                  x20=xdata[2]-xdata[0];x21=xdata[2]-xdata[1];
                  x30=xdata[3]-xdata[0];x32=xdata[3]-xdata[2];

               If Range minimum < Point < Range maximum
                  Cubic interpolation
                  Know (x0,y0), (x1,y1), (x2,y2) and (x3,y3) then at given x:
                  y=y0*ce(1-a)+y1*ad(1-b)+y2*bf(1-c)+y3*(1-d)(1-e)(1-f)
                  a=(x-x0)/(x1-x0)
                  b=(x-x1)/(x2-x1)
                  c=(x-x2)/(x0-x2)
                  d=(x-x3)/(x1-x3)
                  e=(x-x3)/(x0-x3)
                  f=(x-x3)/(x2-x3)

               If Point >= Range maximum
                  Linear extrapolation
                  Create a tangent line at the end of the known data and extend it beyond that limit
                  y=y0+m(x-x0) where m is the slope i.e. first derivative

                  first_derivative=F(size-1,size-2,size-3,size-4)=F(0,1,2,3)=[f1 + f2 + f3]/[f4]

                  f1=(y1-y0)*(x2-x0)^{2}*(x3-x0)^{2}*(x3-x2)=(ydata[size-2]-ydata[size-1])*x20*x20*x30*x30*x32
                  f2=(y2-y0)*(x3-x0)^{2}*(x1-x0)^{2}*(x1-x3)=(ydata[size-3]-ydata[size-1])*x30*x30*x10*x10*x13
                  f3=(y3-y0)*(x1-x0)^{2}*(x2-x0)^{2}*(x2-x1)=(ydata[size-4]-ydata[size-1])*x10*x10*x20*x20*x21
                  f4=(x1-x0)*(x2-x0)*(x3-x0)*(x2-x1)*(x3-x2)*(x3-x1)=x10*x20*x30*x21*x32*-x13
                  x10=xdata[size-2]-xdata[size-1];x13=xdata[size-2]-xdata[size-4];
                  x20=xdata[size-3]-xdata[size-1];x21=xdata[size-3]-xdata[size-2];
                  x30=xdata[size-4]-xdata[size-1];x32=xdata[size-4]-xdata[size-3];

   denom,temp: Temporary variables */

/* Sanity checks */
   if(size<2){ /* Number of x OR y must be >= 2 */
      print_error("Interpolation error","Improved Akima Error","The number of data points must be >= 2");
   }
   if(points<1){ /* Number of points for which interpolation is desired must be > 0 */
      print_error("Interpolation error","Improved Akima Error","Number of points for which interpolation is desired must be > 0");
   }
   for(i=0;i<(size-1);i++){ /* Check to ensure that the x data is sequentially increasing and that xdata[i+1]!=xdata[i] */
      if((xdata[i]>xdata[i+1])||((xdata[i+1]-xdata[i])<LDBL_EPSILON)){
         sprintf(string,"Improved Akima Error: %.*Lg>=%.*Lg",LDBLP,xdata[i],LDBLP,xdata[i+1]);
         print_error("Interpolation error",string,"The 'x data' must be sequentially increasing");
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* General case - 5 points or more */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
   if(size>4){ /* General case - 5 points or more */
      if(degree<3){degree=3;} /* 3 is the minimum value for the degree of the polynomials for the interpolating function */
      left=0; /* Acceleration: Only need to zero here as start point will move up through the array */
      for(i=0;i<points;i++){ /* Cycle through the points for which interpolation is desired */
         right=size-1;
         if((xinterp[i]-xdata[0])<LDBL_EPSILON){ /* Point <= Range minimum */
            interval=0;
         }
         else if(xinterp[i]<xdata[right]){ /* Range minimum < Point < Range maximum */
            /* Locate the interval that includes the desired point by binary search */
            mid=(int)((long double)(left+right)/2.0L);
            while(mid>left){
               if(xinterp[i]<xdata[mid]){right=mid;}
               else{left=mid;}
               mid=(int)((long double)(left+right)/2.0L);
            }
            interval=++mid;
         }
         else{ /* Point >= Range maximum */
            interval=size;
         }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Subcase 1: Point <= minimum of data range: Extrapolation */
         if(interval==0){ /* Point <= Range minimum */
/* [Equation 8] first_derivative=F(0,1,2,3)=[f1 + f2 + f3]/[f4]

   f1=(y1-y0)*(x2-x0)^{2}*(x3-x0)^{2}*(x3-x2)
   f2=(y2-y0)*(x3-x0)^{2}*(x1-x0)^{2}*(x1-x3)
   f3=(y3-y0)*(x1-x0)^{2}*(x2-x0)^{2}*(x2-x1)
   f4=(x1-x0)*(x2-x0)*(x3-x0)*(x2-x1)*(x3-x2)*(x3-x1)

   f1=(ydata[1]-ydata[0])*x20*x20*x30*x30*x32;
   f2=(ydata[2]-ydata[0])*x30*x30*x10*x10*x13;
   f3=(ydata[3]-ydata[0])*x10*x10*x20*x20*x21;
   f4=x10*x20*x30*x21*x32*-x13;
   Note: Only required when the interval is not the same as the one for the previous point */
            if(interval!=previous_interval){
               previous_interval=interval;
               x10=xdata[1]-xdata[0];x13=xdata[1]-xdata[3];
               x20=xdata[2]-xdata[0];x21=xdata[2]-xdata[1];
               x30=xdata[3]-xdata[0];x32=xdata[3]-xdata[2];
               first_derivative=(((ydata[1].dat[part]-ydata[0].dat[part])*x20*x20*x30*x30*x32)+
                                 ((ydata[2].dat[part]-ydata[0].dat[part])*x30*x30*x10*x10*x13)+
                              ((ydata[3].dat[part]-ydata[0].dat[part])*x10*x10*x20*x20*x21))/(x10*x20*x30*x21*x32*-x13);
            }
            /* Evaluate yinterp[i] */
            yinterp[i].dat[part]=ydata[0].dat[part]+first_derivative*(xinterp[i]-xdata[0]);
         }
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Subcase 2: Point >= maximum of data range: Extrapolation */
         else if(interval==size){ /* Point >= Range maximum */
/* [Equation 8] first_derivative=F(size-1,size-2,size-3,size-4)=F(0,1,2,3)=[f1 + f2 + f3]/[f4]

   f1=(y1-y0)*(x2-x0)^{2}*(x3-x0)^{2}*(x3-x2)
   f2=(y2-y0)*(x3-x0)^{2}*(x1-x0)^{2}*(x1-x3)
   f3=(y3-y0)*(x1-x0)^{2}*(x2-x0)^{2}*(x2-x1)
   f4=(x1-x0)*(x2-x0)*(x3-x0)*(x2-x1)*(x3-x2)*(x3-x1)

   f1=(ydata[size-2]-ydata[size-1])*x20*x20*x30*x30*x32;
   f2=(ydata[size-3]-ydata[size-1])*x30*x30*x10*x10*x13;
   f3=(ydata[size-4]-ydata[size-1])*x10*x10*x20*x20*x21;
   f4=x10*x20*x30*x21*x32*-x13;
   Note: Only required when the interval is not the same as the one for the previous point */
            if(interval!=previous_interval){
               previous_interval=interval;
               x10=xdata[size-2]-xdata[size-1];x13=xdata[size-2]-xdata[size-4];
               x20=xdata[size-3]-xdata[size-1];x21=xdata[size-3]-xdata[size-2];
               x30=xdata[size-4]-xdata[size-1];x32=xdata[size-4]-xdata[size-3];
               first_derivative=(((ydata[size-2].dat[part]-ydata[size-1].dat[part])*x20*x20*x30*x30*x32)+
                                 ((ydata[size-3].dat[part]-ydata[size-1].dat[part])*x30*x30*x10*x10*x13)+
                                 ((ydata[size-4].dat[part]-ydata[size-1].dat[part])*x10*x10*x20*x20*x21))/(x10*x20*x30*x21*x32*-x13);
            }
            /* Evaluate yinterp[i] */
            yinterp[i].dat[part]=ydata[size-1].dat[part]+first_derivative*(xinterp[i]-xdata[size-1]);
         }
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Subcase 3: Range minimum < Point < Range maximum: Interpolation */
         else{ /* Range minimum < Point < Range maximum */
/* Calculates the coefficients of the third-degree polynomial (for
   degree<=3) or the factors for the higher-degree polynomials (for
   degree>3), when the interval is not the same as the one for the
   previous desired point. */
            if(interval!=previous_interval){
               previous_interval=interval;
               for(endpoint=0;endpoint<2;endpoint++){ /* Loop over endpoints */
                  /* Calculates the estimate of the first derivative at an endpoints */
                  id0=interval+endpoint-1; /* i */
/* The final estimate of the first derivative of the interpolating function is

   y' is the first derivative of the function i.e. the slope of the curve
   y'_{i}=[y'_{imm}w_{imm}+y'_{im}w_{im}+y'_{ip}w_{ip}+y'_{ipp}w_{ipp}]/[w_{imm}+w_{im}+w_{ip}+w_{ipp}]

   pe=1 sum_primary_estimate=y'_{imm}w_{imm}
        sum_weighting_factor=w_{imm}
   pe=2 sum_primary_estimate=y'_{imm}w_{imm}+y'_{im}w_{im}
        sum_weighting_factor=w_{imm}+w_{im}
   pe=3 sum_primary_estimate=y'_{imm}w_{imm}+y'_{im}w_{im}+y'_{ip}w_{ip}
        sum_weighting_factor=w_{imm}+w_{im}+w_{ip}
   pe=4 sum_primary_estimate=y'_{imm}w_{imm}+y'_{im}w_{im}+y'_{ip}w_{ip}+y'_{ipp}w_{ipp}
        sum_weighting_factor=w_{imm}+w_{im}+w_{ip}+w_{ipp}
   Note that the difference between finite and infinite is whether or not the volatility factor
   is essentially zero
   weighting_factor=1/[V*D] if V~0 then the weight is infinite otherwise it is finite */
                  sum_primary_estimate_finite=0.0L; /* Initalise value to zero */
                  sum_weighting_factor_finite=0.0L; /* Initalise value to zero */
                  sum_primary_estimate_infinite=0.0L; /* Initalise value to zero */
                  sum_weighting_factor_infinite=0.0L; /* Initalise value to zero */
/* There are 4 sets of 4 data points that include the data point P_{i}:
   P_{i-3} to P_{i}   given by [Equation 12] y'_{imm}=F(i,i-3,i-2,i-1)
   P_{i-2} to P_{i+1} given by [Equation 12] y'_{im}=F(i,i-2,i-1,i+1)
   P_{i-1} to P_{i+2} given by [Equation 12] y'_{ip}=F(i,i-1,i+1,i+2)
   P_{i} to P_{i+3}   given by [Equation 12] y'_{ipp}=F(i,i+1,i+2,i+3)
   Calculate the 4 primary estimates as the first derivative of a 3rd-degree polynomial fitted to each set
   Note that for F(i,j,k,l) the first index must be id0=i, which is the point in question, however, the
   remaining indices can be in arbitrary order i.e. id1,id2 & id3 can be in any order

   pe=1 i-3,i-2,i-1,i
        id0=i
        id1=i-3
        id2=i-2
        id3=i-1
   pe=2 i-2,i-1,i,i+1
        Only need to change id1=i-3 to id1=i+1
        id0=i
        id2=i-2
        id3=i-1
        id1=i+1
   pe=3 i-1,i,i+1,i+2
        Only need to change id2=i-2 to id2=i+2
        id0=i
        id3=i-1
        id1=i+1
        id2=i+2
   pe=4 i,i+1,i+2,i+3
        Only need to change id3=i-1 to id3=i+3
        id0=i
        id1=i+1
        id2=i+2
        id3=i+3 */
                  for(pe=1;pe<5;pe++){ /* Loop of the primary estimates */
                     if(pe==1){ /* i-3,i-2,i-1,i */
                        id1=id0-3;/* i-3 */id2=id0-2;/* i-2 */id3=id0-1;/* i-1 */
                     }
                     else if(pe==2){ /* i-2,i-1,i,i+1 */
                        id1=id0+1; /* id1=i-3 to id1=i+1 */
                     }
                     else if(pe==3){ /* i-1,i,i+1,i+2 */
                        id2=id0+2; /* id2=i-2 to id2=i+2 */
                     }
                     else{ /* i,i+1,i+2,i+3 */
                        id3=id0+3; /* id3=i-1 to id3=i+3 */
                     }
/* Only perform calculation of the primary estimates if all 4 data points are within the data range 0 to size-1 */
                     if((id1>=0)&&(id2>=0)&&(id3>=0)&&(id1<size)&&(id2<size)&&(id3<size)){
/* [Equation 8] first_derivative=F(0,1,2,3)=F(id0,id1,id2,id3)=[f1 + f2 + f3]/[f4]

   f1=(y1-y0)*(x2-x0)^{2}*(x3-x0)^{2}*(x3-x2)
   f2=(y2-y0)*(x3-x0)^{2}*(x1-x0)^{2}*(x1-x3)
   f3=(y3-y0)*(x1-x0)^{2}*(x2-x0)^{2}*(x2-x1)
   f4=(x1-x0)*(x2-x0)*(x3-x0)*(x2-x1)*(x3-x2)*(x3-x1)

   f1=(ydata[id1]-ydata[id0])*x20*x20*x30*x30*x32;
   f2=(ydata[id2]-ydata[id0])*x30*x30*x10*x10*x13;
   f3=(ydata[3]-ydata[id0])*x10*x10*x20*x20*x21;
   f4=x10*x20*x30*x21*x32*-x13; */
                        x10=xdata[id1]-xdata[id0];x13=xdata[id1]-xdata[id3];
                        x20=xdata[id2]-xdata[id0];x21=xdata[id2]-xdata[id1];
                        x30=xdata[id3]-xdata[id0];x32=xdata[id3]-xdata[id2];
                        primary_estimate=(((ydata[id1].dat[part]-ydata[id0].dat[part])*x20*x20*x30*x30*x32)+
                                          ((ydata[id2].dat[part]-ydata[id0].dat[part])*x30*x30*x10*x10*x13)+
                                          ((ydata[id3].dat[part]-ydata[id0].dat[part])*x10*x10*x20*x20*x21))/(x10*x20*x30*x21*x32*-x13);
/* Calculate the volatility factor V, the distance factor D
   [Equation 9] Volatility factor V=sum[y-(b0+b1x)]^{2}

   [Equation 10] b0=(sum[x^{2}]sum[y]-sum[x]sum[y])/(4*sum[x^{2}]-sum[x]^{2})
   [Equation 10] b1=(4*sum[xy]-sum[x]sum[y])/(4*sum[x^{2}]-sum[x]^{2})

   x0=xdata[id0]  y0=ydata[id0]
   x1=xdata[id1]  y0=ydata[id0]
   x2=xdata[id2]  y2=ydata[id2]
   x3=xdata[id3]  y3=ydata[id3]

   sx=(x0+x1+x2+x3)
   sy=(y0+y1+y2+y3)
   sxy=(x0y0+x1y1+x2y2+x3y3)
   sxs=(x0x0+x1x1+x2x2+x3x3)
   sx2=sx*sx
   denom=4*sxs-sx2
   b0=[sxs*sy-sx*sxy]/denom
   b1=[4*sxy-sx*sy]/denom 

   v0=y0-(b0+b1*x0)
   v1=y1-(b0+b1*x1)
   v2=y2-(b0+b1*x2)
   v3=y3-(b0+b1*x3)
   volatility_factor=v0*v0+v1*v1+v2*v2+v3*v3

   [Equation 11] Distance factor D(0,1,2,3)=(x1-x0)^{2}+(x2-x0)^{2}+(x3-x0)^{2}
   d0=x1-x0
   d1=x2-x0
   d2=x3-x0
   distance_factor=d0*d0+d1*d1+d2*d2 */
                        /* Volatility factor V(0,1,2,3) */
                        sx=xdata[id0]+xdata[id1]+xdata[id2]+xdata[id3]; /* sum[x] */
                        sy=ydata[id0].dat[part]+ydata[id1].dat[part]+ydata[id2].dat[part]+ydata[id3].dat[part]; /* sum[y] */ 
                        sxy=xdata[id0]*ydata[id0].dat[part]+xdata[id1]*ydata[id1].dat[part]+
                              xdata[id2]*ydata[id2].dat[part]+xdata[id3]*ydata[id3].dat[part]; /* sum[xy] */
                        sxs=xdata[id0]*xdata[id0]+xdata[id1]*xdata[id1]+xdata[id2]*xdata[id2]+xdata[id3]*xdata[id3]; /* sum[x^{2}] */
                        sx2=sx*sx; /* (sum[x])^{2} */
                        denom=4.0L*sxs-sx2; /* 4*sum[x^{2}]-(sum[x])^{2} */
                        b0=(sxs*sy-sx*sxy)/denom; /* b0=(sum[x^{2}]sum[y]-sum[x]sum[y])/(4*sum[x^{2}]-sum[x]^{2}) */
                        b1=(4.0L*sxy-sx*sy)/denom; /* b1=(4*sum[xy]-sum[x]sum[y])/(4*sum[x^{2}]-sum[x]^{2}) */
                        v0=ydata[id0].dat[part]-(b0+b1*xdata[id0]); /* v0=y0-(b0+b1*x0) */
                        v1=ydata[id1].dat[part]-(b0+b1*xdata[id1]); /* v1=y1-(b0+b1*x1) */
                        v2=ydata[id2].dat[part]-(b0+b1*xdata[id2]); /* v2=y2-(b0+b1*x2) */
                        v3=ydata[id3].dat[part]-(b0+b1*xdata[id3]); /* v3=y3-(b0+b1*x3) */
                        volatility_factor=v0*v0+v1*v1+v2*v2+v3*v3; /* V=sum[y-(b0+b1x)]^{2} */
                        /* Distance factor D(0,1,2,3) */
                        d0=xdata[id1]-xdata[id0]; /* d0=x1-x0 */
                        d1=xdata[id2]-xdata[id0]; /* d1=x2-x0 */
                        d2=xdata[id3]-xdata[id0]; /* d2=x3-x0 */
                        distance_factor=d0*d0+d1*d1+d2*d2; /* D(0,1,2,3)=(x1-x0)^{2}+(x2-x0)^{2}+(x3-x0)^{2} */
/* The method weights the 4 primary estimates using the reciprocal of the product
   of the volatility factor and the distance factor:

   P_{i-3} to P_{i}     y'_{imm}=F(i,i-3,i-2,i-1) corresponding weight [Equation 13] w_{imm}=1/[V(i,i-3,i-2,i-1)D(i,i-3,i-2,i-1)]
   P_{i-2} to P_{i+1}   y'_{im}=F(i,i-2,i-1,i+1)  corresponding weight [Equation 13] w_{im}=1/[V(i,i-2,i-1,i+1)D(i,i-2,i-1,i+1)]
   P_{i-1} to P_{i+2}   y'_{ip}=F(i,i-1,i+1,i+2)  corresponding weight [Equation 13] w_{ip}=1/[V(i,i-1,i+1,i+2)D(i,i-1,i+1,i+2)]
   P_{i} to P_{i+3}     y'_{ipp}=F(i,i+1,i+2,i+3) corresponding weight [Equation 13] w_{ipp}=1/[V(i,i+1,i+2,i+3)D(i,i+1,i+2,i+3)] */
/* Check whether the volatility factor is essentially zero */
                        if(volatility_factor>((ydata[id0].dat[part]*ydata[id0].dat[part]+
                                                ydata[id1].dat[part]*ydata[id1].dat[part]+
                                                ydata[id2].dat[part]*ydata[id2].dat[part]+
                                                ydata[id3].dat[part]*ydata[id3].dat[part])*LDBL_EPSILON)){
                           /* Finite weight */
                           weighting_factor=1.0L/(volatility_factor*distance_factor);
                           sum_primary_estimate_finite+=primary_estimate*weighting_factor;
                           sum_weighting_factor_finite+=weighting_factor;
                        }
                        else{
                           /* Infinite weight */
                           sum_primary_estimate_infinite+=primary_estimate;
                           sum_weighting_factor_infinite+=1.0L;
                        }
                     } /* End check 4 data points are within data range */
/* Calculate the final estimate of the first derivative of the interpolating function is

   [Equation 14] y'_{i}=[y'_{imm}w_{imm}+y'_{im}w_{im}+y'_{ip}w_{ip}+y'_{ipp}w_{ipp}]/[w_{imm}+w_{im}+w_{ip}+w_{ipp}] */
                     if(sum_weighting_factor_infinite<0.5L){ /* If no infinite weights exist */
                        final_estimate=sum_primary_estimate_finite/sum_weighting_factor_finite;
                     }
                     else{ /* If infinite weights exist */
                        final_estimate=sum_primary_estimate_infinite/sum_weighting_factor_infinite;
                     }
                  } /* End for primary estimates and calculation of the final estimate */
                  if(endpoint==0){
                     final_estimate0=final_estimate; /* Store final_estimate for endpoint=0 */
                  }
               } /* End for endpoints */

               if(degree==3){
/* Calculate the coefficients of the 3rd-degree polynomial
   The 3rd-degree polynomial for the y value in the interval x_{i} and x_{i+1} is

   [Equation 3] y=a0+a1(x-x_{i})+a2(x-x_{i})^{2}+a3(x-x_{i})^{3}

   [Equation 4] a0=y_{i}
   [Equation 4] a1=y'_{i}
   [Equation 4] a2=-[2(y'_{i}-m_{i})+(y'_{i+1}-m_{i})]/[x_{i+1}-x_{i}]
   [Equation 4] a3=[(y'_{i}-m_{i})+(y'_{i+1}-m_{i})]/[(x_{i+1}-x_{i})^{2}]

   where m_{i} is the slope of the line segment connecting P_{i} and P_{i+1}

   [Equation 5] m_{i}=(y_{i+1}-y_{i})/(x_{i+1}-x_{i}) */
                  denom=xdata[interval]-xdata[interval-1];
                  slope=(ydata[interval].dat[part]-ydata[interval-1].dat[part])/denom;
                  a0=ydata[interval-1].dat[part];
                  a1=final_estimate0;
                  a2=-((2.0L*(final_estimate0-slope)+(final_estimate-slope))/denom);
                  a3=((final_estimate0-slope)+(final_estimate-slope))/(denom*denom);
               }
               else{ /* nth-degree polynomial n>3 */
/* The nth-degree polynomial is represented by

   [Equation 19] v(u)=A0[u^{n}-u]+A1[(1-u)^{n}-(1-u)]
   [Equation 20] A0=[v'_{0}+(n-1)v'_{1}]/[n(n-2)]
   [Equation 20] A1=-[(n-1)v'_{0}+v'_{1}]/[n(n-2)]

   u=0,v=0,v'_{0}=(y'_{i}-m_{i})/(x_{i+1}-x_{i}) @ P_{i}
   u=1,v=0,v'_{1}=(y'_{i+1}-m_{i})/(x_{i+1}-x_{i}) @ P_{i}

   where m_{i} is the slope of the line segment connecting P_{i} and P_{i+1}

   [Equation 5] m_{i}=(y_{i+1}-y_{i})/(x_{i+1}-x_{i}) */
                  temp=xdata[interval]-xdata[interval-1];
                  slope=(ydata[interval].dat[part]-ydata[interval-1].dat[part])/temp;
                  v0=(final_estimate0-slope)*temp; /* v'_{0} */
                  v1=(final_estimate-slope)*temp; /* v'_{1} */
                  denom=(long double)(degree*(degree-2)); /* n(n-2) */
                  a0=(v0+(degree-1)*v1)/denom; /* A0 */
                  a1=-(((degree-1)*v0+v1)/denom); /* A1 */
               }
            } /* End if interval */

            if(degree==3){ /* 3rd-degree polynomial */
/* Evaluate [Equation 3] */
               temp=xinterp[i]-xdata[interval-1];
               yinterp[i].dat[part]=a0+temp*(a1+temp*(a2+a3*temp));
            }
            else{ /* nth-degree polynomial */
/* Evaluate [Equation 19] */
               u=(xinterp[i]-xdata[interval-1])/(xdata[interval]-xdata[interval-1]);
               temp=1.0L-u;
               v=a0*(power(u,degree)-u)+a1*(power(temp,degree)-temp);
               /* Evaluate [Equation 15] y-y_{i}=(y_{i+1}-y_{i})*u+v */
               yinterp[i].dat[part]=ydata[interval-1].dat[part]+(ydata[interval].dat[part]-ydata[interval-1].dat[part])*u+v;
            }
         } /* End for 3 subcases */
      } /* End cycle through the points for which interpolation is desired */
   } /* End general case - 5 points or more */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Special cases - 4 points or less */
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
   else if(size==2){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Subcase 1: 2 data points - Linear interpolation and extrapolation
   Know (x0,y0) and (x1,y1) then at given x:
   y=y0(1-a)+a*y1
   a=(x-x0)/(x1-x0) */
      for(i=0;i<points;i++){ /* Cycle through the points for which interpolation is desired */
         a=(xinterp[i]-xdata[0])/(xdata[1]-xdata[0]);
         yinterp[i].dat[part]=(1.0L-a)*ydata[0].dat[part]+a*ydata[1].dat[part];
      }
   }
   else if(size==3){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Subcase 2: 3 data points - Quadratic interpolation and linear extrapolation
   Quadratic interpolation - Know (x0,y0), (x1,y1) and (x2,y2) then at given x:
   y=y0*c(1-a)+y1*a(1-b)+y2*b(1-c)
   a=(x-x0)/(x1-x0)
   b=(x-x1)/(x2-x1)
   c=(x-x2)/(x0-x2) 

   Linear extrapolation
   Create a tangent line at the end of the known data and extend it beyond that limit
   y=y0+m(x-x0) where m is the slope i.e. first derivative */
      for(i=0;i<points;i++){ /* Cycle through the points for which interpolation is desired */
         a=(xinterp[i]-xdata[0])/(xdata[1]-xdata[0]);
         b=(xinterp[i]-xdata[1])/(xdata[2]-xdata[1]);
         c=(xinterp[i]-xdata[2])/(xdata[0]-xdata[2]);
         if(xinterp[i]<(xdata[0]+LDBL_EPSILON)){ /* Point <= Range minimum - Linear extrapolation */
/* first_derivative=F(0,1,2)=[f1 + f2]/[f3 + f4]

   f1=(y1-y0)*(x0-x2)^{2}
   f2=(y0-y2)*(x1-x0)^{2}
   f3=(x0-x2)^{2}(x1-x0)
   f4=(x1-x0)^{2}(x0-x2)

   f1=(ydata[1]-ydata[0])*x20*x20;
   f2=(ydata[0]-ydata[2])*x10*x10;
   f3=x20*x20*x10;
   f4=x10*x10*x20; */
            x10=xdata[1]-xdata[0];
            x20=xdata[0]-xdata[2];
            first_derivative=(((ydata[1].dat[part]-ydata[0].dat[part])*x20*x20)+
                              ((ydata[0].dat[part]-ydata[2].dat[part])*x10*x10))/(x20*x20*x10+x10*x10*x20);
            yinterp[i].dat[part]=ydata[0].dat[part]+first_derivative*(xinterp[i]-xdata[0]);
         }
         else if(xinterp[i]<(xdata[size-1]-LDBL_EPSILON)){ /* Range minimum < Point < Range maximum - Quadratic interpolation */
            yinterp[i].dat[part]=ydata[0].dat[part]*c*(1.0L-a)+ydata[1].dat[part]*a*(1.0L-b)+ydata[2].dat[part]*b*(1.0L-c);
         }
         else{ /* Point >= Range maximum - Linear extrapolation */
/* first_derivative=F(size-1,size-2,size-3)=F(0,1,2)=[f1 + f2]/[f3 + f4]

   f1=(y1-y0)*(x0-x2)^{2}
   f2=(y0-y2)*(x1-x0)^{2}
   f3=(x0-x2)^{2}(x1-x0)
   f4=(x1-x0)^{2}(x0-x2)

   f1=(ydata[size-2]-ydata[size-1])*x20*x20;
   f2=(ydata[size-1]-ydata[2])*x10*x10;
   f3=x20*x20*x10;
   f4=x10*x10*x20; */
            x10=xdata[size-2]-xdata[size-1];
            x20=xdata[size-1]-xdata[size-3];
            first_derivative=(((ydata[size-2].dat[part]-ydata[size-1].dat[part])*x20*x20)+
                              ((ydata[size-1].dat[part]-ydata[size-3].dat[part])*x10*x10))/(x20*x20*x10+x10*x10*x20);
            yinterp[i].dat[part]=ydata[size-1].dat[part]+first_derivative*(xinterp[i]-xdata[size-1]);
         }
      }
   }
   else if(size==4){
/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Subcase 3: 4 data points - Cubic interpolation and linear extrapolation
   Cubic interpolation - Know (x0,y0), (x1,y1), (x2,y2) and (x3,y3) then at given x:
   y=y0*ce(1-a)+y1*ad(1-b)+y2*bf(1-c)+y3*(1-d)(1-e)(1-f)
   a=(x-x0)/(x1-x0)
   b=(x-x1)/(x2-x1)
   c=(x-x2)/(x0-x2)
   d=(x-x3)/(x1-x3)
   e=(x-x3)/(x0-x3)
   f=(x-x3)/(x2-x3)

   Linear extrapolation
   Create a tangent line at the end of the known data and extend it beyond that limit
   y=y0+m(x-x0) where m is the slope i.e. first derivative */
      for(i=0;i<points;i++){ /* Cycle through the points for which interpolation is desired */
         a=(xinterp[i]-xdata[0])/(xdata[1]-xdata[0]);
         b=(xinterp[i]-xdata[1])/(xdata[2]-xdata[1]);
         c=(xinterp[i]-xdata[2])/(xdata[0]-xdata[2]);
         d=(xinterp[i]-xdata[3])/(xdata[1]-xdata[3]);
         e=(xinterp[i]-xdata[3])/(xdata[0]-xdata[3]);
         f=(xinterp[i]-xdata[3])/(xdata[2]-xdata[3]);
         if(xinterp[i]<(xdata[0]+LDBL_EPSILON)){ /* Point <= Range minimum - Linear extrapolation */
/* first_derivative=F(0,1,2,3)=[f1 + f2 + f3]/[f4]

   f1=(y1-y0)*(x2-x0)^{2}*(x3-x0)^{2}*(x3-x2)
   f2=(y2-y0)*(x3-x0)^{2}*(x1-x0)^{2}*(x1-x3)
   f3=(y3-y0)*(x1-x0)^{2}*(x2-x0)^{2}*(x2-x1)
   f4=(x1-x0)*(x2-x0)*(x3-x0)*(x2-x1)*(x3-x2)*(x3-x1)

   f1=(ydata[1]-ydata[0])*x20*x20*x30*x30*x32;
   f2=(ydata[2]-ydata[0])*x30*x30*x10*x10*x13;
   f3=(ydata[3]-ydata[0])*x10*x10*x20*x20*x21;
   f4=x10*x20*x30*x21*x32*-x13; */
            x10=xdata[1]-xdata[0];x13=xdata[1]-xdata[3];
            x20=xdata[2]-xdata[0];x21=xdata[2]-xdata[1];
            x30=xdata[3]-xdata[0];x32=xdata[3]-xdata[2];
            first_derivative=(((ydata[1].dat[part]-ydata[0].dat[part])*x20*x20*x30*x30*x32)+
                              ((ydata[2].dat[part]-ydata[0].dat[part])*x30*x30*x10*x10*x13)+
                              ((ydata[3].dat[part]-ydata[0].dat[part])*x10*x10*x20*x20*x21))/(x10*x20*x30*x21*x32*-x13);
            yinterp[i].dat[part]=ydata[0].dat[part]+first_derivative*(xinterp[i]-xdata[0]);
         }
         else if(xinterp[i]<(xdata[size-1]-LDBL_EPSILON)){ /* Range minimum < Point < Range maximum - Cubic interpolation */
            yinterp[i].dat[part]=ydata[0].dat[part]*c*e*(1.0L-a)+ydata[1].dat[part]*a*d*(1.0L-b)+
                        ydata[2].dat[part]*b*f*(1.0L-c)+ydata[3].dat[part]*(1.0L-d)*(1.0L-e)*(1.0L-f);
         }
         else{ /* Point >= Range maximum - Linear extrapolation */
/* first_derivative=F(size-1,size-2,size-3,size-4)=F(0,1,2,3)=[f1 + f2 + f3]/[f4]

   f1=(y1-y0)*(x2-x0)^{2}*(x3-x0)^{2}*(x3-x2)
   f2=(y2-y0)*(x3-x0)^{2}*(x1-x0)^{2}*(x1-x3)
   f3=(y3-y0)*(x1-x0)^{2}*(x2-x0)^{2}*(x2-x1)
   f4=(x1-x0)*(x2-x0)*(x3-x0)*(x2-x1)*(x3-x2)*(x3-x1)

   f1=(ydata[size-2]-ydata[size-1])*x20*x20*x30*x30*x32;
   f2=(ydata[size-3]-ydata[size-1])*x30*x30*x10*x10*x13;
   f3=(ydata[size-4]-ydata[size-1])*x10*x10*x20*x20*x21;
   f4=x10*x20*x30*x21*x32*-x13; */
            x10=xdata[size-2]-xdata[size-1];x13=xdata[size-2]-xdata[size-4];
            x20=xdata[size-3]-xdata[size-1];x21=xdata[size-3]-xdata[size-2];
            x30=xdata[size-4]-xdata[size-1];x32=xdata[size-4]-xdata[size-3];
            first_derivative=(((ydata[size-2].dat[part]-ydata[size-1].dat[part])*x20*x20*x30*x30*x32)+
                              ((ydata[size-3].dat[part]-ydata[size-1].dat[part])*x30*x30*x10*x10*x13)+
                              ((ydata[size-4].dat[part]-ydata[size-1].dat[part])*x10*x10*x20*x20*x21))/(x10*x20*x30*x21*x32*-x13);
            yinterp[i].dat[part]=ydata[size-1].dat[part]+first_derivative*(xinterp[i]-xdata[size-1]);
         }
      }
   }
}

long double power(long double x,int n){
/* Instead of using pow fn from math.h
   Returns x^{n} */

   int k;
   long double operand=1.0L;

   for(k=0;k<n;k++){
      operand*=x;
   }

   return operand;
}
