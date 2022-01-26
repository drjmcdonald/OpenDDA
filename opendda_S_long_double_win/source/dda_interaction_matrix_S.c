/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Computes the DDA equation i.e. calculates the 1st column of the
   interaction matrix with the diagonal equal to zero to preserve
   Toeplitz structure

   tc: Controls the tensor component 0,1,2,3,4,5=xx,xy,xz,yy,yz,zz
   a,b,i,j,k,index0*,
   array0: Loop and array index control variables
   interaction_matrix[6][Ka*Ja*Pa]: Stores the 6 indepedent interaction matrix tensor components
   r: r_{i)=(1,1,1),r_{j}=(i,j,k}
      r=|r_{ij}|=|r_{i}-r_{j}|=|(x[0],x[1],x[2])|
   r2: |r_{ij}|^{2}=x[0]^{2}+x[1]^{2}+x[2]^{2}
   r3: |r_{ij}|^{3}
   ikr: ikr_{ij}, k is the wavenumber
   eikr: e^{ikr_{ij}}, k is the wavenumber

   interaction_matrix_{ij}=[interaction_matrix[0]   interaction_matrix[1]   interaction_matrix[2]]   
                           [interaction_matrix[1]   interaction_matrix[3]   interaction_matrix[4]]
                           [interaction_matrix[2]   interaction_matrix[4]   interaction_matrix[5]]

   interaction_matrix_{ij}=[e^{ikr_{ij}}/r_{ij}^{3}]*[((ikr_{ij}-1)/r_{ij}^{2})(3r_{ij}*r_{ij}-
                            r_{ij}^{2}I)+k^{2}(r_{ij}*r_{ij}-r_{ij}^{2}I)]

   where

   r_{ij}*r_{ij}=[r_{x}^{2}   r_{x}r_{y}   r_{x}r_{z}] i.e. dyadic product   
                 [r_{y}r_{x}  r_{y}^{2}    r_{y}r_{z}]
                 [r_{z}r_{x}  r_{z}r_{y}   r_{z}^{2} ]

   interaction_matrix_{ij}=frac0*[frac1*index1+index2]

   frac0: [e^{ikr}/r^{3}]
   frac1: ((ikr-1)/r^{2})
   index1: (3r*r-r^{2}) for diagonal entries
           (3r*r) for off-diagonal entries
   index2: k^{2}(r*r-r^{2}) for diagonal entries
           k^{2}(r*r) for off-diagonal entries

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dda_interaction_matrix_S.h"

void dda_interaction_matrix_1st_column_zero_diagonal(void){

   int a,b,i,j,k,tc;
   long int array0,index0k,index0j,index0i;
   long double r,r2,r3,k2,x0[3],index1,index2;
   ldcomplex frac0,frac1,ikr,eikr;

   k2=wavenumber*wavenumber; /* k^{2} */

   /* Only need to compute if 1st time or wavenumber has changed */
   if((wavenumber_previous<0.0L)||fabsl(wavenumber_previous-wavenumber)>LDBL_EPSILON){
      wavenumber_previous=wavenumber; /* Store previous wavenumber for comparison */
      dft_interaction_matrix_flag=0; /* Every time the interaction matrix is changed the
      DFT of the interaction matrix must be recalculated */
      for(k=0;k<target.P;k++){
         index0k=k*target.Ma;
         x0[2]=((long double)k)*target.dipole_spacing; /* z */
         for(j=0;j<target.J;j++){
            index0j=index0k+j*target.Ka;
            x0[1]=((long double)j)*target.dipole_spacing; /* y */
            for(i=0;i<target.K;i++){
               if(i==0&&j==0&&k==0){continue;}
               index0i=index0j+i;
               x0[0]=((long double)i)*target.dipole_spacing; /* x */
               r2=(x0[0]*x0[0]+x0[1]*x0[1]+x0[2]*x0[2]); /* r^2 */
               r=sqrtl(r2); /* r */
               r3=r2*r; /* r^3 */
               ikr.dat[0]=0.0L;ikr.dat[1]=wavenumber*r; /* ikr=k*r*i */
               /* Complex exponential where real part of 'in' is zero
                  in=a+ib=0.0+ib
                  out=exp(in)=exp(a)*[cos(b)+isin(b)]=exp(0.0)*[cos(b)+isin(b)]
                       =1.0*[cos(b)+isin(b)]=cos(b)+isin(b) */
               eikr.dat[0]=cosl(ikr.dat[1]);eikr.dat[1]=sinl(ikr.dat[1]); /* e^{ikr} */
               frac0.dat[0]=eikr.dat[0]*(1.0L/r3);frac0.dat[1]=eikr.dat[1]*(1.0L/r3); /* frac0=eikr/r3 */
               frac1.dat[0]=(-1.0L/r2);frac1.dat[1]=(ikr.dat[1])*(1.0L/r2); /* frac1=(ikr-1)/r2 */
               tc=0;
               for(a=0;a<3;a++){
                  for(b=a;b<3;b++){
                     array0=tc*target.Na+index0i;
                     if(a==b){ /* Diagonal entries */
                        index1=3.0L*x0[a]*x0[b]-r2; /* (3r*r-r^{2}) */
                        index2=k2*(x0[a]*x0[b]-r2); /* k^{2}(r*r-r^{2}) */
                        interaction_matrix[array0].dat[0]=frac0.dat[0]*index2+(frac0.dat[0]*(frac1.dat[0]*index1)-
                                             frac0.dat[1]*(frac1.dat[1]*index1)); /* interaction_matrix_{ii}=frac0*[frac1*index1+index2] Real component */
                        interaction_matrix[array0].dat[1]=frac0.dat[1]*index2+(frac0.dat[0]*(frac1.dat[1]*index1)+
                                             frac0.dat[1]*(frac1.dat[0]*index1)); /* interaction_matrix_{ii}=frac0*[frac1*index1+index2] Imaginary component */
                     }
                     else{ /* Off-diagonal entries */
                        index1=3.0L*x0[a]*x0[b]; /* For the off-diagonal entries (r^{2}-3r*r) becomes (3r*r) */
                        index2=k2*(x0[a]*x0[b]); /* For the off-diagonal entries k^{2}(r*r-r^{2}) becomes k^{2}(r*r)  */
                        interaction_matrix[array0].dat[0]=frac0.dat[0]*index2+(frac0.dat[0]*(frac1.dat[0]*index1)-
                                             frac0.dat[1]*(frac1.dat[1]*index1)); /* interaction_matrix_{ij}=frac0*[frac1*index1+index2] Real component */
                        interaction_matrix[array0].dat[1]=frac0.dat[1]*index2+(frac0.dat[0]*(frac1.dat[1]*index1)+
                                             frac0.dat[1]*(frac1.dat[0]*index1)); /* interaction_matrix_{ij}=frac0*[frac1*index1+index2] Imaginary component */
                     }
                     tc++;
                  }
               }
            }
         }
      }

      for(tc=0;tc<6;tc++){ /* Zero the zero displacement elements */
         interaction_matrix[tc*target.Na]=zero;
      }
   }
}
