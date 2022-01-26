/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Multiplies the interaction matrix by the given vector using DFTs
   Since the fundamental elements in the symmetric interaction matrix are tensors,
   there are 6 independent components arrays interaction_matrix[6] and three vector
   components arrays vector[3]

   xsign,ysign,zsign: In order to facilitate the use of DFTs the Toeplitz blocks
                      must be converted to circulant blocks. The dyadic products
                      in the DDA equation
                       _                              _
                      | (rx)^{2}   (rx)(ry)   (rx)(rz) |   
                      | (ry)(rx)   (ry)^{2}   (ry)(rz) |
                      |_(rz)(rx)   (rz)(ry)   (rz)^{2}_|
   
                      contain vector sign information. For the diagonal elements,
                      this sign information is irrelevant as the terms are squared,
                      however, for the off-diagonal elements, the sign information
                      creates antisymmetries, i.e. for the off-diagonal component
                      (rx)(ry) the dipole level (dipoles in x-direction lines
                      interacting, x direction) and the line level (lines in each 
                      plane interacting, y direction) are antisymmetric
                      [interaction_matrix=-interaction_matrix^{T}],
                      whereas the plane level (xy planes interacting, z direction)
                      is symmetric 
                        
                      ANTISYMMETRIC TOEPLITZ          CIRCULANT
                         _            _        _                      _
                        |  0   -B   +C |      |  0   -B   +C   -C   +B |
                        | +B    0   -B | -->  | +B    0   -B   +C   -C |
                        |_-C   +B    0_|      | -C   +B    0   -B   +C |
                                              | +C   -C   +B    0   -B |
                                              |_-B   +C   -C   +B    0_|
                        i.e. mirror the 1st col excluding 1st entry and multiply by -1
                             at same time to create circulant from antisymmetric Toeplitz

                        COMPONENT     x     y     z
                           xx         -1    +1    +1
                           xy         -1    -1    +1        
                           xz         -1    +1    -1
                           yy         +1    -1    +1
                           yz         +1    -1    -1
                           zz         +1    +1    -1

   isign,jsign,ksign: Since only the [0->K-1][0->J-1][0->P-1] quadrant of the six tensor
                      components is stored explicitly, in order to extract the pertinent
                      data for the Fourier Domain multiplications, the sign disparities
                      between the stored quadrant and the remaining 7 quadrants for the
                      off-diagonal components, due to the antisymmetries in the off-diagonal
                      components, can be taken into account using

                             _
                      isign=| +1 if 0<=i<Ka
                            |_-1 if Ka<=i<Kpp

                             _
                      jsign=| +1 if 0<=j<Ja
                            |_-1 if Ja<=j<Jpp

                             _
                      ksign=| +1 if 0<=k<Pa
                            |_-1 if Pa<=k<Ppp

   tc,vc: tc controls the tensor component 0,1,2,3,4,5=xx,xy,xz,yy,yz,zz
          vc controls the vector component 0,1,2=x,y,z

   grid_points: FFTW does not normalise the transforms and grid_points is used for DFT normalisation

   i,j,k,ii,jj,kk,index0,index1,element: Loop and array index control variables

   xdft[Kpp],
   ydft[Jpp],
   zdft[Ppp]: Scratch arrays used for the DFT of the six tensor components

   dft_interaction_matrix_flag: Flag for the DFT of the 6 independent tensor components.
   interaction_matrix[6][Ka*Ja*Pa]: Stores the 6 indepedent interaction matrix tensor components
   vector[3][K*J*P]: Stores the 3 vector components of the vector to be multiplied by interaction_matrix
   zero_padded_vector[3][K*Jpp*P]: Stores the zero-padded vector components for the multiplication
                           and acts as a scratch space for the FD multiplications

BASIC FLOWCHART:
   
   1. Zero-pad the 3 vector components
   2. DFT of the 6 independent tensor components
   3. DFT of the 3 vector components
   4. FD element-wise tensor-vector multiplication
   5. iDFT of the resultant
   6. Normalisation
   7. Take diagonal contributions, that were ignored to obtain the requisite Toeplitz structure
      to permit the use of DFTs, into account

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dftmatvec_S.h"

void dftmatvec_S(fcomplex *vector){

   int i,j,k,ii,jj,kk,tc,vc;
   long int index0k,index0j,index0i,index1k,index1j,index1i;
   long int array0,array1,array2,array3,element;
   float isign,jsign,ksign,grid_points;
   float xsign[6]={1.0F,-1.0F,-1.0F,1.0F,1.0F,1.0F};
   float ysign[6]={1.0F,-1.0F,1.0F,1.0F,-1.0F,1.0F};
   float zsign[6]={1.0F,1.0F,-1.0F,1.0F,-1.0F,1.0F};
   size_t size0,size1,size2;
   double start,finish;
   fcomplex m0,m1,m2,m3,m4,m5,v0,v1,v2;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
   if(timing.enabled){ /* Increment the number of DFT multiplications and start timing */
      timing.dft_engine_count++;
      start=walltime();
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Set appropriate section of zero_padded_vector to zero */
   size0=(target.Jpp-target.J)*target.K*sizeof(fcomplex); /* Zero-padding size */
   for(vc=0;vc<3;vc++){
      array0=vc*target.Nv+target.M;
      for(k=0;k<target.P;k++){
         memset(&zero_padded_vector[array0+k*target.Mv],0,size0);
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Embed vector into zero_padded_vector */
   for(vc=0;vc<3;vc++){
      element=0;
      array0=vc*target.Nd;
      array1=vc*target.Nv;
      for(k=0;k<target.P;k++){
         index0k=k*target.M;index1k=k*target.Mv+array1;
         for(j=0;j<target.J;j++){
            index0j=index0k+j*target.K;index1j=index1k+j*target.K;
            for(i=0;i<target.K;i++){
               index0i=index0j+i;
               if(target.occupied[(int)((double)index0i/32.0)]&(1<<(index0i%32))){ /* Lattice site is occupied */
                  zero_padded_vector[index1j+i]=vector[array0+element];
                  element++;
               }
               else{
                  zero_padded_vector[index1j+i]=zero;
               }
            }
         }
      }
   }
   if(timing.enabled){ /* Finish timing and add contribution */
      finish=walltime();
      timing.zero_padded_vector+=(finish-start);
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Perform the DFT of the 6 independent tensor components of the interaction matrix */
   if(dft_interaction_matrix_flag==0){
      dft_interaction_matrix_flag=1;
      
      if(timing.enabled){ /* Increment number of tensor components DFTs and start timing */
         timing.dft_tensor_components_count++;
         start=walltime();
      }
      size0=target.kf*sizeof(fcomplex); /* Zero-padding size */
      for(tc=0;tc<6;tc++){ /* xDFT */
         array0=tc*target.Na;
         for(k=0;k<target.P;k++){ /* Extract the J*P sections into xdft one at a time, make them circulant, DFT
                                     and place the first K element of the result back into interaction_matrix */
            index0k=k*target.Ma+array0;
            for(j=0;j<target.J;j++){
               index0j=index0k+j*target.Ka;
               for(i=0;i<target.K;i++){ 
                  index0i=index0j+i;
                  xdft[i]=interaction_matrix[index0i];
                  if(i>0){ /* Make circulant */
                     xdft[target.Kpp-i]=fcomplex_scale(interaction_matrix[index0i],xsign[tc]); 
                  }
               }
               memset(&xdft[target.K],0,size0); /* Set the zero-fill K to K+kf-1 to zero */
               fftwf_execute(plan_tcxdft); /* Execute xDFT */
               for(i=0;i<target.Ka;i++){ /* Copy 1st Ka elements of xdft back into interaction_matrix */
                  interaction_matrix[index0j+i]=xdft[i];
               }
            }
         }
      }

      size0=target.jf*sizeof(fcomplex); /* Zero-padding size */
      for(tc=0;tc<6;tc++){ /* yDFT */
         array0=tc*target.Na;
         for(k=0;k<target.P;k++){ /* Extract the Ka*P sections into ydft one at a time, make them circulant, DFT
                                     and place the first J element of the result back into interaction_matrix */
            index0k=k*target.Ma+array0;
            for(i=0;i<target.Ka;i++){
               index0i=index0k+i;
               for(j=0;j<target.J;j++){
                  index0j=index0i+j*target.Ka;
                  ydft[j]=interaction_matrix[index0j];
                  if(j>0){ /* Make circulant */
                     ydft[target.Jpp-j]=fcomplex_scale(interaction_matrix[index0j],ysign[tc]); 
                  }
               }
               memset(&ydft[target.J],0,size0); /* Set the zero-fill J to J+jf-1 to zero */
               fftwf_execute(plan_tcydft); /* Execute yDFT */
               for(j=0;j<target.Ja;j++){ /* Copy 1st Ja elements of ydft back into interaction_matrix */
                  interaction_matrix[index0i+j*target.Ka]=ydft[j];
               }
            }
         }
      }

      size0=target.pf*sizeof(fcomplex); /* Zero-padding size */
      for(tc=0;tc<6;tc++){ /* zDFT */
         array0=tc*target.Na;
         for(j=0;j<target.Ja;j++){ /* Extract the Ja*P sections into zdft one at a time, make them circulant, DFT
                                      and place the first P element of the result back into interaction_matrix */
            index0j=j*target.Ka+array0;
            for(i=0;i<target.Ka;i++){
               index0i=index0j+i;
               for(k=0;k<target.P;k++){
                  index0k=index0i+k*target.Ma;
                  zdft[k]=interaction_matrix[index0k];
                  if(k>0){ /* Make circulant */
                     zdft[target.Ppp-k]=fcomplex_scale(interaction_matrix[index0k],zsign[tc]);
                  }
               }
               memset(&zdft[target.P],0,size0); /* Set the zero-fill P to P+pf-1 to zero */
               fftwf_execute(plan_tczdft); /* Execute zDFT */
               for(k=0;k<target.Pa;k++){ /* Copy 1st Pa elements of zdft back into interaction_matrix */
                  interaction_matrix[index0i+k*target.Ma]=zdft[k];
               }
            }
         }
      }
      if(timing.enabled){ /* Finish timing and add contribution */
         finish=walltime();
         timing.dft_tensor_components+=(finish-start);
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Perform the DFT of the 3 vector components of zero_padded_vector, the
   subsequent Fourier domain multplication with the 6 independent tensor
   components and the iDFT of the resultant */
   if(timing.enabled){start=walltime();} /* Start timing */
   for(vc=0;vc<3;vc++){ /* yDFT */
      array0=vc*target.Nv;
      for(k=0;k<target.P;k++){
         index0k=k*target.Mv+array0;
         fftwf_execute_dft(plan_vcydft,(fftwf_complex *)(&zero_padded_vector[index0k]),(fftwf_complex *)(&zero_padded_vector[index0k]));
      }
   }

   array0=target.Kpp*target.P;
   size0=target.K*sizeof(fcomplex); /* Copy size */
   size1=target.Kpp*(target.Ppp-target.P)*sizeof(fcomplex); /* Zero-padding size */
   size2=(target.Kpp-target.K)*sizeof(fcomplex); /* Zero-padding size */
   for(j=0;j<target.Jpp;j++){
      jj=j>=target.Ja?(target.Jpp-j):j; /* Extract pertinent data from the interaction_matrix */
      jsign=j>=target.Ja?-1.0F:1.0F; /* Sign information */
      index0j=j*target.K;index1j=jj*target.Ka;
      for(vc=0;vc<3;vc++){  /* Extract the 3 component sections into xzscratch and zero-pad them */
         array1=vc*target.Nv+index0j;
         array2=vc*target.Qpp;
         for(k=0;k<target.P;k++){
            memcpy(&xzscratch[array2+k*target.Kpp],&zero_padded_vector[array1+k*target.Mv],size0); /* Copy data */
            memset(&xzscratch[array2+k*target.Kpp+target.K],0,size2); /* Zero-padding */
         }
         memset(&xzscratch[array2+array0],0,size1); /* Zero-padding */
      } /* End vc loop */
      for(vc=0;vc<3;vc++){ /* xDFT */
         fftwf_execute_dft(plan_vcxdft,(fftwf_complex *)(&xzscratch[vc*target.Qpp]),(fftwf_complex *)(&xzscratch[vc*target.Qpp]));
      }
      for(vc=0;vc<3;vc++){ /* zDFT */
         fftwf_execute_dft(plan_vczdft,(fftwf_complex *)(&xzscratch[vc*target.Qpp]),(fftwf_complex *)(&xzscratch[vc*target.Qpp]));
      }
      /* Fourier Domain multiplication */
      for(k=0;k<target.Ppp;k++){
         kk=k>=target.Pa?(target.Ppp-k):k; /* Extract pertinent data from the interaction_matrix */
         ksign=k>=target.Pa?-1.0F:1.0F; /* Sign information */
         index0k=k*target.Kpp;index1k=index1j+kk*target.Ma;
         for(i=0;i<target.Kpp;i++){
            ii=i>=target.Ka?(target.Kpp-i):i; /* Extract pertinent data from the interaction_matrix */
            isign=i>=target.Ka?-1.0F:1.0F; /* Sign information */
            index0i=index0k+i;index1i=index1k+ii;
            array1=target.Qpp+index0i;
            array2=target.Qpp+array1; /* 2*target.Qpp+index0i */
            m0=interaction_matrix[index1i];
            m1=interaction_matrix[target.Na+index1i];
            m2=interaction_matrix[2*target.Na+index1i];
            m3=interaction_matrix[3*target.Na+index1i];
            m4=interaction_matrix[4*target.Na+index1i];
            m5=interaction_matrix[5*target.Na+index1i];
            v0=xzscratch[index0i];
            v1=xzscratch[array1];
            v2=xzscratch[array2];
            xzscratch[index0i].dat[0]=(m0.dat[0]*v0.dat[0]-m0.dat[1]*v0.dat[1])+/*[0]*/
                           (m1.dat[0]*v1.dat[0]-m1.dat[1]*v1.dat[1])*isign*jsign+/*[1]isign*jsign*/
                           (m2.dat[0]*v2.dat[0]-m2.dat[1]*v2.dat[1])*isign*ksign;/*[2]isign*ksign*/
            xzscratch[index0i].dat[1]=(m0.dat[0]*v0.dat[1]+m0.dat[1]*v0.dat[0])+/*[0]*/
                           (m1.dat[0]*v1.dat[1]+m1.dat[1]*v1.dat[0])*isign*jsign+/*[1]isign*jsign*/
                           (m2.dat[0]*v2.dat[1]+m2.dat[1]*v2.dat[0])*isign*ksign;/*[2]isign*ksign*/
            xzscratch[array1].dat[0]=(m1.dat[0]*v0.dat[0]-m1.dat[1]*v0.dat[1])*jsign*isign+/*[1]jsign*isign*/
                           (m3.dat[0]*v1.dat[0]-m3.dat[1]*v1.dat[1])+/*[3]*/
                           (m4.dat[0]*v2.dat[0]-m4.dat[1]*v2.dat[1])*jsign*ksign;/*[4]jsign*ksign*/
            xzscratch[array1].dat[1]=(m1.dat[0]*v0.dat[1]+m1.dat[1]*v0.dat[0])*jsign*isign+/*[1]jsign*isign*/
                           (m3.dat[0]*v1.dat[1]+m3.dat[1]*v1.dat[0])+/*[3]*/
                           (m4.dat[0]*v2.dat[1]+m4.dat[1]*v2.dat[0])*jsign*ksign;/*[4]jsign*ksign*/
            xzscratch[array2].dat[0]=(m2.dat[0]*v0.dat[0]-m2.dat[1]*v0.dat[1])*ksign*isign+/*[2]ksign*isign*/
                           (m4.dat[0]*v1.dat[0]-m4.dat[1]*v1.dat[1])*ksign*jsign+/*[4]ksign*jsign*/
                           (m5.dat[0]*v2.dat[0]-m5.dat[1]*v2.dat[1]);/*[5]*/
            xzscratch[array2].dat[1]=(m2.dat[0]*v0.dat[1]+m2.dat[1]*v0.dat[0])*ksign*isign+/*[2]ksign*isign*/
                           (m4.dat[0]*v1.dat[1]+m4.dat[1]*v1.dat[0])*ksign*jsign+/*[4]ksign*jsign*/
                           (m5.dat[0]*v2.dat[1]+m5.dat[1]*v2.dat[0]);/*[5]*/
         } /* End i loop */
      } /* End k loop */
      for(vc=0;vc<3;vc++){ /* ziDFT */
         fftwf_execute_dft(plan_vczdfti,(fftwf_complex *)(&xzscratch[vc*target.Qpp]),(fftwf_complex *)(&xzscratch[vc*target.Qpp]));
      }
      for(vc=0;vc<3;vc++){ /* xiDFT */
         fftwf_execute_dft(plan_vcxdfti,(fftwf_complex *)(&xzscratch[vc*target.Qpp]),(fftwf_complex *)(&xzscratch[vc*target.Qpp]));
      }
      for(vc=0;vc<3;vc++){ /* Copy 1st P elements back into the zero_padded_vector */
         array1=vc*target.Nv+index0j;
         array2=vc*target.Qpp;
         for(k=0;k<target.P;k++){
            memcpy(&zero_padded_vector[array1+k*target.Mv],&xzscratch[array2+k*target.Kpp],size0); /* Copy data */
         }
      } /* End vc loop */
   } /* End j loop */

   for(vc=0;vc<3;vc++){ /* yiDFT */
      array0=vc*target.Nv;
      for(k=0;k<target.P;k++){
         index0k=k*target.Mv+array0;
         fftwf_execute_dft(plan_vcydfti,(fftwf_complex *)(&zero_padded_vector[index0k]),(fftwf_complex *)(&zero_padded_vector[index0k]));
      }
   }
   if(timing.enabled){ /* Finish timing and add contribution */
      finish=walltime();
      timing.zdft_zidft_and_fd_mult+=(finish-start);
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Normalisation + Take diagonal contributions into account + Compress data */
   if(timing.enabled){start=walltime();} /* Start timing */
   grid_points=1.0F/(float)target.Npp; /* Normalisation value */
   for(vc=0;vc<3;vc++){
      element=0;
      array0=vc*target.Nd;
      array1=vc*target.Nv;
      for(k=0;k<target.P;k++){
         index0k=k*target.M;index1k=k*target.Mv+array1;
         for(j=0;j<target.J;j++){
            index0j=index0k+j*target.K;index1j=index1k+j*target.K;
            for(i=0;i<target.K;i++){
               index0i=index0j+i;index1i=index1j+i;
               if(target.occupied[(int)((double)index0i/32.0)]&(1<<(index0i%32))){ /* Lattice site is occupied */
                  /*interaction_matrix_diagonal(vc)*x(vc)*/
                  array2=array0+element;
                  array3=array1+element; /* Index for compression */
                  m0=interaction_matrix_diagonal[array2];
                  v0=vector[array2];
                  zero_padded_vector[array3].dat[0]=(zero_padded_vector[index1i].dat[0]*grid_points)+(m0.dat[0]*v0.dat[0]-m0.dat[1]*v0.dat[1]);
                  zero_padded_vector[array3].dat[1]=(zero_padded_vector[index1i].dat[1]*grid_points)+(m0.dat[0]*v0.dat[1]+m0.dat[1]*v0.dat[0]);
                  element++;
               }
            }
         }
      }
   }
   if(timing.enabled){ /* Finish timing and add contribution */
      finish=walltime();
      timing.norm_diag_compress+=(finish-start);
   }
}
