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

   isign,ksign: Since only the [0->K-1][0->J-1][0->P-1] quadrant of the six tensor
                      components is stored expficitly, in order to extract the pertinent
                      data for the Fourier Domain multiplications, the sign disparities
                      between the stored quadrant and the remaining 7 quadrants for the
                      off-diagonal components, due to the antisymmetries in the off-diagonal
                      components, can be taken into account using

                             _
                      isign=| +1 if 0<=i<Ka
                            |_-1 if Ka<=i<Kpp

                             _
                      ksign=| +1 if 0<=k<Pa
                            |_-1 if Pa<=k<Ppp

   tc,vc: tc controls the tensor component 0,1,2,3,4,5=xx,xy,xz,yy,yz,zz
          vc controls the vector component 0,1,2=x,y,z

   grid_points: FFTW does not normalise the transforms and grid_points is used for DFT normalisation

   i,j,k,ii,kk,index0*,index1*: Loop and array index control variables

   xdft[Kpp],
   Note full y direction in stored in the MPI case so a scratch array is not required,
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
#include "dftmatvec_MPI.h"

void dftmatvec_MPI(fcomplex *vector){

   int i,j,k,ii,kk,tc,vc;
   long int index0i,index0j,index0k,index1k,index1j,index1i;
   long int array0,array1,array2,array3,element;
   float isign,ksign,grid_points;
   float xsign[6]={1.0F,-1.0F,-1.0F,1.0F,1.0F,1.0F};
   float ysign[6]={1.0F,-1.0F,1.0F,1.0F,-1.0F,1.0F};
   float zsign[6]={1.0F,1.0F,-1.0F,1.0F,-1.0F,1.0F};
   double start,finish;
   size_t size0,size1,size2;
   fcomplex m0,m1,m2,m3,m4,m5,v0,v1,v2;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
   if(timing.enabled){ /* Increment the number of DFT multiplications and start timing */
      timing.dft_engine_count++;
      start=MPI_Wtime();
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Set appropriate section of zero_padded_vector to zero */
   size0=(target.Jpp-target.J)*target.K*sizeof(fcomplex); /* Zero-padding size */
   for(vc=0;vc<3;vc++){
      array0=vc*parallel.alloc_padded+target.M;
      for(k=0;k<target.P;k++){
         if(k%parallel.np==parallel.myrank){
            parallel.plane=((int)((double)k/parallel.dnp));
            memset(&zero_padded_vector[array0+parallel.plane*target.Mv],0,size0);
         }
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Embed vector into zero_padded_vector */
   for(vc=0;vc<3;vc++){
      array0=vc*parallel.alloc_vector;
      array1=vc*parallel.alloc_padded;
      for(k=0;k<target.P;k++){
         if(k%parallel.np==parallel.myrank){ /* Plane k is on proc (k%np) */
            parallel.plane=((int)((double)k/parallel.dnp)); /* Plane k is plane ((int)(k/np)) on that proc */
            index0k=k*target.M;index1k=parallel.plane*target.Mv+array1;
            element=target.populated[parallel.plane]; /* Start index for plane k for the stored non-zero lattice sites */
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
   }
   if(timing.enabled){ /* Finish timing and add contribution */
      finish=MPI_Wtime();
      timing.zero_padded_vector+=(finish-start);
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Perform the DFT of the 6 independent tensor components of the interaction matrix */
   if(dft_interaction_matrix_flag==0){
      dft_interaction_matrix_flag=1;
      
      if(timing.enabled){ /* Increment number of tensor components DFTs and start timing */
         timing.dft_tensor_components_count++;
         start=MPI_Wtime();
      }
      size0=target.kf*sizeof(fcomplex); /* Zero-padding size */
      for(tc=0;tc<6;tc++){ /* xDFT */
         array0=tc*parallel.alloc_tensor;
         for(k=0;k<target.P;k++){ /* Extract the J*P sections into xdft one at a time, make them circulant, DFT
                                     and place the first K element of the result back into interaction_matrix */
            if(k%parallel.np==parallel.myrank){ /* Plane k is on proc (k%np) */
               parallel.plane=((int)((double)k/parallel.dnp)); /* Plane k is plane ((int)(k/np)) on that proc */
               index0k=parallel.plane*target.KaJpp+array0;
               for(j=0;j<target.J;j++){
                  index0j=index0k+j*target.Ka;
                  for(i=0;i<target.K;i++){ 
                     index0i=index0j+i;
                     xdft[i]=interaction_matrix[index0i];
                     if(i>0&&i<target.K){ /* Make circulant */
                        xdft[target.Kpp-i].dat[0]=interaction_matrix[index0i].dat[0]*xsign[tc];
                        xdft[target.Kpp-i].dat[1]=interaction_matrix[index0i].dat[1]*xsign[tc]; 
                     }
                  }
                  memset(&xdft[target.K],0,size0); /* Set the zero-fill K to K+kf-1 to zero */
                  fftwf_execute(plan_tcxdft); /* Execute xDFT */
                  for(i=0;i<target.Ka;i++){ /* Copy 1st Ka elements of xdft back into interaction_matrix */
                     index0i=index0j+i;
                     interaction_matrix[index0i]=xdft[i];
                  }
               }
            }
         }
      }

/* Note: For MPI version, have full y direction data so do not need to copy
   into external scratch vector ydft like for OMP and serial versions */
      for(tc=0;tc<6;tc++){ /* yDFT */
         array0=tc*parallel.alloc_tensor;
         for(k=0;k<target.P;k++){
            if(k%parallel.np==parallel.myrank){ /* Plane k is on proc (k%np) */
               parallel.plane=((int)((double)k/parallel.dnp)); /* Plane k is plane ((int)(k/np)) on that proc */
               index0k=parallel.plane*target.KaJpp+array0;
               for(i=0;i<target.Ka;i++){ /* Make circulant */
                  index0i=index0k+i;
                  for(j=1;j<target.J;j++){ 
                     index0j=index0i+j*target.Ka;
                     index1j=index0i+(target.Jpp-j)*target.Ka;
                     interaction_matrix[index1j].dat[0]=interaction_matrix[index0j].dat[0]*ysign[tc];
                     interaction_matrix[index1j].dat[1]=interaction_matrix[index0j].dat[1]*ysign[tc];
                  }
                  for(j=0;j<target.jf;j++){ /* Set the zero-fill J to J+jf-1 to zero */
                     interaction_matrix[index0k+(target.J+j)*target.Ka+i]=zero;
                  }
               }
                /* Execute yDFT */
               fftwf_execute_dft(plan_tcydft,(fftwf_complex *)(&interaction_matrix[index0k]),(fftwf_complex *)(&interaction_matrix[index0k]));
            }
         }
      }
      if(timing.enabled){ /* Finish timing and add contribution */
         finish=MPI_Wtime();
         timing.dft_tensor_components+=(finish-start);
      }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Transpose the 6 independent tensor components of the interaction
   matrix from xy-planes to xz-planes (Distributed transpose) */
/* Transposes the 6 tensor components separately. These transposes come from an earlier development
   code and are essentially a safety net for `transpose_tensor_component_6in1_2'. Since this transpose
   algorithm sends full planes the message sizes are large so it is possible that the message size would
   be too large for all 6 tensor components together.

   transpose_tensor_component_0: Transposes interaction_matrix[6][Ka*Jpp*Pa]
                                 line-by-line (lines in the x-direction)
   transpose_tensor_component_1: Transposes interaction_matrix[6][Ka*Jpp*Pa]
                                 plane-by-plane (includes a local transpose
                                 and only works if message size is below the
                                 MPI EAGER limit
   transpose_tensor_component_2: Transposes interaction_matrix[6][Ka*Jpp*Pa]
                                 plane-by-plane (includes a local transpose) 

   Note: Due to the EAGER limit, transposes

         transpose_tensor_component_0()
         transpose_tensor_component_1()
         transpose_tensor_component_6in1_0()
         transpose_tensor_component_6in1_1()

   are only compliant with MPICH1 and NOT MPICH2 and OPENMPI in their default
   configurations. Transposes

         transpose_tensor_component_2()
         transpose_tensor_component_6in1_2()

   are fully compliant and are the recommended choice */
      if(timing.enabled){start=MPI_Wtime();} /* Start timing */
      if(parallel.tensor_transpose==0){
         for(tc=0;tc<6;tc++){
            transpose_tensor_component_0(target.Jpp,target.Pa,parallel.limitY_tensor,parallel.limitZ_tensor,tc);
         }
      }
      else if(parallel.tensor_transpose==1){
         for(tc=0;tc<6;tc++){
            transpose_tensor_component_1(target.Jpp,target.Pa,parallel.planesY_tensor,parallel.planesZ_tensor,parallel.Znp_tensor,tc);
         }
      }
      else if(parallel.tensor_transpose==2){
         for(tc=0;tc<6;tc++){
            transpose_tensor_component_2(target.Jpp,target.Pa,parallel.planesY_tensor,parallel.planesZ_tensor,parallel.Znp_tensor,tc);
         }
      }
/* Transposes the 6 tensor components together */
      else if(parallel.tensor_transpose==3){
         transpose_tensor_component_6in1_0(target.Jpp,target.Pa,parallel.limitY_tensor,parallel.limitZ_tensor);
      }
      else if(parallel.tensor_transpose==4){
         transpose_tensor_component_6in1_1(target.Jpp,target.Pa,parallel.planesY_tensor,parallel.planesZ_tensor,parallel.Znp_tensor);
      }
      else{
         transpose_tensor_component_6in1_2(target.Jpp,target.Pa,parallel.planesY_tensor,parallel.planesZ_tensor,parallel.Znp_tensor);
      }
      if(timing.enabled){ /* Finish timing and add contribution */
         finish=MPI_Wtime();
         timing.transpose_tensor_forward+=(finish-start);
      }

      if(timing.enabled){start=MPI_Wtime();} /* Start timing */
      size0=target.pf*sizeof(fcomplex); /* Zero-padding size */
      for(tc=0;tc<6;tc++){ /* zDFT */
         array0=tc*parallel.alloc_tensor;
         for(j=0;j<target.Jpp;j++){ /* Extract the Jpp*P sections into zdft one at a time, make them circulant, DFT
                                       and place the first P element of the result back into interaction_matrix */
            if(j%parallel.np==parallel.myrank){ /* Plane j is on proc (j%np) */
               parallel.plane=((int)((double)j/parallel.dnp)); /* Plane j is plane ((int)(j/np)) on that proc */
               index0j=parallel.plane*target.Qa+array0;
               for(i=0;i<target.Ka;i++){
                  index0i=index0j+i;
                  for(k=0;k<target.P;k++){
                     index0k=index0i+k*target.Ka;
                     zdft[k]=interaction_matrix[index0k];
                     if(k>0&&k<target.P){ /* Make circulant */
                        zdft[target.Ppp-k].dat[0]=interaction_matrix[index0k].dat[0]*zsign[tc];
                        zdft[target.Ppp-k].dat[1]=interaction_matrix[index0k].dat[1]*zsign[tc];
                     }
                  }
                  memset(&zdft[target.P],0,size0); /* Set the zero-fill P to P+pf-1 to zero */
                  fftwf_execute(plan_tczdft); /* Execute zDFT */
                  for(k=0;k<target.Pa;k++){ /* Copy 1st Pa elements of zdft back into interaction_matrix */
                     index0k=index0i+k*target.Ka;
                     interaction_matrix[index0k]=zdft[k];
                  }
               }
            }
         }
      }
      if(timing.enabled){ /* Finish timing and add contribution */
         finish=MPI_Wtime();
         timing.dft_tensor_components+=(finish-start);
      }
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Perform the DFT of the 3 vector components of zero_padded_vector, the
   subsequent Fourier domain multplication with the 6 independent tensor
   components and the iDFT of the resultant */
   if(timing.enabled){start=MPI_Wtime();} /* Start timing */
   for(vc=0;vc<3;vc++){ /* yDFT */
      array0=vc*parallel.alloc_padded;
      for(k=0;k<target.P;k++){
         if(k%parallel.np==parallel.myrank){ /* Plane k is on proc (k%np) */
            parallel.plane=((int)((double)k/parallel.dnp)); /* Plane k is plane ((int)(k/np)) on that proc */
            index0k=parallel.plane*target.Mv+array0;
            fftwf_execute_dft(plan_vcydft,(fftwf_complex *)(&zero_padded_vector[index0k]),(fftwf_complex *)(&zero_padded_vector[index0k]));
         }
      }
   }
   if(timing.enabled){ /* Finish timing and add contribution */
      finish=MPI_Wtime();
      timing.zdft_zidft_and_fd_mult+=(finish-start);
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Transpose the 3 zero-padded vector components from xy-planes to xz-planes
   (Distributed transpose) */
/* Transposes the 3 zero-padded vector components separately. These transposes come
   from an earlier development code and are essentially a safety net for
   `transpose_vector_component_3in1_2'. Since this transpose algorithm sends full planes the
   message sizes are large so it is possible that the message size would
   be too large for all 3 vector components together. 

   transpose_vector_component_0: Transposes zero_padded_vector[3][K*Jpp*P]
                                 line-by-line (lines in the x-direction)
   transpose_vector_component_1: Transposes zero_padded_vector[3][K*Jpp*P]
                                 plane-by-plane (includes a local transpose
                                 and only works if message size is below the
                                 MPI EAGER limit)
   transpose_vector_component_2: Transposes zero_padded_vector[3][K*Jpp*P]
                                 plane-by-plane (includes a local transpose) 

   Note: Due to the EAGER limit, transposes

         transpose_vector_component_0()
         transpose_vector_component_1()
         transpose_vector_components_3in1_0()
         transpose_vector_components_3in1_1()

   are only compliant with MPICH1 and NOT MPICH2 and OPENMPI in their default
   configurations. Transposes

         transpose_vector_component_2()
         transpose_vector_component_3in1_2()

   are fully compliant and are the recommended choice */
   if(timing.enabled){start=MPI_Wtime();} /* Start timing */
   if(parallel.padded_transpose==0){
      for(vc=0;vc<3;vc++){
         transpose_vector_component_0(target.Jpp,target.P,parallel.limitY_padded,parallel.limitZ_padded,vc);
      }
   }
   else if(parallel.padded_transpose==1){
      for(vc=0;vc<3;vc++){
         transpose_vector_component_1(target.Jpp,target.P,parallel.planesY_padded,parallel.planesZ_padded,parallel.Znp_padded,vc);
      }
   }
   else if(parallel.padded_transpose==2){
      for(vc=0;vc<3;vc++){
         transpose_vector_component_2(target.Jpp,target.P,parallel.planesY_padded,parallel.planesZ_padded,parallel.Znp_padded,vc);
      }
   }
/* Transposes the 3 zero-padded vector components together */
   else if(parallel.padded_transpose==3){
      transpose_vector_component_3in1_0(target.Jpp,target.P,parallel.limitY_padded,parallel.limitZ_padded);
   }
   else if(parallel.padded_transpose==4){
      transpose_vector_component_3in1_1(target.Jpp,target.P,parallel.planesY_padded,parallel.planesZ_padded,parallel.Znp_padded);
   }
   else{
      transpose_vector_component_3in1_2(target.Jpp,target.P,parallel.planesY_padded,parallel.planesZ_padded,parallel.Znp_padded);
   }
   if(timing.enabled){ /* Finish timing and add contribution */
      finish=MPI_Wtime();
      timing.transpose_vector_forward+=(finish-start);
   }

   if(timing.enabled){start=MPI_Wtime();} /* Start timing */
   size0=target.K*sizeof(fcomplex); /* Copy size */
   size1=target.Kpp*(target.Ppp-target.P)*sizeof(fcomplex); /* Zero-padding size */
   size2=(target.Kpp-target.K)*sizeof(fcomplex); /* Zero-padding size */
   for(j=0;j<target.Jpp;j++){
      /* No need for jsign as for the MPI version are storing full j direction data */
      if(j%parallel.np==parallel.myrank){ /* Plane j is on proc (j%np) */
         parallel.plane=((int)((double)j/parallel.dnp)); /* Plane j is plane ((int)(j/np)) on that proc */
         index0j=parallel.plane*target.Q;index1j=parallel.plane*target.Qa;
         /* Now using j instead of jj as for the MPI version are storing full j direction data */
         for(vc=0;vc<3;vc++){  /* Extract the 3 component sections into xzscratch and zero-pad them */
            array1=vc*parallel.alloc_padded+index0j;
            array2=vc*target.Qpp;
            for(k=0;k<target.P;k++){
               memcpy(&xzscratch[array2+k*target.Kpp],&zero_padded_vector[array1+k*target.K],size0); /* Copy data */
               memset(&xzscratch[array2+k*target.Kpp+target.K],0,size2); /* Zero-padding */
            }
            memset(&xzscratch[array2+target.PKpp],0,size1); /* Zero-padding */
         }
         for(vc=0;vc<3;vc++){ /* xDFT */
            fftwf_execute_dft(plan_vcxdft,(fftwf_complex *)(&xzscratch[vc*target.Qpp]),(fftwf_complex *)(&xzscratch[vc*target.Qpp]));
         }
         for(vc=0;vc<3;vc++){ /* zDFT */
            fftwf_execute_dft(plan_vczdft,(fftwf_complex *)(&xzscratch[vc*target.Qpp]),(fftwf_complex *)(&xzscratch[vc*target.Qpp]));
         }
         /* Fourier Domain multiplication */
         for(k=0;k<target.Ppp;k++){
            kk=k>=target.Pa?(target.Ppp-k):k; /* Extract pertinent data from stored section of the interaction matrix */
            ksign=k>=target.Pa?-1.0F:1.0F; /* Sign information */
            index0k=k*target.Kpp;index1k=index1j+kk*target.Ka;
            for(i=0;i<target.Kpp;i++){
               ii=i>=target.Ka?(target.Kpp-i):i; /* Extract pertinent data from stored section of the interaction matrix */
               isign=i>=target.Ka?-1.0F:1.0F; /* Sign information */
               index0i=index0k+i;index1i=index1k+ii;
               array1=target.Qpp+index0i;
               array2=target.Qpp+array1; /* 2*target.Qpp+index0i */
               m0=interaction_matrix[index1i]; /* xx tensor component */
               m1=interaction_matrix[parallel.alloc_tensor+index1i]; /* xy=yx tensor component */
               m2=interaction_matrix[2*parallel.alloc_tensor+index1i]; /* xz=xz tensor component */
               m3=interaction_matrix[3*parallel.alloc_tensor+index1i]; /* yy tensor component */
               m4=interaction_matrix[4*parallel.alloc_tensor+index1i]; /* yz=zy tensor component */
               m5=interaction_matrix[5*parallel.alloc_tensor+index1i]; /* zz tensor component */
               v0=xzscratch[index0i]; /* x vector component */
               v1=xzscratch[array1]; /* y vector component */
               v2=xzscratch[array2]; /* z vector component */
               xzscratch[index0i].dat[0]=(m0.dat[0]*v0.dat[0]-m0.dat[1]*v0.dat[1])+/*[0]*/
                              (m1.dat[0]*v1.dat[0]-m1.dat[1]*v1.dat[1])*isign+/*[1]isign*/
                              (m2.dat[0]*v2.dat[0]-m2.dat[1]*v2.dat[1])*isign*ksign;/*[2]isign*ksign*/
               xzscratch[index0i].dat[1]=(m0.dat[0]*v0.dat[1]+m0.dat[1]*v0.dat[0])+/*[0]*/
                              (m1.dat[0]*v1.dat[1]+m1.dat[1]*v1.dat[0])*isign+/*[1]isign*/
                              (m2.dat[0]*v2.dat[1]+m2.dat[1]*v2.dat[0])*isign*ksign;/*[2]isign*ksign*/
               xzscratch[array1].dat[0]=(m1.dat[0]*v0.dat[0]-m1.dat[1]*v0.dat[1])*isign+/*[1]isign*/
                              (m3.dat[0]*v1.dat[0]-m3.dat[1]*v1.dat[1])+/*[3]*/
                              (m4.dat[0]*v2.dat[0]-m4.dat[1]*v2.dat[1])*ksign;/*[4]ksign*/
               xzscratch[array1].dat[1]=(m1.dat[0]*v0.dat[1]+m1.dat[1]*v0.dat[0])*isign+/*[1]isign*/
                              (m3.dat[0]*v1.dat[1]+m3.dat[1]*v1.dat[0])+/*[3]*/
                              (m4.dat[0]*v2.dat[1]+m4.dat[1]*v2.dat[0])*ksign;/*[4]ksign*/
               xzscratch[array2].dat[0]=(m2.dat[0]*v0.dat[0]-m2.dat[1]*v0.dat[1])*ksign*isign+/*[2]ksign*isign*/
                              (m4.dat[0]*v1.dat[0]-m4.dat[1]*v1.dat[1])*ksign+/*[4]ksign*/
                              (m5.dat[0]*v2.dat[0]-m5.dat[1]*v2.dat[1]);/*[5]*/
               xzscratch[array2].dat[1]=(m2.dat[0]*v0.dat[1]+m2.dat[1]*v0.dat[0])*ksign*isign+/*[2]ksign*isign*/
                              (m4.dat[0]*v1.dat[1]+m4.dat[1]*v1.dat[0])*ksign+/*[4]ksign*/
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
            array1=vc*parallel.alloc_padded+index0j;
            array2=vc*target.Qpp;
            for(k=0;k<target.P;k++){
               memcpy(&zero_padded_vector[array1+k*target.K],&xzscratch[array2+k*target.Kpp],size0); /* Copy data */
            }
         } /* End vc loop */
      } /* End if */
   } /* End j loop */
   if(timing.enabled){ /* Finish timing and add contribution */
      finish=MPI_Wtime();
      timing.zdft_zidft_and_fd_mult+=(finish-start);
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Transpose the 3 zero-padded vector components from xz-planes to xy-planes
   (Distributed transpose) */
/* Transposes the 3 zero-padded vector components separately. These transposes come
   from an earlier development code and are essentially a safety net for
   `transpose_vector_component_3in1_2'. Since this transpose algorithm sends full planes the
   message sizes are large so it is possible that the message size would
   be too large for all 3 vector components together. 

   transpose_vector_component_0: Transposes zero_padded_vector[3][K*Jpp*P]
                                 line-by-line (lines in the x-direction)
   transpose_vector_component_1: Transposes zero_padded_vector[3][K*Jpp*P]
                                 plane-by-plane (includes a local transpose
                                 and only works if message size is below the
                                 MPI EAGER limit
   transpose_vector_component_2: Transposes zero_padded_vector[3][K*Jpp*P]
                                 plane-by-plane (includes a local transpose) 

   Note: Due to the EAGER limit, transposes

         transpose_vector_component_0()
         transpose_vector_component_1()
         transpose_vector_component_3in1_0()
         transpose_vector_component_3in1_1()

   are only compliant with MPICH1 and NOT MPICH2 and OPENMPI in their default
   configurations. Transposes

         transpose_vector_component_2()
         transpose_vector_component_3in1_2()

   are fully compliant and are the recommended choice */
   if(timing.enabled){start=MPI_Wtime();} /* Start timing */
   if(parallel.padded_transpose==0){
      for(vc=0;vc<3;vc++){transpose_vector_component_0(target.P,target.Jpp,parallel.limitZ_padded,parallel.limitY_padded,vc);}
   }
   else if(parallel.padded_transpose==1){
      for(vc=0;vc<3;vc++){transpose_vector_component_1(target.P,target.Jpp,parallel.planesZ_padded,parallel.planesY_padded,parallel.Ynp,vc);}
   }
   else if(parallel.padded_transpose==2){
      for(vc=0;vc<3;vc++){transpose_vector_component_2(target.P,target.Jpp,parallel.planesZ_padded,parallel.planesY_padded,parallel.Ynp,vc);}
   }
/* Transposes the 3 zero-padded vector components together */
   else if(parallel.padded_transpose==3){
      transpose_vector_component_3in1_0(target.P,target.Jpp,parallel.limitZ_padded,parallel.limitY_padded);
   }      
   else if(parallel.padded_transpose==4){
      transpose_vector_component_3in1_1(target.P,target.Jpp,parallel.planesZ_padded,parallel.planesY_padded,parallel.Ynp);
   }
   else{
      transpose_vector_component_3in1_2(target.P,target.Jpp,parallel.planesZ_padded,parallel.planesY_padded,parallel.Ynp);
   }
   if(timing.enabled){ /* Finish timing and add contribution */
      finish=MPI_Wtime();
      timing.transpose_vector_reverse+=(finish-start);
   }

   if(timing.enabled){start=MPI_Wtime();} /* Start timing */
   for(vc=0;vc<3;vc++){ /* yiDFT */
      array0=vc*parallel.alloc_padded;
      for(k=0;k<target.P;k++){
         if(k%parallel.np==parallel.myrank){ /* Plane k is on proc (k%np) */
            parallel.plane=((int)((double)k/parallel.dnp)); /* Plane k is plane ((int)(k/np)) on that proc */
            index0k=parallel.plane*target.Mv+array0;
            fftwf_execute_dft(plan_vcydfti,(fftwf_complex *)(&zero_padded_vector[index0k]),(fftwf_complex *)(&zero_padded_vector[index0k]));
         }
      }
   }

   if(timing.enabled){ /* Finish timing and add contribution */
      finish=MPI_Wtime();
      timing.zdft_zidft_and_fd_mult+=(finish-start);
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Normalisation + Take diagonal contributions into account + Compress data */
   if(timing.enabled){start=MPI_Wtime();} /* Start timing */
   grid_points=1.0F/(float)target.Npp; /* Normalisation value */
   for(vc=0;vc<3;vc++){
      array0=vc*parallel.alloc_vector;
      array1=vc*parallel.alloc_padded;
      for(k=0;k<target.P;k++){
         if(k%parallel.np==parallel.myrank){ /* Plane k is on proc (k%np) */
            parallel.plane=((int)((double)k/parallel.dnp)); /* Plane k is plane ((int)(k/np)) on that proc */
            index0k=k*target.M;index1k=parallel.plane*target.Mv+array1;
            element=target.populated[parallel.plane]; /* Start index for plane k for the stored non-zero lattice sites */
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
   }
   if(timing.enabled){ /* Finish timing and add contribution */
      finish=MPI_Wtime();
      timing.norm_diag_compress+=(finish-start);
   }
}
