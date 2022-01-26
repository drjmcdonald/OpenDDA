/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Calculates an estmation of the memory requirements
   per node

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "guesstimate_MPI.h"

attributes_iterative iterative={0};
attributes_target target={0};
attributes_parallel parallel={0};
MPI_Op add_dcomplex; /* User-defined MPI function to add complex_numbers */

int main(int argc,char **argv){

   int i,j,k,m,error,shape,file,index0i,index0j,index0k,position[3];
   int maxtensor=0,maxpadded=0,maxvector=0;
	int index0,buffind,bufsize,maxbuffer=0,tensor_mpi_buffer=0;
   double sumtensor=0.0,sumpadded=0.0,sumvector=0.0,sumbuffer=0.0;
	double tmptensor,tmppadded,tmpvector,tmpbuffer;
   char in[50],string[200],compare[200],filename[100];
   double base_memory=0.0,sum_memory=0.0,iter_sum_memory,iter_memory,onegb=1024.0*1024.0*1024.0;
   double ii,jj,kk,i0,j0,k0,isqrd,jsqrd,ksqrd,Ksqd,Jsqd,Psqd;
   double sizedc=(double)sizeof(dcomplex);
   FILE *data;

   mpi_initialise(&argc,&argv); /* Initialise the MPI environment */

   if(parallel.myrank==0){ /* Restrict to master */
/* ------------------------------------------------------------------------- */
/* K */
      fprintf(stdout,"Please enter K [the integer x dimension of the dipole array]\n");fflush(stdout);
      target.K=-1;
      while(target.K<1){
         error=0; /* Initialise error */
         reset_string(in); /* Reset string */
         i=0; /* Initialise i */
         scanf("%s",in); /* Read string in */
         do{ /* Check for invalid input */
            if((isdigit(in[i++]))==0){error=1;}
         } while((in[i]!='\n')&&(in[i]!='\0'));
         if(error==1){ /* If any character is not a digit then print an error */
            target.K=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nArgument for K contained an invalid character!\n\nPlease re-enter\n");fflush(stdout);
         }
         else{
            target.K=atoi(in); /* Convert to an integer value */
            if(target.K<1){ /* Check for any an invalid value */
               target.K=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nK MUST be an integer > 0\n\nPlease re-enter\n");fflush(stdout);
            }
         }
      }

/* ------------------------------------------------------------------------- */
/* J */
      fprintf(stdout,"Please enter J [the integer y dimension of the dipole array]\n");fflush(stdout);
      target.J=-1;
      while(target.J<1){
         error=0; /* Initialise error */
         reset_string(in); /* Reset string */
         i=0; /* Initialise i */
         scanf("%s",in); /* Read string in */
         do{ /* Check for invalid input */
            if((isdigit(in[i++]))==0){error=1;}
         } while((in[i]!='\n')&&(in[i]!='\0'));
         if(error==1){ /* If any character is not a digit then print an error */
            target.J=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nArgument for J contained an invalid character!\n\nPlease re-enter\n");fflush(stdout);
         }
         else{
            target.J=atoi(in); /* Convert to an integer value */
            if(target.J<1){ /* Check for any an invalid value */
               target.J=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nJ MUST be an integer > 0\n\nPlease re-enter\n");fflush(stdout);
            }
         }
      }

/* ------------------------------------------------------------------------- */
/* P */
      fprintf(stdout,"Please enter P [the integer z dimension of the dipole array]\n");fflush(stdout);
      target.P=-1;
      while(target.P<1){
         error=0; /* Initialise error */
         reset_string(in); /* Reset string */
         i=0; /* Initialise i */
         scanf("%s",in); /* Read string in */
         do{ /* Check for invalid input */
            if((isdigit(in[i++]))==0){error=1;}
         } while((in[i]!='\n')&&(in[i]!='\0'));
         if(error==1){ /* If any character is not a digit then print an error */
            target.P=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nArgument for P contained an invalid character!\n\nPlease re-enter\n");fflush(stdout);
         }
         else{
            target.P=atoi(in); /* Convert to an integer value */
            if(target.P<1){ /* Check for any an invalid value */
               target.P=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nP MUST be an integer > 0\n\nPlease re-enter\n");fflush(stdout);
            }
         }
      }
   } /* End if(parallel.myrank==0){ Restrict to master */

   MPI_Bcast(&target.K,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&target.J,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&target.P,1,MPI_INT,0,MPI_COMM_WORLD);

/* FFTW is best at handling sizes of the form 2^{a}3^{b}5^{c}7^{d}11^{e}13^{f},
   where (e + f) is either 0 or 1, and the other exponents are arbitrary. */
   find_dft_size(target.K,&target.Kp,&target.Kpp,&target.kf);
   find_dft_size(target.J,&target.Jp,&target.Jpp,&target.jf);
   find_dft_size(target.P,&target.Pp,&target.Ppp,&target.pf);

   target.M=target.K*target.J; /* M=K*J */
   target.N=target.K*target.J*target.P; /* N=K*J*P */
   target.N3=3*target.N; /* N3=N*3 */
   target.Ka=target.K+target.kf; /* Ka=K+kf */
   target.Ja=target.J+target.jf; /* Ja=J+jf */
   target.Pa=target.P+target.pf; /* Pa=P+pf */
   target.Na=target.Ka*target.Ja*target.Pa; /* Na=Ka*Ja*Pa */
   target.Mpp=target.Kpp*target.Jpp; /* Mpp=Kpp*Jpp */
   target.Qpp=target.Kpp*target.Ppp; /* Qpp=Kpp*Ppp */
   target.KaJpp=target.Ka*target.Jpp; /* KaJpp=Ka*Jpp */
   target.Npp=target.Kpp*target.Jpp*target.Ppp; /* Npp=Kpp*Jpp*Ppp */

/* ------------------------------------------------------------------------- */
/* Target shape and number of occupied dipoles */
   if(parallel.myrank==0){ /* Restrict to master */
      fprintf(stdout,"Please select the required shape using the appropriate integer\n");fflush(stdout);
      fprintf(stdout,"   0. ellipsoid\n");fflush(stdout);
      fprintf(stdout,"   1. cuboid\n");fflush(stdout);
      fprintf(stdout,"   2. custom from file\n");fflush(stdout);
      shape=-1;
      while(shape<0){
         error=0; /* Initialise error */
         reset_string(in); /* Reset string */
         i=0; /* Initialise i */
         scanf("%s",in); /* Read string in */
         do{ /* Check for invalid input */
            if((isdigit(in[i++]))==0){error=1;}
         } while((in[i]!='\n')&&(in[i]!='\0'));
         if(error==1){ /* If any character is not a digit then print an error */
            shape=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nArgument for shape contained an invalid character!\n\nPlease re-enter\n");fflush(stdout);
         }
         else{
            shape=atoi(in); /* Convert to an integer value */
            if((shape<0)||(shape>2)){ /* Check for any an invalid value */
               shape=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nThe shape identifier must be an integer in the range 0..2\n\nPlease re-enter\n");fflush(stdout);
            }
         }
      }
   } /* End if(parallel.myrank==0){ Restrict to master */

   /* Setup bitfielded occupied/unoccupied flag array */
   target.occ_ctrl=(int)((double)target.N/32.0+1.0);
   uint_calloc(&target.occupied,(size_t)target.occ_ctrl);

   if(parallel.myrank==0){ /* Restrict to master */
      if(shape!=2){ /* Built-in target construction */
         if(shape==0){ /* Ellipsoid */
            Ksqd=(double)(target.K*target.K);
            Jsqd=(double)(target.J*target.J);
            Psqd=(double)(target.P*target.P);
   
/* Origin */
            i0=((double)target.K-1.0)/2.0;
            j0=((double)target.J-1.0)/2.0;
            k0=((double)target.P-1.0)/2.0;
   
            target.Nd=0; /* Initialise Nd */
            for(k=0;k<target.P;k++){
               index0k=k*target.M;
               kk=(double)k-k0;
               ksqrd=(kk*kk)/Psqd;
               if(ksqrd<0.25){
                  for(j=0;j<target.J;j++){
                     index0j=index0k+j*target.K;
                     jj=(double)j-j0;
                     jsqrd=(jj*jj)/Jsqd;
                     if((jsqrd+ksqrd)<0.25){
                        for(i=0;i<target.K;i++){
                           ii=(double)i-i0;
                           isqrd=(ii*ii)/Ksqd;
                           if((isqrd+jsqrd+ksqrd)<0.25){
                              index0i=index0j+i;
                              target.occupied[(int)((double)index0i/32.0)]|=(1<<(index0i%32)); /* Set occupied flag */
                              target.Nd++; /* Increment Nd */
                           }
                        }
                     }
                  }
               }
            }
         }
         else if(shape==1){ /* Cuboid */
            for(k=0;k<target.P;k++){
               index0k=k*target.M;
               for(index0j=0;index0j<target.M;index0j++){
                  index0i=index0k+index0j;
                  target.occupied[(int)((double)index0i/32.0)]|=(1<<(index0i%32)); /* Set occupied flag */
                  /* Homogeneous & isotropic cuboid. Set composition data array target.composition[0..2]=0 already */
                  target.Nd++; /* Increment Nd */
               }
            }
         }
      }
      else{ /* Custom target from file */
         fprintf(stdout,"Please enter a valid file name in this directory\n");fflush(stdout);
         file=-1;
         while(file<0){
            error=0; /* Initialise error */
            reset_string(filename); /* Reset string */
            i=0; /* Initialise i */
            scanf("%s",filename); /* Read string in */
            if((data=fopen(filename,"r"))==NULL){error=1;} /* Check if file exists */
            if(error==1){ /* If file does not exist then print an error */
               file=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nThe file \"%s\" could not be opened \n\nPlease re-enter\n",filename);fflush(stdout);
            }
         }
         target.Nd=0; /* Initialise Nd */
         reset_string(string);
         while((fgets(string,sizeof(string),data))!=NULL){
            if(string[0]=='#'){ /* Ignore comments */
               continue;
            }
            else{
               i=-1;
               /* Read in the x, y, and z coordinate */
               for(m=0;m<3;m++){
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
                  k=0;reset_string(compare);
                  do{ /* Copy everything after '=' OR ',' to ',' OR '\n' to compare string,
                     ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
                     if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
                     else{compare[k++]=string[i];}
                  } while((string[i+1]!=',')&&(string[i+1]!='\n'));
                  position[m]=atoi(compare); /* string to int conversion */
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               }
               /* Check the x, y & z data */
               if(position[0]<0){ /* x data MUST obey 0<x<K */
                  reset_string(string);
                  sprintf(string,"An invalid x coordinate was encountered in the target input file \"%s\"",target.shape);
                  print_error("Target construction error",string,"x data MUST be >= 0. Please check the target data input file");
               }
               if(position[0]>=target.K){ /* x data MUST obey 0<x<K */
                  reset_string(string);
                  sprintf(string,"An invalid x coordinate was encountered in the target input file \"%s\"",target.shape);
                  print_error("Target construction error",string,"x data MUST be < K. Please check the target data input file");
               }
               if(position[1]<0){ /* y data MUST obey 0<y<J */
                  reset_string(string);
                  sprintf(string,"An invalid y coordinate was encountered in the target input file \"%s\"",target.shape);
                  print_error("Target construction error",string,"y data MUST be >= 0. Please check the target data input file");
               }
               if(position[1]>=target.J){ /* y data MUST obey 0<y<J */
                  reset_string(string);
                  sprintf(string,"An invalid y coordinate was encountered in the target input file \"%s\"",target.shape);
                  print_error("Target construction error",string,"y data MUST be < J. Please check the target data input file");
               }
               if(position[2]<0){ /* z data MUST obey 0<z<P */
                  reset_string(string);
                  sprintf(string,"An invalid z coordinate was encountered in the target input file \"%s\"",target.shape);
                  print_error("Target construction error",string,"z data MUST be >= 0. Please check the target data input file");
               }
               if(position[2]>=target.P){ /* z data MUST obey 0<z<P */
                  reset_string(string);
                  sprintf(string,"An invalid z coordinate was encountered in the target input file \"%s\"",target.shape);
                  print_error("Target construction error",string,"z data MUST be < P. Please check the target data input file");
               }
               /* Set the occupied flag and increment Nd */
               index0i=position[2]*target.M+position[1]*target.K+position[0];
               target.occupied[(int)((double)index0i/32.0)]|=(1<<(index0i%32)); /* Set occupied flag */
               target.Nd++; /* Increment Nd */
            }
         }
         fclose(data);
      }
   } /* End if(parallel.myrank==0){ Restrict to master */

   MPI_Bcast(&target.Nd,1,MPI_LONG,0,MPI_COMM_WORLD);
   MPI_Bcast(target.occupied,target.occ_ctrl,MPI_UNSIGNED,0,MPI_COMM_WORLD);
      
/* ------------------------------------------------------------------------- */
/* Number of MPI nodes */
   parallel.dnp=(double)parallel.np;

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Determine memory allocation sizes for each node */
   memory_allocation_sizes();

/* Count the number of occupied lattice sites for each processor */
   for(k=0;k<target.P;k++){
      if(k%parallel.np==parallel.myrank){
         parallel.plane=((int)((double)k/parallel.dnp));
         index0k=k*target.M;
         for(j=0;j<target.M;j++){
            index0i=index0k+j;
            if(target.occupied[(int)((double)index0i/32.0)]&(1<<(index0i%32))){ /* Lattice site is occupied */
               parallel.alloc_vector++; /* Increment the number of occupied lattice sites on this processor */
            }
         }
      }
   }
   uint_free(&target.occupied);

   MPI_Reduce(&parallel.alloc_tensor,&maxtensor,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD); /* Find max */
   MPI_Reduce(&parallel.alloc_padded,&maxpadded,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD); /* Find max */
   MPI_Reduce(&parallel.alloc_vector,&maxvector,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD); /* Find max */

   tmptensor=(double)parallel.alloc_tensor;
   tmppadded=(double)parallel.alloc_padded;
   tmpvector=(double)parallel.alloc_vector;
   MPI_Reduce(&tmptensor,&sumtensor,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); /* Find sum */
   MPI_Reduce(&tmppadded,&sumpadded,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); /* Find sum */
   MPI_Reduce(&tmpvector,&sumvector,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); /* Find sum */

/* ------------------------------------------------------------------------- */
/* Transpose algorithm */
   if(parallel.myrank==0){ /* Restrict to master */
      fprintf(stdout,"Transpose algorithm? 0..5\n");fflush(stdout);
      parallel.tensor_transpose=-1;
      while((parallel.tensor_transpose<0)||(parallel.tensor_transpose>5)){
         error=0; /* Initialise error */
         reset_string(in); /* Reset string */
         i=0; /* Initialise i */
         scanf("%s",in); /* Read string in */
         do{ /* Check for invalid input */
            if((isdigit(in[i++]))==0){error=1;}
         } while((in[i]!='\n')&&(in[i]!='\0'));
         if(error==1){ /* If any character is not a digit then print an error */
            parallel.tensor_transpose=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nArgument for transpose algorithm contained an invalid character!\n\nPlease re-enter\n");fflush(stdout);
         }
         else{
            parallel.tensor_transpose=atoi(in); /* Convert to an integer value */
            if((parallel.tensor_transpose<0)||(parallel.tensor_transpose>5)){ /* Check for any an invalid value */
               parallel.tensor_transpose=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nTranspose algorithm must be in the range 0..5\n\nPlease re-enter\n");fflush(stdout);
            }
         }
      }
   } /* End if(parallel.myrank==0){ Restrict to master */
   MPI_Bcast(&parallel.tensor_transpose,1,MPI_INT,0,MPI_COMM_WORLD);

   if(parallel.tensor_transpose==0||parallel.tensor_transpose==1||parallel.tensor_transpose==3||parallel.tensor_transpose==4){
     tensor_mpi_buffer=0.0;
   }
   else{
      index0=0;i=0;i0=0;buffind=0;bufsize=0;
      while(i<target.Jpp||i0<target.Jpp||index0<parallel.planesZ_tensor){
         while((i<((index0+1)*parallel.np)||(i*parallel.planesY_tensor)<((index0+1)*target.Pa))&&i<target.Jpp){
            buffind++;
            i++;
         }
         for(j=0;j<parallel.np;j++){
            if(i0<target.Jpp){
               buffind--;
               i0++;
            }
         }
         index0++;
         bufsize=buffind>bufsize?buffind:bufsize;
      }
      bufsize+=parallel.np;
      if(parallel.tensor_transpose==2){
         tensor_mpi_buffer=(bufsize*parallel.planesY_tensor*target.Ka);
      }
      else if(parallel.tensor_transpose==5){
         tensor_mpi_buffer=(bufsize*parallel.planesY_tensor*target.Ka*6);
      }
   }
   MPI_Reduce(&tensor_mpi_buffer,&maxbuffer,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD); /* Find max */
   tmpbuffer=(double)tensor_mpi_buffer;
   MPI_Reduce(&tmpbuffer,&sumbuffer,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); /* Find sum */

   if(parallel.myrank==0){ /* Restrict to master */
      parallel.alloc_vector3=maxvector*3;

      /* Main arrays */
      base_memory+=((double)parallel.alloc_vector3); /* incident_E_field[3][Nd] */
      sum_memory+=(sumvector*3.0);
      base_memory+=((double)maxpadded*3.0); /* zero_padded_vector[3][P][Jpp][K] */
      sum_memory+=(sumpadded*3.0);
      base_memory+=((double)maxtensor*6.0); /* interaction_matrix[6][Pa][Ja][Ka] */
      sum_memory+=(sumtensor*6.0);
      base_memory+=((double)parallel.alloc_vector3); /* dipole_polarisation[3][Nd] */
      sum_memory+=(sumvector*3.0);

      /* MPI buffer */
      base_memory+=tensor_mpi_buffer; /* mpi buffer */
      sum_memory+=(sumbuffer);

      /* DFT scratch arrays */
      base_memory+=((double)target.Kpp); /* xdft[Kpp] */
      sum_memory+=(parallel.dnp*(double)target.Kpp);
      /* no y scratch array in MPI version */
      base_memory+=((double)target.Ppp); /* zdft[Ppp] */
      sum_memory+=(parallel.dnp*(double)target.Ppp);

      /* DFT and FD multiplication scratch arrays */
      base_memory+=((double)target.Qpp*3.0); /* xzscratch[3*Qpp] */
      sum_memory+=(parallel.dnp*3.0*(double)target.Qpp);

/* ------------------------------------------------------------------------- */
/* Preconditioning */
      fprintf(stdout,"Preconditioning enabled OR disabled? 0 OR 1\n");fflush(stdout);
      iterative.precond=-1;
      while((iterative.precond<0)||(iterative.precond>1)){
         error=0; /* Initialise error */
         reset_string(in); /* Reset string */
         i=0; /* Initialise i */
         scanf("%s",in); /* Read string in */
         do{ /* Check for invalid input */
            if((isdigit(in[i++]))==0){error=1;}
         } while((in[i]!='\n')&&(in[i]!='\0'));
         if(error==1){ /* If any character is not a digit then print an error */
            iterative.precond=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nArgument for preconditioning contained an invalid character!\n\nPlease re-enter\n");fflush(stdout);
         }
         else{
            iterative.precond=atoi(in); /* Convert to an integer value */
            if((iterative.precond<0)||(iterative.precond>1)){ /* Check for any an invalid value */
               iterative.precond=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nArgument MUST be an either 0 OR 1\n\nPlease re-enter\n");fflush(stdout);
            }
         }
      }

      /* Preconditioning and interaction matrix diagonal */
      if(iterative.precond==1){
         base_memory+=((double)parallel.alloc_vector3); /* point_jacobi[3][Nd] */
         sum_memory+=(sumvector*3.0);
      }
      base_memory+=((double)parallel.alloc_vector3); /* interaction_matrix_diagonal[3][Nd] */
      sum_memory+=(sumvector*3.0);

/* ------------------------------------------------------------------------- */
/* Number of starting vectors */
      fprintf(stdout,"Number of starting vectors for the iterative schemes:\nRBiCGSTAB & ML(n)BiCGSTAB*? 1..50\n");fflush(stdout);
      iterative.vec=-1;
      while((iterative.vec<1)||(iterative.vec>50)){
         error=0; /* Initialise error */
         reset_string(in); /* Reset string */
         i=0; /* Initialise i */
         scanf("%s",in); /* Read string in */
         do{ /* Check for invalid input */
            if((isdigit(in[i++]))==0){error=1;}
         } while((in[i]!='\n')&&(in[i]!='\0'));
         if(error==1){ /* If any character is not a digit then print an error */
            iterative.vec=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nArgument for the number of starting vectors contained an invalid character!\n\nPlease re-enter\n");fflush(stdout);
         }
         else{
            iterative.vec=atoi(in); /* Convert to an integer value */
            if((iterative.vec<1)||(iterative.vec>50)){ /* Check for any an invalid value */
               iterative.vec=-1;fprintf(stdout,"\nMemory 'guesstimation' error\nArgument MUST be > 0 AND < 50\n\nPlease re-enter\n");fflush(stdout);
            }
         }
      }

/* ------------------------------------------------------------------------- */
      print_section_title(stdout,"Maximum allocated memory per node");
/* ------------------------------------------------------------------------- */
      fprintf(stdout,"\nNumber of occupied sites in the %dx%dx%d array = %ld\n",target.K,target.J,target.P,target.Nd);fflush(stdout);

      fprintf(stdout,"\nIterative Scheme    Total (GB)   Maximum per node (GB)\n");fflush(stdout);

      /* bicg */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)parallel.alloc_vector3); /* r */
      iter_memory+=((double)parallel.alloc_vector3); /* rtilde */
      iter_memory+=((double)parallel.alloc_vector3); /* p */
      iter_memory+=((double)parallel.alloc_vector3); /* ptilde */
      iter_memory+=((double)parallel.alloc_vector3); /* q */
      iter_sum_memory+=(sumvector*3.0*5.0);
      fprintf(stdout,"\nBiCG                %.3lfGB    %.3lfGB\n",(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
/* ------------------------------------------------------------------------- */
      /* bicg_sym */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)parallel.alloc_vector3); /* r */
      iter_memory+=((double)parallel.alloc_vector3); /* p */
      iter_sum_memory+=(sumvector*3.0*2.0);
      fprintf(stdout,"BiCG_sym            %.3lfGB    %.3lfGB\n",(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
/* ------------------------------------------------------------------------- */
      /* bicgstab */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)parallel.alloc_vector3); /* r */
      iter_memory+=((double)parallel.alloc_vector3); /* rtilde */
      iter_memory+=((double)parallel.alloc_vector3); /* p */
      iter_memory+=((double)parallel.alloc_vector3); /* s */
      iter_memory+=((double)parallel.alloc_vector3); /* v */
      iter_sum_memory+=(sumvector*3.0*5.0);
      fprintf(stdout,"BiCGSTAB            %.3lfGB    %.3lfGB\n",(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
/* ------------------------------------------------------------------------- */
      /* cg */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)parallel.alloc_vector3); /* r */
      iter_memory+=((double)parallel.alloc_vector3); /* p */
      iter_sum_memory+=(sumvector*3.0*2.0);
      fprintf(stdout,"CG                  %.3lfGB    %.3lfGB\n",(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
/* ------------------------------------------------------------------------- */
      /* cgs */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)parallel.alloc_vector3); /* r */
      iter_memory+=((double)parallel.alloc_vector3); /* rtilde */
      iter_memory+=((double)parallel.alloc_vector3); /* p */
      iter_memory+=((double)parallel.alloc_vector3); /* u */
      iter_memory+=((double)parallel.alloc_vector3); /* q */
      iter_sum_memory+=(sumvector*3.0*5.0);
      fprintf(stdout,"CGS                 %.3lfGB    %.3lfGB\n",(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
/* ------------------------------------------------------------------------- */
      /* mlbicgstab_orig */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)(iterative.vec*parallel.alloc_vector3)); /* g */
      iter_memory+=((double)(iterative.vec*parallel.alloc_vector3)); /* omega */
      iter_memory+=((double)(iterative.vec*parallel.alloc_vector3)); /* q */
      iter_sum_memory+=(sumvector*3.0*3.0*(double)iterative.vec);
      if(iterative.vec>1){
         iter_memory+=((double)((iterative.vec-1)*parallel.alloc_vector3)); /* d */
         iter_sum_memory+=(sumvector*3.0*(double)(iterative.vec-1));
      }
      iter_memory+=((double)parallel.alloc_vector3); /* r */
      iter_memory+=((double)parallel.alloc_vector3); /* u */
      iter_memory+=((double)parallel.alloc_vector3); /* zd */
      iter_memory+=((double)parallel.alloc_vector3); /* zg */
      iter_memory+=((double)parallel.alloc_vector3); /* zomega */
      iter_sum_memory+=(sumvector*3.0*5.0);
      if(iterative.precond!=0){ /* Preconditioning enabled */
         iter_memory+=((double)parallel.alloc_vector3); /* gtilde */
         iter_memory+=((double)parallel.alloc_vector3); /* utilde */
         iter_sum_memory+=(sumvector*3.0*2.0);
      }
      fprintf(stdout,"ML(%d)BiCGSTAB_orig  %.3lfGB    %.3lfGB\n",iterative.vec,(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
/*   ------------------------------------------------------------------------- */
      /* mlbicgstab */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)(iterative.vec*parallel.alloc_vector3)); /* g */
      iter_memory+=((double)(iterative.vec*parallel.alloc_vector3)); /* omega */
      iter_memory+=((double)(iterative.vec*parallel.alloc_vector3)); /* q */
      iter_sum_memory+=(sumvector*3.0*3.0*(double)iterative.vec);
      if(iterative.vec>1){
         iter_memory+=((double)((iterative.vec-1)*parallel.alloc_vector3)); /* d */
         iter_sum_memory+=(sumvector*3.0*(double)(iterative.vec-1));
      }
      iter_memory+=((double)parallel.alloc_vector3); /* r */
      iter_memory+=((double)parallel.alloc_vector3); /* u */
      iter_memory+=((double)parallel.alloc_vector3); /* zd */
      iter_memory+=((double)parallel.alloc_vector3); /* zg */
      iter_memory+=((double)parallel.alloc_vector3); /* zomega */
      iter_sum_memory+=(sumvector*3.0*5.0);
      if(iterative.precond!=0){ /* Preconditioning enabled */
         iter_memory+=((double)parallel.alloc_vector3); /* gtilde */
         iter_memory+=((double)parallel.alloc_vector3); /* utilde */
         iter_sum_memory+=(sumvector*3.0*2.0);
      }
      fprintf(stdout,"ML(%d)BiCGSTAB       %.3lfGB    %.3lfGB\n",iterative.vec,(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
/* ------------------------------------------------------------------------- */
      /* mlbicgstab_ss */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)(iterative.vec*parallel.alloc_vector3)); /* g */
      iter_memory+=((double)(iterative.vec*parallel.alloc_vector3)); /* omega */
      iter_memory+=((double)(iterative.vec*parallel.alloc_vector3)); /* q */
      iter_sum_memory+=(sumvector*3.0*3.0*(double)iterative.vec);
      iter_memory+=((double)parallel.alloc_vector3); /* r */
      iter_memory+=((double)parallel.alloc_vector3); /* zg */
      iter_memory+=((double)parallel.alloc_vector3); /* zomega */
      iter_sum_memory+=(sumvector*3.0*3.0);
      fprintf(stdout,"ML(%d)BiCGSTAB_ss    %.3lfGB    %.3lfGB\n",iterative.vec,(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
/* ------------------------------------------------------------------------- */
      /* qmr */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)parallel.alloc_vector3); /* r */
      iter_memory+=((double)parallel.alloc_vector3); /* p */
      iter_memory+=((double)parallel.alloc_vector3); /* q */
      iter_memory+=((double)parallel.alloc_vector3); /* d */
      iter_memory+=((double)parallel.alloc_vector3); /* s */
      iter_memory+=((double)parallel.alloc_vector3); /* omegatilde */
      iter_memory+=((double)parallel.alloc_vector3); /* vtilde */
      iter_sum_memory+=(sumvector*3.0*7.0);
      fprintf(stdout,"QMR                 %.3lfGB    %.3lfGB\n",(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
/* ------------------------------------------------------------------------- */
      /* qmr_sym */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)parallel.alloc_vector3); /* r */
      iter_memory+=((double)parallel.alloc_vector3); /* p */
      iter_memory+=((double)parallel.alloc_vector3); /* pold */
      iter_memory+=((double)parallel.alloc_vector3); /* v */
      iter_memory+=((double)parallel.alloc_vector3); /* vtilde */
      iter_sum_memory+=(sumvector*3.0*5.0);
      fprintf(stdout,"QMR_sym             %.3lfGB    %.3lfGB\n",(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
/* ------------------------------------------------------------------------- */
      /* rbicgstab */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)((iterative.vec+1)*parallel.alloc_vector3)); /* u */
      iter_memory+=((double)((iterative.vec+1)*parallel.alloc_vector3)); /* r */
      iter_sum_memory+=(sumvector*2.0*3.0*(double)(iterative.vec+1));
      iter_memory+=((double)parallel.alloc_vector3); /* rtilde */
      iter_sum_memory+=(sumvector*3.0);
      fprintf(stdout,"RBiCGSTAB(%d)        %.3lfGB    %.3lfGB\n",iterative.vec,(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
/* ------------------------------------------------------------------------- */
      /* tfqmr */
      iter_memory=0.0;
      iter_sum_memory=0.0;
      iter_memory+=((double)parallel.alloc_vector3); /* r */
      iter_memory+=((double)parallel.alloc_vector3); /* rtilde */
      iter_memory+=((double)parallel.alloc_vector3); /* d */
      iter_memory+=((double)parallel.alloc_vector3); /* v */
      iter_memory+=((double)parallel.alloc_vector3); /* ya */
      iter_memory+=((double)parallel.alloc_vector3); /* yb */
      iter_memory+=((double)parallel.alloc_vector3); /* omega */
      iter_sum_memory+=(sumvector*3.0*7.0);
      fprintf(stdout,"TFQMR               %.3lfGB    %.3lfGB\n",(sum_memory+iter_sum_memory)*sizedc/onegb,(base_memory+iter_memory)*sizedc/onegb);fflush(stdout);
   } /* End if(parallel.myrank==0){ Restrict to master */

   mpi_finalise();
   return 0;
}
