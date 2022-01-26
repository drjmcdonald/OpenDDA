/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Calculates an estmation of the memory requirements
   per node

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "guesstimate_S.h"

attributes_iterative iterative={0};
attributes_target target={0};

int main(void){

   int i,j,k,error,shape,file;
   char in[50],string[200],filename[100];
   double base_memory=0.0,iter_memory,onegb=1024.0*1024.0*1024.0;
   double ii,jj,kk,i0,j0,k0,isqrd,jsqrd,ksqrd,Ksqd,Jsqd,Psqd;
   double sizedc=(double)sizeof(fcomplex);
   FILE *data;

/* ------------------------------------------------------------------------- */
/* K */
   printf("Please enter K [the integer x dimension of the dipole array]\n");
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
         target.K=-1;printf("\nMemory 'guesstimation' error\nArgument for K contained an invalid character!\n\nPlease re-enter\n");
      }
      else{
         target.K=atoi(in); /* Convert to an integer value */
         if(target.K<1){ /* Check for any an invalid value */
            target.K=-1;printf("\nMemory 'guesstimation' error\nK MUST be an integer > 0\n\nPlease re-enter\n");
         }
      }
   }

/* ------------------------------------------------------------------------- */
/* J */
   printf("Please enter J [the integer y dimension of the dipole array]\n");
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
         target.J=-1;printf("\nMemory 'guesstimation' error\nArgument for J contained an invalid character!\n\nPlease re-enter\n");
      }
      else{
         target.J=atoi(in); /* Convert to an integer value */
         if(target.J<1){ /* Check for any an invalid value */
            target.J=-1;printf("\nMemory 'guesstimation' error\nJ MUST be an integer > 0\n\nPlease re-enter\n");
         }
      }
   }

/* ------------------------------------------------------------------------- */
/* P */
   printf("Please enter P [the integer z dimension of the dipole array]\n");
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
         target.P=-1;printf("\nMemory 'guesstimation' error\nArgument for P contained an invalid character!\n\nPlease re-enter\n");
      }
      else{
         target.P=atoi(in); /* Convert to an integer value */
         if(target.P<1){ /* Check for any an invalid value */
            target.P=-1;printf("\nMemory 'guesstimation' error\nP MUST be an integer > 0\n\nPlease re-enter\n");
         }
      }
   }

/* ------------------------------------------------------------------------- */
/* Target shape and number of occupied dipoles */
   printf("Please select the required shape using the appropriate integer\n");
   printf("   0. ellipsoid\n");
   printf("   1. cuboid\n");
   printf("   2. custom from file\n");
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
         shape=-1;printf("\nMemory 'guesstimation' error\nArgument for shape contained an invalid character!\n\nPlease re-enter\n");
      }
      else{
         shape=atoi(in); /* Convert to an integer value */
         if((shape<0)||(shape>2)){ /* Check for any an invalid value */
            shape=-1;printf("\nMemory 'guesstimation' error\nThe shape identifier must be an integer in the range 0..2\n\nPlease re-enter\n");
         }
      }
   }

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
            kk=(double)k-k0;
            ksqrd=(kk*kk)/Psqd;
            if(ksqrd<0.25){
               for(j=0;j<target.J;j++){
                  jj=(double)j-j0;
                  jsqrd=(jj*jj)/Jsqd;
                  if((jsqrd+ksqrd)<0.25){
                     for(i=0;i<target.K;i++){
                        ii=(double)i-i0;
                        isqrd=(ii*ii)/Ksqd;
                        if((isqrd+jsqrd+ksqrd)<0.25){
                           target.Nd++; /* Increment Nd */
                        }
                     }
                  }
               }
            }
         }
      }
      else if(shape==1){ /* Cuboid */
         target.Nd=target.K*target.J*target.P;
      }
   }
   else{ /* Custom target from file */
      printf("Please enter a valid file name in this directory\n");
      file=-1;
      while(file<0){
         error=0; /* Initialise error */
         reset_string(filename); /* Reset string */
         i=0; /* Initialise i */
         scanf("%s",filename); /* Read string in */
         if((data=fopen(filename,"r"))==NULL){error=1;} /* Check if file exists */
         if(error==1){ /* If file does not exist then print an error */
            file=-1;printf("\nMemory 'guesstimation' error\nThe file \"%s\" could not be opened \n\nPlease re-enter\n",filename);
         }
      }
      target.Nd=0; /* Initialise Nd */
      reset_string(string);
      while((fgets(string,sizeof(string),data))!=NULL){
         if(string[0]=='#'){ /* Ignore comments */
            continue;
         }
         else{
            target.Nd++; /* Increment Nd */
         }
      }
      fclose(data);
   }

   find_dft_size(target.K,&target.Kp,&target.Kpp,&target.kf);
   find_dft_size(target.J,&target.Jp,&target.Jpp,&target.jf);
   find_dft_size(target.P,&target.Pp,&target.Ppp,&target.pf);

   target.N=target.K*target.J*target.P; /* N=K*J*P */
   target.N3=3*target.N; /* N3=N*3 */
   target.Ka=target.K+target.kf; /* Ka=K+kf */
   target.Ja=target.J+target.jf; /* Ja=J+jf */
   target.Pa=target.P+target.pf; /* Pa=P+pf */
   target.Na=target.Ka*target.Ja*target.Pa; /* Na=Ka*Ja*Pa */
   target.Qpp=target.Kpp*target.Ppp; /* Qpp=Kpp*Ppp */
   target.Npp=target.Kpp*target.Jpp*target.Ppp; /* Npp=Kpp*Jpp*Ppp */
   target.Nd3=target.Nd*3; /* Nd3=Nd*3 */
   target.Nv=target.K*target.Jpp*target.P; /* Nv=K*Jpp*P */

   /* Main arrays */
   base_memory+=((double)target.Nd3); /* incident_E_field[3][Nd] */
   base_memory+=((double)target.Nv*3.0); /* zero_padded_vector[3][P][Jpp][K] */
   base_memory+=((double)target.Na*6.0); /* interaction_matrix[6][Pa][Ja][Ka] */
   base_memory+=((double)target.Nd3); /* dipole_polarisation[3][Nd] */

   /* DFT scratch arrays */
   base_memory+=((double)target.Kpp); /* xdft[Kpp] */
   base_memory+=((double)target.Jpp); /* ydft[Jpp] */
   base_memory+=((double)target.Ppp); /* zdft[Ppp] */

   /* DFT and FD multiplication scratch arrays */
   base_memory+=((double)target.Qpp*3.0); /* xzscratch[3*Qpp] */

/* ------------------------------------------------------------------------- */
/* Preconditioning */
   printf("Preconditioning enabled OR disabled? 0 OR 1\n");
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
         iterative.precond=-1;printf("\nMemory 'guesstimation' error\nArgument for preconditioning contained an invalid character!\n\nPlease re-enter\n");
      }
      else{
         iterative.precond=atoi(in); /* Convert to an integer value */
         if((iterative.precond<0)||(iterative.precond>1)){ /* Check for any an invalid value */
            iterative.precond=-1;printf("\nMemory 'guesstimation' error\nArgument MUST be an either 0 OR 1\n\nPlease re-enter\n");
         }
      }
   }

   /* Preconditioning and interaction matrix diagonal */
   if(iterative.precond==1){
      base_memory+=((double)target.Nd3); /* point_jacobi[3][Nd] */
   }
   base_memory+=((double)target.Nd3); /* interaction_matrix_diagonal[3][Nd] */

/* ------------------------------------------------------------------------- */
/* Number of starting vectors */
   printf("Number of starting vectors for the iterative schemes:\nRBiCGSTAB & ML(n)BiCGSTAB*? 1..50\n");
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
         iterative.vec=-1;printf("\nMemory 'guesstimation' error\nArgument for the number of starting vectors contained an invalid character!\n\nPlease re-enter\n");
      }
      else{
         iterative.vec=atoi(in); /* Convert to an integer value */
         if((iterative.vec<1)||(iterative.vec>50)){ /* Check for any an invalid value */
            iterative.vec=-1;printf("\nMemory 'guesstimation' error\nArgument MUST be > 0 AND < 50\n\nPlease re-enter\n");
         }
      }
   }

/* ------------------------------------------------------------------------- */
   print_section_title(stderr,"Approximate memory usage");
/* ------------------------------------------------------------------------- */
   printf("\nNumber of occupied sites in the %dx%dx%d array = %ld\n",target.K,target.J,target.P,target.Nd);
   /* bicg */
   iter_memory=0.0;
   iter_memory+=((double)target.Nd3); /* r */
   iter_memory+=((double)target.Nd3); /* rtilde */
   iter_memory+=((double)target.Nd3); /* p */
   iter_memory+=((double)target.Nd3); /* ptilde */
   iter_memory+=((double)target.Nd3); /* q */
   printf("\nBiCG                %.3lfGB\n",(base_memory+iter_memory)*sizedc/onegb);
/* ------------------------------------------------------------------------- */
   /* bicg_sym */
   iter_memory=0.0;
   iter_memory+=((double)target.Nd3); /* r */
   iter_memory+=((double)target.Nd3); /* p */
   printf("BiCG_sym            %.3lfGB\n",(base_memory+iter_memory)*sizedc/onegb);
/* ------------------------------------------------------------------------- */
   /* bicgstab */
   iter_memory=0.0;
   iter_memory+=((double)target.Nd3); /* r */
   iter_memory+=((double)target.Nd3); /* rtilde */
   iter_memory+=((double)target.Nd3); /* p */
   iter_memory+=((double)target.Nd3); /* s */
   iter_memory+=((double)target.Nd3); /* v */
   printf("BiCGSTAB            %.3lfGB\n",(base_memory+iter_memory)*sizedc/onegb);
/* ------------------------------------------------------------------------- */
   /* cg */
   iter_memory=0.0;
   iter_memory+=((double)target.Nd3); /* r */
   iter_memory+=((double)target.Nd3); /* p */
   printf("CG                  %.3lfGB\n",(base_memory+iter_memory)*sizedc/onegb);
/* ------------------------------------------------------------------------- */
   /* cgs */
   iter_memory=0.0;
   iter_memory+=((double)target.Nd3); /* r */
   iter_memory+=((double)target.Nd3); /* rtilde */
   iter_memory+=((double)target.Nd3); /* p */
   iter_memory+=((double)target.Nd3); /* u */
   iter_memory+=((double)target.Nd3); /* q */
   printf("CGS                 %.3lfGB\n",(base_memory+iter_memory)*sizedc/onegb);
/* ------------------------------------------------------------------------- */
   /* mlbicgstab_orig */
   iter_memory=0.0;
   iter_memory+=((double)(iterative.vec*target.Nd3)); /* g */
   iter_memory+=((double)(iterative.vec*target.Nd3)); /* omega */
   iter_memory+=((double)(iterative.vec*target.Nd3)); /* q */
   if(iterative.vec>1){
      iter_memory+=((double)((iterative.vec-1)*target.Nd3)); /* d */
   }
   iter_memory+=((double)target.Nd3); /* r */
   iter_memory+=((double)target.Nd3); /* u */
   iter_memory+=((double)target.Nd3); /* zd */
   iter_memory+=((double)target.Nd3); /* zg */
   iter_memory+=((double)target.Nd3); /* zomega */
   if(iterative.precond!=0){ /* Preconditioning enabled */
      iter_memory+=((double)target.Nd3); /* gtilde */
      iter_memory+=((double)target.Nd3); /* utilde */
   }
   printf("ML(%d)BiCGSTAB_orig  %.3lfGB\n",iterative.vec,(base_memory+iter_memory)*sizedc/onegb);
/* ------------------------------------------------------------------------- */
   /* mlbicgstab */
   iter_memory=0.0;
   iter_memory+=((double)(iterative.vec*target.Nd3)); /* g */
   iter_memory+=((double)(iterative.vec*target.Nd3)); /* omega */
   iter_memory+=((double)(iterative.vec*target.Nd3)); /* q */
   if(iterative.vec>1){
      iter_memory+=((double)((iterative.vec-1)*target.Nd3)); /* d */
   }
   iter_memory+=((double)target.Nd3); /* r */
   iter_memory+=((double)target.Nd3); /* u */
   iter_memory+=((double)target.Nd3); /* zd */
   iter_memory+=((double)target.Nd3); /* zg */
   iter_memory+=((double)target.Nd3); /* zomega */
   if(iterative.precond!=0){ /* Preconditioning enabled */
      iter_memory+=((double)target.Nd3); /* gtilde */
      iter_memory+=((double)target.Nd3); /* utilde */
   }
   printf("ML(%d)BiCGSTAB       %.3lfGB\n",iterative.vec,(base_memory+iter_memory)*sizedc/onegb);
/* ------------------------------------------------------------------------- */
   /* mlbicgstab_ss */
   iter_memory=0.0;
   iter_memory+=((double)(iterative.vec*target.Nd3)); /* g */
   iter_memory+=((double)(iterative.vec*target.Nd3)); /* omega */
   iter_memory+=((double)(iterative.vec*target.Nd3)); /* q */
   iter_memory+=((double)target.Nd3); /* r */
   iter_memory+=((double)target.Nd3); /* zg */
   iter_memory+=((double)target.Nd3); /* zomega */
   printf("ML(%d)BiCGSTAB_ss    %.3lfGB\n",iterative.vec,(base_memory+iter_memory)*sizedc/onegb);
/* ------------------------------------------------------------------------- */
   /* qmr */
   iter_memory=0.0;
   iter_memory+=((double)target.Nd3); /* r */
   iter_memory+=((double)target.Nd3); /* p */
   iter_memory+=((double)target.Nd3); /* q */
   iter_memory+=((double)target.Nd3); /* d */
   iter_memory+=((double)target.Nd3); /* s */
   iter_memory+=((double)target.Nd3); /* omegatilde */
   iter_memory+=((double)target.Nd3); /* vtilde */
   printf("QMR                 %.3lfGB\n",(base_memory+iter_memory)*sizedc/onegb);
/* ------------------------------------------------------------------------- */
   /* qmr_sym */
   iter_memory=0.0;
   iter_memory+=((double)target.Nd3); /* r */
   iter_memory+=((double)target.Nd3); /* p */
   iter_memory+=((double)target.Nd3); /* pold */
   iter_memory+=((double)target.Nd3); /* v */
   iter_memory+=((double)target.Nd3); /* vtilde */
   printf("QMR_sym             %.3lfGB\n",(base_memory+iter_memory)*sizedc/onegb);
/* ------------------------------------------------------------------------- */
   /* rbicgstab */
   iter_memory=0.0;
   iter_memory+=((double)((iterative.vec+1)*target.Nd3)); /* u */
   iter_memory+=((double)((iterative.vec+1)*target.Nd3)); /* r */
   iter_memory+=((double)target.Nd3); /* rtilde */
   printf("RBiCGSTAB(%d)        %.3lfGB\n",iterative.vec,(base_memory+iter_memory)*sizedc/onegb);
/* ------------------------------------------------------------------------- */
   /* tfqmr */
   iter_memory=0.0;
   iter_memory+=((double)target.Nd3); /* r */
   iter_memory+=((double)target.Nd3); /* rtilde */
   iter_memory+=((double)target.Nd3); /* d */
   iter_memory+=((double)target.Nd3); /* v */
   iter_memory+=((double)target.Nd3); /* ya */
   iter_memory+=((double)target.Nd3); /* yb */
   iter_memory+=((double)target.Nd3); /* omega */
   printf("TFQMR               %.3lfGB\n",(base_memory+iter_memory)*sizedc/onegb);

   return 0;
}
