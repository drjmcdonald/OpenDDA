/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Read target data from file

   i,j,k,index0*: Loop and array index control variables

   Nd: Nd is the number of lattice sites actually occupied by dipoles

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "target_from_file.h"

void target_from_file(void){

   int i,k,m,position[3],value,shift;
   long int index0i;
   char string[200],compare[200];
   FILE *data;

   if((data=fopen(target.shape,"r"))==NULL){
      reset_string(string);
      sprintf(string,"Could not create the output file \"%s\"",target.shape);
      print_error("File IO error",string,NULL);
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
         /* Read in the x, y, and z material properties */
         for(m=0;m<3;m++){
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' OR ',' to ',' OR '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
               else{compare[k++]=string[i];}
            } while((string[i+1]!=',')&&(string[i+1]!='\n'));
            /* Set composition data */
            value=atoi(compare); /* string to int conversion */
            /* Check the composition data */
            if((value<0)||(value>=refractive_index.number_of_materials)){ /* Composition value MUST be > 0 and < number of materials */
               reset_string(string);
               sprintf(string,"An invalid composition value was encountered in the target input file \"%s\"",target.shape);
               print_error("Target construction error",string,"Composition value >= 0 and < number of materials. Please check the target data input file");
            }
            if(refractive_index.number_of_materials>1){ /* target.composition only exists when the number of constituent materials > 1 */
               shift=(index0i%target.masks_per_integer)*target.mask_size;
               target.composition[m][(int)((double)index0i/(double)target.masks_per_integer)]|=(value<<shift);
            }
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
         }
      }
   }
   fclose(data);
   if(target.output_xyz_vtk==1){
      /* target xyz coordinates are output to a vtk file for subsequent 3D visualisation */
      if(parallel.myrank==0){ /* Restrict to master */
         target_xyz_vtk();
      }
   }
}
