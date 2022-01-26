/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Creates a VTK file containing the target xyz's 

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "target_xyz_vtk.h"

void target_xyz_vtk(void){

   char filename[100],string[200];
   FILE *vtkout;
   int i,j,k;
   long int index0i,index0j,index0k;

   reset_string(filename);
   sprintf(filename,"target_%s_xyz.vtk",target.shape); /* Output file for the target xyz's */
   if((vtkout=fopen(filename,"w"))==NULL){
      reset_string(string);
      sprintf(string,"Could not create the output file \"%s\"",filename);
      print_error("File IO error",string,NULL);
   }

   fprintf(vtkout,"# vtk DataFile Version 2.0\n");
   fprintf(vtkout,"Shape: %s Information: xyz of occupied lattice sites\n",target.shape);
   fprintf(vtkout,"ASCII\n");
   fprintf(vtkout,"DATASET POLYDATA\n");
   fprintf(vtkout,"POINTS %ld double\n",target.Nd);

   for(k=0;k<target.P;k++){
      index0k=k*target.M;
      for(j=0;j<target.J;j++){
         index0j=index0k+j*target.K;
         for(i=0;i<target.K;i++){
            index0i=index0j+i;
            if(target.occupied[(int)((double)index0i/32.0)]&(1<<(index0i%32))){ /* Lattice site is occupied */
               fprintf(vtkout,"%.*g %.*g %.*g\n",DBLP,(double)i,DBLP,(double)j,DBLP,(double)k);
            }
         }
      }
   }

   fprintf(vtkout,"VERTICES %d %ld\n",1,target.Nd+1);
   fprintf(vtkout,"%ld\t",target.Nd);
   for(index0i=0;index0i<target.Nd;index0i++){
      fprintf(vtkout,"%ld\t",index0i);
   }

   fprintf(vtkout,"\nPOINT_DATA %ld\n",target.Nd);
   fprintf(vtkout,"SCALARS xyz_position double %d\n",1);
   fprintf(vtkout,"LOOKUP_TABLE default\n");
   for(k=0;k<target.P;k++){
      index0k=k*target.M;
      for(j=0;j<target.J;j++){
         index0j=index0k+j*target.K;
         for(i=0;i<target.K;i++){
            index0i=index0j+i;
            if(target.occupied[(int)((double)index0i/32.0)]&(1<<(index0i%32))){ /* Lattice site is occupied */
               fprintf(vtkout,"%.*g\n",DBLP,(double)k);
            }
         }
      }
   }

   fclose(vtkout);
}
