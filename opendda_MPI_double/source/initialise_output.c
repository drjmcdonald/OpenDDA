/* Written by James Mc Donald 2006
   drjmcdonald@gmail.com

   Initialise program output

   Copyright (C) 2006 James Mc Donald
   Computational Astrophysics Laboratory
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "initialise_output.h"

void initialise_output(void){

   int i,count;
   char string[200],title0[32]="Written by James Mc Donald 2006",title1[35]="Copyright (C) 2006 James Mc Donald",title2[55]="drjmcdonald@gmail.com",title3[38]="Computational Astrophysics Laboratory",title4[39]="National University of Ireland, Galway",title5[55]="This code is covered by the GNU General Public License";

/* ->->->->->->->->->->->-> Incident wave ->->->->->->->->->->->-> */
   count=1;
   reset_string(string);
   sprintf(string,"OUTPUT_incident_wave.info"); /* Output file for the incident wave properties */
   if((output.incident_wave=fopen(string,"w"))==NULL){
      reset_string(string);
      sprintf(string,"Could not create the output file \"OUTPUT_incident_wave.info\"");
      print_error("File IO error",string,NULL);
   }
   print_out(output.incident_wave,title0,title1,title2);
   print_out(output.incident_wave,title3,title4,title5);
   print_section_title(output.incident_wave,"Incident polarisation state [Normalised]");

   /* Incident polarisations */
   fprintf(output.incident_wave,"\n%d. Incident polarisation e0=[%.*g%+.*gi]x+[%.*g%+.*gi]y\n",count++,DBLP,incident_LF.polarisation[0].dat[0],DBLP,incident_LF.polarisation[0].dat[1],DBLP,incident_LF.polarisation[1].dat[0],DBLP,incident_LF.polarisation[1].dat[1]);
   if(polarisations==1){
      fprintf(output.incident_wave,"\n%d. Calculations for the orthonormal polarisation state e1 were NOT included\n",count++);
   }
   else{
      fprintf(output.incident_wave,"\n%d. Calculations for the orthonormal polarisation state\n   e1=[%.*g%+.*gi]x+[%.*g%+.*gi]y were included\n",count++,DBLP,incident_LF.orthonormal[0].dat[0],DBLP,incident_LF.orthonormal[0].dat[1],DBLP,incident_LF.orthonormal[1].dat[0],DBLP,incident_LF.orthonormal[1].dat[1]);
   }

   /* Incident wavelengths */
   print_section_title(output.incident_wave,"Incident wavelength");
   if(wavelength.flag==0){ /* Use minimum,increment,maximum from control.input */
      fprintf(output.incident_wave,"\n%d. Incident wavelength [micrometers]:\n\n   Used minimum,increment,maximum from control.input\n   minimum=%.*g\n   increment=%.*g\n   maximum=%.*g\n\n",count++,DBLP,wavelength.limits[0],DBLP,wavelength.limits[1],DBLP,wavelength.limits[2]); /* Print description */
      fprintf(output.incident_wave,"   Number of incident wavelength values: %d\n\n",wavelength.list_length); /* Print number of values */
      wavelength.value=wavelength.limits[0];
      fprintf(output.incident_wave,"   0. \t%.*g\n",DBLP,wavelength.value); /* Print values */
      for(i=1;i<wavelength.list_length;i++){
         wavelength.value+=wavelength.limits[1];
         fprintf(output.incident_wave,"   %d. \t%.*g\n",i,DBLP,wavelength.value);
      }
   }
   else if(wavelength.flag==1){ /* Data was read from the file INPUT_incident_wavelength.input */
      fprintf(output.incident_wave,"\n%d. Incident wavelength [micrometers]:\n\n   Data was read in from the file \"INPUT_incident_wavelength.input\"\n\n",count++); /* Print description */
      fprintf(output.incident_wave,"   Number of incident wavelength values: %d\n\n",wavelength.list_length); /* Print number of values */
      for(i=0;i<wavelength.list_length;i++){ /* Print values */
         fprintf(output.incident_wave,"   %d. \t%.*g\n",i,DBLP,wavelength.list[i]);
      }
   }
   fclose(output.incident_wave);

/* ->->->->->->->->->->->-> Effective radius ->->->->->->->->->->->-> */
   count=1;
   reset_string(string);
   sprintf(string,"OUTPUT_target_effective_radius.info"); /* Output file for the effective radius */
   if((output.eff_radius=fopen(string,"w"))==NULL){
      reset_string(string);
      sprintf(string,"Could not create the output file \"OUTPUT_target_effective_radius.info\"");
      print_error("File IO error",string,NULL);
   }
   print_out(output.eff_radius,title0,title1,title2);
   print_out(output.eff_radius,title3,title4,title5);
   print_section_title(output.eff_radius,"Effective radius [micrometers]");

   if(radius.flag==0){ /* Use minimum,increment,maximum from control.input */
      fprintf(output.eff_radius,"\n%d. Effective radius [micrometers]:\n\n   Used minimum,increment,maximum from control.input\n   minimum=%.*g\n   increment=%.*g\n   maximum=%.*g\n\n",count++,DBLP,radius.limits[0],DBLP,radius.limits[1],DBLP,radius.limits[2]); /* Print description */
      fprintf(output.eff_radius,"   Number of effective radius values: %d\n\n",radius.list_length); /* Print number of values */
      radius.value=radius.limits[0];
      fprintf(output.eff_radius,"   0. \t%.*g\n",DBLP,radius.value); /* Print values */
      for(i=1;i<radius.list_length;i++){
         radius.value+=radius.limits[1];
         fprintf(output.eff_radius,"   %d. \t%.*g\n",i,DBLP,radius.value);
      }
   }
   else if(radius.flag==1){ /* Data was read from the file INPUT_target_effective_radius.input */
      fprintf(output.eff_radius,"\n%d. Effective radius [micrometers]:\n\n   Data was read in from the file\n   \"INPUT_target_effective_radius.input\"\n\n",count++); /* Print description */
      fprintf(output.eff_radius,"   Number of effective radius values: %d\n\n",radius.list_length); /* Print number of values */
      for(i=0;i<radius.list_length;i++){ /* Print values */
         fprintf(output.eff_radius,"   %d. \t%.*g\n",i,DBLP,radius.list[i]);
      }
   }
   fclose(output.eff_radius);

/* ->->->->->->->->->->->-> Complex refractive index vs Wavelength ->->->->->->->->->->->-> */
   for(i=0;i<refractive_index.number_of_materials;i++){
      reset_string(string);
      /* Output file for complex refractive index vs wavelength for material i */
      sprintf(string,"OUTPUT_refractive_index_vs_incident_wavelength_%d.info",i);
      if((output.refindex_vs_wave[i]=fopen(string,"w"))==NULL){
         reset_string(string);
         sprintf(string,"Could not create the output file \"OUTPUT_refractive_index_vs_incident_wavelength_%d.info\"",i);
         print_error("File IO error",string,NULL);
      }
      print_out(output.refindex_vs_wave[i],title0,title1,title2);
      print_out(output.refindex_vs_wave[i],title3,title4,title5);
      reset_string(string);
      sprintf(string,"Complex refractive index vs Wavelength [micrometers] for material %d",i);
      print_section_title(output.refindex_vs_wave[i],string);
      fprintf(output.refindex_vs_wave[i],"\nWavelength\t\tComplex refractive index\n");
   }

/* ->->->->->->->->->->->-> Target description ->->->->->->->->->->->-> */
   count=1;
   reset_string(string);
   sprintf(string,"OUTPUT_target_description.info"); /* Output file for the target description */
   if((output.target_descrip=fopen(string,"w"))==NULL){
      reset_string(string);
      sprintf(string,"Could not create the output file \"OUTPUT_target_description.info\"");
      print_error("File IO error",string,NULL);
   }
   print_out(output.target_descrip,title0,title1,title2);
   print_out(output.target_descrip,title3,title4,title5);
   print_section_title(output.target_descrip,"Target description");

   /* Target shape */
   fprintf(output.target_descrip,"\n%d. Target shape description: %s\n",count++,target.shape);

   /* Dipole array dimensions */
   fprintf(output.target_descrip,"\n%d. Number of dipoles in the x direction, K: %d\n",count++,target.K);
   fprintf(output.target_descrip,"\n%d. Number of dipoles in the y direction, J: %d\n",count++,target.J);
   fprintf(output.target_descrip,"\n%d. Number of dipoles in the z direction, P: %d\n",count++,target.P);

   /* Number of dipoles */
   fprintf(output.target_descrip,"\n%d. Number of dipoles in extended rectangular array, N: %ld\n",count++,target.N);
   fprintf(output.target_descrip,"\n%d. Number of dipoles in actual target, Nd: %ld\n",count++,target.Nd);

   /* Material information */   
   fprintf(output.target_descrip,"\n%d. Number of different materials in the target: %d\n",count++,refractive_index.number_of_materials);

   /* Target description data */
   print_section_title(output.target_descrip,"Target parameterisation");
   fprintf(output.target_descrip,"\n Column 0: Wavelength [Micrometers]\n");
   fprintf(output.target_descrip," Column 1: Effective radius [Micrometers]\n");
   fprintf(output.target_descrip," Column 2: Size parameter\n");
   fprintf(output.target_descrip," Column 3: Dipole spacing [Micrometers]\n");
   fprintf(output.target_descrip," Column 4: Dipoles per wavelength\n\n");

/* ->->->->->->->->->->->-> Target orientations ->->->->->->->->->->->-> */
   count=1;
   reset_string(string);
   sprintf(string,"OUTPUT_target_orientations.info"); /* Output file for the target orientations */
   if((output.target_orien=fopen(string,"w"))==NULL){
      reset_string(string);
      sprintf(string,"Could not create the output file \"OUTPUT_target_orientations.info\"");
      print_error("File IO error",string,NULL);
   }
   print_out(output.target_orien,title0,title1,title2);
   print_out(output.target_orien,title3,title4,title5);
   print_section_title(output.target_orien,"Target orientations [Degrees]");
   /* Euler phi */
   if(euler_phi.flag==0){ /* Used minimum,increment,maximum from control.input */
      fprintf(output.target_orien,"\n%d. 1st Euler angle - phi - Azimuthal angle:\n\n   Used minimum,increment,maximum from control.input\n   minimum=%.*g\n   increment=%.*g\n   maximum=%.*g\n\n",count++,DBLP,euler_phi.limits[0],DBLP,euler_phi.limits[1],DBLP,euler_phi.limits[2]); /* Print description */
      fprintf(output.target_orien,"   Number of Euler phi values: %d\n\n",euler_phi.list_length); /* Print number of values */
      euler_phi.value=euler_phi.limits[0];
      fprintf(output.target_orien,"   0. \t%.*g\n",DBLP,euler_phi.value); /* Print values */
      for(i=1;i<euler_phi.list_length;i++){
         euler_phi.value+=euler_phi.limits[1];
         fprintf(output.target_orien,"   %d. \t%.*g\n",i,DBLP,euler_phi.value);
      }
   }
   else if(euler_phi.flag==1){ /* Data was read from the file INPUT_target_orientation_euler_phi.input */
      fprintf(output.target_orien,"\n%d. 1st Euler angle - phi - Azimuthal angle:\n\n   Data was read in from the file\n   \"INPUT_target_orientation_euler_phi.input\"\n\n",count++); /* Print description */
      fprintf(output.target_orien,"   Number of Euler phi values: %d\n\n",euler_phi.list_length); /* Print number of values */
      for(i=0;i<euler_phi.list_length;i++){ /* Print values */
         fprintf(output.target_orien,"   %d. \t%.*g\n",i,DBLP,euler_phi.list[i]);
      }
   }
   /* Euler theta */
   if(euler_theta.flag==0){ /* Used minimum,increment,maximum from control.input */
      fprintf(output.target_orien,"\n%d. 2st Euler angle - theta - Polar/zenith angle:\n\n   N.B. Distributed via cos(theta)\n   Used minimum,increment,maximum from control.input\n   minimum=%.*g\n   increment=%.*g\n   maximum=%.*g\n\n",count++,DBLP,euler_theta.limits[0],DBLP,euler_theta.limits[1],DBLP,euler_theta.limits[2]); /* Print description */
      fprintf(output.target_orien,"   Number of Euler theta values: %d\n\n",euler_theta.list_length); /* Print number of values */
      euler_theta.value=euler_theta.limits[0];
      fprintf(output.target_orien,"   0. \t%+.*g\t\t%.*g\n",DBLP,euler_theta.value,DBLP,radians_to_degrees(acos(euler_theta.value))); /* Print values */
      for(i=1;i<euler_theta.list_length;i++){
         euler_theta.value+=euler_theta.limits[1];
         fprintf(output.target_orien,"   %d. \t%+.*g\t\t%.*g\n",i,DBLP,euler_theta.value,DBLP,radians_to_degrees(acos(euler_theta.value)));
      }
   }
   else if(euler_theta.flag==1){ /* Data was read from the file INPUT_target_orientation_euler_theta.input */
      fprintf(output.target_orien,"\n%d. 2st Euler angle - theta - Polar/zenith angle:\n\n   N.B. Distributed via cos(theta)\n   Data was read in from the file\n   \"INPUT_target_orientation_euler_theta.input\"\n\n",count++); /* Print description */
      fprintf(output.target_orien,"   Number of Euler theta values: %d\n\n",euler_theta.list_length); /* Print number of values */
      for(i=0;i<euler_theta.list_length;i++){ /* Print values */
         fprintf(output.target_orien,"   %d. \t%+.*g\t\t%.*g\n",i,DBLP,euler_theta.list[i],DBLP,radians_to_degrees(acos(euler_theta.list[i])));
      }
   }   
   /* Euler psi */
   if(euler_psi.flag==0){ /* Used minimum,increment,maximum from control.input */
      fprintf(output.target_orien,"\n%d. 3rd Euler angle - psi - Azimuthal angle:\n\n   Used minimum,increment,maximum from control.input\n   minimum=%.*g\n   increment=%.*g\n   maximum=%.*g\n\n",count++,DBLP,euler_psi.limits[0],DBLP,euler_psi.limits[1],DBLP,euler_psi.limits[2]); /* Print description */
      euler_psi.value=euler_psi.limits[0];
      fprintf(output.target_orien,"   Number of Euler psi values: %d\n\n",euler_psi.list_length); /* Print number of values */
      fprintf(output.target_orien,"   0. \t%.*g\n",DBLP,euler_psi.value); /* Print values */
      for(i=1;i<euler_psi.list_length;i++){
         euler_psi.value+=euler_psi.limits[1];
         fprintf(output.target_orien,"   %d. \t%.*g\n",i,DBLP,euler_psi.value);
      }
   }
   else if(euler_psi.flag==1){ /* Data was read from the file INPUT_target_orientation_euler_psi.input */
      fprintf(output.target_orien,"\n%d. 3rd Euler angle - psi - Azimuthal angle:\n\n   Data was read in from the file\n   \"INPUT_target_orientation_euler_psi.input\"\n\n",count++); /* Print description */
      fprintf(output.target_orien,"   Number of Euler psi values: %d\n\n",euler_psi.list_length); /* Print number of values */
      for(i=0;i<euler_psi.list_length;i++){ /* Print values */
         fprintf(output.target_orien,"   %d. \t%.*g\n",i,DBLP,euler_psi.list[i]);
      }
   }
   fclose(output.target_orien);

/* ->->->->->->->->->->->-> Scattering angles ->->->->->->->->->->->-> */
   count=1;
   reset_string(string);
   sprintf(string,"OUTPUT_scattering_angles.info"); /* Output file for the scattering angles */
   if((output.scatt_angle=fopen(string,"w"))==NULL){
      reset_string(string);
      sprintf(string,"Could not create the output file \"OUTPUT_scattering_angles.info\"");
      print_error("File IO error",string,NULL);
   }
   print_out(output.scatt_angle,title0,title1,title2);
   print_out(output.scatt_angle,title3,title4,title5);
   print_section_title(output.scatt_angle,"Scattering angles [Degrees]");
   
   /* phi */
   if(phi.flag==0){ /* Used minimum,increment,maximum from control.input */
      fprintf(output.scatt_angle,"\n%d. Scattering azimuthal angle phi:\n\n   Used minimum,increment,maximum from control.input\n   minimum=%.*g\n   increment=%.*g\n   maximum=%.*g\n\n",count++,DBLP,phi.limits[0],DBLP,phi.limits[1],DBLP,phi.limits[2]); /* Print description */
      phi.value=phi.limits[0];
      fprintf(output.scatt_angle,"   Number of phi values: %d\n\n",phi.list_length); /* Print number of values */
      fprintf(output.scatt_angle,"   0. \t%.*g\n",DBLP,phi.value); /* Print values */
      for(i=1;i<phi.list_length;i++){
         phi.value+=phi.limits[1];
         fprintf(output.scatt_angle,"   %d. \t%.*g\n",i,DBLP,phi.value);
      }
   }
   else if(phi.flag==1){ /* Data was read from the file INPUT_scattering_direction_phi.input */
      fprintf(output.scatt_angle,"\n%d. Scattering azimuthal angle phi:\n\n   Data was read in from the file\n   \"INPUT_scattering_direction_phi.input\"\n\n",count++); /* Print description */
      fprintf(output.scatt_angle,"   Number of phi values: %d\n\n",phi.list_length); /* Print number of values */
      for(i=0;i<phi.list_length;i++){ /* Print values */
         fprintf(output.scatt_angle,"   %d. \t%.*g\n",i,DBLP,phi.list[i]);
      }
   }
   /* theta */
   if(theta.flag==0){ /* Use minimum,increment,maximum from control.input */
      fprintf(output.scatt_angle,"\n%d. Scattering polar/zenith angle theta:\n\n   Used minimum,increment,maximum from control.input\n   minimum=%.*g\n   increment=%.*g\n   maximum=%.*g\n\n",count++,DBLP,theta.limits[0],DBLP,theta.limits[1],DBLP,theta.limits[2]); /* Print description */
      fprintf(output.scatt_angle,"   Number of theta values: %d\n\n",theta.list_length); /* Print number of values */
      theta.value=theta.limits[0];
      fprintf(output.scatt_angle,"   0. \t%.*g\n",DBLP,theta.value); /* Print values */
      for(i=1;i<theta.list_length;i++){
         theta.value+=theta.limits[1];
         fprintf(output.scatt_angle,"   %d. \t%.*g\n",i,DBLP,theta.value);
      }
   }
   else if(theta.flag==1){ /* Data was read from the file INPUT_scattering_direction_theta.input */
      fprintf(output.scatt_angle,"\n%d. Scattering polar/zenith angle theta:\n\n   Data was read in from the file\n   \"INPUT_scattering_direction_theta.input\"\n\n",count++); /* Print description */
      fprintf(output.scatt_angle,"   Number of theta values: %d\n\n",theta.list_length); /* Print number of values */
      for(i=0;i<theta.list_length;i++){ /* Print values */
         fprintf(output.scatt_angle,"   %d. \t%.*g\n",i,DBLP,theta.list[i]);
      }
   }
   fclose(output.scatt_angle);

/* ->->->->->->->->->->->-> Iterative parameters ->->->->->->->->->->->-> */
   count=1;
   reset_string(string);
   sprintf(string,"OUTPUT_iterative.info"); /* Output file for the iterative parameters */
   if((output.iter=fopen(string,"w"))==NULL){
      reset_string(string);
      sprintf(string,"Could not create the output file \"OUTPUT_iterative.info\"");
      print_error("File IO error",string,NULL);
   }
   print_out(output.iter,title0,title1,title2);
   print_out(output.iter,title3,title4,title5);
   print_section_title(output.iter,"Iterative parameters");
   /* Iterative scheme - Choice */
   reset_string(string);
   if((strcmp(iterative.scheme,"bicg")==0)){
      sprintf(string,"[bicg]\n   BiConjugate-Gradients");
   }
   if((strcmp(iterative.scheme,"bicg_sym")==0)){
      sprintf(string,"[bicg_sym]\n   BiConjugate-Gradients for symmetric systems");
   }
   if((strcmp(iterative.scheme,"bicgstab")==0)){
      sprintf(string,"[bicgstab]\n   Stabilised version of the BiConjugate-Gradients");
   }
   if((strcmp(iterative.scheme,"cg")==0)){
      sprintf(string,"[cg]\n   Conjugate-Gradients");
   }
   if((strcmp(iterative.scheme,"cgs")==0)){
      sprintf(string,"[cgs]\n   Conjugate-Gradients Squared");
   }
   if((strcmp(iterative.scheme,"mlbicgstab_orig")==0)){
      sprintf(string,"[mlbicgstab_orig]\n   A variant of BiCGSTAB based on multiple Lanczos starting\n   vectors (Author's original algorithm)");
   }
   if((strcmp(iterative.scheme,"mlbicgstab")==0)){
      sprintf(string,"[mlbicgstab]\n   A variant of BiCGSTAB based on multiple Lanczos starting\n   vectors (Author's reformulated algorithm)");
   }
   if((strcmp(iterative.scheme,"mlbicgstab_ss")==0)){
      sprintf(string,"[mlbicgstab_ss]\n   A variant of BiCGSTAB based on multiple Lanczos starting\n   vectors (Author's space saving algorithm)");
   }
   if((strcmp(iterative.scheme,"qmr")==0)){
      sprintf(string,"[qmr]\n   Quasi-minimal residual with coupled two-term recurrences");
   }
   if((strcmp(iterative.scheme,"qmr_sym")==0)){
      sprintf(string,"[qmr_sym]\n   Quasi-minimal residual for symmetric systems");
   }
   if((strcmp(iterative.scheme,"rbicgstab")==0)){
      sprintf(string,"[rbicgstab]\n   Restarted, stabilised version of the BiConjugate-Gradients");
   }
   if((strcmp(iterative.scheme,"tfqmr")==0)){
      sprintf(string,"[tfqmr]\n   Transpose-free quasi-minimal residual");
   }
   fprintf(output.iter,"\n%d. Iterative scheme: %s\n",count++,string);
   /* Iterative scheme - Number of starting vectors */
   if(((strcmp(iterative.scheme,"mlbicgstab_orig")==0))||
      ((strcmp(iterative.scheme,"mlbicgstab")==0))||
      ((strcmp(iterative.scheme,"mlbicgstab_ss")==0))||
      ((strcmp(iterative.scheme,"rbicgstab")==0))){
      fprintf(output.iter,"\n%d. Number of starting vectors: %d\n",count++,iterative.vec);
   }
   /* Iterative scheme - Convergence tolerance */
   fprintf(output.iter,"\n%d. Convergence tolerance: %.*g\n",count++,DBLP,iterative.tolerance);
   /* Iterative scheme - Breakdown tolerance */
   fprintf(output.iter,"\n%d. Breakdown tolerance: %.*g\n",count++,DBLP,iterative.breakdown);
   /* Iterative scheme - Maximum number of iterations */
   fprintf(output.iter,"\n%d. Maximum number of iterations: %d\n",count++,iterative.maximum);
   /* Iterative scheme - Initial guess */
   if(iterative.initial_guess==0){
      fprintf(output.iter,"\n%d. Initial guess for the unknown polarisations x=A^{-1}b:\n   Set x=0\n",count++);
   }
   else if(iterative.initial_guess==1){
      fprintf(output.iter,"\n%d. Initial guess for the unknown polarisations x=A^{-1}b:\n   Set x=1\n",count++);
   }
   else if(iterative.initial_guess==2){
      fprintf(output.iter,"\n%d. Initial guess for the unknown polarisations x=A^{-1}b:\n   Set x=b i.e. incident electric-field\n",count++);
   }
   else if(iterative.initial_guess==3){
      fprintf(output.iter,"\n%d. Initial guess for the unknown polarisations x=A^{-1}b:\n   Set x=(1/alpha) i.e. (1/polarisability)\n",count++);
   }
   /* Iterative scheme - Preconditioning */
   if(iterative.precond==0){
      fprintf(output.iter,"\n%d. Preconditioning: No preconditioning was used\n",count++);
   }
   else if(iterative.precond==1){
      fprintf(output.iter,"\n%d. Preconditioning: Point-Jacobi preconditioning was used\n",count++);
   }

/* ->->->->->->->->->->->-> Program execution information ->->->->->->->->->->->-> */
   if(timing.enabled){
      reset_string(string);
      sprintf(string,"OUTPUT_program_execution.info"); /* Output file for the timing information */
      if((output.time=fopen(string,"w"))==NULL){
         reset_string(string);
         sprintf(string,"Could not create the output file \"OUTPUT_program_execution.info\"");
         print_error("File IO error",string,NULL);
      }
      print_out(output.time,title0,title1,title2);
      print_out(output.time,title3,title4,title5);
      print_section_title(output.time,"MPI information");
      fprintf(output.time,"\n  Number of processors: %d\n",parallel.np);
      fprintf(output.time,"  Distributed tensor component transpose algorithm: %d\n",parallel.tensor_transpose);
      fprintf(output.time,"  Distributed vector component transpose algorithm: %d\n",parallel.padded_transpose);
      print_section_title(output.time,"Runtime, iterative kernel & DFT engine timing information");
      fprintf(output.time,"\n  Note: The number in square brackets after each time specifies\n        the number of executions in the average per processor\n"); 
   }
}
