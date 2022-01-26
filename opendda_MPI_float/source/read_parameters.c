/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Reads in the user defined parameters using control.input

   Copyright (C) 2006 James Mc Donald
   Computational Astrophysics Laboratory
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "read_parameters.h"

void read_parameters(void){

   int i,j,k,m,a,b;
   fcomplex normalise;
   char string[200]={0},compare[200]={0},temp[60]={0};
   FILE *input,*data;

   /* Set incident direction in the lab frame to +z */
   incident_LF.n_inc[0]=0.0F;
   incident_LF.n_inc[1]=0.0F;   
   incident_LF.n_inc[2]=1.0F;

   if((input=fopen("control.input","r"))==NULL){
      print_error("File IO error","Could not open the file \"control.input\"",NULL);
   }

   while((fgets(string,sizeof(string),input))!=NULL){
      if(string[0]=='#'){ /* Ignore comments */
         continue;
      }
      else{
         i=-1;reset_string(compare);
         do{ /* Copy everything up to and including the '=' to compare string */
            i++;compare[i]=string[i];
         } while((compare[i]!='=')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
/* *********************************************************************************************************************************************
Wavelengths - [micrometers]
If "wavelength.input" exists read in data else use minimum,increment,maximum in "control.input" */
         if(strcmp(compare,"Wavelengths=")==0){
            if((data=fopen("INPUT_incident_wavelength.input","r"))!=NULL){
               /* Read data from file into the structure wavelength */
               read_from_file(&data,&wavelength,"INPUT_incident_wavelength.input","Wavelength");
            }
            else{ /* Use minimum,increment,maximum in "control.input" */
               for(m=0;m<3;m++){
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
                  k=0;reset_string(compare);
                  do{ /* Copy everything after '=' OR ',' to ',' OR '\n' to compare string,
                     ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
                     if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
                     else{compare[k++]=string[i];}
                  } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
                  wavelength.limits[m]=(float)atof(compare); /* string to float conversion */
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               }
               /* Check wavelength values */
               if((wavelength.limits[0]<0.0F)||(wavelength.limits[0]<FLT_EPSILON)){
                  print_error("Wavelength description error","Wavelength minimum MUST be > 0","Please check the \'control.input\' file");
               }
               if(wavelength.limits[1]<0.0F){
                  print_error("Wavelength description error","Wavelength increment MUST be >= 0","Please check the \'control.input\' file");
               }
               if((wavelength.limits[2]<0.0F)||(wavelength.limits[2]<FLT_EPSILON)){
                  print_error("Wavelength description error","Wavelength maximum MUST be > 0","Please check the \'control.input\' file");
               }
               if(wavelength.limits[0]>wavelength.limits[2]){
                  print_error("Wavelength description error","Minimum wavelength > Maximum wavelength","Please check the \'control.input\' file");
               }
            }
         }
/* *********************************************************************************************************************************************
Incident polarisation - Normalised Jones vector
e0=[x_{real}+x_{imag}i,y_{real}+y_{imag}i,0+0i]
e1=[0,0,1] CROSS conj(e0) */
         else if(strcmp(compare,"Incident polarisation e0=")==0){
            /* Read in complex data i.e. x, y and z components of e0 */
            for(a=0;a<3;a++){ /* a=0 x,a=1 y,a=2 z */
               for(b=0;b<2;b++){ /* b=0 real component,b=1 imaginary component */
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
                  k=0;reset_string(compare);
                  do{ /* Copy everything after '=' to '+' OR '-' OR ',' OR '\n' to compare string,
                     ignoring punctuation (except fullstop and signs) and spaces, tabs, newlines etc.. */
                     if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-')&&(string[i]!='+'))){continue;}
                     else{compare[k++]=string[i];}
                  } while((string[i+1]!='+')&&(string[i+1]!='-')&&(string[i+1]!=',')&&(string[i+1]!='I')&&(string[i+1]!='i')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
                  incident_LF.polarisation[a].dat[b]=(float)atof(compare); /* string to float conversion */
                  while(((ispunct(string[i+1])!=0)&&(string[i+1]!='-'))||(string[i+1]=='i')||(string[i+1]=='I')){i++;} /* Skip 'i', 'I' and all punctuation (except '-') */
               }
            }
            /* Check the z component of the incident polarisation vector [SHOULD BE ZERO] */
            if(fcomplex_abs(incident_LF.polarisation[2])>FLT_EPSILON){
               print_error("Incident polarisation description error","Incident polarisation vector should not have a `z' component in the lab frame","Please check the \'control.input\' file");
            }
            /* Normalise the incident polarisation vector i.e. 1/sqrt(x^2+y^2+z^2) */
            normalise=fcomplex_add(fcomplex_mul(incident_LF.polarisation[0],fcomplex_conj(incident_LF.polarisation[0])),
                                    fcomplex_mul(incident_LF.polarisation[1],fcomplex_conj(incident_LF.polarisation[1])));
            /* Check normalisation value */
            if(normalise.dat[0]<FLT_EPSILON){
               print_error("Incident polarisation description error","The incident polarisation vector e0 CANNOT be the zero vector","Please check the \'control.input\' file");
            }
            normalise.dat[0]=1.0F/sqrtf(normalise.dat[0]);
            incident_LF.polarisation[0]=fcomplex_scale(incident_LF.polarisation[0],normalise.dat[0]);
            incident_LF.polarisation[1]=fcomplex_scale(incident_LF.polarisation[1],normalise.dat[0]);
            /* Read in whether to include calculations for the orthonormal polarisation state */
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
               else{compare[k++]=(char)tolower(string[i]);}
            } while(string[i+1]!='\n');
            if(strcmp(compare,"yes")==0){ /* Include the orthonormal polarisation state */
               polarisations=2; /* 2 polarisation states, incident & orthonormal */
               /* Calculate the orthonormal polarisation state */
               /* =======================================================================
                  a=[x0,y0,z0],b=[x1,y1,z1]
                  a CROSS b = [y0z1-y1z0,x1z0-x0z1,x0y1-x1y0]

                  e1=z CROSS conj(e0)
                  a=[0,0,z0]=[0,0,1]
                  b=[conj(x1),conj(y1),0]
                  e1 = a CROSS b = [-conj(y1)z0,conj(x1)z0,0]=[-conj(y1),conj(x1),0] */ 
               incident_LF.orthonormal[0]=fcomplex_mul(fcomplex_conj(incident_LF.polarisation[1]),minusone); /* x-component */
               incident_LF.orthonormal[1]=fcomplex_conj(incident_LF.polarisation[0]); /* y-component */
            }
            else if(strcmp(compare,"no")==0){ /* Exclude the orthonormal polarisation state */
               polarisations=1; /* Only 1 incident polarisation state */
            }
            else{ /* Invalid argument for "Orthonormal polarisation e1=" */
               reset_string(string);
               sprintf(string,"Invalid argument \"%s\" for \"%s\" in \'control.input\'",compare,temp);
               reset_string(compare);
               sprintf(compare,"The argument MUST be \'yes\' OR \'no\'");
               print_error("Incident polarisation description error",string,compare);
            }
         }
/* *********************************************************************************************************************************************
Effective radii - [micrometers]
If "effective_radius.input" exists read in data else use minimum,increment,maximum in "control.input" */
         else if(strcmp(compare,"Effective radii=")==0){
            if((data=fopen("INPUT_target_effective_radius.input","r"))!=NULL){
               /* Read data from file into the structure radius */
               read_from_file(&data,&radius,"INPUT_target_effective_radius.input","Effective radius");
            }
            else{ /* Use minimum,increment,maximum in "control.input" */
               for(m=0;m<3;m++){
                  k=0;reset_string(compare);
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
                  do{ /* Copy everything after '=' OR ',' to ',' OR '\n' to compare string,
                     ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
                     if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
                     else{compare[k++]=string[i];}
                  } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
                  radius.limits[m]=(float)atof(compare); /* string to float conversion */
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               }
               /* Check effective radii values */
               if((radius.limits[0]<0.0F)||(radius.limits[0]<FLT_EPSILON)){
                  print_error("Effective radius description error","Effective radius minimum MUST be > 0","Please check the \'control.input\' file");
               }
               if(radius.limits[1]<0.0F){
                  print_error("Effective radius description error","Effective radius increment MUST be >= 0","Please check the \'control.input\' file");
               }
               if((radius.limits[2]<0.0F)||(radius.limits[2]<FLT_EPSILON)){
                  print_error("Effective radius description error","Effective radius maximum MUST be > 0","Please check the \'control.input\' file");
               }
               if(radius.limits[0]>radius.limits[2]){
                  print_error("Effective radius description error","Minimum effective radius > Maximum effective radius","Please check the \'control.input\' file");
               }
            }
         }
/* *********************************************************************************************************************************************
Target shape and target construction parameters */
         else if(strcmp(compare,"Shape=")==0){
            /* Extract the shape information */
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-')&&(string[i]!='_'))){continue;}
               else{compare[k++]=(char)tolower(string[i]);}
            } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
            reset_string(target.shape);
            strncpy(target.shape,compare,strlen(compare));
            /* Check to see if the shape identifier is recognised OR is a valid file name */
            if((strcmp(target.shape,"ellipsoid")!=0)&&
               (strcmp(target.shape,"cuboid")!=0)){ /* If not a recognised shape identifier */
               if((data=fopen(target.shape,"r"))!=NULL){ /* Target data file exists */
                  target.from_file=1; /* Set flag to indicate that the target data is to be read from a file */
                  fclose(data);
               }
               else{ /* The specified identifier is NOT recognised AND is NOT a valid file name */
                  reset_string(string);
                  sprintf(string,"\'%s\' is not a recognised shape and is not the name of an accessible data file",target.shape);
                  print_error("Target shape description error",string,"Please check the \'control.input\' file");
               }
            }
            /* Read in the 6 target construction parameters from the control.input file */
            j=0;
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            for(m=0;m<6;m++){ /* Read in the 6 target construction parameters */
               k=0;reset_string(compare);
               do{ /* Copy everything after '=' OR ',' to ',' OR '\n' to compare string,
                  ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
                  if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
                  else{compare[k++]=(char)toupper(string[i]);}
               } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
               if(strcmp(compare,"NULL")!=0){
                  j++;
                  target.construction_parameters[m]=atoi(compare); /* string to int conversion */
                  if(target.construction_parameters[m]<=0){
                     reset_string(string);
                     sprintf(string,"Target construction parameter %d MUST be > 0",m);
                     print_error("Target shape description error",string,"Please check the \'control.input\' file");
                  }
               }
               while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            }
            /* Check that the shape and construction parameters are consistent */
            if((strcmp(target.shape,"ellipsoid")==0)||(strcmp(target.shape,"cuboid")==0)){
               if(j<3){ /* The shape must be accompanied by the first 3 target construction parameters */
                  reset_string(string);
                  sprintf(string,"The \"%s\" shape MUST be accompanied by the first 3 target construction parameters",target.shape);
                  print_error("Target construction parameter error",string,"Please check the \'control.input\' file");
               }
            }
         }
/* *********************************************************************************************************************************************
The number of dielectric materials for the target [MAX=16] and the corresponding complex
refractive indices, m=n+ki, for the target, 0 to (number_of_materials-1) */
         else if(strcmp(compare,"Number of dielectric materials for the target=")==0){
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
               else{compare[k++]=string[i];}
            } while(string[i+1]!='\n');
            refractive_index.number_of_materials=atoi(compare); /* string to int conversion */
            /* Check number_of_materials value */
            if(refractive_index.number_of_materials<1){
               print_error("Target description error","The number of dielectric materials for the target MUST be > 0","Please check the \'control.input\' file");
            }

            /* Allocate memory to store various parameters relevant to the specification of the refractive index */
            /* Allocate memory to store the current complex refractive index values for materials 0 to (number_of_materials-1) */
            fcomplex_malloc(&refractive_index.value,(size_t)refractive_index.number_of_materials);
             /* Allocate memory to store the flag which indicates whether the data for material 'j' has been read from file (flag=1) or the value from the file "contol.input" is being used (flag=0) */
            int_calloc(&refractive_index.flag,(size_t)refractive_index.number_of_materials);
            /* Allocate memory to store the length of the list that has been read in from file for the material 'j' */
            int_calloc(&refractive_index.list_length,(size_t)refractive_index.number_of_materials);
            /* Allocate memory to store the output file pointers for the various materials */
            file_malloc_ptr(&output.refindex_vs_wave,(size_t)refractive_index.number_of_materials);
            /* Allocate memory for the array pointers to the wavelength and complex refractive index arrays */
            float_malloc_ptr(&refractive_index.wavelengths,(size_t)refractive_index.number_of_materials);
            fcomplex_malloc_ptr(&refractive_index.refractive_indices,(size_t)refractive_index.number_of_materials);

            for(j=0;j<refractive_index.number_of_materials;j++){ /* Read in the j different complex refractive indices */
               reset_string(string);
               reset_string(compare);
               fgets(string,sizeof(string),input);
               i=-1;
               do{ /* Copy everything up to and including the '=' to compare string */
                  i++;compare[i]=string[i];
               } while((compare[i]!='=')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
               reset_string(temp);
               sprintf(temp,"Material %d=",j);
               if(strcmp(compare,temp)==0){
                  reset_string(temp);
                  sprintf(temp,"INPUT_refractive_index_vs_wavelength_material_%d.input",j);
                  if((data=fopen(temp,"r"))!=NULL){ /* Check if the appropriate input file exists */
                     /* Read data for material 'j' from file into the structure refractive_index */
                     read_from_file_refractive_index_vs_wavelength(&data,&refractive_index,temp,j);
                  }
                  else{ /* Read complex refractive index for material 'j' from "control.input" */
                     while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
                     /* Read in complex data */
                     for(a=0;a<2;a++){ /* a=0 real component,a=1 imaginary component */
                        k=0;reset_string(compare);
                        do{ /* Copy everything after '=' to '+' OR '-' OR ',' OR '\n' to compare string,
                           ignoring punctuation (except fullstop and signs) and spaces, tabs, newlines etc.. */
                           if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-')&&(string[i]!='+'))){continue;}
                           else{compare[k++]=(char)toupper(string[i]);}
                        } while((string[i+1]!='+')&&(string[i+1]!='-')&&(string[i+1]!=',')&&(string[i+1]!='I')&&(string[i+1]!='i')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
                        while(((ispunct(string[i+1])!=0)&&(string[i+1]!='-'))||(string[i+1]=='i')||(string[i+1]=='I')){i++;} /* Skip 'i', 'I' and all punctuation (except '-') */
                        if(strcmp(compare,"NULL")==0){ /* Print error: No valid input file or value in 'control.input' was given */
                           reset_string(string);
                           reset_string(compare);
                           sprintf(string,"No input file for material '%d' AND the value in \'control.input\' is set to \'NULL\'",j);
                           sprintf(compare,"Please check the \'control.input\' file OR create the appropriate input file");
                           print_error("Complex refractive index description error",string,compare);
                        }
                        refractive_index.value[j].dat[a]=(float)atof(compare); /* string to float conversion */
                     }
                  }
               } /* End strcmp(compare,"Material j=") */
               else{
                  reset_string(string);
                  reset_string(compare);
                  sprintf(string,"No input file for material '%d' AND the line for \"Material %d=\" in \'control.input\' is corrupt or missing",j,j);
                  sprintf(compare,"Please check the \'control.input\' file OR create the appropriate input file");
                  print_error("Complex refractive index description error",string,compare);
               }
            } /* End for j=0 to (number_of_materials-1) */
            do{ /* Skip the rest of the "Material j=" input lines */
               fgets(string,sizeof(string),input);
            } while((char)toupper(string[0])=='M');
         }
/* *********************************************************************************************************************************************
Target orientations */
/* Euler angle - phi
If "target_orientation_euler_phi.input" exists read in data else use minimum,increment,maximum in "control.input" */
         else if(strcmp(compare,"Euler phi=")==0){
            if((data=fopen("INPUT_target_orientation_euler_phi.input","r"))!=NULL){
               /* Read data from file into the structure euler_phi */
               read_from_file(&data,&euler_phi,"INPUT_target_orientation_euler_phi.input","Target orientation Euler phi");
            }
            else{ /* Use minimum,increment,maximum in "control.input" */
               while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               for(m=0;m<3;m++){ /* Use minimum,increment,maximum in "control.input" */
                  k=0;reset_string(compare);
                  do{ /* Copy everything after '=' OR ',' to ',' OR '\n' to compare string,
                     ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
                     if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
                     else{compare[k++]=string[i];}
                  } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
                  euler_phi.limits[m]=(float)atof(compare); /* string to float conversion */
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               }
               /* Check phi values */
               if(euler_phi.limits[0]<0.0F){
                  print_error("Target orientation description error","Euler phi minimum MUST be >= 0","Please check the \'control.input\' file");
               }
               if(euler_phi.limits[1]<0.0F){
                  print_error("Target orientation description error","Euler phi increment MUST be >= 0","Please check the \'control.input\' file");
               }
               if(euler_phi.limits[2]<0.0F){
                  print_error("Target orientation description error","Euler phi maximum MUST be >= 0","Please check the \'control.input\' file");
               }
               if(euler_phi.limits[0]>euler_phi.limits[2]){
                  print_error("Target orientation description error","Minimum Euler phi > Maximum Euler phi","Please check the \'control.input\' file");
               }
            }
         }
/* Euler angle - theta
If "target_orientation_euler_theta.input" exists read in data else use minimum,increment,maximum in "control.input" */
         else if(strcmp(compare,"Euler theta=")==0){
            if((data=fopen("INPUT_target_orientation_euler_theta.input","r"))!=NULL){
               /* Read data from file into the structure euler_theta */
               read_from_file(&data,&euler_theta,"INPUT_target_orientation_euler_theta.input","Target orientation Euler theta");
            }
            else{ /* Use minimum,increment,maximum in "control.input" */
               while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               for(m=0;m<3;m++){ /* Read in minimum,increment,maximum */
                  k=0;reset_string(compare);
                  do{ /* Copy everything after '=' OR ',' to ',' OR '\n' to compare string,
                     ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
                     if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
                     else{compare[k++]=string[i];}
                  } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
                  euler_theta.limits[m]=(float)atof(compare); /* string to float conversion */
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               }
               /* Check theta values */
               if(euler_theta.limits[1]<0.0F){
                  print_error("Target orientation description error","Euler theta increment MUST be >= 0","Please check the \'control.input\' file");
               }
               if(euler_theta.limits[0]>euler_theta.limits[2]){
                  print_error("Target orientation description error","Minimum Euler theta > Maximum Euler theta","Please check the \'control.input\' file");
               }
            }
         }
/* Euler angle - psi
If "target_orientation_euler_psi.input" exists read in data else use minimum,increment,maximum in "control.input" */
         else if(strcmp(compare,"Euler psi=")==0){
            if((data=fopen("INPUT_target_orientation_euler_psi.input","r"))!=NULL){
               /* Read data from file into the structure euler_psi */
               read_from_file(&data,&euler_psi,"INPUT_target_orientation_euler_psi.input","Target orientation Euler psi");
            }
            else{ /* Use minimum,increment,maximum in "control.input" */
               while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               for(m=0;m<3;m++){ /* Read in minimum,increment,maximum */
                  k=0;reset_string(compare);
                  do{ /* Copy everything after '=' OR ',' to ',' OR '\n' to compare string,
                     ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
                     if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
                     else{compare[k++]=string[i];}
                  } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
                  euler_psi.limits[m]=(float)atof(compare); /* string to float conversion */
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               }
               /* Check psi values */
               if(euler_psi.limits[0]<0.0F){
                  print_error("Target orientation description error","Euler psi minimum MUST be >= 0","Please check the \'control.input\' file");
               }
               if(euler_psi.limits[1]<0.0F){
                  print_error("Target orientation description error","Euler psi increment MUST be >= 0","Please check the \'control.input\' file");
               }
               if(euler_psi.limits[2]<0.0F){
                  print_error("Target orientation description error","Euler psi maximum MUST be >= 0","Please check the \'control.input\' file");
               }
               if(euler_psi.limits[0]>euler_psi.limits[2]){
                  print_error("Target orientation description error","Minimum Euler psi > Maximum Euler psi","Please check the \'control.input\' file");
               }
            }
         }
/* *********************************************************************************************************************************************
Scattering directions */
/* Spherical coordinate azimuth angle - phi
If "scattering_direction_phi.input" exists read in data else use minimum,increment,maximum in "control.input" */
         else if(strcmp(compare,"Scattering angle phi=")==0){
            if((data=fopen("INPUT_scattering_direction_phi.input","r"))!=NULL){
               /* Read data from file into the structure phi */
               read_from_file(&data,&phi,"INPUT_scattering_direction_phi.input","Scattering angle phi");
            }
            else{ /* Use minimum,increment,maximum in "control.input" */
               while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               for(m=0;m<3;m++){ /* Use minimum,increment,maximum in "control.input" */
                  k=0;reset_string(compare);
                  do{ /* Copy everything after '=' OR ',' to ',' OR '\n' to compare string,
                     ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
                     if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
                     else{compare[k++]=string[i];}
                  } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
                  phi.limits[m]=(float)atof(compare); /* string to float conversion */
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               }
               /* Check phi values */
               if(phi.limits[0]<0.0F){
                  print_error("Scattering direction description error","phi minimum MUST be >= 0","Please check the \'control.input\' file");
               }
               if(phi.limits[1]<0.0F){
                  print_error("Scattering direction description error","phi increment MUST be >= 0","Please check the \'control.input\' file");
               }
               if(phi.limits[2]<0.0F){
                  print_error("Scattering direction description error","phi maximum MUST be >= 0","Please check the \'control.input\' file");
               }
               if(phi.limits[0]>phi.limits[2]){
                  print_error("Scattering direction description error","Minimum phi > Maximum phi","Please check the \'control.input\' file");
               }
            }
         }
/* Spherical coordinate polar/zenith angle - theta
If "scattering_direction_theta.input" exists read in data else use minimum,increment,maximum in "control.input" */
         else if(strcmp(compare,"Scattering angle theta=")==0){
            if((data=fopen("INPUT_scattering_direction_theta.input","r"))!=NULL){
               /* Read data from file into the structure theta */
               read_from_file(&data,&theta,"INPUT_scattering_direction_theta.input","Scattering angle theta");
            }
            else{ /* Use minimum,increment,maximum in "control.input" */
               while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               for(m=0;m<3;m++){ /* Use minimum,increment,maximum in "control.input" */
                  k=0;reset_string(compare);
                  do{ /* Copy everything after '=' OR ',' to ',' OR '\n' to compare string,
                     ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
                     if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
                     else{compare[k++]=string[i];}
                  } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
                  theta.limits[m]=(float)atof(compare); /* string to float conversion */
                  while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
               }
               /* Check theta values */
               if(theta.limits[0]<0.0F){
                  print_error("Scattering direction description error","theta minimum MUST be >= 0","Please check the \'control.input\' file");
               }
               if(theta.limits[1]<0.0F){
                  print_error("Scattering direction description error","theta increment MUST be >= 0","Please check the \'control.input\' file");
               }
               if(theta.limits[2]<0.0F){
                  print_error("Scattering direction description error","theta maximum MUST be >= 0","Please check the \'control.input\' file");
               }
               if(theta.limits[0]>theta.limits[2]){
                  print_error("Scattering direction description error","Minimum theta > Maximum theta","Please check the \'control.input\' file");
               }
            }
         }
/* *********************************************************************************************************************************************
Iterative scheme - Choice and Number of starting vectors which is only relevant for 'mlbicgstab_orig','mlbicgstab','mlbicgstab_ss' & 'rbicgstab' */
         else if(strcmp(compare,"Iterative scheme=")==0){
            /* Extract the iterative scheme */
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and underscore) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='_'))){continue;}
               else{compare[k++]=(char)tolower(string[i]);}
            } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
            reset_string(iterative.scheme);
            strncpy(iterative.scheme,compare,strlen(compare));
            /* Check the iterative scheme */
            if((strcmp(iterative.scheme,"bicg")!=0)&&
               (strcmp(iterative.scheme,"bicg_sym")!=0)&&
               (strcmp(iterative.scheme,"bicgstab")!=0)&&
               (strcmp(iterative.scheme,"cg")!=0)&&
               (strcmp(iterative.scheme,"cgs")!=0)&&
               (strcmp(iterative.scheme,"mlbicgstab_orig")!=0)&&
               (strcmp(iterative.scheme,"mlbicgstab")!=0)&&
               (strcmp(iterative.scheme,"mlbicgstab_ss")!=0)&&
               (strcmp(iterative.scheme,"qmr")!=0)&&
               (strcmp(iterative.scheme,"qmr_sym")!=0)&&
               (strcmp(iterative.scheme,"rbicgstab")!=0)&&
               (strcmp(iterative.scheme,"tfqmr")!=0)){
               reset_string(string);
               sprintf(string,"Scheme identifier \'%s\' not recognised",iterative.scheme);
               print_error("Iterative scheme description error",string,"Please check the \'control.input\' file");
            }
            /* Extract the number of starting vectors */
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and underscore) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='_'))){continue;}
               else{compare[k++]=(char)toupper(string[i]);}
            } while(string[i+1]!='\n');
            /* Check for illegal NULL */
            if(((strcmp(compare,"NULL")==0)&&(strcmp(iterative.scheme,"mlbicgstab_orig")==0))||
               ((strcmp(compare,"NULL")==0)&&(strcmp(iterative.scheme,"mlbicgstab")==0))||
               ((strcmp(compare,"NULL")==0)&&(strcmp(iterative.scheme,"mlbicgstab_ss")==0))||
               ((strcmp(compare,"NULL")==0)&&(strcmp(iterative.scheme,"rbicgstab")==0))){
               reset_string(string);
               sprintf(string,"The iterative scheme \'%s\' requires a non-NULL value for the number of starting vectors",iterative.scheme);
               print_error("Iterative scheme description error",string,"Please check the \'control.input\' file");
            }
            else if(strcmp(compare,"NULL")==0){ /* If NULL set to ZERO */
               iterative.vec=0;
            }
            else{
               iterative.vec=atoi(compare); /* string to int conversion */
            }
            /* Check vec value */
            if((((iterative.vec<1)||(iterative.vec>50))&&(strcmp(iterative.scheme,"mlbicgstab_orig")==0))||
               (((iterative.vec<1)||(iterative.vec>50))&&(strcmp(iterative.scheme,"mlbicgstab")==0))||
               (((iterative.vec<1)||(iterative.vec>50))&&(strcmp(iterative.scheme,"mlbicgstab_ss")==0))||
               (((iterative.vec<1)||(iterative.vec>50))&&(strcmp(iterative.scheme,"rbicgstab")==0))){
               reset_string(string);
               sprintf(string,"The number of starting vectors for \'%s\' MUST be > 0 and < 51",iterative.scheme);
               print_error("Iterative scheme description error",string,"Please check the \'control.input\' file");
            }
         }
/* *********************************************************************************************************************************************
Iterative scheme - Convergence tolerance */
         else if(strcmp(compare,"Convergence tolerance=")==0){
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
               else{compare[k++]=string[i];}
            } while(string[i+1]!='\n');
            iterative.tolerance=(float)atof(compare); /* string to float conversion */
            /* Check tolerance value */
            if(iterative.tolerance<0.0F){
               print_error("Iterative scheme description error","The convergence tolerance MUST be > 0","Please check the \'control.input\' file");
            }
         }
/* *********************************************************************************************************************************************
Iterative scheme - Breakdown tolerance */
         else if(strcmp(compare,"Breakdown tolerance=")==0){
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
               else{compare[k++]=string[i];}
            } while(string[i+1]!='\n');
            iterative.breakdown=(float)atof(compare); /* string to float conversion */
            /* Check breakdown value */
            if(iterative.breakdown<0.0F){
               print_error("Iterative scheme description error","The breakdown tolerance MUST be > 0","Please check the \'control.input\' file");
            }
         }
/* *********************************************************************************************************************************************
Iterative scheme - Maximum number of iterations */
         else if(strcmp(compare,"Maximum number of iterations=")==0){
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
               else{compare[k++]=string[i];}
            } while(string[i+1]!='\n');
            iterative.maximum=(int)atof(compare); /* string to int conversion */
            /* Check breakdown value */
            if(iterative.maximum<0){
               print_error("Iterative scheme description error","The maximum number of iterations MUST be > 0","Please check the \'control.input\' file");
            }
         }
/* *********************************************************************************************************************************************
Iterative scheme - Initial guess */
         else if(strcmp(compare,"Initial guess=")==0){
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
               else{compare[k++]=string[i];}
            } while(string[i+1]!='\n');
            iterative.initial_guess=atoi(compare); /* string to int conversion */
            /* Check initialguess value */
            if((iterative.initial_guess<0)||(iterative.initial_guess>3)){
               print_error("Iterative scheme description error","The initial guess flag MUST be in the range 0..3","Please check the \'control.input\' file");
            }
         }
/* *********************************************************************************************************************************************
Iterative scheme - Preconditioning */
         else if(strcmp(compare,"Preconditioning=")==0){
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
               else{compare[k++]=string[i];}
            } while(string[i+1]!='\n');
            iterative.precond=atoi(compare); /* string to int conversion */
            /* Check precond value */
            if((iterative.precond!=0)&&(iterative.precond!=1)){
               print_error("Iterative scheme description error","The preconditioning flag MUST be set to either 0 or 1","Please check the \'control.input\' file");
            }
         }
/* *********************************************************************************************************************************************
Degree of the interpolating polynomial for the improved Akima method */
         else if(strcmp(compare,"Polynomial degree=")==0){
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
               else{compare[k++]=string[i];}
            } while(string[i+1]!='\n');
            degree=atoi(compare); /* string to int conversion */
            /* Check degree value */
            if(degree<3){
               print_error("Interpolation description error","The degree of the interpolating polynomial MUST be >= 3 [Default=3]","Please check the \'control.input\' file");
            }
         }
/* *************************************************************************************
Print out the target xyz coordinate data to a vtk file for 3D visualisation using a
GUI like Paraview [www.paraview.org] OR Mayavi [mayavi.sourceforge.net] */
         else if(strcmp(compare,"VTK target xyz coordinate data=")==0){
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
               else{compare[k++]=string[i];}
            } while(string[i+1]!='\n');
            target.output_xyz_vtk=atoi(compare); /* string to int conversion */
            /* Check output_xyz_vtk value */
            if((target.output_xyz_vtk<0)||(target.output_xyz_vtk>1)){
               print_error("Target description error","The VTK target xyz coordinate output flag MUST be set to either 0 or 1","Please check the \'control.input\' file");
            }
         }
/* *********************************************************************************************************************************************
Timing */
         else if(strcmp(compare,"Timing=")==0){
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
               else{compare[k++]=string[i];}
            } while(string[i+1]!='\n');
            timing.enabled=atoi(compare); /* string to int conversion */
            /* Check value */
            if((timing.enabled!=0)&&(timing.enabled!=1)){
               print_error("Timing description error","The timing flag MUST be set to either 0 or 1","Please check the \'control.input\' file");
            }
         }
/* *********************************************************************************************************************************************
Output precision */
         else if(strcmp(compare,"Output precision [FLOAT,DOUBLE]=")==0){
            /* Extract the output precision for the float data type */
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and underscore) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='_'))){continue;}
               else{compare[k++]=(char)tolower(string[i]);}
            } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
            FLTP=atoi(compare); /* string to int conversion */
            /* Check FLTP value */
            if((FLTP<0)||(FLTP>FLT_DIG)){
               reset_string(string);
               sprintf(string,"The float output precision MUST be in the range 0 to %d",FLT_DIG);
               print_error("Output precision description error",string,"Please check the \'control.input\' file");
            }
            if(FLTP==0){ /* if set to ZERO then set equal to the maximum */
               FLTP=FLT_DIG;
            }
            /* Extract the output precision for the double data type */
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and underscore) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='_'))){continue;}
               else{compare[k++]=string[i];}
            } while(string[i+1]!='\n');
            DBLP=atoi(compare); /* string to int conversion */
            /* Check DBLP value */
            if((DBLP<0)||(DBLP>DBL_DIG)){
               reset_string(string);
               sprintf(string,"The double output precision MUST be in the range 0 to %d",DBL_DIG);
               print_error("Output precision description error",string,"Please check the \'control.input\' file");
            }
            if(DBLP==0){ /* if set to ZERO then set equal to the maximum */
               DBLP=DBL_DIG;
            }
         }
/* *********************************************************************************************************************************************
Distributed transpose algorithm */
         else if(strcmp(compare,"Distributed transpose algorithm=")==0){
            /* Extract the distributed transpose algorithm for the 6 independent tensor components of the interaction matrix */
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and underscore) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='_'))){continue;}
               else{compare[k++]=(char)tolower(string[i]);}
            } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
            parallel.tensor_transpose=atoi(compare); /* string to int conversion */
            /* Check tensor_transpose value */
            if((parallel.tensor_transpose<0)||(parallel.tensor_transpose>5)){
               print_error("Distributed tarnspose description error","The distributed tensor transpose algorithm MUST be in the range 0..5","Please check the \'control.input\' file");
            }
            /* Extract the distributed transpose algorithm for the 3 zero-padded vector components */
            while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
            k=0;reset_string(compare);
            do{ /* Copy everything after '=' to '\n' to compare string,
               ignoring punctuation (except fullstop and underscore) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='_'))){continue;}
               else{compare[k++]=string[i];}
            } while(string[i+1]!='\n');
            parallel.padded_transpose=atoi(compare); /* string to int conversion */
            /* Check padded_transpose value */
            if((parallel.padded_transpose<0)||(parallel.padded_transpose>5)){
               print_error("Distributed tarnspose description error","The distributed vector transpose algorithm MUST be in the range 0..5","Please check the \'control.input\' file");
            }
         }
         else{ /* Print error */
            j=0;
            k=(int)strlen(compare);
            for(i=0;i<k;i++){
               if((isspace(compare[i]))==0){j=1;} /* Character is not a space, tab, newline etc.. */
            }
            if(j==0){ /* Input line is blank */
               print_error("Parameter file input error","A forbidden blank line was encountered","Please check the \'control.input\' file");
            }
            else{
               reset_string(string);
               for(i=0;i<k;i++){
                  if(compare[i]=='\t'){ /* Replace tabs with spaces for correct error printing */
                     compare[i]=' ';
                  }
               }
               sprintf(string,"The input description string \"%s\" is invalid",compare);
               print_error("Parameter file input error",string,"Please check the \'control.input\' file");
            }
         }
      }
      reset_string(string);
   }
   fclose(input);
}

void read_from_file(FILE **file_ptr,attributes_general *structure_ptr,char *file_string,char *data_string){
/* Reads in the float data from the file pointed to by *file_ptr (ignoring comments) to the array
   pointed to by *array_ptr, sets the flag to indicate that the data has been read in from a
   file and sets the length of the array */

   int k,count=0;
   float previous;
   char string[200],string1[200],string2[200];

   reset_string(string);

   (*structure_ptr).flag=1; /* Set the flag to indicate that data has been read from a file */
   while((fgets(string,sizeof(string),*file_ptr))!=NULL){ /* For memory allocation, count the number of data entries in the file */
      /* Count the number of data entries in the file */
      if(string[0]=='#'){continue;} /* Ignore comments */
      else{count++;}
   }
   fclose(*file_ptr); /* Close file */
   float_malloc(&(*structure_ptr).list,(size_t)count); /* Allocate memory to store the data */
   (*structure_ptr).list_length=count; /* Set the length of the array to count */

   *file_ptr=fopen(file_string,"r"); /* Reopen the data file for reading */
   reset_string(string);
   count=0; /* Reset counter */
   /* Read in the data */
   while((fgets(string,sizeof(string),*file_ptr))!=NULL){
      if(string[0]=='#'){continue;} /* Ignore comments */
      else{
         k=(int)strlen(string);
         string[k-1]='\0'; /* Strip the '\n' */
         ((*structure_ptr).list)[count]=(float)atof(string); /* string to long float conversion */
         if((strcmp(data_string,"Wavelength")==0)||(strcmp(data_string,"Effective radius")==0)){
            /* Wavelength and effective radius values must be > 0 */
            if(((*structure_ptr).list)[count]<FLT_EPSILON){ /* Input data MUST be > 0 */
               reset_string(string);
               reset_string(string1);
               reset_string(string2);
               sprintf(string,"Please check the file %s",file_string);
               sprintf(string1,"%s input error",data_string);
               sprintf(string2,"([%d:%.*g] <= 0) %s data MUST be > 0",count,LDBL_DIG,((*structure_ptr).list)[count],data_string);
               print_error(string1,string2,string);
            }
         }
         else if(strcmp(data_string,"Target orientation Euler theta")!=0){
            /* If not Wavelength, Effective radius OR Euler theta i.e.
               Euler phi
               Euler psi
               phi
               theta
               the data MUST be >= 0
               Note: Euler theta is distributed via cosine [-1,1] so it can be negative */
            if(((*structure_ptr).list)[count]<0.0F){ /* Input data MUST be positive */
               reset_string(string);
               reset_string(string1);
               reset_string(string2);
               sprintf(string,"Please check the file %s",file_string);
               sprintf(string1,"%s input error",data_string);
               sprintf(string2,"([%d:%.*g] < 0) %s data MUST be >= 0",count,LDBL_DIG,((*structure_ptr).list)[count],data_string);
               print_error(string1,string2,string);
            }
         }
         if(count==0){ /* Set initial value of previous */
            previous=((*structure_ptr).list)[count];
         }
         else{ /* Check that the data is monotonically increasing */
            if((previous>((*structure_ptr).list)[count])||((((*structure_ptr).list)[count]-previous)<FLT_EPSILON)){
               /* If data is not monotonically increasing */
               reset_string(string);
               reset_string(string1);
               reset_string(string2);
               sprintf(string,"Please check the file %s",file_string);
               sprintf(string1,"%s input error",data_string);
               sprintf(string2,"([%d:%.*g] >= [%d:%.*g]) %s data MUST be monotonically increasing",count-1,LDBL_DIG,((*structure_ptr).list)[count-1],count,LDBL_DIG,((*structure_ptr).list)[count],data_string);
               print_error(string1,string2,string);
            }
            previous=((*structure_ptr).list)[count]; /* Update previous */
         }
         count++; /* Next entry */
      }
      reset_string(string);
   }
   fclose(*file_ptr);
}

void read_from_file_refractive_index_vs_wavelength(FILE **file_ptr,attributes_refractive_index *structure_ptr,char *file_string,int material){
/* Reads in the wavelength and complex refractive index data from the file
   "INPUT_refractive_index_vs_wavelength_material_`material'.input" (ignoring comments) into the
   structure refractive_index, sets the flag to indicate that
   the data has been read in from a file and sets the length of the array */

   int a,i,k,count=0;
   float previous; /* Check if wavelength data is monotonically increasing */
   char string[200],compare[200];

   reset_string(string);
   reset_string(compare);

   (*structure_ptr).flag[material]=1; /* Set the flag to indicate that data has been read from a file */

   while((fgets(string,sizeof(string),*file_ptr))!=NULL){ /* For memory allocation, count the number of data entries in the file */
      /* Count the number of data entries in the file */
      if(string[0]=='#'){continue;} /* Ignore comments */
      else{count++;}
   }
   fclose(*file_ptr); /* Close file */   

   float_malloc(&(*structure_ptr).wavelengths[material],(size_t)count); /* Allocate memory to store the wavelength data */
   fcomplex_malloc(&(*structure_ptr).refractive_indices[material],(size_t)count); /* Allocate memory to store the complex refractive index data */
   (*structure_ptr).list_length[material]=count; /* Set the length of the array to count */

   *file_ptr=fopen(file_string,"r"); /* Reopen the data file for reading */
   reset_string(string);
   count=0; /* Reset counter */
   /* Read in the data */
   while((fgets(string,sizeof(string),*file_ptr))!=NULL){
      if(string[0]=='#'){continue;} /* Ignore comments */
      else{
         i=-1;
         while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
         k=0;reset_string(compare);
         do{ /* Extract the wavelength data
            Copy everything to ',' to compare string,
            ignoring punctuation (except fullstop and minus) and spaces, tabs, newlines etc.. */
            if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-'))){continue;}
            else{compare[k++]=string[i];}
         } while((string[i+1]!=',')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
         ((*structure_ptr).wavelengths[material])[count]=(float)atof(compare); /* string to long float conversion */
         if(((*structure_ptr).wavelengths[material])[count]<FLT_EPSILON){ /* Wavelength data MUST be positive */
            reset_string(string);
            reset_string(compare);
            sprintf(compare,"([%d:%.*g] <= 0) The wavelength data MUST be > 0",count,LDBL_DIG,((*structure_ptr).wavelengths[material])[count]);
            sprintf(string,"Please check the file INPUT_ref_index_vs_wavelength_material_%d.input",material);
            print_error("Complex refractive index vs wavelength input error",compare,string);
         }
         if(count==0){ /* Set initial value of previous */
            previous=((*structure_ptr).wavelengths[material])[count];
         }
         else{ /* Check that the wavelength data is monotonically increasing */
            if((previous>((*structure_ptr).wavelengths[material])[count])||((((*structure_ptr).wavelengths[material])[count]-previous)<FLT_EPSILON)){
               /* If data is not monotonically increasing */
               reset_string(string);
               reset_string(compare);
               sprintf(string,"Please check the file INPUT_ref_index_vs_wavelength_material_%d.input",material);
               sprintf(compare,"([%d:%.*g] >= [%d:%.*g]) Wavelength data MUST be monotonically increasing",count-1,LDBL_DIG,((*structure_ptr).wavelengths[material])[count-1],count,LDBL_DIG,((*structure_ptr).wavelengths[material])[count]);
               print_error("Complex refractive index vs wavelength input error",compare,string);
            }
            previous=((*structure_ptr).wavelengths[material])[count]; /* Update previous */
         }
         while(((ispunct(string[i+1]))!=0)&&(string[i+1]!='-')){i++;} /* Skip all punctuation (except '-') */
         /* Extract the complex refractive index data */
         for(a=0;a<2;a++){ /* a=0 real component,a=1 imaginary component */
            k=0;reset_string(compare);
            do{ /* Copy everything after ',' to '+' OR '-' OR ',' OR '\n' to compare string,
               ignoring punctuation (except fullstop and signs) and spaces, tabs, newlines etc.. */
               if(((isspace(string[++i]))!=0)||(((ispunct(string[i]))!=0)&&(string[i]!='.')&&(string[i]!='-')&&(string[i]!='+'))){continue;}
               else{compare[k++]=string[i];}
            } while((string[i+1]!='+')&&(string[i+1]!='-')&&(string[i+1]!=',')&&(string[i+1]!='I')&&(string[i+1]!='i')&&(string[i+1]!='\0')&&(string[i+1]!='\n'));
            while(((ispunct(string[i+1])!=0)&&(string[i+1]!='-'))||(string[i+1]=='i')||(string[i+1]=='I')){i++;} /* Skip 'i', 'I' and all punctuation (except '-') */
            ((*structure_ptr).refractive_indices[material])[count].dat[a]=(float)atof(compare); /* string to float conversion */
         }
         count++; /* Next entry */
      }
      reset_string(string);
   }
   fclose(*file_ptr);
}
