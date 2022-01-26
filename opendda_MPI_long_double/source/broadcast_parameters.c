/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Broadcasts parameters to the slave nodes

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "broadcast_parameters.h"

void broadcast_control_parameters(void){

   long begin;
   int i,count=26;
   int blocklens[26]={6,1,50,1,3,3,15,1,1,1,1,1,1,1,3,3,3,3,3,3,3,1,1,1,1,1};
   MPI_Datatype types[26]={MPI_INT,MPI_INT,MPI_CHAR,MPI_INT,ldcomplex_type,ldcomplex_type,MPI_CHAR,
                  MPI_INT,MPI_LONG_DOUBLE,MPI_LONG_DOUBLE,MPI_INT,MPI_INT,MPI_INT,
                  MPI_INT,MPI_LONG_DOUBLE,MPI_LONG_DOUBLE,MPI_LONG_DOUBLE,
                  MPI_LONG_DOUBLE,MPI_LONG_DOUBLE,MPI_LONG_DOUBLE,MPI_LONG_DOUBLE,
                  MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT};
   MPI_Aint displacements[26];
   MPI_Datatype broadcast_control_type;

   /* Define a structure for broadcasting the control parameters */
   struct broadcast_control{
      int construction_parameters[6]; /* target construction parameters */
      int from_file; /* target data read from file flag */
      char shape[50]; /* target description string */
      int polarisations; /* number of incident polarisations */
      ldcomplex polarisation[3]; /* Incident polarisation state */
      ldcomplex orthonormal[3]; /* Orthonormal polarisation state */
      char iterative_scheme[15]; /* iterative scheme */
      int starting_vector; /* number of starting vectors */
      long double iterative_tolerance; /* convergence tolerance */
      long double iterative_breakdown; /* breakdown tolerance */
      int iterative_maximum; /* maximum number of iterations */
      int iterative_initial_guess; /* initial guess */
      int iterative_precond; /* preconditioning */
      int number_of_materials; /* number of materials in the target */
      long double wavelength_limits[3]; /* wavelength limits */
      long double radius_limits[3]; /* radius limits */
      long double euler_phi_limits[3]; /* Euler phi limits */
      long double euler_theta_limits[3]; /* Euler theta limits */
      long double euler_psi_limits[3]; /* Euler psi limits */
      long double theta_limits[3]; /* theta limits */
      long double phi_limits[3]; /* phi limits */
      int timing_enabled; /* timing flag */
      int DBLP; /* double output precision */
      int LDBLP; /* long double output precision */
      int tensor_transpose; /* tensor transpose algorithm */
      int padded_transpose; /* vector transpose algorithm */
   };
   struct broadcast_control control;

   /* Calculate the data offsets fot the MPI_type_struct datatype */
   MPI_Address(&control,&displacements[0]);
   MPI_Address(&control.from_file,&displacements[1]);
   MPI_Address(control.shape,&displacements[2]);
   MPI_Address(&control.polarisations,&displacements[3]);
   MPI_Address(control.polarisation,&displacements[4]);
   MPI_Address(control.orthonormal,&displacements[5]);
   MPI_Address(control.iterative_scheme,&displacements[6]);
   MPI_Address(&control.starting_vector,&displacements[7]);
   MPI_Address(&control.iterative_tolerance,&displacements[8]);
   MPI_Address(&control.iterative_breakdown,&displacements[9]);
   MPI_Address(&control.iterative_maximum,&displacements[10]);
   MPI_Address(&control.iterative_initial_guess,&displacements[11]);
   MPI_Address(&control.iterative_precond,&displacements[12]);
   MPI_Address(&control.number_of_materials,&displacements[13]);
   MPI_Address(control.wavelength_limits,&displacements[14]);
   MPI_Address(control.radius_limits,&displacements[15]);
   MPI_Address(control.euler_phi_limits,&displacements[16]);
   MPI_Address(control.euler_theta_limits,&displacements[17]);
   MPI_Address(control.euler_psi_limits,&displacements[18]);
   MPI_Address(control.theta_limits,&displacements[19]);
   MPI_Address(control.phi_limits,&displacements[20]);
   MPI_Address(&control.timing_enabled,&displacements[21]);
   MPI_Address(&control.DBLP,&displacements[22]);
   MPI_Address(&control.LDBLP,&displacements[23]);
   MPI_Address(&control.tensor_transpose,&displacements[24]);
   MPI_Address(&control.padded_transpose,&displacements[25]);
   begin=displacements[0];
   for(i=0;i<count;i++){
      displacements[i]-=begin;
   }

   if(parallel.myrank==0){ /* Restrict to master */
      /* Cram all pertinent data into a single structure for a single broadcast */
      for(i=0;i<6;i++){ /* target description string */
         control.construction_parameters[i]=target.construction_parameters[i];
      }
      control.from_file=target.from_file; /* target data read from file flag */
      for(i=0;i<50;i++){ /* target description string */
         control.shape[i]=target.shape[i];
      }
      control.polarisations=polarisations; /* number of incident polarisations */
      for(i=0;i<3;i++){ /* Incident polarisation state */
         control.polarisation[i]=incident_LF.polarisation[i];
      }
      for(i=0;i<3;i++){ /* Incident polarisation state */
         control.orthonormal[i]=incident_LF.orthonormal[i];
      }
      for(i=0;i<15;i++){ /* iterative scheme */
         control.iterative_scheme[i]=iterative.scheme[i];
      }
      control.starting_vector=iterative.vec; /* number of starting vectors */
      control.iterative_tolerance=iterative.tolerance; /* convergence tolerance */
      control.iterative_breakdown=iterative.breakdown; /* breakdown tolerance */
      control.iterative_maximum=iterative.maximum; /* maximum number of iterations */
      control.iterative_initial_guess=iterative.initial_guess; /* initial guess */
      control.iterative_precond=iterative.precond; /* preconditioning */
      control.number_of_materials=refractive_index.number_of_materials; /* number of materials in the target */
      for(i=0;i<3;i++){ /* wavelength limits */
         control.wavelength_limits[i]=wavelength.limits[i];
      }
      for(i=0;i<3;i++){ /* radius limits */
         control.radius_limits[i]=radius.limits[i];
      }
      for(i=0;i<3;i++){ /* Euler phi limits */
         control.euler_phi_limits[i]=euler_phi.limits[i];
      }
      for(i=0;i<3;i++){ /* Euler theta limits */
         control.euler_theta_limits[i]=euler_theta.limits[i];
      }
      for(i=0;i<3;i++){ /* Euler psi limits */
         control.euler_psi_limits[i]=euler_psi.limits[i];
      }
      for(i=0;i<3;i++){ /* theta limits */
         control.theta_limits[i]=theta.limits[i];
      }
      for(i=0;i<3;i++){ /* phi limits */
         control.phi_limits[i]=phi.limits[i];
      }
      control.timing_enabled=timing.enabled; /* timing flag */
      control.DBLP=DBLP; /* double output precision */
      control.LDBLP=LDBLP; /* long double output precision */
      control.tensor_transpose=parallel.tensor_transpose; /* tensor transpose algorithm */
      control.padded_transpose=parallel.padded_transpose; /* vector transpose algorithm */
   }

   /* Create and commit the MPI Datatype */
   MPI_Type_struct(count,blocklens,displacements,types,&broadcast_control_type);   
   MPI_Type_commit(&broadcast_control_type);

   /* Broadcast the control parameter structure */
   MPI_Bcast(&control,1,broadcast_control_type,0,MPI_COMM_WORLD);

   if(parallel.myrank!=0){ /* Restrict to slaves */
      /* Extract the data */
      for(i=0;i<6;i++){ /* target description string */
         target.construction_parameters[i]=control.construction_parameters[i];
      }
      target.from_file=control.from_file; /* target data read from file flag */
      for(i=0;i<50;i++){ /* target description string */
         target.shape[i]=control.shape[i];
      }
      polarisations=control.polarisations; /* number of incident polarisations */
      for(i=0;i<3;i++){ /* Incident polarisation state */
         incident_LF.polarisation[i]=control.polarisation[i];
      }
      for(i=0;i<3;i++){ /* Incident polarisation state */
         incident_LF.orthonormal[i]=control.orthonormal[i];
      }
      for(i=0;i<15;i++){ /* iterative scheme */
         iterative.scheme[i]=control.iterative_scheme[i];
      }
      iterative.vec=control.starting_vector; /* number of starting vectors */
      iterative.tolerance=control.iterative_tolerance; /* convergence tolerance */
      iterative.breakdown=control.iterative_breakdown; /* breakdown tolerance */
      iterative.maximum=control.iterative_maximum; /* maximum number of iterations */
      iterative.initial_guess=control.iterative_initial_guess; /* initial guess */
      iterative.precond=control.iterative_precond; /* preconditioning */
      refractive_index.number_of_materials=control.number_of_materials; /* number of materials in the target */
      for(i=0;i<3;i++){ /* wavelength limits */
         wavelength.limits[i]=control.wavelength_limits[i];
      }
      for(i=0;i<3;i++){ /* radius limits */
         radius.limits[i]=control.radius_limits[i];
      }
      for(i=0;i<3;i++){ /* Euler phi limits */
         euler_phi.limits[i]=control.euler_phi_limits[i];
      }
      for(i=0;i<3;i++){ /* Euler theta limits */
         euler_theta.limits[i]=control.euler_theta_limits[i];
      }
      for(i=0;i<3;i++){ /* Euler psi limits */
         euler_psi.limits[i]=control.euler_psi_limits[i];
      }
      for(i=0;i<3;i++){ /* theta limits */
         theta.limits[i]=control.theta_limits[i];
      }
      for(i=0;i<3;i++){ /* phi limits */
         phi.limits[i]=control.phi_limits[i];
      }
      timing.enabled=control.timing_enabled; /* timing flag */
      DBLP=control.DBLP; /* double output precision */
      LDBLP=control.LDBLP; /* long double output precision */
      parallel.tensor_transpose=control.tensor_transpose; /* tensor transpose algorithm */
      parallel.padded_transpose=control.padded_transpose; /* vector transpose algorithm */
   }

   /* Free the MPI Datatype */
   MPI_Type_free(&broadcast_control_type);
}

void broadcast_length_parameters(void){

   int count=7,*length;

   /* Allocate memory to store the length integers */
   int_malloc(&length,(size_t)count);

   if(parallel.myrank==0){ /* Restrict to master */
      /* Cram all pertinent data into a single array for a single broadcast */
      length[0]=wavelength.list_length;
      length[1]=radius.list_length;
      length[2]=euler_phi.list_length;
      length[3]=euler_theta.list_length;
      length[4]=euler_psi.list_length;
      length[5]=theta.list_length;
      length[6]=phi.list_length;
   }

   /* Broadcast the list length array */
   MPI_Bcast(length,count,MPI_INT,0,MPI_COMM_WORLD);

   if(parallel.myrank!=0){ /* Restrict to slaves */
      /* Extract the data */
      wavelength.list_length=length[0];
      radius.list_length=length[1];
      euler_phi.list_length=length[2];
      euler_theta.list_length=length[3];
      euler_psi.list_length=length[4];
      theta.list_length=length[5];
      phi.list_length=length[6];
   }

   /* Free allocated memory */
   int_free(&length);
}

void broadcast_physical_parameters(void){

   long begin;
   int i,count=8;
   int blocklens[8];
   MPI_Datatype types[8]={MPI_LONG_DOUBLE,MPI_LONG_DOUBLE,MPI_LONG_DOUBLE,
               MPI_LONG_DOUBLE,MPI_LONG_DOUBLE,MPI_LONG_DOUBLE,
               MPI_LONG_DOUBLE,ldcomplex_type};
   MPI_Aint displacements[8];
   MPI_Datatype broadcast_physical_type;

   /* Define a structure for broadcasting the physical parameters */
   struct broadcast_physical{
      long double wavelength_value;
      long double wavenumber;
      long double radius_value;
      long double target_dipole_spacing;
      long double euler_phi_value;
      long double euler_theta_value;
      long double euler_psi_value;
      ldcomplex *refractive_index_value;
   };
   struct broadcast_physical physical;

   /* Allocate memory to store the current refractive index values */
   ldcomplex_malloc(&physical.refractive_index_value,(size_t)refractive_index.number_of_materials);

   /* Set the block lengths for the MPI_type_struct datatype */
   blocklens[0]=1;
   blocklens[1]=1;
   blocklens[2]=1;
   blocklens[3]=1;
   blocklens[4]=1;
   blocklens[5]=1;
   blocklens[6]=1;
   blocklens[7]=refractive_index.number_of_materials;

   /* Calculate the data offsets fot the MPI_type_struct datatype */
   MPI_Address(&physical,&displacements[0]);
   MPI_Address(&physical.wavenumber,&displacements[1]);
   MPI_Address(&physical.radius_value,&displacements[2]);
   MPI_Address(&physical.target_dipole_spacing,&displacements[3]);
   MPI_Address(&physical.euler_phi_value,&displacements[4]);
   MPI_Address(&physical.euler_theta_value,&displacements[5]);
   MPI_Address(&physical.euler_psi_value,&displacements[6]);
   MPI_Address(physical.refractive_index_value,&displacements[7]);
   begin=displacements[0];
   for(i=0;i<count;i++){
      displacements[i]-=begin;
   }

   if(parallel.myrank==0){ /* Restrict to master */
      /* Cram all pertinent data into a single structure for a single broadcast */
      physical.wavelength_value=wavelength.value;
      physical.wavenumber=wavenumber;
      physical.radius_value=radius.value;
      physical.target_dipole_spacing=target.dipole_spacing;
      physical.euler_phi_value=euler_phi.value;
      physical.euler_theta_value=euler_theta.value;
      physical.euler_psi_value=euler_psi.value;
      for(i=0;i<refractive_index.number_of_materials;i++){
         physical.refractive_index_value[i]=refractive_index.value[i];
      } 
   }

   /* Create and commit the MPI Datatype */
   MPI_Type_struct(count,blocklens,displacements,types,&broadcast_physical_type);   
   MPI_Type_commit(&broadcast_physical_type);

   /* Broadcast the physical parameter structure */
   MPI_Bcast(&physical,1,broadcast_physical_type,0,MPI_COMM_WORLD);

   if(parallel.myrank!=0){ /* Restrict to slaves */
      /* Extract the data */
      wavelength.value=physical.wavelength_value;
      wavenumber=physical.wavenumber;
      radius.value=physical.radius_value;
      target.dipole_spacing=physical.target_dipole_spacing;
      euler_phi.value=physical.euler_phi_value;
      euler_theta.value=physical.euler_theta_value;
      euler_psi.value=physical.euler_psi_value;
      for(i=0;i<refractive_index.number_of_materials;i++){
         refractive_index.value[i]=physical.refractive_index_value[i];
      }
   }

   /* Free the MPI Datatype */
   MPI_Type_free(&broadcast_physical_type);

   /* Free memory that stored the current refractive index values */
   ldcomplex_free(&physical.refractive_index_value);
}

void broadcast_target_from_file(void){

   int i;
   unsigned int *file_target;

   if(refractive_index.number_of_materials==1){
      /* Broadcast the bitfielded array of 0/1 flags for whether a site is occupied by a dipole or not */
      MPI_Bcast(target.occupied,target.occ_ctrl,MPI_UNSIGNED,0,MPI_COMM_WORLD);
   }
   else{
      /* Composition array is only required if more than one material is used for target construction */
      /* Allocate memory to store the bitfielded arrays */
      uint_malloc(&file_target,(size_t)(target.occ_ctrl+3*target.mat_ctrl));

      if(parallel.myrank==0){ /* Restrict to master */
         /* Cram all pertinent data into a single array for a single broadcast */
         memcpy(file_target,target.occupied,(size_t)target.occ_ctrl*sizeof(unsigned int));
         for(i=0;i<3;i++){
            memcpy(&file_target[target.occ_ctrl+i*target.mat_ctrl],target.composition[i],(size_t)target.mat_ctrl*sizeof(unsigned int));
         }
      }
 
      /* Broadcast the bitfielded target data */
      MPI_Bcast(file_target,target.occ_ctrl+3*target.mat_ctrl,MPI_UNSIGNED,0,MPI_COMM_WORLD);

      if(parallel.myrank!=0){ /* Restrict to slaves */
         /* Extract the data */
         memcpy(target.occupied,file_target,(size_t)target.occ_ctrl*sizeof(unsigned int));
         for(i=0;i<3;i++){
            memcpy(target.composition[i],&file_target[target.occ_ctrl+i*target.mat_ctrl],(size_t)target.mat_ctrl*sizeof(unsigned int));
         }
      }

      /* Free allocated memory */
      uint_free(&file_target);
   }
}
