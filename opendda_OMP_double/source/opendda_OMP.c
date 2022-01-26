/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Master control file for the OpenDDA
   See files definitions.h OR attributes_type.h for
   explanation of the variables

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "opendda_OMP.h"

int main(void){

   int j,k;
   long int index0i,index0k;

   read_parameters(); /* Read in user-defined control parameters */

#pragma omp parallel
   #pragma omp master
      { /* Get the number of threads available */
         threads_available=omp_get_num_threads();
      }
      /* If valid, set the number of threads to use for parallel regions */
      if((openmp_threads>0)&&(openmp_threads<threads_available)){
         omp_set_num_threads(openmp_threads);
      }
      else{
         openmp_threads=threads_available;
      }

   if(timing.enabled){timing.overall[0]=walltime();} /* Initialise timing */

   set_domain(); /* Set K, J and P for the target shape */

   target.M=target.K*target.J; /* M=K*J */
   target.N=target.K*target.J*target.P; /* N=K*J*P */
   target.N3=3*target.N; /* N*3 */

/* FFTW is best at handling sizes of the form 2^{a}3^{b}5^{c}7^{d}11^{e}13^{f},
   where (e + f) is either 0 or 1, and the other exponents are arbitrary. */
   find_dft_size(target.K,&target.Kp,&target.Kpp,&target.kf);
   find_dft_size(target.J,&target.Jp,&target.Jpp,&target.jf);
   find_dft_size(target.P,&target.Pp,&target.Ppp,&target.pf);

   target.Ka=target.K+target.kf; /* Ka=K+kf */
   target.Ja=target.J+target.jf; /* Ja=J+jf */
   target.Pa=target.P+target.pf; /* Pa=P+pf */
   target.Ma=target.Ja*target.Ka; /* Ma=Ka*Ja */
   target.Na=target.Ka*target.Ja*target.Pa; /* Na=Ka*Ja*Pa */
   target.Mpp=target.Kpp*target.Jpp; /* Mpp=Kpp*Jpp */
   target.Mv=target.K*target.Jpp; /* Mv=K*Jpp */
   target.Qpp=target.Kpp*target.Ppp; /* Qpp=Kpp*Ppp */
   target.Npp=target.Kpp*target.Jpp*target.Ppp; /* Npp=Kpp*Jpp*Ppp */
   target.Nv=target.K*target.Jpp*target.P; /* Nv=Kpp*Jpp*P */

/* Initialise bitfielding for target parameterisation */
   bitfielding();

/* Build the target */
   build_target();
   target.Nd3=3*target.Nd; /* Nd*3 */

/* If the corresponding flag=0 then the code is using minimum,increment,maximum
   from the "control.input" file and list_length=[(maximum-minimum)/increment]+1 */
   initialise(wavelength.flag,&wavelength.list_length,wavelength.limits); /* Wavelength */
   initialise(radius.flag,&radius.list_length,radius.limits); /* Effective radius */
   initialise(euler_phi.flag,&euler_phi.list_length,euler_phi.limits); /* Target orientation - Euler angle phi */
   initialise(euler_theta.flag,&euler_theta.list_length,euler_theta.limits); /* Target orientation - Euler angle theta */
   initialise(euler_psi.flag,&euler_psi.list_length,euler_psi.limits); /* Target orientation - Euler angle psi */
   initialise(phi.flag,&phi.list_length,phi.limits); /* Scattering direction - phi */
   initialise(theta.flag,&theta.list_length,theta.limits); /* Scattering direction - theta */

   initialise_output(); /* Initialise program output */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Allocate memory */
   memory_allocate_OMP();

   /* To conserve memory array elements that correspond to vacant lattice sites
      and are thus zero are not stored. The populated array stores the starting
      index for the non-zero elements for each of the P xy-planes */
   for(k=0;k<target.P;k++){target.populated[k]=0;}
#pragma omp parallel for private(index0k,j,index0i) schedule(dynamic)
   for(k=0;k<(target.P-1);k++){
      index0k=k*target.M;
      for(j=0;j<target.M;j++){
         index0i=index0k+j;
         if(target.occupied[(int)((double)index0i/32.0)]&(1<<(index0i%32))){ /* Lattice site is occupied */
            target.populated[k+1]++; /* Increment the number of dipoles stored in plane k */
         }
      }
   }
   for(k=0;k<(target.P-1);k++){
      target.populated[k+1]+=target.populated[k];
   }

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Create the FFTW DFT plans */
   dft_plan(1);

   /* -------------- Cycle through the incident wavelengths -------------- */
   wavelength.value=wavelength.limits[0]; /* Set wavelength to minimum from the "control.input" file.
   If flag=1 then data was read form file and this value will be overwritten */
   for(wl=0;wl<wavelength.list_length;wl++){
      if(wavelength.flag==1){wavelength.value=wavelength.list[wl];} /* Data was read from file into 'wavelength.list' */
      wavenumber=(2.0*dpi)/wavelength.value; /* k=(2*PI)/wavelength */

      /* Dielectric properties: calculate the complex refractive index for each of the dielectric materials in the target */
      for(j=0;j<refractive_index.number_of_materials;j++){ /* Loop over the constituent materials */
         /* Note: there are 2 possible scenarios for each material `j'

         (1) refractive_index.flag[j]=1 [Interpolation table read in from file]
            A file INPUT_ref_index_vs_wavelength_material_`j'.input exists and the complex refractive index vs
            wavelength data, of length refractive_index.list_length[j], was read in the arrays
            refractive_index.wavelengths[j] & refractive_index.refractive_indices[j] for subsequent interpolation
            of the appropriate complex refractive index for a given wavelength.
            In this case, refractive_index.value[j]=0.0+0.0i awaiting interpolation

         (2) refractive_index.flag[j]=0 [Constant value from control.input]
            If the file INPUT_ref_index_vs_wavelength_material_`j'.input does not exist, the constant complex refractive
            index was given in the control.input file using the line "Material `j'=a+bi".
            In this case, refractive_index.value[j]=a+bi, refractive_index.list_length[j]=1 and the array pointers
            refractive_index.wavelengths[j]=refractive_index.refractive_indices[j]=NULL */
         if(refractive_index.flag[j]==1){
         /* Complex refractive index vs wavelength interpolation table was read in from file for material `j'
            Interpolate the real and imaginary components of the refractive index */
            for(k=0;k<2;k++){
               interpolation(degree,refractive_index.list_length[j],refractive_index.wavelengths[j],
               refractive_index.refractive_indices[j],1,&wavelength.value,&refractive_index.value[j],k);
            }
         }
         /* Output current wavelength to OUTPUT_OpenDDA_refractive_index_vs_wavelength_`j'.info */
         fprintf(output.refindex_vs_wave[j],"\n%.*g",DBLP,wavelength.value);
          /* Output refractive index for material `j' for the current material to OUTPUT_OpenDDA_refractive_index_vs_wavelength_`j'.info */
         fprintf(output.refindex_vs_wave[j],"\t\t\t%.*g%+.*gi",DBLP,refractive_index.value[j].dat[0],DBLP,refractive_index.value[j].dat[1]);
      }

      /* -------------- Cycle through the effective radii -------------- */
      radius.value=radius.limits[0]; /* Set radius to minimum from the "control.input" file.
      If flag=1 then data was read form file and this value will be overwritten */
      for(er=0;er<radius.list_length;er++){
         if(radius.flag==1){radius.value=radius.list[er];} /* Data was read from file into 'radius.list' */

         /* Calculate the size parameter and the wavenumber
         size_parameter~wavenumber*effective_radius=[2*PI*effective_radius]/wavelength
         where the effective radius is the radius of a sphere of equal volume

         volume=(4/3)*PI*(effective_radius)^{3}=Nd*(dipole_spacing)^{3}
         where Nd=number of lattice sites occupied by dipoles

         dipole_spacing=effective_radius*[(4*PI)/(3*Nd)]^{1/3}
         dipoles_per_wavelength=wavelength/dipole_spacing */
         size_parameter=wavenumber*radius.value;
         target.dipole_spacing=radius.value*pow(((4.0*dpi)/(3.0*(double)target.Nd)),(1.0/3.0));
         target.dipoles_per_wavelength=wavelength.value/target.dipole_spacing;

         /* Output target description data */
         fprintf(output.target_descrip," %.*g\t\t%.*g\t\t%.*g\t\t%.*g\t\t%.*g\n",DBLP,wavelength.value,DBLP,radius.value,DBLP,size_parameter,DBLP,target.dipole_spacing,DBLP,target.dipoles_per_wavelength);
         
         /* -------------- Cycle through the target orientations -------------- */
         /* -------------- Cycle through the Euler theta -------------- */
         euler_theta.value=euler_theta.limits[0]; /* Set euler_theta to minimum from the "control.input" file.
         If flag=1 then data was read form file and this value will be overwritten */
         for(t0=0;t0<euler_theta.list_length;t0++){
            if(euler_theta.flag==1){euler_theta.value=euler_theta.list[t0];} /* Data was read from file into 'euler_theta.list' */

            /* -------------- Cycle through the Euler psi -------------- */
            euler_psi.value=euler_psi.limits[0]; /* Set euler_psi to minimum from the "control.input" file.
            If flag=1 then data was read form file and this value will be overwritten */
            for(p0=0;p0<euler_psi.list_length;p0++){
               if(euler_psi.flag==1){euler_psi.value=euler_psi.list[p0];} /* Data was read from file into 'euler_psi.list' */

               /* -------------- Cycle through the Euler phi -------------- */
               euler_phi.value=euler_phi.limits[0]; /* Set euler_phi to minimum from the "control.input" file.
               If flag=1 then data was read form file and this value will be overwritten */
               for(f0=0;f0<euler_phi.list_length;f0++){
                  if(euler_phi.flag==1){euler_phi.value=euler_phi.list[f0];} /* Data was read from file into 'euler_phi.list' */

                  /* Calculate the scattering vectors */
                  scattering_vectors(phi.value,theta.value);

                  /* Rotate the incident polarisation vectors and scattering vectors LF to TF */
                  rotation_euler(euler_phi.value,euler_theta.value,euler_psi.value);

                  /* Loop over the incident polarisations */
                  for(polarisation_state=0;polarisation_state<polarisations;polarisation_state++){

                      /* Calculate the incident electric field at the lattice sites */
                     incident_electric_field_OMP();

                      /* Calculate the dipole polarisabilities */
                     dipole_polarisabilities_OMP();

                     /* Create the first column of the DDA matrix excluding the diagonal */
                     dda_interaction_matrix_1st_column_zero_diagonal();

                     /* Set the initial guess for the iterative solver */
                     set_initial_guess();

                     /* Solve the linear system
                        interaction_matrix*polarisations=incident_E_field */
                     iterative_solution();

                     /* Calculate the cross sections */
                     cross_sections_OMP();

                     /* Calculate the efficiencies */
                     efficiencies();
                  } /* End loop over the polarisation states */
                  if(euler_phi.flag==0){euler_phi.value+=euler_phi.limits[1];} /* Increment the Euler angle phi */
               } /* End loop over Euler phi, f0 */
               if(euler_psi.flag==0){euler_psi.value+=euler_psi.limits[1];} /* Increment the Euler angle psi */
            } /* End loop over Euler psi, p0 */
            if(euler_theta.flag==0){euler_theta.value+=euler_theta.limits[1];} /* Increment the Euler angle theta */
         } /* End loop over Euler theta, t0 */
         if(radius.flag==0){radius.value+=radius.limits[1];} /* Increment the effective radius */
      } /* End loop over effective radii, er */
      if(wavelength.flag==0){wavelength.value+=wavelength.limits[1];} /* Increment the wavelength */
   } /* End loop over incident wavelengths, wl */
   for(j=0;j<refractive_index.number_of_materials;j++){
      fclose(output.refindex_vs_wave[j]); /* Close complex refractive index vs wavelength output files */
   }
   fclose(output.iter); /* Close iterative output file */
   fclose(output.target_descrip); /* Close target description output file */

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Destroy all FFTW DFT plans */
   dft_plan(0);

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Free allocated memory back to the heap */
   memory_free_OMP();

/* ->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->-> */
/* Finalise timing information */
   if(timing.enabled){timing_info();}

   return 0;
}

void initialise(int flag,int *list_length,double *limits){
/* If flag=0 then the code is using minimum,increment,maximum from the
   "control.input" file and list_length=[(maximum-minimum)/increment]+1
   N.B. If increment < breakdown then minimum=maximum and list_length=1*/

   if(flag==0){
      if(limits[1]<DBL_EPSILON){
         *list_length=1;
      }
      else{
         *list_length=(int)((limits[2]-limits[0])/limits[1]+1.0);
      }
   }
}
