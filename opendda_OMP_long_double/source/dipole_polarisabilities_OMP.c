/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Compute the dipole polarisabilities and invert them
   as the diagonal of the interaction matrix is set to
   1/polarisability

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "dipole_polarisabilities_OMP.h"

void dipole_polarisabilities_OMP(void){

   int j,k,vc,material,shift;
   long int array0,index0i,index0k,element;
   long double c1,c2,c3,b1,b2,b3,kd,d3,S;
   ldcomplex m2,current_polarisation[3],term0,term1,alpha_cm,temp;

   c1=-5.942421945054091L;
   c2=0.517881857406567L;
   c3=4.006974687341422L;
   b1=c1/ldpi;
   b2=c2/ldpi;
   b3=-(3.0L*c2+c3)/ldpi;
/* b1=-1.891531652986228,b2=0.164846915087734,b3=-1.770000401932182 */
/* b1=-1.8915316,b2=0.1648469,b3=-1.7700004
   Values from Beyond "Clausius-Mossotti - Wave propagation on a polarizable point lattice
   and the Discrete Dipole Approximation",
   B. T. Draine and Jeremy Goodman, 
   The Astrophysical Journal, Volume 405, Number 2, Pages 685-697, 1993 */
/* b1=-1.8915316L;
   b2=0.1648469L;
   b3=-1.7700004L; */

   if(polarisation_state==0){ /* Incident polarisation e0 */ 
      current_polarisation[0]=incident_TF.polarisation[0]; /* e0x in the target frame */
      current_polarisation[1]=incident_TF.polarisation[1]; /* e0y in the target frame */
      current_polarisation[2]=incident_TF.polarisation[2]; /* e0z in the target frame */
   }
   else if(polarisation_state==1){ /* Orthonormal polarisation e1 */
      current_polarisation[0]=incident_TF.orthonormal[0]; /* e1x in the target frame */
      current_polarisation[1]=incident_TF.orthonormal[1]; /* e1y in the target frame */
      current_polarisation[2]=incident_TF.orthonormal[2]; /* e1z in the target frame */
   }
   else{
      print_error("Dipole polarisability calculation error","Incident polarisation error","The current polarisation state MUST be 0 OR 1");
   }

   kd=wavenumber*target.dipole_spacing; /* kd */
   d3=target.dipole_spacing*target.dipole_spacing*target.dipole_spacing; /* d3=d^{3} */

   /* Calculate S=sum_{i}(incident_TF.n_inc[i]current_polarisation[i])^{2} */
   S=0.0L;
   for(vc=0;vc<3;vc++){
      S+=(incident_TF.n_inc[vc]*ldcomplex_abs(current_polarisation[vc]))*
         (incident_TF.n_inc[vc]*ldcomplex_abs(current_polarisation[vc]));
   }

   if(refractive_index.number_of_materials>1){ /* Composition array is only used
      if more than one material is used for target construction */
      /* alpha_LDR=alpha_cm/(1+[alpha_cm/d^{3}][(b1+m^{2}b2+m^{2}b3S)(kd)^{2}-(2/3)i(kd)^{3}) */
      term1=ldcomplex_scale(onei,((2.0L/3.0L)*kd*kd*kd)); /* (2/3)i(kd)^{3} */
      for(vc=0;vc<3;vc++){
         array0=vc*target.Nd;
#pragma omp parallel for private(index0k,element,j,index0i,shift,material,m2,alpha_cm,term0) schedule(dynamic)
         for(k=0;k<target.P;k++){
            index0k=k*target.M;
            element=target.populated[k]; /* Start index for plane k for the stored non-zero lattice sites */
            for(j=0;j<target.M;j++){
               index0i=index0k+j;
               if(target.occupied[(int)((double)index0i/32.0)]&(1<<(index0i%32))){ /* Lattice site is occupied */
                  /* Calculate the Clausius-Mossotti relation, alpha_cm=(3d^{3}(m^{2}-1))/(4pi(m^{2}+2)), where
                     m is the complex refractive index */
               /* Extract the material composition from the dynamically bitfielded composition array */
                  shift=(index0i%target.masks_per_integer)*target.mask_size;
                  material=(target.composition[vc][(int)((double)index0i/target.masks_per_integer)]>>shift)&target.mask;
                  m2=ldcomplex_mul(refractive_index.value[material],refractive_index.value[material]); /* m^{2} */
                  alpha_cm=ldcomplex_scale(ldcomplex_div(ldcomplex_sub(m2,one),
                           ldcomplex_add(m2,two)),((3.0L*d3)/(4.0L*ldpi))); /* Clausius-Mossotti polarisability */
                  term0=ldcomplex_scale(ldcomplex_add(ldcomplex_scale(one,b1),
                           ldcomplex_scale(m2,(b2+b3*S))),(kd*kd)); /* [b1+(b2+b3*S)m^{2}](kd)^{2} */
                  /* alpha_LDR=alpha_cm/(1+[alpha_cm/d^{3}][(b1+m^{2}b2+m^{2}b3S)(kd)^{2}-(2/3)i(kd)^{3})
                     Since the inverse of a diagonal matrix is
                                _     _          _                    _ 
                               | a 0 0 |        | a^{-1}   0      0    |
                           A = | 0 b 0 | A^{-1}=|   0    b^{-1}   0    |
                               |_0 0 c_|        |_  0      0    c^{-1}_|
   
                     1/alpha_LDR=(1+[alpha_cm/d^{3}][(b1+m^{2}b2+m^{2}b3S)(kd)^{2}-(2/3)i(kd)^{3})/alpha_cm */
                  interaction_matrix_diagonal[array0+element]=ldcomplex_div(ldcomplex_add(one,
                        ldcomplex_mul(ldcomplex_scale(alpha_cm,(1.0L/d3)),ldcomplex_sub(term0,term1))),alpha_cm);
                  element++;
               }
            }
         }
      }
   }
   else{ /* Target is homogeneous and isotropic */
      /* alpha_LDR=alpha_cm/(1+[alpha_cm/d^{3}][(b1+m^{2}b2+m^{2}b3S)(kd)^{2}-(2/3)i(kd)^{3}) */
      term1=ldcomplex_scale(onei,((2.0L/3.0L)*kd*kd*kd)); /* (2/3)i(kd)^{3} */
      /* Calculate the Clausius-Mossotti relation, alpha_cm=(3d^{3}(m^{2}-1))/(4pi(m^{2}+2)), where
         m is the complex refractive index */
      m2=ldcomplex_mul(refractive_index.value[0],refractive_index.value[0]); /* m^{2} */
      alpha_cm=ldcomplex_scale(ldcomplex_div(ldcomplex_sub(m2,one),
               ldcomplex_add(m2,two)),((3.0L*d3)/(4.0L*ldpi))); /* Clausius-Mossotti polarisability */
      term0=ldcomplex_scale(ldcomplex_add(ldcomplex_scale(one,b1),
               ldcomplex_scale(m2,(b2+b3*S))),(kd*kd)); /* [b1+(b2+b3*S)m^{2}](kd)^{2} */
      /* alpha_LDR=alpha_cm/(1+[alpha_cm/d^{3}][(b1+m^{2}b2+m^{2}b3S)(kd)^{2}-(2/3)i(kd)^{3})
         Since the inverse of a diagonal matrix is
                    _     _          _                    _ 
                   | a 0 0 |        | a^{-1}   0      0    |
               A = | 0 b 0 | A^{-1}=|   0    b^{-1}   0    |
                   |_0 0 c_|        |_  0      0    c^{-1}_|
   
         1/alpha_LDR=(1+[alpha_cm/d^{3}][(b1+m^{2}b2+m^{2}b3S)(kd)^{2}-(2/3)i(kd)^{3})/alpha_cm */
      temp=ldcomplex_div(ldcomplex_add(one,ldcomplex_mul(ldcomplex_scale(alpha_cm,(1.0L/d3)),
                        ldcomplex_sub(term0,term1))),alpha_cm);
#pragma omp parallel for schedule(dynamic,target.M)
      for(index0i=0;index0i<target.Nd3;index0i++){
         interaction_matrix_diagonal[index0i]=temp;
      }
   }
}
