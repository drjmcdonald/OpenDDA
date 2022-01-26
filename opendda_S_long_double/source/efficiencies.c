/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Calculates the various efficiencies

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "efficiencies.h"

void efficiencies(void){

   /* Geometrical cross section pi*radius^{2} */
   cross_section.geometrical=ldpi*radius.value*radius.value;

   /* Q_{ext} Extinction efficiency */
   efficiency.extinction=cross_section.extinction/cross_section.geometrical;
   printf("Qext=%+.*Lg\n",LDBLP,efficiency.extinction);

   /* Q_{abs} Absorption efficiency */
   efficiency.absorption=cross_section.absorption/cross_section.geometrical;
   printf("Qabs=%+.*Lg\n",LDBLP,efficiency.absorption);

   /* Q_{sca} Scattering efficiency */
   efficiency.scattering=cross_section.scattering/cross_section.geometrical;
   printf("Qsca=%+.*Lg\n",LDBLP,efficiency.scattering);
}
