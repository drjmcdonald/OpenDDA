/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Reset a string by setting all used elements to '\0'

   Copyright (C) 2006 James Mc Donald
   Computational Astrophysics Laboratory
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "reset_string.h"

void reset_string(char *string){

   int i,j;

   j=(int)strlen(string); /* Length of string */
   for(i=0;i<j;i++){ /* Reset string */
      string[i]='\0';
   }
}
