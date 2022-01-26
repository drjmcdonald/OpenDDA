/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function to print details to stderr or file in a consistent fashion

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "print_details.h"

void print_error(char *error_string0,char *error_string1,char *error_string2){
/* Print error message and exit */

   int i,length=0;
   char tilteerror[15]="OpenDDA Error:\0",symbol='O';

/* Find the length of the longest component of the error message */
   length=(int)strlen(tilteerror);
   length=(int)strlen(error_string0)>length?(int)strlen(error_string0):length;
   if(error_string1!=NULL){
      length=(int)strlen(error_string1)>length?(int)strlen(error_string1):length;
   }
   if(error_string2!=NULL){
      length=(int)strlen(error_string2)>length?(int)strlen(error_string2):length;
   }

   length+=4;
   fprintf(stderr,"\n\n");
   for(i=0;i<length;i++){ /* Print top line */
      fprintf(stderr,"%c",symbol);
   }
   fprintf(stderr,"\n");

   fprintf(stderr,"%c",symbol); /* Print empty line */
   for(i=0;i<(length-2);i++){
      fprintf(stderr," ");
   }
   fprintf(stderr,"%c\n",symbol);

   fprintf(stderr,"%c",symbol); /* Print title line */
   fprintf(stderr," %s",tilteerror);
   for(i=0;i<(length-3-(int)strlen(tilteerror));i++){
      fprintf(stderr," ");
   }
   fprintf(stderr,"%c\n",symbol);

   fprintf(stderr,"%c",symbol); /* Print empty line */
   for(i=0;i<(length-2);i++){
      fprintf(stderr," ");
   }
   fprintf(stderr,"%c\n",symbol);

   fprintf(stderr,"%c %s",symbol,error_string0); /* Print error string 0 */
   for(i=0;i<(length-3-(int)strlen(error_string0));i++){
      fprintf(stderr," ");
   }
   fprintf(stderr,"%c\n",symbol);

   fprintf(stderr,"%c",symbol); /* Print empty line */
   for(i=0;i<(length-2);i++){
      fprintf(stderr," ");
   }
   fprintf(stderr,"%c\n",symbol);

   if(error_string1!=NULL){ /* Print error string 1 if it exists */
      fprintf(stderr,"%c %s",symbol,error_string1);
      for(i=0;i<(length-3-(int)strlen(error_string1));i++){
         fprintf(stderr," ");
      }
      fprintf(stderr,"%c\n",symbol);
   }

   if(error_string2!=NULL){ /* Print error string 2 if it exists */
      fprintf(stderr,"%c %s",symbol,error_string2);
      for(i=0;i<(length-3-(int)strlen(error_string2));i++){
         fprintf(stderr," ");
      }
      fprintf(stderr,"%c\n",symbol);
   }

   fprintf(stderr,"%c",symbol); /* Print empty line */
   for(i=0;i<(length-2);i++){
      fprintf(stderr," ");
   }
   fprintf(stderr,"%c\n",symbol);

   for(i=0;i<length;i++){ /* Print bottom line */
      fprintf(stderr,"%c",symbol);
   }
   fprintf(stderr,"\n\n\n");

   exit(EXIT_FAILURE);
}

void print_out(FILE *file_ptr,char *print_string0,char *print_string1,char *print_string2){
/* Print header information for output files */
   int i,length=0;
   char title[22]="OpenDDA [Version 1.0]",symbol='O';

/* Find the length of the longest component of the message */
   length=(int)strlen(title);
   if(print_string0!=NULL){
      length=(int)strlen(print_string0)>length?(int)strlen(print_string0):length;
   }
   if(print_string1!=NULL){
      length=(int)strlen(print_string1)>length?(int)strlen(print_string1):length;
   }
   if(print_string2!=NULL){
      length=(int)strlen(print_string2)>length?(int)strlen(print_string2):length;
   }

   length+=4;
   fprintf(file_ptr,"\n ");
   for(i=0;i<length;i++){ /* Print top line */
      fprintf(file_ptr,"%c",symbol);
   }
   fprintf(file_ptr,"\n");

   fprintf(file_ptr," %c",symbol); /* Print empty line */
   for(i=0;i<(length-2);i++){
      fprintf(file_ptr," ");
   }
   fprintf(file_ptr,"%c\n",symbol);

   fprintf(file_ptr," %c",symbol); /* Print title line */
   fprintf(file_ptr," %s",title);
   for(i=0;i<(length-3-(int)strlen(title));i++){
      fprintf(file_ptr," ");
   }
   fprintf(file_ptr,"%c\n",symbol);

   fprintf(file_ptr," %c",symbol); /* Print empty line */
   for(i=0;i<(length-2);i++){
      fprintf(file_ptr," ");
   }
   fprintf(file_ptr,"%c\n",symbol);

   if(print_string0!=NULL){ /* Print string 0 if it exists */
      fprintf(file_ptr," %c %s",symbol,print_string0); /* Print string 0 */
      for(i=0;i<(length-3-(int)strlen(print_string0));i++){
         fprintf(file_ptr," ");
      }
      fprintf(file_ptr,"%c\n",symbol);
   }

   if(print_string1!=NULL){ /* Print string 1 if it exists */
      fprintf(file_ptr," %c %s",symbol,print_string1); /* Print string 1 */
      for(i=0;i<(length-3-(int)strlen(print_string1));i++){
         fprintf(file_ptr," ");
      }
      fprintf(file_ptr,"%c\n",symbol);
   }

   if(print_string2!=NULL){ /* Print string 2 if it exists */
      fprintf(file_ptr," %c %s",symbol,print_string2); /* Print string 2 */
      for(i=0;i<(length-3-(int)strlen(print_string2));i++){
         fprintf(file_ptr," ");
      }
      fprintf(file_ptr,"%c\n",symbol);
   }

   fprintf(file_ptr," %c",symbol); /* Print empty line */
   for(i=0;i<(length-2);i++){
      fprintf(file_ptr," ");
   }
   fprintf(file_ptr,"%c\n ",symbol);

   for(i=0;i<length;i++){ /* Print bottom line */
      fprintf(file_ptr,"%c",symbol);
   }
   fprintf(file_ptr,"\n");
}

void print_section_title(FILE *file_ptr,char *print_string0){
/* Print section title for output files */

   int i,length=0;
   char symbol='O';

/* Find the length of the longest component of the message */
   if(print_string0!=NULL){
      length=(int)strlen(print_string0);
   }

   length+=4;
   fprintf(file_ptr,"\n ");
   for(i=0;i<length;i++){ /* Print top line */
      fprintf(file_ptr,"%c",symbol);
   }

   if(print_string0!=NULL){ /* Print string 0 if it exists */
      fprintf(file_ptr,"\n %c %s",symbol,print_string0); /* Print string 0 */
      for(i=0;i<(length-3-(int)strlen(print_string0));i++){
         fprintf(file_ptr," ");
      }
      fprintf(file_ptr,"%c\n ",symbol);
   }

   for(i=0;i<length;i++){ /* Print bottom line */
      fprintf(file_ptr,"%c",symbol);
   }
   fprintf(file_ptr,"\n");
}
