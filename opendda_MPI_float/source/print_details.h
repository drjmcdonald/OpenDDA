/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for print_details

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void print_error(char *error_string0,char *error_string1,char *error_string2);
void print_out(FILE *file_ptr,char *print_string0,char *print_string1,char *print_string2);
void print_section_title(FILE *file_ptr,char *print_string0);
