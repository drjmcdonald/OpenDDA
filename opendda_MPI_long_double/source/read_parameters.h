/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for read_parameters

   Copyright (C) 2006 James Mc Donald
   Computational Astrophysics Laboratory
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "ldcomplex_alloc.h"
#include "ldcomplex_math_MPI.h"
#include "definitions.h"
#include "long_double_alloc.h"
#include "file_alloc.h"
#include "int_alloc.h"
#include "print_details.h"
#include "reset_string.h"

void read_parameters(void);
void read_from_file(FILE **file_ptr,attributes_general *structure_ptr,char *file_string,char *data_string);
void read_from_file_refractive_index_vs_wavelength(FILE **file_ptr,attributes_refractive_index *structure_ptr,char *file_string,int material);
void get_domain_size(FILE **file_ptr);
