/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Initialise the bitfielded arrays for specifying whether the
   the lattice site is occupied OR unoccupied and the x, y and
   z composition

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License */
#include "bitfielding.h"

void bitfielding(void){

   int vc;

   /* Setup bitfielded occupied/unoccupied flag array */ 
   target.occ_ctrl=(int)((double)(target.N-1)/32.0+1.0);
   uint_calloc(&target.occupied,(size_t)target.occ_ctrl);

   if(refractive_index.number_of_materials>1){ /* Composition array is only required
      if more than one material is used for target construction */
      /* Setup dynamic bitfielded array for material composition */
      /* Find mask size */
      target.mask_size=1;
      while((1<<target.mask_size)<refractive_index.number_of_materials){
         target.mask_size++;
      }
      target.mask=(1<<target.mask_size)-1; /* Actual bitfielding extraction mask */

      /* Find masks per integer i.e. how many numbers we are storing in the integer */
      target.masks_per_integer=(int)(32.0/(double)target.mask_size);

      /* Calculate how many integers are required to store the N data values */
      target.mat_ctrl=(int)(((double)(target.N-1)/(double)target.masks_per_integer)+1.0);

      /* Allocate the memory for the x, y & z composition data */
      for(vc=0;vc<3;vc++){
         uint_calloc(&target.composition[vc],(size_t)target.mat_ctrl);
      }
   }

/* NOTES on dynamic bitfielding:

   (1) To set the x value for dipole `i' to `x_value' left-shift `x_value'
       by the appropriate amount

          shift=(index0i%masks_per_integer)*mask_size
          x_value<<shift

       and then use LOGICAL OR to set

          composition[0][(int)((double)i/(double)masks_per_integer)]|=(x_value<<shift)

       where

          (int)((double)i/(double)masks_per_integer) is the integer in which the value is to be bitfielded

   (2) To extract the x value for dipole `i' right-shift the mask by the appropriate amount

          shift=(index0i%masks_per_integer)*mask_size
          composition[0][(int)((double)index0i/masks_per_integer)]>>shift

       and the use LOGICAL AND with the mask to extract

          (composition[0][(int)((double)index0i/masks_per_integer)]>>shift)&mask

       where

          (int)((double)i/(double)masks_per_integer) is the integer in which the value is to be bitfielded

   Example: 13 dipoles 0..12 with 8 possible material values 0..7
            Want to store the following

            Dipole       x material       y material       z material
              0                 3                1                3
              1                1                4                4
              2                7                3                5
              3                3                7                1
            * 4               5                1                2
              5                6                0                6
              6                2                2                7
              7                0                5                0
              8                4                6                2
              9                3                6                4
              10                7                3                3
              11                5                4                1
              12                6                0                5

      (a) There are 8 possible material values 0..7 which can be stored in 3bits 000...111 in binary
      (b) There are 13 dipoles which means that (13)(3bits)=39bits is required for each of the x,y & z
          directions i.e. 2 integers:
          Integer to store x, y & z data 32bits:[integer_1][integer_0]
Dipole -------------------------------- 12  11  10       9   8   7   6   5   4   3   2   1   0 
 x dirn [XX XXX XXX XXX XXX XXX XXX XXX 000 000 000][XX 000 000 000 000 000 000 000 000 000 000]
 y dirn [XX XXX XXX XXX XXX XXX XXX XXX 000 000 000][XX 000 000 000 000 000 000 000 000 000 000]
 z dirn [XX XXX XXX XXX XXX XXX XXX XXX 000 000 000][XX 000 000 000 000 000 000 000 000 000 000]

      mask_size=3
      mask=(1<<3)-1=2^{3}-1=8-1=7=111 in binary
      masks_per_integer=(int)(32.0/(double)mask_size)=(int)(32.0/3.0)=10

      (c) To set x, y & z data for dipole i=4 i.e. 5, 1, 2, see * above

         5=101
         1=001
         2=010

         x_value<<((i%masks_per_integer)*mask_size)=5<<((4%10)*3)=5<<12=101 000 000 000 000
            [XX 000 000 000 000 000 000 000 000 000 000]|[101 000 000 000 000]=
            [XX 000 000 000 000 000 101 000 000 000 000]

         y_value<<((i%masks_per_integer)*mask_size)=1<<((4%10)*3)=1<<12=001 000 000 000 000
            [XX 000 000 000 000 000 000 000 000 000 000]|[001 000 000 000 000]=
            [XX 000 000 000 000 000 001 000 000 000 000]

         z_value<<((i%masks_per_integer)*mask_size)=2<<((4%10)*3)=2<<12=010 000 000 000 000
            [XX 000 000 000 000 000 000 000 000 000 000]|[010 000 000 000 000]=
            [XX 000 000 000 000 000 010 000 000 000 000]

Dipole -------------------------------- 12  11  10       9   8   7   6   5   4   3   2   1   0 
 x dirn [XX XXX XXX XXX XXX XXX XXX XXX 000 000 000][XX 000 000 000 000 000 101 000 000 000 000]
 y dirn [XX XXX XXX XXX XXX XXX XXX XXX 000 000 000][XX 000 000 000 000 000 001 000 000 000 000]
 z dirn [XX XXX XXX XXX XXX XXX XXX XXX 000 000 000][XX 000 000 000 000 000 010 000 000 000 000]

      (d) To extract x, y & z data for dipole i=4

          Right-shift the target array and USE LOGICAL AND with extraction mask

          x: {[XX 000 000 000 000 000 101 000 000 000 000]>>((4%10)*3)}&(111)=
             {[XX 000 000 000 000 000 101 000 000 000 000]>>12}&(111)=101=5

          y: {[XX 000 000 000 000 000 001 000 000 000 000]>>((4%10)*3)}&(111)=
             {[XX 000 000 000 000 000 001 000 000 000 000]>>12}&(111)=001=1

          z: {[XX 000 000 000 000 000 010 000 000 000 000]>>((4%10)*3)}&(111)=
             {[XX 000 000 000 000 000 010 000 000 000 000]>>12}&(111)=010=2 */
}
