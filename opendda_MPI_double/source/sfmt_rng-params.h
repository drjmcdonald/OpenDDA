#ifndef SFMT_PARAMS_H
#define SFMT_PARAMS_H

/* Added by James Mc Donald 13-Aud-2007 */
#define MEXP 19937

#if !defined(MEXP)
#ifdef __GNUC__
  #warning "MEXP is not defined. I assume MEXP is 19937."
#endif
  #define MEXP 19937
#endif
/*-----------------
  BASIC DEFINITIONS
  -----------------*/
/** Mersenne Exponent. The period of the sequence 
 *  is a multiple of 2^MEXP-1.
 * #define MEXP 19937 */
/** SFMT generator has an internal state array of 128-bit integers,
 * and N is its size. */
#define N (MEXP / 128 + 1)
/** N32 is the size of internal state array when regarded as an array
 * of 32-bit integers.*/
#define N32 (N * 4)
/** N64 is the size of internal state array when regarded as an array
 * of 64-bit integers.*/
#define N64 (N * 2)

/*----------------------
  the parameters of SFMT
  following definitions are in paramsXXXX.h file.
  ----------------------*/
/** the pick up position of the array.
#define POS1 122 
*/

/** the parameter of shift left as four 32-bit registers.
#define SL1 18
 */

/** the parameter of shift left as one 128-bit register. 
 * The 128-bit integer is shifted by (SL2 * 8) bits. 
#define SL2 1 
*/

/** the parameter of shift right as four 32-bit registers.
#define SR1 11
*/

/** the parameter of shift right as one 128-bit register. 
 * The 128-bit integer is shifted by (SL2 * 8) bits. 
#define SR2 1 
*/

/** A bitmask, used in the recursion.  These parameters are introduced
 * to break symmetry of SIMD.
#define MSK1 0xdfffffefU
#define MSK2 0xddfecb7fU
#define MSK3 0xbffaffffU
#define MSK4 0xbffffff6U 
*/

/** These definitions are part of a 128-bit period certification vector.
#define PARITY1	0x00000001U
#define PARITY2	0x00000000U
#define PARITY3	0x00000000U
#define PARITY4	0xc98e126aU
*/

#if MEXP == 607
  #include "sfmt_rng-params607.h"
#elif MEXP == 1279
  #include "sfmt_rng-params1279.h"
#elif MEXP == 2281
  #include "sfmt_rng-params2281.h"
#elif MEXP == 4253
  #include "sfmt_rng-params4253.h"
#elif MEXP == 11213
  #include "sfmt_rng-params11213.h"
#elif MEXP == 19937
  #include "sfmt_rng-params19937.h"
#elif MEXP == 44497
  #include "sfmt_rng-params44497.h"
#elif MEXP == 86243
  #include "sfmt_rng-params86243.h"
#elif MEXP == 132049
  #include "sfmt_rng-params132049.h"
#elif MEXP == 216091
  #include "sfmt_rng-params216091.h"
#else
#ifdef __GNUC__
  #error "MEXP is not valid."
  #undef MEXP
#else
  #undef MEXP
#endif

#endif

#endif /* SFMT_PARAMS_H */
