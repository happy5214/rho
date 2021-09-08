/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#ifndef _FACTOR_H_
#define _FACTOR_H_

//support libraries
#include <gmp.h>
#include "types.h"

#define NUM_POLYS 3

typedef struct
{
	mpz_t gmp_n;
	mpz_t gmp_f;
	uint32 iterations;
	uint32 num_poly;
	uint32 *polynomials;
	uint32 curr_poly;			//current polynomial in the list of polynomials
	double ttime;

} rho_obj_t;

typedef struct
{
	rho_obj_t rho_obj;			//info for any rho work
} fact_obj_t;

#endif //_FACTOR_H
