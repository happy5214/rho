/******************************************************************************
 * Common fact_obj_t code.
 *
 * Most of the code in this file is from yafu, which has been placed into the
 * public domain by its author, Ben Buhrow. The original comment from the yafu
 * version is reproduced below. The remainder has been placed under the CC0 1.0
 * Universal license. To view a copy of this license, visit
 * http://creativecommons.org/publicdomain/zero/1.0.
 *****************************************************************************/

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
       				   --bbuhrow@gmail.com 3/26/10
----------------------------------------------------------------------*/

#include "factor.h"

void init_factobj(fact_obj_t *fobj) {
	// initialize stuff for rho
	fobj->rho_obj.num_poly = NUM_POLYS;
	fobj->rho_obj.polynomials = (uint32 *)malloc(fobj->rho_obj.num_poly * sizeof(uint32));
	for (int i = 0; i < NUM_POLYS; i++) {
		fobj->rho_obj.polynomials[i] = NUM_POLYS - i;
	}
	mpz_init(fobj->rho_obj.gmp_n);
	mpz_init(fobj->rho_obj.gmp_f);
	fobj->rho_obj.iterations = MAX_ITERATIONS;
	fobj->rho_obj.curr_poly = 0;
}

void free_factobj(fact_obj_t *fobj)
{
	// free stuff in rho
	free(fobj->rho_obj.polynomials);
	mpz_clear(fobj->rho_obj.gmp_n);
	mpz_clear(fobj->rho_obj.gmp_f);
}
