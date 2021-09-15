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
#include "rhoTypes.h"

void init_factobj(fact_obj_t *fobj);
void free_factobj(fact_obj_t *fobj);
void alloc_factobj(fact_obj_t *fobj);

/*--------------DECLARATIONS FOR MANAGING FACTORS FOUND -----------------*/

//yafu
void add_to_factor_list(fact_obj_t *fobj, mpz_t n, FinishingState finishingState);
void print_factors(fact_obj_t *fobj);
void clear_factor_list(fact_obj_t *fobj);
void delete_from_factor_list(fact_obj_t *fobj, mpz_t n);

int is_mpz_prp(mpz_t n);

#endif //_FACTOR_H
