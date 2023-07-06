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

#include <stdlib.h>

#include "factor.h"

void init_factobj(fact_obj_t *fobj)
{
	// get space for everything
	alloc_factobj(fobj);

	// initialize stuff for rho
	fobj->rho_obj.iterations = MAX_ITERATIONS;
	fobj->rho_obj.curr_poly = 0;
}

void free_factobj(fact_obj_t *fobj)
{
	// free stuff in rho
	free(fobj->rho_obj.polynomials);
	mpz_clear(fobj->rho_obj.gmp_n);
	mpz_clear(fobj->rho_obj.gmp_f);

	clear_factor_list(fobj);
	free(fobj->fobj_factors);
}

void alloc_factobj(fact_obj_t *fobj)
{
	int i;

	fobj->rho_obj.num_poly = NUM_POLYS;
	fobj->rho_obj.polynomials = (uint32 *)malloc(fobj->rho_obj.num_poly * sizeof(uint32));
	for (i = 0; i < NUM_POLYS; i++) {
		fobj->rho_obj.polynomials[i] = NUM_POLYS - i;
	}
	mpz_init(fobj->rho_obj.gmp_n);
	mpz_init(fobj->rho_obj.gmp_f);

	fobj->allocated_factors = 8;
	fobj->fobj_factors = (factor_t *)malloc(8 * sizeof(factor_t));
	for (i = 0; i < fobj->allocated_factors; i++)
	{
		fobj->fobj_factors[i].type = UNKNOWN;
		fobj->fobj_factors[i].count = 0;
	}

	fobj->num_factors = 0;

	return;
}

/*
 * This function is from arith3.c.
 */
int is_mpz_prp(mpz_t n)
{
	return mpz_probab_prime_p(n, 25);
}

void add_to_factor_list(fact_obj_t *fobj, mpz_t n, FinishingState finishingState)
{
	//stick the number n into the global factor list
	uint32 i;
	int found = 0, v = 0;

	if (fobj->num_factors >= fobj->allocated_factors)
	{
		fobj->allocated_factors *= 2;
		fobj->fobj_factors = (factor_t *)realloc(fobj->fobj_factors,
			fobj->allocated_factors * sizeof(factor_t));
	}

	//look to see if this factor is already in the list
	for (i=0; i < fobj->num_factors && !found; i++)
	{
		if (mpz_cmp(n, fobj->fobj_factors[i].factor) == 0)
		{
			found = 1;
			fobj->fobj_factors[i].count++;
			return;
		}
	}

	//else, put it in the list
	mpz_init(fobj->fobj_factors[fobj->num_factors].factor);
	mpz_set(fobj->fobj_factors[fobj->num_factors].factor, n);
	fobj->fobj_factors[fobj->num_factors].count = 1;
	if (is_mpz_prp(n))
	{
		if (mpz_cmp_ui(n, 100000000) < 0)
			fobj->fobj_factors[fobj->num_factors].type = PRIME;
		else
			fobj->fobj_factors[fobj->num_factors].type = PRP;
	}
	else
		fobj->fobj_factors[fobj->num_factors].type = COMPOSITE;

	fobj->fobj_factors[fobj->num_factors].finishingState = finishingState;
	fobj->fobj_factors[fobj->num_factors].polynomial = fobj->rho_obj.curr_poly;
	fobj->num_factors++;

	return;
}

void delete_from_factor_list(fact_obj_t *fobj, mpz_t n)
{
	//remove the number n from the global factor list
	uint32 i;

	//find the factor
	for (i=0;i<fobj->num_factors; i++)
	{
		if (mpz_cmp(n,fobj->fobj_factors[i].factor) == 0)
		{
			int j;
			// copy everything above this in the list back one position
			for (j=i; j<fobj->num_factors-1; j++)
			{
				mpz_set(fobj->fobj_factors[j].factor, fobj->fobj_factors[j+1].factor);
				fobj->fobj_factors[j].count = fobj->fobj_factors[j+1].count;
			}
			// remove the last one in the list
			fobj->fobj_factors[j].count = 0;
			mpz_clear(fobj->fobj_factors[j].factor);

			fobj->num_factors--;
			break;
		}
	}

	return;
}

void clear_factor_list(fact_obj_t *fobj)
{
	uint32 i;

	//clear this info
	for (i=0; i<fobj->num_factors; i++)
	{
		fobj->fobj_factors[i].count = 0;
		mpz_clear(fobj->fobj_factors[i].factor);
	}
	fobj->num_factors = 0;

	return;
}

static void print_factor(fact_obj_t *fobj, factor_t factor);

void print_factors(fact_obj_t *fobj)
{
	uint32 i, j;

    for (i = 0; i < fobj->num_factors; i++)
	{
		for (j = 0; j < fobj->fobj_factors[i].count;j++)
		{
			print_factor(fobj, fobj->fobj_factors[i]);
		}
	}

	if (mpz_cmp_ui(fobj->rho_obj.gmp_n, 1) > 0)
	{
#if DEBUG
		gmp_printf("Cofactor: %Zd\n", fobj->rho_obj.gmp_n);
#else
		gmp_printf("%Zd\n", fobj->rho_obj.gmp_n);
#endif
	}
}

static void print_factor(fact_obj_t *fobj, factor_t factor) {
#if DEBUG
	gmp_printf("Factor: %Zd\n", factor.factor);			// Print the factor
	printf("Polynomial: x^2+%d\n", fobj->rho_obj.polynomials[factor.polynomial]);	// Print the polynomial used
	printf("Ending index: %d\n", factor.finishingState.final_index);		// Print the ending index
	printf("Function calls: %d\n", factor.finishingState.function_calls);		// Print the number of function calls
#else
	gmp_printf("%Zd\n", factor.factor);			// Print the factor
#endif
}
