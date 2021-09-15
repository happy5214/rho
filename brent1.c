/******************************************************************************
 * Pollard's rho algorithm using Brent's cycle-finding algorithm.
 *
 * This version does not skip indices k for 2^i < k < 3*2^(i-1), as suggested
 * by Brent in his paper. This is an oversight corrected in brent2.c.
 *
 * Copyright 2017, 2021, Alexander Jones.
 *
 * Based on code from yafu, which has been placed into the public domain by its
 * author, Ben Buhrow.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; see the file LICENSE.  If not, see http://www.gnu.org/licenses/
 * or write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 ******************************************************************************/

#include "rho.h"

FinishingState run_rho(fact_obj_t *fobj) {
	uint32 c = fobj->rho_obj.curr_poly;
	mpz_t x, y, product, curr_gcd, temp, f;

	uint32_t i, skip_counter, power;
	int iterations;
	FinishingState finishingState = {0, 0, 0};

	// Initialize local bigints
	mpz_init(x);				// "Tortoise"
	mpz_init_set_ui(y, X_0);		// "Hare"
	mpz_init(temp);				// Temporary storage
	mpz_init(f);				// Found factor
	mpz_init_set_ui(constant, polys[c]);	// Constant in polynomial
	mpz_init_set_ui(curr_gcd, 1);		// Current GCD

	// Starting state of algorithm
	power = 1;				// Current power of two
	i = 0;					// Loop counter
	iterations = 0;				// Rho iteration count

	do {
		mpz_set(x, y);

		skip_counter = 0;
		do {
			square(y, y);

			mpz_sub(temp, x, y); //q = q*abs(x-y) mod n
			mpz_abs(temp, temp);
			mpz_gcd(curr_gcd, temp, fobj->rho_obj.gmp_n);
			iterations++;
			skip_counter++;
		} while (skip_counter < power && mpz_get_ui(curr_gcd) == 1 && iterations < max_iterations);
		power *= 2;
	} while (mpz_get_ui(curr_gcd) == 1 && iterations < max_iterations);

	if (mpz_cmp(curr_gcd, fobj->rho_obj.gmp_n) == 0 || mpz_get_ui(curr_gcd) == 1) {
		mpz_set_ui(f, 0);
	} else {
		mpz_set(f, curr_gcd);
	}
#if DEBUG
	finishingState.final_index = iterations;
#endif

free:
	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(constant);
	mpz_clear(temp);
	mpz_clear(curr_gcd);
	mpz_set(fobj->rho_obj.gmp_f, f);
	mpz_clear(f);

	return finishingState;
}
