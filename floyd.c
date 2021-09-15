/******************************************************************************
 * Pollard's rho algorithm using Floyd's cycle-finding algorithm.
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
	mpz_t x, y, curr_gcd, temp, f;

	uint32_t i, skip_counter, power;
	int iterations;
	int function_calls = 0;

	// Initialize local bigints
	mpz_init_set_ui(x, X_0);		// "Tortoise"
	mpz_init_set_ui(y, X_0);		// "Hare"
	mpz_init(temp);				// Temporary storage
	mpz_init(f);				// Found factor
	mpz_init_set_ui(constant, polys[c]);	// Constant in polynomial
	mpz_init_set_ui(curr_gcd, 1);		// Current GCD

	// Starting state of algorithm
	iterations = 0;				// Rho iteration count

	do {
		g(x, x, fobj->rho_obj.gmp_n, temp, &function_calls);

		for (i = 0; i < 2; i++) {
			g(y, y, fobj->rho_obj.gmp_n, temp, &function_calls);
		}

		mpz_sub(temp, x, y);
		mpz_abs(temp, temp);
		mpz_gcd(curr_gcd, temp, fobj->rho_obj.gmp_n);
		iterations++;
	} while (mpz_get_ui(curr_gcd) == 1 && iterations < max_iterations);
#if DEBUG
	int final_index = iterations * 2;
#endif

	if (mpz_cmp(curr_gcd, fobj->rho_obj.gmp_n) == 0 || mpz_get_ui(curr_gcd) == 1) {
		mpz_set_ui(f, 0);
	} else {
		mpz_set(f, curr_gcd);
	}

free:
	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(constant);
	mpz_clear(temp);
	mpz_clear(curr_gcd);
	mpz_set(fobj->rho_obj.gmp_f, f);
	mpz_clear(f);

#if DEBUG
	FinishingState returnValue = {function_calls, final_index};
#else
	FinishingState returnValue = {0, 0};
#endif
	return returnValue;
}
