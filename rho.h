/******************************************************************************
 * Main header file.
 *
 * Copyright 2017, 2019, 2021, Alexander Jones.
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

#ifndef RHO_H
#define RHO_H 1

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include "factor.h"

#include "rhoTypes.h"

#define MIN(a,b) ((a) < (b)? (a) : (b))
#define MAX(a,b) ((a) > (b)? (a) : (b))

#define C_MAX 10
#define X_0 0
#define G_CONSTANT 1
#define GCD_STEP 10

extern int gcd_step, max_iterations;
extern uint32 *polys;
extern mpz_t constant;

#if DEBUG
#define square(out,in) (g((out), (in), fobj->rho_obj.gmp_n, temp, &finishingState))
#else
#define square(out,in) (g((out), (in), fobj->rho_obj.gmp_n, temp))
#endif

#if DEBUG
static inline void g(mpz_t output, mpz_t input, mpz_t n, mpz_t temp, FinishingState *finishingState) {
#else
static inline void g(mpz_t output, mpz_t input, mpz_t n, mpz_t temp) {
#endif
	mpz_mul(temp, input, input);
	mpz_add(temp, temp, constant);
	mpz_tdiv_r(output, temp, n);
#if DEBUG
	finishingState->function_calls++;
#endif
}

FinishingState run_rho(fact_obj_t *fobj);

#endif // RHO_H
