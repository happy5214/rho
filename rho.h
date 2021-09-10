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

#define MIN(a,b) ((a) < (b)? (a) : (b))
#define MAX(a,b) ((a) > (b)? (a) : (b))

#define MAX_NUM_SIZE 200

#define C_MAX 10
#define X_0 0
#define G_CONSTANT 1
#define GCD_STEP 10

typedef enum { false = 0, true = 1 } bool;

typedef struct {
	int floyd;
	int brent;
	int yafu;
	int montgomery;
} MaxLength;

typedef struct {
	char number[MAX_NUM_SIZE];
	int aliquotSeq;
	short aliquotTerm;
	MaxLength maxLengths[NUM_POLYS];
} Composite;

int gcd_step, max_iterations;
uint32_t *polys;

#if DEBUG
int function_calls;
int final_index;
#endif

mpz_t constant;

static inline void g(mpz_t output, mpz_t input, mpz_t n, mpz_t temp) {
	mpz_mul(temp, input, input);
	mpz_add(temp, temp, constant);
	mpz_tdiv_r(output, temp, n);
#if DEBUG
	function_calls++;
#endif
}

void run_rho(fact_obj_t *fobj);

int readComposites(Composite **composites, char *filename);

void init_factobj(fact_obj_t *fobj);
void free_factobj(fact_obj_t *fobj);

#endif // RHO_H
