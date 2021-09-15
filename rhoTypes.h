/******************************************************************************
 * Yafu and custom typedefs and enums.
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

#ifndef RHOTYPES_H
#define RHOTYPES_H 1

#include <gmp.h>

#include "types.h"

/*-------------------------------CUSTOM----------------------------------*/

#define MAX_NUM_SIZE 200
#define NUM_POLYS 3

typedef enum { false = 0, true = 1 } bool;

typedef enum { PRIME = 0, PRP = 1, COMPOSITE = 2, UNKNOWN = 3 } FactorType;

typedef struct {
	int final_index;
	int function_calls;
	int gcd_calls;
} FinishingState;

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

/*--------------------------FROM YAFU.H----------------------------------*/

typedef struct
{
	mpz_t factor;
	int count;
	FactorType type;
	FinishingState finishingState;
	uint32 polynomial;
} factor_t;

/*-------------------------FROM FACTOR.H---------------------------------*/

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

	//global storage for a list of factors
	factor_t *fobj_factors;
	uint32 num_factors;
	uint32 allocated_factors;
} fact_obj_t;

#endif // RHOTYPES_H
