/******************************************************************************
 * Main module.
 *
 * Copyright 2017, 2019, 2021, Alexander Jones.
 *
 * Based on code from yafu, which has been placed into the public domain by its
 * author, Ben Buhrow.
 *
 * Portions of the main function are based on code from Arg_parser, licensed
 * under the 2-clause BSD license by its author, Antonio Diaz Diaz. Arg_parser
 * is included as carg_parser.c and carg_parser.h in this distribution, and the
 * full license can be found in those files.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; see the file LICENSE.  If not, see http://www.gnu.org/licenses/
 * or write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 ******************************************************************************/

#include <stdio.h>
#include <string.h>

#include <gmp.h>

#include "carg_parser.h"
#include "rho.h"
#include "types.h"

static int loop_count = LOOP_COUNT;

static bool only_one_poly = false;
static int single_poly;

int gcd_step = GCD_STEP;
int max_iterations = MAX_ITERATIONS;

uint32 *polys;
mpz_t constant;

static void rho_loop(fact_obj_t *fobj);
static bool rho_inner(fact_obj_t *fobj);

/**
 * Run rho algorithm.
 *
 * @param compositeNumber: The number to factor.
 * @return 0 on success
 */
static int rho(char *composite) {
	fact_obj_t fobj;
	init_factobj(&fobj);
	mpz_set_str(fobj.rho_obj.gmp_n, composite, 10);
	polys = fobj.rho_obj.polynomials;
	rho_loop(&fobj);
	print_factors(&fobj);
	free_factobj(&fobj);
	return 0;				// Always return 0 if there's no error
}

static void rho_loop(fact_obj_t *fobj) {
	if ((mpz_cmp_ui(fobj->rho_obj.gmp_n, 1) == 0) || (mpz_cmp_ui(fobj->rho_obj.gmp_n, 0) == 0))
		return;

	if (mpz_cmp_ui(fobj->rho_obj.gmp_n, 2) == 0)
		return;

	bool fullyFactored = false;
	if (only_one_poly) {
		fobj->rho_obj.curr_poly = NUM_POLYS - single_poly;
		rho_inner(fobj);		// Actually run rho algorithm (dependent on algorithm linked)
	} else {
		fobj->rho_obj.curr_poly = 0;		// An index for polys
		while (fobj->rho_obj.curr_poly < NUM_POLYS) {	// Loop through polynomials
			fullyFactored = rho_inner(fobj);		// Actually run rho algorithm (dependent on algorithm linked)
			if (fullyFactored) {		// We found a factor!
				break;
			}
			fobj->rho_obj.curr_poly++;
		}
	}
}

static bool rho_inner(fact_obj_t *fobj) {
	//for each different constant, first check primalty because each
	//time around the number may be different
	if (is_mpz_prp(fobj->rho_obj.gmp_n)) {
		FinishingState dummyState = {-1, -1};
		add_to_factor_list(fobj, fobj->rho_obj.gmp_n, dummyState);
		mpz_set_ui(fobj->rho_obj.gmp_n, 1);
		return true;
	}

	//call rho algorithm
	FinishingState finishingState = run_rho(fobj);

	//check to see if 'f' is non-trivial
	if ((mpz_cmp_ui(fobj->rho_obj.gmp_f, 1) > 0)
		&& (mpz_cmp(fobj->rho_obj.gmp_f, fobj->rho_obj.gmp_n) < 0)) {
		//non-trivial factor found

		add_to_factor_list(fobj, fobj->rho_obj.gmp_f, finishingState);

		//reduce input
		mpz_tdiv_q(fobj->rho_obj.gmp_n, fobj->rho_obj.gmp_n, fobj->rho_obj.gmp_f);
	}
	return false;
}

static const char * const program_year = "2023";

static void show_version() {
	printf("Rho Calculator\n");
	printf("Copyright (C) %s Alexander Jones.\n", program_year);
	printf("Based on yafu, which has been released into the public domain by Ben Buhrow.\n");
	printf("License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>\n"
			"This is free software: you are free to change and redistribute it.\n"
			"There is NO WARRANTY, to the extent permitted by law.\n");
}

/**
 * Start a rho run.
 *
 * @param argc: The number of arguments.
 * @param argv: The array of arguments.
 * @return Value returned by rho() (currently always 0)
 */
int main(const int argc, const char * const argv[]) {
	char composite[MAX_NUM_SIZE];
	*composite = '\0';

	// Legal command-line arguments
	const struct ap_Option options[] =
		{
		{ 'V', "version",    ap_no    },	// Display the version information
		{ 'p', "polynomial", ap_yes   },	// Use a specific polynomial
		{ 'g', "gcd-step",   ap_yes   },	// How many GCD calculations to merge at once
		{ 'i', "iterations", ap_yes   },	// The iterations limit for the algorithm
		{ 'l', "loop",       ap_yes   } };	// An optional number of times to repeat the whole thing

	// Parse the arguments
	struct Arg_parser parser;
	int argind;
	ap_init( &parser, argc, argv, options, false );

	for( argind = 0; argind < ap_arguments( &parser ); ++argind ) {
		const int code = ap_code( &parser, argind );
		const char * const arg = ap_argument( &parser, argind );
		switch (code) {
			case 'V': show_version(); return 0;
			case 'p': only_one_poly = true; single_poly = strtol(arg, NULL, 10); break;
			case 'g': gcd_step = strtol(arg, NULL, 10); break;
			case 'i': max_iterations = strtol(arg, NULL, 10); break;
			case 'l': loop_count = strtol(arg, NULL, 10); break;
			case '\0':
				strncpy(composite, arg, MAX_NUM_SIZE - 1);
				composite[MAX_NUM_SIZE - 1] = '\0';
				break;
		}
		if (!code) {
			break;
		}
	}

	if (!(*composite)) {
		fprintf(stderr, "No composite provided.\n");
		return 1;
	}

	// Tail-call the rho starter
	return rho(composite);
}
