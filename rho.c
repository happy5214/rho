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

#include "carg_parser.h"
#include "rho.h"

// The array of composites used as test data.
char *compositeNums[] = {"108279571507399440659395651945074688076362401221378706580648545958263360543467201098809802671699082208985556979828647789",
// Factor: 70779587 (yafu: 7345, floyd: 6492, mcf: 7346)
"4606167572662656995739741462385108329495475747757164040978611808463116059373031233073985683826684323406449922365172139669",
// Factor: 1502687 (yafu: 3337, floyd: 2576, mcf: 3338)
"31495110750042400109854812254550375129055121194217936834507224183974466104881629618957416613947399087099978098510182327",
// Factor: None
"8998945954365211872301965528648321447489390217286011554329542928047815512233706999895496926187837189129631709094897681394591",
// Factor: None
"8889499927206867186591563004160081593161829693237104303925676323878532410062555813650853069890284719",
// Factor: None
"5003154769843638353453740563696190734503519409422149272769455573992592523175954553518500340172963423158925596801711",
// 236840:763 / Factor: 8420273 (yafu: 7995, floyd: 3220, mcf: 3662)
"594179638812617875151285541893498077141147253708062585710636172246742180826673262674321882458319750815552607",
// 236840:763 / Factor: 310429271 (yafu: 33067, floyd: 38016, mcf: 25250)
"447992281473092503403639"
// Factor: 89392903
};

static int loop_count = LOOP_COUNT;

static bool only_one_poly = false;
static int single_poly;

int gcd_step = GCD_STEP;
int max_iterations = MAX_ITERATIONS;

/**
 * Run rho loop.
 *
 * @param comp_index: The index in composites chosen for the test.
 * @param import: Whether to import composites from a file.
 * @param importFilename: The file to import from.
 * @return 0 on success
 */
static int rho(int comp_index, bool import, char *importFilename) {
	fact_obj_t fobj;
	init_factobj(&fobj);
	polys = fobj.rho_obj.polynomials;
	char *composite;
	if (import) {
		Composite *composites;
		int compositeCount = readComposites(&composites, importFilename);
		if (comp_index >= compositeCount)
			return 1;
		composite = composites[comp_index].number;
		free(composites);
	} else {
		composite = compositeNums[comp_index];
	}
	mpz_set_str(fobj.rho_obj.gmp_n, composite, 10);
	uint64_t factor;			// The factor found by rho
	for (int i = 0; i < loop_count; i++) {	// Loop m times, where m is given on command line or 1
		if (only_one_poly) {
			fobj.rho_obj.curr_poly = NUM_POLYS - single_poly;
			run_rho(&fobj);		// Actually run rho algorithm (dependent on algorithm linked)
			factor = mpz_get_ui(fobj.rho_obj.gmp_f);
		} else {
			fobj.rho_obj.curr_poly = 0;		// An index for polys
			while (fobj.rho_obj.curr_poly < NUM_POLYS) {	// Loop through polynomials
				run_rho(&fobj);		// Actually run rho algorithm (dependent on algorithm linked)
				factor = mpz_get_ui(fobj.rho_obj.gmp_f);
				if (factor) {		// We found a factor!
					break;
				}
				fobj.rho_obj.curr_poly++;
			}
		}
	}
#if DEBUG
	printf("Factor: %ld\n", factor);			// Print the factor
	if (fobj.rho_obj.curr_poly != NUM_POLYS) {
		printf("Polynomial: x^2+%d\n", fobj.rho_obj.polynomials[fobj.rho_obj.curr_poly]);	// Print the polynomial used
	}
	printf("Ending index: %d\n", final_index);		// Print the ending index
	printf("Function calls: %d\n", function_calls);		// Print the ending index
#else
	printf("%ld\n", factor);			// Print the factor
#endif
	free_factobj(&fobj);
	return 0;				// Always return 0 if there's no error
}

/**
 * Start a rho run.
 *
 * @param argc: The number of arguments.
 * @param argv: The array of arguments.
 * @return Value returned by rho() (currently always 0)
 */
int main(const int argc, const char * const argv[]) {
	int comp_index = 0;				// A composite index
	bool import = false;
	char importFilename[30];

	// Legal command-line arguments
	const struct ap_Option options[] =
		{
		{ 'I', "import",     ap_maybe },	// Whether to import the composites from a file
		{ 'p', "polynomial", ap_yes   },	// Use a specific polynomial
		{ 'c', "composite",  ap_yes   },	// An index for the composite array
		{ 'g', "gcd-step",   ap_yes   },	// How many GCD calculations to merge at once
		{ 'i', "iterations", ap_yes   },	// The iterations limit for the algorithm
		{ 'l', "loop",       ap_yes   } };	// An optional number of times to repeat the whole thing

	// Parse the arguments
	struct Arg_parser parser;
	int argind;
	ap_init( &parser, argc, argv, options, 0 );

	for( argind = 0; argind < ap_arguments( &parser ); ++argind ) {
		const int code = ap_code( &parser, argind );
		const char * const arg = ap_argument( &parser, argind );
		switch (code) {
			case 'I':
				import = true;
				if (arg[0]) {
					strncpy(importFilename, arg, 30);
				} else {
					strcpy(importFilename, "composites.txt");
				}
				break;
			case 'p': only_one_poly = true; single_poly = strtol(arg, NULL, 10); break;
			case 'c': comp_index = strtol(arg, NULL, 10); break;
			case 'g': gcd_step = strtol(arg, NULL, 10); break;
			case 'i': max_iterations = strtol(arg, NULL, 10); break;
			case 'l': loop_count = strtol(arg, NULL, 10); break;
		}
	}

	// Tail-call the rho starter
	return rho(comp_index, import, importFilename);
}
