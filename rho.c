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
 * Load the composite number to test.
 *
 * @param composite: The buffer to store the composite number in.
 * @param comp_index: The index in composites chosen for the test.
 * @param import: Whether to import composites from a file.
 * @param importFilename: The file to import from.
 * @param compositeNumber: A provided composite number.
 * @return 0 on success
 */
static bool loadComposite(char *composite, int comp_index, bool import, char *importFilename, char *compositeNumber) {
	if (compositeNumber[0]) {
		strncpy(composite, compositeNumber, MAX_NUM_SIZE - 1);
	} else if (import) {
		Composite *composites;
		int compositeCount = readComposites(&composites, importFilename);
		if (comp_index >= compositeCount)
			return false;
		strncpy(composite, composites[comp_index].number, MAX_NUM_SIZE - 1);
		free(composites);
	} else {
		strncpy(composite, compositeNums[comp_index], MAX_NUM_SIZE - 1);
	}
	composite[MAX_NUM_SIZE - 1] = '\0';
	return true;
}

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
	char composite[MAX_NUM_SIZE];
	char compositeNumber[MAX_NUM_SIZE];
	compositeNumber[0] = '\0';
	char importFilename[30];

	// Legal command-line arguments
	const struct ap_Option options[] =
		{
		{ 'I', "import",     ap_maybe },	// Whether to import the composites from a file
		{ 'p', "polynomial", ap_yes   },	// Use a specific polynomial
		{ 'c', "composite",  ap_yes   },	// An index for the composite array
		{ 'n', "number",     ap_yes   },	// A specific number to factor
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
					strncpy(importFilename, arg, 30 - 1);
					importFilename[29] = '\0';
				} else {
					strcpy(importFilename, "composites.txt");
				}
				break;
			case 'n':
				strncpy(compositeNumber, arg, MAX_NUM_SIZE - 1);
				compositeNumber[MAX_NUM_SIZE - 1] = '\0';
				break;
			case 'p': only_one_poly = true; single_poly = strtol(arg, NULL, 10); break;
			case 'c': comp_index = strtol(arg, NULL, 10); break;
			case 'g': gcd_step = strtol(arg, NULL, 10); break;
			case 'i': max_iterations = strtol(arg, NULL, 10); break;
			case 'l': loop_count = strtol(arg, NULL, 10); break;
		}
	}

	if (!loadComposite(composite, comp_index, import, importFilename, compositeNumber)) {
        return 1;
    }
	// Tail-call the rho starter
	return rho(composite);
}
