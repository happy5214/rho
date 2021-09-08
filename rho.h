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

int final_index;
int gcd_step, max_iterations;
int function_calls;
uint32_t *polys;

mpz_t constant;

static inline void g(mpz_t output, mpz_t input, mpz_t n, mpz_t temp) {
	mpz_mul(temp, input, input);
	mpz_add(temp, temp, constant);
	mpz_tdiv_r(output, temp, n);
    function_calls++;
}

void run_rho(fact_obj_t *fobj);

int readComposites(Composite **composites, char *filename);

void init_factobj(fact_obj_t *fobj);
void free_factobj(fact_obj_t *fobj);

#endif // RHO_H
