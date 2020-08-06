/* 
 * 
 * 
 */

#ifndef ENCODING_H
#define ENCODING_H

#include <gmp.h>

// Defaults for encoding

#define ENCODING_MOD_BITS_DEFAULT 128
#define ENCODING_FRAC_BITS_DEFAULT 32

typedef struct EncodingParamsTag {
    unsigned int mod_bits;
    unsigned int frac_bits;
} encoding_params_t;

// Encoding for Paillier with single multiplication (assumes mpz already initialised)

void encode_from_dbl(mpz_t res, double x, unsigned int mults, encoding_params_t *encoding_params);
double decode_to_dbl(mpz_t e, unsigned int mults, encoding_params_t *encoding_params);

#endif