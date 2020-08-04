/* 
 * 
 * 
 */

#ifndef ENCODING_H
#define ENCODING_H

#include <gmp.h>

// Defaults for encoding

#define MOD_BITS 256
#define FRAC_BITS 64

// Encoding for Paillier with single multiplication (assumes mpz already initialised)

void encode_from_dbl(mpz_t res, double x, unsigned int mults, unsigned int mod_bits, unsigned int frac_bits);
double decode_to_dbl(mpz_t e, unsigned int mults, unsigned int mod_bits, unsigned int frac_bits);

#endif