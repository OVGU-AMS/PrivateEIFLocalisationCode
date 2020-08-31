/* 
 * 
 */

#include "encoding.h"


// n+2m - mod bits
// m - frac bits



// 8888888888
// 888
// 888
// 8888888    88888b.   .d8888b
// 888        888 "88b d88P"
// 888        888  888 888
// 888        888  888 Y88b.
// 8888888888 888  888  "Y8888P




void encode_from_dbl(mpz_t res, double x, unsigned int mults, mpz_t modulus, encoding_params_t *encoding_params){
    int sign = 0;
    unsigned int adj_frac_bits;
    mpz_t frac_factor;
    mpf_t scaled_x, frac_factor_f;

    // Inits
    mpz_init(frac_factor);
    mpf_init(scaled_x);
    mpf_init(frac_factor_f);

    // Make x positive
    if (x<0){
        sign = 1;
        x = -x;
    }

    // Compute scaling factor frac_factor for number of mults
    if (mults > 0){
        adj_frac_bits = (encoding_params->frac_bits)*(mults+1);
    } else {
        adj_frac_bits = encoding_params->frac_bits;
    }
    mpz_ui_pow_ui(frac_factor, 2, adj_frac_bits);

    // Scale x (as float and then convert to int)
    mpf_set_z(frac_factor_f, frac_factor);
    mpf_set_d(scaled_x, x);
    mpf_mul(scaled_x, frac_factor_f, scaled_x);
    mpz_set_f(res, scaled_x);

    // res should already be less than the modulus (otherwise it's an encoding overflow)
    mpz_mod(res, res, modulus);

    // Account for sign
    if (sign){
        mpz_sub(res, modulus, res);
    }

    // Frees
    mpz_clear(frac_factor);
    mpf_clear(scaled_x);
    mpf_clear(frac_factor_f);

    return;
}

// 8888888b.
// 888  "Y88b
// 888    888
// 888    888  .d88b.   .d8888b
// 888    888 d8P  Y8b d88P"
// 888    888 88888888 888
// 888  .d88P Y8b.     Y88b.
// 8888888P"   "Y8888   "Y8888P




double decode_to_dbl(mpz_t e, unsigned int mults, mpz_t modulus, encoding_params_t *encoding_params){
    int sign = 0;
    unsigned int adj_frac_bits;
    double res;
    mpz_t e_mod, modulus_hlf, frac_factor;
    mpf_t e_f, frac_factor_f;

    // Inits
    mpz_init(e_mod);
    mpz_init(modulus_hlf);
    mpz_init(frac_factor);
    mpf_init(e_f);
    mpf_init(frac_factor_f);

    // Save half modulus
    mpz_fdiv_q_ui(modulus_hlf, modulus, 2);

    // Compute scaling factor frac_factor for number of mults
    if (mults > 0){
        adj_frac_bits = (encoding_params->frac_bits)*(mults+1);
    } else {
        adj_frac_bits = encoding_params->frac_bits;
    }
    mpz_ui_pow_ui(frac_factor, 2, adj_frac_bits);

    // e should already be less than the modulus
    mpz_mod(e_mod, e, modulus);

    // Remove sign from encoding
    if (mpz_cmp(e_mod, modulus_hlf) >= 0){
        sign = 1;
        mpz_sub(e_mod, modulus, e_mod);
    }

    // Compute value by converting to float and computing division
    mpf_set_z(frac_factor_f, frac_factor);
    mpf_set_z(e_f, e_mod);
    mpf_div(e_f, e_f, frac_factor_f);
    res = mpf_get_d(e_f);
    
    // Reapply sign to result
    if (sign){
        res = -res;
    }

    // Frees
    mpz_clear(e_mod);
    mpz_clear(modulus_hlf);
    mpz_clear(frac_factor);
    mpf_clear(e_f);
    mpf_clear(frac_factor_f);

    return res;
}
