/* 
 * 
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




void encode_from_dbl(mpz_t res, double x, unsigned int mod_bits, unsigned int frac_bits){
    int sign = 0;
    mpz_t mx_enc, frac_factor;
    mpf_t scaled_x, mx_enc_f, frac_factor_f;

    mpz_init(mx_enc);
    mpz_init(frac_factor);

    mpf_init(scaled_x);
    mpf_init(mx_enc_f);
    mpf_init(frac_factor_f);

    if (x<0){
        sign = 1;
        x = -x;
    }

    mpz_ui_pow_ui(mx_enc, 2, mod_bits);
    mpz_ui_pow_ui(frac_factor, 2, frac_bits);

    
    mpf_set_z(frac_factor_f, frac_factor);
    mpf_set_d(scaled_x, x);
    mpf_mul(scaled_x, frac_factor_f, scaled_x);

    mpz_set_f(res, scaled_x);

    if (sign){
        mpz_sub(res, mx_enc, res);
    }

    mpz_clear(mx_enc);
    mpz_clear(frac_factor);

    mpf_clear(scaled_x);
    mpf_clear(mx_enc_f);
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




double decode_to_dbl(mpz_t e, int mult, unsigned int mod_bits, unsigned int frac_bits){
    int sign = 0;
    double res;
    mpz_t e_copy, mx_enc, mx_enc_hlf, frac_factor;
    mpf_t e_f, frac_factor_f;

    mpz_init(e_copy);
    mpz_init(mx_enc);
    mpz_init(mx_enc_hlf);
    mpz_init(frac_factor);

    mpz_ui_pow_ui(mx_enc, 2, mod_bits);
    mpz_ui_pow_ui(mx_enc_hlf, 2, mod_bits-1);
    if (mult){
        mpz_ui_pow_ui(frac_factor, 2, 2*frac_bits);
    } else {
        mpz_ui_pow_ui(frac_factor, 2, frac_bits);
    }
    
    
    mpf_init(e_f);
    mpf_init(frac_factor_f);

    mpz_powm_ui(e_copy, e, 1, mx_enc);

    

    if (mpz_cmp(e_copy, mx_enc_hlf) >= 0){
        sign = 1;
        mpz_sub(e_copy, mx_enc, e_copy);
    }

    mpf_set_z(frac_factor_f, frac_factor);
    mpf_set_z(e_f, e_copy);
    mpf_div(e_f, e_f, frac_factor_f);

    res = mpf_get_d(e_f);
    
    if (sign){
        res = -res;
    }

    mpz_clear(e_copy);
    mpz_clear(mx_enc);
    mpz_clear(mx_enc_hlf);
    mpz_clear(frac_factor);

    mpf_clear(e_f);
    mpf_clear(frac_factor_f);

    return res;
}

// 8888888888                         .d88888b.  8888888b.  8888888888
// 888                               d88P" "Y88b 888   Y88b 888
// 888                               888     888 888    888 888
// 8888888    88888b.   .d8888b      888     888 888   d88P 8888888
// 888        888 "88b d88P"         888     888 8888888P"  888
// 888        888  888 888           888     888 888 T88b   888
// 888        888  888 Y88b.         Y88b. .d88P 888  T88b  888
// 8888888888 888  888  "Y8888P       "Y88888P"  888   T88b 8888888888




unsigned long int encode_dbl_cmp(double x, unsigned int mod_bits, unsigned int frac_bits){
    int sign = 0;
    unsigned long int enc;
    mpz_t res, mx_enc_hlf, frac_factor;
    mpf_t scaled_x, frac_factor_f;

    mpz_init(res);
    mpz_init(mx_enc_hlf);
    mpz_init(frac_factor);

    mpf_init(scaled_x);
    mpf_init(frac_factor_f);

    if (x<0){
        sign = 1;
        x = -x;
    }
    mpz_ui_pow_ui(mx_enc_hlf, 2, mod_bits-1);
    mpz_ui_pow_ui(frac_factor, 2, frac_bits);

    
    mpf_set_z(frac_factor_f, frac_factor);
    mpf_set_d(scaled_x, x);
    mpf_mul(scaled_x, frac_factor_f, scaled_x);

    mpz_set_f(res, scaled_x);

    if (!sign){
        mpz_add(res, mx_enc_hlf, res);
    }

    enc = mpz_get_ui(res);

    mpz_clear(res);
    mpz_clear(mx_enc_hlf);
    mpz_clear(frac_factor);

    mpf_clear(scaled_x);
    mpf_clear(frac_factor_f);

    return enc;
}
