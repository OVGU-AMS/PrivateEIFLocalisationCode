/*
 * 
 */

#include "encoded_paillier_agg.h"
#include "libpaillier-0.8/paillier.h"

// Local
void init_rand_agg(gmp_randstate_t rnd, paillier_get_rand_t get_rand, int bytes);
int PKCS1_MGF1(unsigned char *mask, long len, const unsigned char *seed, long seedlen, const EVP_MD *dgst);


// 888    d8P
// 888   d8P
// 888  d8P
// 888d88K      .d88b.  888  888  .d88b.   .d88b.  88888b.
// 8888888b    d8P  Y8b 888  888 d88P"88b d8P  Y8b 888 "88b
// 888  Y88b   88888888 888  888 888  888 88888888 888  888
// 888   Y88b  Y8b.     Y88b 888 Y88b 888 Y8b.     888  888
// 888    Y88b  "Y8888   "Y88888  "Y88888  "Y8888  888  888
//                           888      888
//                      Y8b d88P Y8b d88P
//                       "Y88P"   "Y88P"

// Use wrapped paillier library to genenrate keys
void key_gen(int paillier_bitsize, pubkey_t **pubkey, prvkey_t **prvkey){
    paillier_keygen(paillier_bitsize, pubkey, prvkey, paillier_get_rand_devrandom);
}

// Generate aggregation keys given the paillier public key parameters
void agg_key_gen(pubkey_t *pubkey, int participants, aggkey_t *key_array){
    // Init and allocate memory for last key (used during generation of others)
    mpz_init(key_array[participants-1]);

    // Required for large random number generation. Copied from paillier library
    gmp_randstate_t rnd;
    init_rand_agg(rnd, paillier_get_rand_devrandom, (pubkey->bits)*2 / 8 + 1);

    // Generate all but one key as random values less than the Paillier modulus
    for (int s=0; s<participants-1; s++){
        // Init and allocate memory for each key
        mpz_init(key_array[s]);

        // Generate random key
        do
            mpz_urandomb(key_array[s], rnd, (pubkey->bits)*2);
        while(mpz_cmp(key_array[s], pubkey->n_squared) >= 0);

        // Keep track of sum of keys
        mpz_add(key_array[participants-1], key_array[participants-1], key_array[s]);
    }

    // Generate the final key, making them all sum to 0
    mpz_neg(key_array[participants-1], key_array[participants-1]);
    mpz_mod(key_array[participants-1], key_array[participants-1], pubkey->n_squared);
}

// 888b     d888 8888888b. 8888888                                d8b          888 d8b                   888    d8b
// 8888b   d8888 888   Y88b  888                                  Y8P          888 Y8P                   888    Y8P
// 88888b.d88888 888    888  888                                               888                       888
// 888Y88888P888 888   d88P  888        .d8888b   .d88b.  888d888 888  8888b.  888 888 .d8888b   8888b.  888888 888  .d88b.  88888b.
// 888 Y888P 888 8888888P"   888        88K      d8P  Y8b 888P"   888     "88b 888 888 88K          "88b 888    888 d88""88b 888 "88b
// 888  Y8P  888 888         888        "Y8888b. 88888888 888     888 .d888888 888 888 "Y8888b. .d888888 888    888 888  888 888  888
// 888   "   888 888         888             X88 Y8b.     888     888 888  888 888 888      X88 888  888 Y88b.  888 Y88..88P 888  888
// 888       888 888       8888888       88888P'  "Y8888  888     888 "Y888888 888 888  88888P' "Y888888  "Y888 888  "Y88P"  888  888




// Serialisation of paillier public key for MPI communication
// Does not allocate buffer
void serialise_pubkey(pubkey_t *pubkey, char *buffer){
    mpz_get_str(buffer, SERIALISATION_BASE, pubkey->n);
}

// Deserialisation of paillier public key for MPI communcation
// Allocates pubkey
pubkey_t* deserialise_pubkey(char *key_serialisation){
    return paillier_pubkey_from_hex(key_serialisation);
}

// Serialisation of aggregation private key for MPI communication
// Does not allocates buffer
void serialise_aggkey(aggkey_t aggkey, char *buffer){
    mpz_get_str(buffer, SERIALISATION_BASE, aggkey);
}

// Deserialisation of aggregation private key for MPI communication
// Allocates aggkey
void deserialise_aggkey(aggkey_t aggkey, char *aggkey_serialisation){
    mpz_init_set_str(aggkey, aggkey_serialisation, SERIALISATION_BASE);
}

// Serialisation of encryptions for MPI communication
// Does not allocates buffer
void serialise_encryption(ciphertext_t *ct, char *buffer){
    mpz_get_str(buffer, SERIALISATION_BASE, ct->c);
}

// Deserialisation of encryptions for MPI communication
// Allocates ciphertext
ciphertext_t* deserialise_encryption(char *encryption_serialisation){
    ciphertext_t *ct;
    ct = paillier_create_enc_zero();
    mpz_set_str(ct->c, encryption_serialisation, SERIALISATION_BASE);
    return ct;
}

// 8888888          d8b 888
//   888            Y8P 888
//   888                888
//   888   88888b.  888 888888
//   888   888 "88b 888 888
//   888   888  888 888 888
//   888   888  888 888 Y88b.
// 8888888 888  888 888  "Y888




ciphertext_t* init_ciphertext(){
    return paillier_create_enc_zero();
}

void copy_encryption(ciphertext_t *dst, ciphertext_t *src){
    mpz_set(dst->c, src->c);
}

void encrypt_zero(pubkey_t *pubkey, ciphertext_t *ct){
    paillier_plaintext_t *p = paillier_plaintext_from_ui(0);
    paillier_enc(ct, pubkey, p, paillier_get_rand_devurandom);
    paillier_freeplaintext(p);
}

void refresh_encryption(pubkey_t *pubkey, ciphertext_t *dst, ciphertext_t *src){
    // Required for large random number generation. Copied from paillier library
    mpz_t r;
    gmp_randstate_t rnd;
    init_rand_agg(rnd, paillier_get_rand_devrandom, (pubkey->bits)*2 / 8 + 1);

    // Init new noise
    mpz_init(r);

    // Generate random key
    do
        mpz_urandomb(r, rnd, (pubkey->bits));
    while(mpz_cmp(r, pubkey->n) >= 0);

    // Compute r^N
    mpz_powm(r, r, pubkey->n, pubkey->n_squared);

    // Multiply exiting encryption mod N^2 refreshing it
    mpz_mul(dst->c, r, src->c);
    mpz_mod(dst->c, dst->c, pubkey->n_squared);
}

// 888    888                                                                         888      d8b
// 888    888                                                                         888      Y8P
// 888    888                                                                         888
// 8888888888  .d88b.  88888b.d88b.   .d88b.  88888b.d88b.   .d88b.  888d888 88888b.  88888b.  888  .d8888b       .d88b.  88888b.  .d8888b
// 888    888 d88""88b 888 "888 "88b d88""88b 888 "888 "88b d88""88b 888P"   888 "88b 888 "88b 888 d88P"         d88""88b 888 "88b 88K
// 888    888 888  888 888  888  888 888  888 888  888  888 888  888 888     888  888 888  888 888 888           888  888 888  888 "Y8888b.
// 888    888 Y88..88P 888  888  888 Y88..88P 888  888  888 Y88..88P 888     888 d88P 888  888 888 Y88b.         Y88..88P 888 d88P      X88
// 888    888  "Y88P"  888  888  888  "Y88P"  888  888  888  "Y88P"  888     88888P"  888  888 888  "Y8888P       "Y88P"  88888P"   88888P'
//                                                                           888                                          888
//                                                                           888                                          888
//                                                                           888                                          888

// Encode and encrypt a float, given public key and number of multiplications for encoding
void encode_and_enc(pubkey_t *pubkey, ciphertext_t *res, double a, unsigned int mults, encoding_params_t *encoding_params){
    paillier_plaintext_t *p = paillier_plaintext_from_ui(0);
    encode_from_dbl(p->m, a, mults, pubkey->n, encoding_params);
    paillier_enc(res, pubkey, p, paillier_get_rand_devurandom);
    paillier_freeplaintext(p);
}

// Encode and put g to power of encoding mod N^2 (encrypt without noise term)
void encode_and_enc_no_noise(pubkey_t *pubkey, ciphertext_t *res, double a, unsigned int mults, encoding_params_t *encoding_params){
    paillier_plaintext_t *p = paillier_plaintext_from_ui(0);
    encode_from_dbl(p->m, a, mults, pubkey->n, encoding_params);
    mpz_powm(res->c, pubkey->n_plusone, p->m, pubkey->n_squared);
    paillier_freeplaintext(p);
}

// Decrypt and decode encryption to a float given necessary keys and the number of multiplications for decoding
double dec_and_decode(pubkey_t *pubkey, prvkey_t *prvkey, ciphertext_t *ct, unsigned int mults, encoding_params_t *encoding_params){
    paillier_plaintext_t *p = paillier_plaintext_from_ui(0);
    paillier_dec(p, pubkey, prvkey, ct);
    return decode_to_dbl(p->m, mults, pubkey->n, encoding_params);
}

// Encode value and multipy given encryption by it
void encode_and_mult_enc(pubkey_t *pubkey, ciphertext_t *res, ciphertext_t *ct, double a, unsigned int mults, encoding_params_t *encoding_params){
    paillier_plaintext_t *p = paillier_plaintext_from_ui(0);
    encode_from_dbl(p->m, a, mults, pubkey->n, encoding_params);
    paillier_exp(pubkey, res, ct, p);
    paillier_freeplaintext(p);
}

// Add two encrypted values together
void add_encs(pubkey_t *pubkey, ciphertext_t *res, ciphertext_t *ct1, ciphertext_t *ct2){
    paillier_mul(pubkey, res, ct1, ct2);
}

//        d8888                                                     888    d8b
//       d88888                                                     888    Y8P
//      d88P888                                                     888
//     d88P 888  .d88b.   .d88b.  888d888 .d88b.   .d88b.   8888b.  888888 888  .d88b.  88888b.
//    d88P  888 d88P"88b d88P"88b 888P"  d8P  Y8b d88P"88b     "88b 888    888 d88""88b 888 "88b
//   d88P   888 888  888 888  888 888    88888888 888  888 .d888888 888    888 888  888 888  888
//  d8888888888 Y88b 888 Y88b 888 888    Y8b.     Y88b 888 888  888 Y88b.  888 Y88..88P 888  888
// d88P     888  "Y88888  "Y88888 888     "Y8888   "Y88888 "Y888888  "Y888 888  "Y88P"  888  888
//                   888      888                      888
//              Y8b d88P Y8b d88P                 Y8b d88P
//               "Y88P"   "Y88P"                   "Y88P"

// Adding aggregation noise
void add_agg_noise(pubkey_t *pubkey, aggkey_t aggkey, ciphertext_t *ct, char *stamp, int stamp_len){
    mpz_t hash;
    mpz_t noise;
    int hash_word_len = mpz_sizeinbase(pubkey->n_squared, 2)/8;
    unsigned char *hash_val_str = (unsigned char *)malloc(hash_word_len*sizeof(unsigned char));
    mpz_init(hash);
    mpz_init(noise);

    // Compute MGF1 hash and scale value to range [0, N^2)
    PKCS1_MGF1(hash_val_str, hash_word_len, (const unsigned char *) stamp, stamp_len, EVP_sha256());
    mpz_import(hash, hash_word_len, 1, sizeof(unsigned char), 0, 0, hash_val_str);
    mpz_mul(hash, hash, pubkey->n_squared);
    mpz_tdiv_q_2exp(hash, hash, mpz_sizeinbase(pubkey->n_squared, 2));

    // Add aggregation noise to the ciphertext
    mpz_powm(noise, hash, aggkey, pubkey->n_squared);
    mpz_mul(ct->c, ct->c, noise);
    mpz_mod(ct->c, ct->c, pubkey->n_squared);

    free(hash_val_str);
}

// Copied and modified (to stop EVP_MD_CTX incomplete type error) from OpenSSL library crypto/rsa/rsa_oaep.c
int PKCS1_MGF1(unsigned char *mask, long len, const unsigned char *seed, long seedlen, const EVP_MD *dgst){
    long i, outlen = 0;
    unsigned char cnt[4];
    EVP_MD_CTX *c;
    unsigned char md[EVP_MAX_MD_SIZE];
    int mdlen;
    int rv = -1;

    c = EVP_MD_CTX_create();
    EVP_MD_CTX_init(c);
    mdlen = EVP_MD_size(dgst);
    if (mdlen < 0){
        goto err;
    }
    for (i = 0; outlen < len; i++){
        cnt[0] = (unsigned char)((i >> 24) & 255);
        cnt[1] = (unsigned char)((i >> 16) & 255);
        cnt[2] = (unsigned char)((i >> 8)) & 255;
        cnt[3] = (unsigned char)(i & 255);
        if (!EVP_DigestInit_ex(c,dgst, NULL)
            || !EVP_DigestUpdate(c, seed, seedlen)
            || !EVP_DigestUpdate(c, cnt, 4)){

            goto err;
        }
        if (outlen + mdlen <= len){
            if (!EVP_DigestFinal_ex(c, mask + outlen, NULL)){
                goto err;
            }
            outlen += mdlen;
        } else {
            if (!EVP_DigestFinal_ex(c, md, NULL)){
                goto err;
            }
            memcpy(mask + outlen, md, len - outlen);
            outlen = len;
        }
    }
    rv = 0;
    err:
    EVP_MD_CTX_destroy(c);
    return rv;
}

// 8888888888
// 888
// 888
// 8888888 888d888 .d88b.   .d88b.       88888b.d88b.   .d88b.  88888b.d88b.   .d88b.  888d888 888  888
// 888     888P"  d8P  Y8b d8P  Y8b      888 "888 "88b d8P  Y8b 888 "888 "88b d88""88b 888P"   888  888
// 888     888    88888888 88888888      888  888  888 88888888 888  888  888 888  888 888     888  888
// 888     888    Y8b.     Y8b.          888  888  888 Y8b.     888  888  888 Y88..88P 888     Y88b 888
// 888     888     "Y8888   "Y8888       888  888  888  "Y8888  888  888  888  "Y88P"  888      "Y88888
//                                                                                                  888
//                                                                                             Y8b d88P
//                                                                                              "Y88P"

void free_pubkey(pubkey_t *pubkey){
    paillier_freepubkey(pubkey);
}

void free_prvkey(prvkey_t *prvkey){
    paillier_freeprvkey(prvkey);
}

void free_aggkey(aggkey_t aggkey){
    mpz_clear(aggkey);
}

void free_ciphertext(ciphertext_t *ct){
    paillier_freeciphertext(ct);
}

// 888                                888
// 888                                888
// 888                                888
// 888      .d88b.   .d8888b  8888b.  888
// 888     d88""88b d88P"        "88b 888
// 888     888  888 888      .d888888 888
// 888     Y88..88P Y88b.    888  888 888
// 88888888 "Y88P"   "Y8888P "Y888888 888




// Copied from libpaillier-X.X library for getting random mpz values for aggregation
void init_rand_agg(gmp_randstate_t rnd, paillier_get_rand_t get_rand, int bytes){
	void* buf;
	mpz_t s;

	buf = malloc(bytes);
	get_rand(buf, bytes);

	gmp_randinit_default(rnd);
	mpz_init(s);
	mpz_import(s, bytes, 1, 1, 0, 0, buf);
	gmp_randseed(rnd, s);
	mpz_clear(s);

	free(buf);
}