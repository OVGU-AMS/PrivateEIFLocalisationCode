/* 
 * 
 * 
 */

#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include <gmp.h>

#include "key_distribution.h"

// Types

typedef struct PlaintextMatrixTag p_mtrx_t;
typedef struct CiphertextMatrixTag c_mtrx_t;

// Memory management

p_mtrx_t* p_mtrx_alloc(int rows, int cols);
c_mtrx_t* c_mtrx_alloc(int rows, int cols);

void p_mtrx_free(p_mtrx_t* m);
void c_mtrx_free(c_mtrx_t* m);

// Getting and setting values

void get_p_mtrx(p_mtrx_t *m, int row, int col, mpz_t out);
void get_c_mtrx(c_mtrx_t *m, int row, int col, mpz_t out);

void set_p_mtrx(p_mtrx_t *m, int row, int col, mpz_t x);

// Encryption and decryption

void encrypt_mtrx(paillier_pubkey_t* pubkey, p_mtrx_t* pm, c_mtrx_t* cm_out);
void decrypt_mtrx(paillier_pubkey_t* pubkey, paillier_prvkey_t* prvkey, c_mtrx_t* cm, p_mtrx_t* pm_out);

// Homomorphic actions

void add_c_mtrx_c_mtrx(paillier_pubkey_t* pubkey, c_mtrx_t* m1, c_mtrx_t* m2, c_mtrx_t* cm_out);
void mult_p_sclr_c_mtrx(paillier_pubkey_t* pubkey, mpz_t s, c_mtrx_t* m, c_mtrx_t* cm_out);

// Misc

void print_p_mtrx(p_mtrx_t* m, double (*decoder)(mpz_t, int, unsigned int, unsigned int), int decode_mult, unsigned int mod_bits, unsigned int frac_bits);


#endif