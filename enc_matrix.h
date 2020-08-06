/* 
 * 
 */

#ifndef ENC_MATRIX_H
#define ENC_MATRIX_H

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include "encoded_paillier_agg.h"

// ====================== // Types // ====================== //

typedef struct CiphertextMatrixTag c_mtrx_t;

// ====================== // Memory management // ====================== //

c_mtrx_t* c_mtrx_alloc(int rows, int cols);
void c_mtrx_free(c_mtrx_t *m);

// ====================== // Init // ====================== //

void init_c_mtrx(c_mtrx_t *m);

// ====================== // Getting and setting values // ====================== //
// Get and set by reference only

ciphertext_t* get_c_mtrx(c_mtrx_t *m, int row, int col);
void set_c_mtrx(c_mtrx_t *m, int row, int col, ciphertext_t *in);

// ====================== // Encryption and decryption // ====================== //
// No compatible size checks between in and out matrices is made

void encrypt_mtrx(pubkey_t *pubkey, gsl_matrix *plain_mat, c_mtrx_t *enc_mat, unsigned int mults, encoding_params_t *encoding_params);
void decrypt_mtrx(pubkey_t *pubkey, prvkey_t *prvkey, c_mtrx_t *enc_mat, gsl_matrix *plain_mat, unsigned int mults, encoding_params_t *encoding_params);
void decrypt_vctr(pubkey_t *pubkey, prvkey_t *prvkey, c_mtrx_t *enc_mat, gsl_vector *plain_vec, unsigned int mults, encoding_params_t *encoding_params);

// ====================== // Homomrphic actions // ====================== //

void add_c_mtrx_c_mtrx(pubkey_t *pubkey, c_mtrx_t *enc_mat1, c_mtrx_t *enc_mat2, c_mtrx_t *out);

// ====================== // Aggregation // ====================== //

void add_agg_noise_c_mtrx(pubkey_t *pubkey, aggkey_t aggkey, c_mtrx_t *m, int timestamp, int identifier);


#endif