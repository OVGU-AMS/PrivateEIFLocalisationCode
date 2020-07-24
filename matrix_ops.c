/* 
 * 
 * 
 */


#include <stdlib.h>
#include <stdio.h>
#include "matrix_ops.h"


// Local functions
unsigned long int ipow(unsigned long int base, unsigned long int exp);


// 88888888888
//     888
//     888
//     888  888  888 88888b.   .d88b.  .d8888b
//     888  888  888 888 "88b d8P  Y8b 88K
//     888  888  888 888  888 88888888 "Y8888b.
//     888  Y88b 888 888 d88P Y8b.          X88
//     888   "Y88888 88888P"   "Y8888   88888P'
//               888 888
//          Y8b d88P 888
//           "Y88P"  888

struct PlaintextMatrixTag {
    int size1;
    int size2;
    paillier_plaintext_t*** data;
};

struct CiphertextMatrixTag {
    int size1;
    int size2;
    paillier_ciphertext_t*** data;
};

// 888b     d888
// 8888b   d8888
// 88888b.d88888
// 888Y88888P888  .d88b.  88888b.d88b.   .d88b.  888d888 888  888
// 888 Y888P 888 d8P  Y8b 888 "888 "88b d88""88b 888P"   888  888
// 888  Y8P  888 88888888 888  888  888 888  888 888     888  888
// 888   "   888 Y8b.     888  888  888 Y88..88P 888     Y88b 888
// 888       888  "Y8888  888  888  888  "Y88P"  888      "Y88888
//                                                            888
//                                                       Y8b d88P
//                                                        "Y88P"

p_mtrx_t* p_mtrx_alloc(int rows, int cols){
    p_mtrx_t *m;
    // Space for matrix
    m = (p_mtrx_t*)malloc(sizeof(p_mtrx_t));
    m->size1 = rows;
    m->size2 = cols;
    // Space for rows
    m->data = (paillier_plaintext_t***)malloc(rows*sizeof(paillier_plaintext_t**));
    // Space for columns
    for(int i=0; i<rows; i++){
        m->data[i] = (paillier_plaintext_t**)malloc(cols*sizeof(paillier_plaintext_t*));
        // Init all values to null
        for(int j=0; j<cols; j++){
            m->data[i][j] = NULL;
        }
    }
    return m;
}

c_mtrx_t* c_mtrx_alloc(int rows, int cols){
    c_mtrx_t *m;
    // Space for matrix
    m = (c_mtrx_t*)malloc(sizeof(c_mtrx_t));
    m->size1 = rows;
    m->size2 = cols;
    // Space for rows
    m->data = (paillier_ciphertext_t***)malloc(rows*sizeof(paillier_ciphertext_t**));
    // Space for columns
    for(int i=0; i<rows; i++){
        m->data[i] = (paillier_ciphertext_t**)malloc(cols*sizeof(paillier_ciphertext_t*));
        // Init all values to null
        for(int j=0; j<cols; j++){
            m->data[i][j] = NULL;
        }
    }
    return m;
}

void p_mtrx_free(p_mtrx_t* m){
    for(int i=0; i<m->size1; i++){
        for(int j=0; j<m->size2; j++){
            // Free plaintext objects if they are present
            if (m->data[i][j] != NULL){
                paillier_freeplaintext(m->data[i][j]);
                m->data[i][j] = NULL;
            }
        }
        // Free columns
        free(m->data[i]);
        m->data[i] = NULL;
    }
    // Free rows
    free(m->data);
    m->data = NULL;
    // Free matrix
    free(m);
    m = NULL;
    return;
}

void c_mtrx_free(c_mtrx_t* m){
    for(int i=0; i<m->size1; i++){
        for(int j=0; j<m->size2; j++){
            // Free plaintext objects if they are present
            if (m->data[i][j] != NULL){
                paillier_freeciphertext(m->data[i][j]);
                m->data[i][j] = NULL;
            }
        }
        // Free columns
        free(m->data[i]);
        m->data[i] = NULL;
    }
    // Free rows
    free(m->data);
    m->data = NULL;
    // Free matrix
    free(m);
    m = NULL;
    return;
}

//  .d8888b.           888           d88P  .d8888b.           888
// d88P  Y88b          888          d88P  d88P  Y88b          888
// 888    888          888         d88P   Y88b.               888
// 888         .d88b.  888888     d88P     "Y888b.    .d88b.  888888
// 888  88888 d8P  Y8b 888       d88P         "Y88b. d8P  Y8b 888
// 888    888 88888888 888      d88P            "888 88888888 888
// Y88b  d88P Y8b.     Y88b.   d88P       Y88b  d88P Y8b.     Y88b.
//  "Y8888P88  "Y8888   "Y888 d88P         "Y8888P"   "Y8888   "Y888




void get_p_mtrx(p_mtrx_t *m, int row, int col, mpz_t out){
    mpz_set(out, (m->data)[row][col]->m);
    return;
}

void get_c_mtrx(c_mtrx_t *m, int row, int col, mpz_t out){
    mpz_set(out, (m->data)[row][col]->c);
    return;
}

void set_p_mtrx(p_mtrx_t *m, int row, int col, mpz_t x){
    if ((m->data)[row][col] == NULL){
        (m->data)[row][col] = paillier_plaintext_from_ui(0);
    }
    mpz_set((m->data)[row][col]->m, x);
    return;
}

// 8888888888                          d88P 8888888b.
// 888                                d88P  888  "Y88b
// 888                               d88P   888    888
// 8888888    88888b.   .d8888b     d88P    888    888  .d88b.   .d8888b
// 888        888 "88b d88P"       d88P     888    888 d8P  Y8b d88P"
// 888        888  888 888        d88P      888    888 88888888 888
// 888        888  888 Y88b.     d88P       888  .d88P Y8b.     Y88b.
// 8888888888 888  888  "Y8888P d88P        8888888P"   "Y8888   "Y8888P




void encrypt_mtrx(paillier_pubkey_t* pubkey, p_mtrx_t* pm, c_mtrx_t* cm_out){
    for (int i=0; i<cm_out->size1;i++){
        for (int j=0; j<cm_out->size2;j++){
            cm_out->data[i][j] = paillier_enc(NULL, pubkey, pm->data[i][j], paillier_get_rand_devurandom);
        }
    }
    return;
}

void decrypt_mtrx(paillier_pubkey_t* pubkey, paillier_prvkey_t* prvkey, c_mtrx_t* cm, p_mtrx_t* pm_out){
    for (int i=0; i<pm_out->size1;i++){
        for (int j=0; j<pm_out->size2;j++){
            pm_out->data[i][j] = paillier_dec(NULL, pubkey, prvkey, cm->data[i][j]);
        }
    }
    return;
}

//  .d88888b.
// d88P" "Y88b
// 888     888
// 888     888 88888b.  .d8888b
// 888     888 888 "88b 88K
// 888     888 888  888 "Y8888b.
// Y88b. .d88P 888 d88P      X88
//  "Y88888P"  88888P"   88888P'
//             888
//             888
//             888

void add_c_mtrx_c_mtrx(paillier_pubkey_t* pubkey, c_mtrx_t* m1, c_mtrx_t* m2, c_mtrx_t* cm_out){
    for (int i=0; i<cm_out->size1;i++){
        for (int j=0; j<cm_out->size2;j++){
            if (cm_out->data[i][j] == NULL){
                cm_out->data[i][j] = paillier_create_enc_zero();
            }
            paillier_mul(pubkey, cm_out->data[i][j], m1->data[i][j], m2->data[i][j]);
        }
    }
    return;
}

void mult_p_sclr_c_mtrx(paillier_pubkey_t* pubkey, mpz_t s, c_mtrx_t* m, c_mtrx_t* cm_out){
    paillier_plaintext_t* sclr;
    sclr = paillier_plaintext_from_ui(0);
    mpz_set(sclr->m, s);

    for (int i=0; i<cm_out->size1;i++){
        for (int j=0; j<cm_out->size2;j++){
            if (cm_out->data[i][j] == NULL){
                cm_out->data[i][j] = paillier_create_enc_zero();
            }
            paillier_exp(pubkey, cm_out->data[i][j], m->data[i][j], sclr);
        }
    }
    paillier_freeplaintext(sclr);
    return;
}

// 888b     d888 d8b
// 8888b   d8888 Y8P
// 88888b.d88888
// 888Y88888P888 888 .d8888b   .d8888b
// 888 Y888P 888 888 88K      d88P"
// 888  Y8P  888 888 "Y8888b. 888
// 888   "   888 888      X88 Y88b.
// 888       888 888  88888P'  "Y8888P




void print_p_mtrx(p_mtrx_t* m, double (*decoder)(mpz_t, int, unsigned int, unsigned int), int decode_mult, unsigned int mod_bits, unsigned int frac_bits){
    fprintf(stderr, "Matrix:\n");
    for(int i=0; i<m->size1; i++){
        fprintf(stderr, "[ ");
        for(int j=0; j<m->size2; j++){
            if (decoder == NULL){
                gmp_fprintf(stderr, "%Zd ", m->data[i][j]->m);   
            } else {
                fprintf(stderr, "%lf ", decoder(m->data[i][j]->m, decode_mult, mod_bits, frac_bits));
            }
        }
        fprintf(stderr, "]\n");
    }
}
