/* 
 * 
 */


#include "enc_matrix.h"



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

struct CiphertextMatrixTag {
    int size1;
    int size2;
    ciphertext_t*** data;
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

c_mtrx_t* c_mtrx_alloc(int rows, int cols){
    c_mtrx_t *m;
    // Space for matrix
    m = (c_mtrx_t*)malloc(sizeof(c_mtrx_t));
    m->size1 = rows;
    m->size2 = cols;
    // Space for rows
    m->data = (ciphertext_t***)malloc(rows*sizeof(ciphertext_t**));
    // Space for columns
    for(int i=0; i<rows; i++){
        m->data[i] = (ciphertext_t**)malloc(cols*sizeof(ciphertext_t*));
        // Init all values to null
        for(int j=0; j<cols; j++){
            m->data[i][j] = NULL;
        }
    }
    return m;
}

void c_mtrx_free(c_mtrx_t *m){
    for(int i=0; i<m->size1; i++){
        for(int j=0; j<m->size2; j++){
            // Free plaintext objects if they are present
            if (m->data[i][j] != NULL){
                free_ciphertext(m->data[i][j]);
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

// 8888888          d8b 888
//   888            Y8P 888
//   888                888
//   888   88888b.  888 888888
//   888   888 "88b 888 888
//   888   888  888 888 888
//   888   888  888 888 Y88b.
// 8888888 888  888 888  "Y888




void init_c_mtrx(c_mtrx_t *m){
    for(int i=0; i<m->size1; i++){
        for(int j=0; j<m->size2; j++){
            if (m->data[i][j] == NULL){
                m->data[i][j] = init_ciphertext();
            }
        }
    }
}

//  .d8888b.           888           d88P  .d8888b.           888
// d88P  Y88b          888          d88P  d88P  Y88b          888
// 888    888          888         d88P   Y88b.               888
// 888         .d88b.  888888     d88P     "Y888b.    .d88b.  888888
// 888  88888 d8P  Y8b 888       d88P         "Y88b. d8P  Y8b 888
// 888    888 88888888 888      d88P            "888 88888888 888
// Y88b  d88P Y8b.     Y88b.   d88P       Y88b  d88P Y8b.     Y88b.
//  "Y8888P88  "Y8888   "Y888 d88P         "Y8888P"   "Y8888   "Y888




ciphertext_t* get_c_mtrx(c_mtrx_t *m, int row, int col){
    return (m->data)[row][col];
}

void set_c_mtrx(c_mtrx_t *m, int row, int col, ciphertext_t *in){
    if ((m->data)[row][col] != NULL){
        free_ciphertext((m->data)[row][col]);
    }
    (m->data)[row][col] = in;
}

// 8888888888                          d88P 8888888b.
// 888                                d88P  888  "Y88b
// 888                               d88P   888    888
// 8888888    88888b.   .d8888b     d88P    888    888  .d88b.   .d8888b
// 888        888 "88b d88P"       d88P     888    888 d8P  Y8b d88P"
// 888        888  888 888        d88P      888    888 88888888 888
// 888        888  888 Y88b.     d88P       888  .d88P Y8b.     Y88b.
// 8888888888 888  888  "Y8888P d88P        8888888P"   "Y8888   "Y8888P




void encrypt_mtrx(pubkey_t *pubkey, gsl_matrix *plain_mat, c_mtrx_t *enc_mat, unsigned int mults, encoding_params_t *encoding_params){
    for (int i=0; i<enc_mat->size1;i++){
        for (int j=0; j<enc_mat->size2;j++){
            if (enc_mat->data[i][j] == NULL){
                enc_mat->data[i][j] = init_ciphertext();
            }
            encode_and_enc(pubkey, enc_mat->data[i][j], gsl_matrix_get(plain_mat, i, j), mults, encoding_params);
        }
    }
    return;
}

void decrypt_mtrx(pubkey_t *pubkey, prvkey_t *prvkey, c_mtrx_t *enc_mat, gsl_matrix *plain_mat, unsigned int mults, encoding_params_t *encoding_params){
    for (int i=0; i<enc_mat->size1;i++){
        for (int j=0; j<enc_mat->size2;j++){
            gsl_matrix_set(plain_mat, i, j, dec_and_decode(pubkey, prvkey, enc_mat->data[i][j], mults, encoding_params));
        }
    }
    return;
}

void decrypt_vctr(pubkey_t *pubkey, prvkey_t *prvkey, c_mtrx_t *enc_mat, gsl_vector *plain_vec, unsigned int mults, encoding_params_t *encoding_params){
    for (int j=0; j<enc_mat->size2;j++){
        gsl_vector_set(plain_vec, j, dec_and_decode(pubkey, prvkey, enc_mat->data[0][j], mults, encoding_params));
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

void add_c_mtrx_c_mtrx(pubkey_t *pubkey, c_mtrx_t *enc_mat1, c_mtrx_t *enc_mat2, c_mtrx_t *out){
    for (int i=0; i<out->size1;i++){
        for (int j=0; j<out->size2;j++){
            if (out->data[i][j] == NULL){
                out->data[i][j] = init_ciphertext();
            }
            add_encs(pubkey, out->data[i][j], enc_mat1->data[i][j], enc_mat2->data[i][j]);
        }
    }
    return;
}

//        d8888                            888b    888          d8b
//       d88888                            8888b   888          Y8P
//      d88P888                            88888b  888
//     d88P 888  .d88b.   .d88b.           888Y88b 888  .d88b.  888 .d8888b   .d88b.
//    d88P  888 d88P"88b d88P"88b          888 Y88b888 d88""88b 888 88K      d8P  Y8b
//   d88P   888 888  888 888  888          888  Y88888 888  888 888 "Y8888b. 88888888
//  d8888888888 Y88b 888 Y88b 888 d8b      888   Y8888 Y88..88P 888      X88 Y8b.
// d88P     888  "Y88888  "Y88888 Y8P      888    Y888  "Y88P"  888  88888P'  "Y8888
//                   888      888
//              Y8b d88P Y8b d88P
//               "Y88P"   "Y88P"

void add_agg_noise_c_mtrx(pubkey_t *pubkey, aggkey_t aggkey, c_mtrx_t *m, int timestamp, int identifier){
    char stamp[4];
    stamp[0] = (char) timestamp;
    stamp[1] = (char) identifier;

    for (int i=0; i<m->size1; i++){
        for (int j=0; j<m->size2; j++){
            stamp[2] = (char) i;
            stamp[3] = (char) j;
            add_agg_noise(pubkey, aggkey, m->data[i][j], stamp, 4);
        }
    }
}