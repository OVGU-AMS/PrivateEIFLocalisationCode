/*
 * 
 */

#ifndef ENC_ENCODED_PAILLIER_AGG_H
#define ENC_ENCODED_PAILLIER_AGG_H

#include <stdlib.h>
#include <gmp.h>
#include "encoding.h"

#define PAILLIER_BITSIZE 1024
#define SERIALISATION_BASE 16
#define MAX_KEY_SERIALISATION_CHARS 1024
#define MAX_ENC_SERIALISATION_CHARS 1024

typedef struct pubkeyTag pubkey_t;
typedef struct prvkeyTag prvkey_t;
typedef struct ciphertextTag ciphertext_t;
typedef mpz_t aggkey_t;

// ====================== // Keygen // ====================== //

void key_gen(int paillier_bitsize, pubkey_t **pubkey, prvkey_t **prvkey);
void agg_key_gen(pubkey_t *pubkey, int participants, aggkey_t *key_array);

// ====================== // MPI serialisation // ====================== //
// Note serialisation does no allocation
// Note deserialisation does allocation

void serialise_pubkey(pubkey_t *pubkey, char *buffer);
pubkey_t* deserialise_pubkey(char *key_serialisation);
void serialise_aggkey(aggkey_t aggkey, char *buffer);
void deserialise_aggkey(aggkey_t aggkey, char *aggkey_serialisation);
void serialise_encryption(ciphertext_t *ct, char *buffer);
ciphertext_t* deserialise_encryption(char *encryption_serialisation);

// ====================== // Init // ====================== //

ciphertext_t* init_ciphertext();

// ====================== // Homomorphic operations // ====================== //

void encode_and_enc(pubkey_t *pubkey, ciphertext_t *res, double a, unsigned int mults);
double dec_and_decode(pubkey_t *pubkey, prvkey_t *prvkey, ciphertext_t *ct, unsigned int mults);
void encode_and_add_enc(pubkey_t *pubkey, ciphertext_t *res, ciphertext_t *ct, double a, unsigned int mults);
void encode_and_mult_enc(pubkey_t *pubkey, ciphertext_t *res, ciphertext_t *ct, double a, unsigned int mults);
void add_encs(pubkey_t *pubkey, ciphertext_t *res, ciphertext_t *ct1, ciphertext_t *ct2);

// ====================== // Aggregation // ====================== //

void add_agg_noise(pubkey_t *pubkey, ciphertext_t *res, ciphertext_t *ct, int stamp);

// ====================== // Free memory // ====================== //

void free_pubkey(pubkey_t *pubkey);
void free_prvkey(prvkey_t *prvkey);
void free_aggkey(aggkey_t aggkey);
void free_ciphertext(ciphertext_t *ct);

// ====================== // Misc // ====================== //


#endif