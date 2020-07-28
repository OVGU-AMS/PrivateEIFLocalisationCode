/*
 *
 */

#ifndef KEY_SETUP_H
#define KEY_SETUP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <gmp.h>

#include "libpaillier-0.8/paillier.h" // Should be the only place this is imported (library doesn't check if header imported more than once)

#define PAILLIER_BITSIZE 1024
#define SERIALISATION_BASE 16
#define MAX_KEY_SERIALISATION_CHARS 1024
#define MAX_ENC_SERIALISATION_CHARS 1024

void gen_phe_keys(int paillier_bitsize, paillier_pubkey_t **pubkey, paillier_prvkey_t **prvkey);
void dist_phe_key(int num_sensors, paillier_pubkey_t *pubkey);
void gen_dist_agg_keys(int num_sensors, paillier_pubkey_t *pubkey);

#endif