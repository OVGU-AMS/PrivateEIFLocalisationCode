/*
 * 
 */

#include "sensor.h"

void run_sensor(int id){
    // Key vars
    char key_str[MAX_KEY_SERIALISATION_CHARS];
    paillier_pubkey_t *pubkey;
    mpz_t agg_key;

    // Get Paillier public key
    MPI_Bcast(key_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, MPI_COMM_WORLD);
    pubkey = paillier_pubkey_from_hex(key_str);

    //printf("id: %d, key: %s\n", id, paillier_pubkey_to_hex(pubkey));

    // Get sensor's private aggregation key
    MPI_Recv(key_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    mpz_init_set_str(agg_key, key_str, SERIALISATION_BASE);

    //printf("id: %d, agg key: %s\n", id, key_str);

    // Open sensor measurements file

    // Start sensor sim
    // Will need to barrier after each send maybe?

}