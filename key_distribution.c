/*
 *
 */

#include "key_distribution.h"


// Broadcasting of the Paillier public key to all sensors (starting from process index 1)
void dist_phe_key(int num_sensors, pubkey_t *pubkey){
    // Serialize key, note character buffer len +1 for null termination
    char key[MAX_KEY_SERIALISATION_CHARS];
    serialise_pubkey(pubkey, key);
    int key_len = strlen(key)+1;

    // Send key to all
    MPI_Bcast(key, key_len, MPI_CHAR, 0, MPI_COMM_WORLD);
}


// Generating and distributing aggregation private keys (starting from process index 1)
void dist_agg_keys(int num_sensors, aggkey_t *aggkeys){
    // Key variables
    char key_str[MAX_KEY_SERIALISATION_CHARS];

    // Serialise and send each key to respective sensor
    for (int s=0; s<num_sensors; s++){
        serialise_aggkey(aggkeys[s], key_str);
        MPI_Send(key_str, strlen(key_str)+1, MPI_CHAR, s+1, 0, MPI_COMM_WORLD);
    }
}
