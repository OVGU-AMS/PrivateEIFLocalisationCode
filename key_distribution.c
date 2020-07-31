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
void dist_agg_keys(int num_sensors, aggkey_t *aggkeys, char *aggkey_strs, MPI_Request *agg_requests){
    // Serialise and send each key to respective sensor
    for (int s=0; s<num_sensors; s++){
        serialise_aggkey(aggkeys[s], aggkey_strs+(s*MAX_KEY_SERIALISATION_CHARS));
        MPI_Isend(aggkey_strs+(s*MAX_KEY_SERIALISATION_CHARS), MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, s+1, 0, MPI_COMM_WORLD, agg_requests+s);
    }
}
