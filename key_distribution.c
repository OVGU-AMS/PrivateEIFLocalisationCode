/*
 *
 */

#include "key_distribution.h"


// Broadcasting of the Paillier public key to all sensors (starting from process index 1)
void dist_phe_key(int num_sensors, pubkey_t *pubkey, paillier_serialisation_params_t *serialisation_params){
    // Serialize key, note character buffer len +1 for null termination
    char *key = (char *)malloc(serialisation_params->paillier_max_key_serialisation_chars * sizeof(char));
    serialise_pubkey(pubkey, key);

    // Send key to all
    MPI_Bcast(key, serialisation_params->paillier_max_key_serialisation_chars, MPI_CHAR, 0, MPI_COMM_WORLD);
    free(key);
}

// Debug - Broadcasting of the Paillier private key to all sensors (starting from process index 1)
void dist_phe_prv_key(int num_sensors, prvkey_t *prvkey, paillier_serialisation_params_t *serialisation_params){
    // Serialize key, note character buffer len +1 for null termination (total overkill for private key - onley serialising lambda, not n, but meh for debugging)
    char *key = (char *)malloc(serialisation_params->paillier_max_key_serialisation_chars * sizeof(char));
    serialise_prvkey(prvkey, key);

    // Send key to all
    MPI_Bcast(key, serialisation_params->paillier_max_key_serialisation_chars, MPI_CHAR, 0, MPI_COMM_WORLD);
    free(key);
}

// Generating and distributing aggregation private keys (starting from process index 1)
void dist_agg_keys(int num_sensors, aggkey_t *aggkeys, char *aggkey_strs, MPI_Request *agg_requests, paillier_serialisation_params_t *serialisation_params){
    // Serialise and send each key to respective sensor
    for (int s=0; s<num_sensors; s++){
        serialise_aggkey(aggkeys[s], aggkey_strs+(s*(serialisation_params->paillier_max_key_serialisation_chars)));
        MPI_Isend(aggkey_strs+(s*(serialisation_params->paillier_max_key_serialisation_chars)), serialisation_params->paillier_max_key_serialisation_chars, MPI_CHAR, s+1, 0, MPI_COMM_WORLD, agg_requests+s);
    }
}
