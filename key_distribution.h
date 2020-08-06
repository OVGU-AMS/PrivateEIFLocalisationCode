/*
 *
 */

#ifndef KEY_DISTRIBUTION_H
#define KEY_DISTRIBUTION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <gmp.h>
#include "encoded_paillier_agg.h"

void dist_phe_key(int num_sensors, pubkey_t *pubkey, paillier_serialisation_params_t *serialisation_params);
void dist_agg_keys(int num_sensors, aggkey_t *aggkeys, char *aggkey_strs, MPI_Request *agg_requests, paillier_serialisation_params_t *serialisation_params);

#endif