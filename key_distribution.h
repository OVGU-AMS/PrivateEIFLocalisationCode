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

void dist_phe_key(int num_sensors, pubkey_t *pubkey);
void dist_agg_keys(int num_sensors, aggkey_t *aggkeys);

#endif