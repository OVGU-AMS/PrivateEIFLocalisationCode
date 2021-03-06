/*
 * 
 */

#ifndef ENC_NAVIGATOR_H
#define ENC_NAVIGATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gmp.h>
#include <gsl/gsl_linalg.h>
#include "encoded_paillier_agg.h"
#include "enc_matrix.h"

void run_navigator(pubkey_t *pubkey, prvkey_t *prvkey, int number_sensors, char *track_filename, char *output_filepath, 
                    encoding_params_t *encoding_params, paillier_serialisation_params_t *serialisation_params);

#endif