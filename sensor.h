/*
 * 
 */

#ifndef ENC_SENSOR_H
#define ENC_SENSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gmp.h>
#include "encoded_paillier_agg.h"
#include "enc_matrix.h"

#define SENSOR_VARIANCE 5.0

void run_sensor(int id, char *sensor_filepath_base, encoding_params_t *encoding_params, paillier_serialisation_params_t *serialisation_params);

#endif