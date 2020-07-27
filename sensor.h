/*
 * 
 */

#ifndef ENC_SENSOR_H
#define ENC_SENSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <gmp.h>
#include "key_setup.h"

#define SENSOR_VARIANCE 5.0

void run_sensor(int id);

#endif