/*
 * 
 */

#ifndef ENC_NAVIGATOR_H
#define ENC_NAVIGATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <gmp.h>
#include <gsl/gsl_linalg.h>

#include "key_setup.h"
#include "encoding.h"

void run_navigator(paillier_pubkey_t *pubkey, paillier_prvkey_t *prvkey);

#endif