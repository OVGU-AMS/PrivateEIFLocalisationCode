/*
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "encoded_paillier_agg.h"
#include "key_distribution.h"
#include "sensor.h"
#include "navigator.h"


// Local functions
void check_mpi_error(int mpi_err, char *message);



int main(int argc, char *argv[]){
    int mpi_err;
    int num_procs;
    int num_sensors;
    int proc_id;

    // Init multi processing
    mpi_err = MPI_Init(&argc, &argv);
    check_mpi_error(mpi_err, "Could not initialise");

    // Get current process id and total number of processes
    mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    check_mpi_error(mpi_err, "Could not check process id");
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    check_mpi_error(mpi_err, "Could not check number of processes");
    num_sensors = num_procs-1;

    // Process 0 acts as key distributer, and then the navigator during the simulation
    if (proc_id == 0){
        pubkey_t *pubkey;
        prvkey_t *prvkey;
        aggkey_t *aggkeys = (aggkey_t*)malloc(num_sensors*sizeof(aggkey_t));

        // Generate Paillier keys
        key_gen(PAILLIER_BITSIZE, &pubkey, &prvkey);

        // Generate aggregation keys
        agg_key_gen(pubkey, num_sensors, aggkeys);

        // Broadcast Paillier key
        dist_phe_key(num_sensors, pubkey);

        // Distributed aggregation keys
        dist_agg_keys(num_sensors, aggkeys);

        // Free aggregation keys
        for (int s=0; s<num_sensors; s++){
            free_aggkey(aggkeys[s]);
        }
        free(aggkeys);
        aggkeys = NULL;

        // Run navigator
        run_navigator(pubkey, prvkey, num_sensors);
        fprintf(stderr, "Navigator finished.\n");

        // Free Paillier keys
        free_pubkey(pubkey);
        free_prvkey(prvkey);

    // Remaining processes are the sensors
    } else {
        run_sensor(proc_id);
        fprintf(stderr, "Sensor %d finished.\n", proc_id);
    }

    // Finish up
    mpi_err = MPI_Finalize();
    check_mpi_error(mpi_err, "Could not finalize");
    return 0;
}

// Return code error printer for MPI calls
void check_mpi_error(int mpi_err, char *message){
    if (mpi_err != 0){
        fprintf(stderr, "MPI Error! Code: %d - %s\n", mpi_err, message);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}