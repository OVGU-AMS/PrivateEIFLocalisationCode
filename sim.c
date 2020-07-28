/*
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "key_setup.h"
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
        
        // Key generation and distribution
        paillier_pubkey_t *pubkey;
        paillier_prvkey_t *prvkey;
        gen_phe_keys(PAILLIER_BITSIZE, &pubkey, &prvkey);
        dist_phe_key(num_sensors, pubkey);
        gen_dist_agg_keys(num_sensors, pubkey);

        // Run navigator
        run_navigator(pubkey, prvkey);


    // Remaining processes are the sensors
    } else {
        run_sensor(proc_id);
    }

    // Finish up
    mpi_err = MPI_Finalize();
    check_mpi_error(mpi_err, "Could not finalize");
    return 0;
}

// Return code error printer for MPI calls
void check_mpi_error(int mpi_err, char *message){
    if (mpi_err != 0){
        fprintf(stderr, "MPI Error! Code: %d - %s", mpi_err, message);
        exit(0);
    }
}