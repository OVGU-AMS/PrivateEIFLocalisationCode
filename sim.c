/*
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include "encoded_paillier_agg.h"
#include "key_distribution.h"
#include "sensor.h"
#include "navigator.h"

// Defaults when no input is given
#define TRACK_FILEPATH "input/debug_track1.txt"
#define NAV_FILEPATH "output/debug_nav_001.txt"
#define SENSOR_FILEPATH_BASE "input/debug_sim_001_sensor%d.txt"


// Local functions
void check_mpi_error(int mpi_err, char *message);


int main(int argc, char *argv[]){
    int mpi_err;
    int num_procs;
    int num_sensors;
    int proc_id;
    char *track_filepath = NULL;
    char *output_filepath = NULL;
    char *sensor_filepath_base = NULL;
    encoding_params_t encoding_params;
    paillier_serialisation_params_t serialisation_params;
    clock_t start_time, end_time;

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
        aggkey_t *aggkeys;
        char *aggkey_strs;
        MPI_Request *agg_requests;

        // If no commandline argument initialise paillier and encoding with defaults
        if (argc == 1){
            track_filepath = TRACK_FILEPATH;
            output_filepath = NAV_FILEPATH;
            encoding_params.frac_bits = ENCODING_FRAC_BITS_DEFAULT;
            serialisation_params.paillier_bitsize = PAILLIER_BITSIZE_DEFAULT;
            serialisation_params.paillier_max_key_serialisation_chars = PAILLIER_MAX_KEY_SERIALISATION_CHARS_DEFAULT;
            serialisation_params.paillier_max_enc_serialisation_chars = PAILLIER_MAX_ENC_SERIALISATION_CHARS_DEFAULT;
        
        // If arguments provided, use them instead
        } else if (argc == 6){
            serialisation_params.paillier_bitsize = atoi(argv[1]);
            encoding_params.frac_bits = atoi(argv[2]);
            output_filepath = argv[3];
            track_filepath = argv[4];

            serialisation_params.paillier_max_key_serialisation_chars = (serialisation_params.paillier_bitsize/2) + 1;
            serialisation_params.paillier_max_enc_serialisation_chars = (serialisation_params.paillier_bitsize/2) + 1;
        
        // Either all or no arguments must be provided, error otherwise
        } else {
            fprintf(stderr, "Incorrect commandline usage!");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Alloc key and key distribution vars
        aggkeys = (aggkey_t*)malloc(num_sensors*sizeof(aggkey_t));
        aggkey_strs = (char *)malloc(num_sensors*(serialisation_params.paillier_max_key_serialisation_chars)*sizeof(char));
        agg_requests = (MPI_Request *)malloc(num_sensors*sizeof(MPI_Request));

        // Generate Paillier keys
        key_gen(serialisation_params.paillier_bitsize, &pubkey, &prvkey);

        // Generate aggregation keys
        agg_key_gen(pubkey, num_sensors, aggkeys);

        // Broadcast Paillier key
        dist_phe_key(num_sensors, pubkey, &serialisation_params);

        // Distributed aggregation keys
        dist_agg_keys(num_sensors, aggkeys, aggkey_strs, agg_requests, &serialisation_params);

        // Free aggregation keys
        for (int s=0; s<num_sensors; s++){
            free_aggkey(aggkeys[s]);
        }
        free(aggkeys);

        // Free aggregation key string buffers
        free(aggkey_strs);

        // 8888888b.                                                    d8b                   888
        // 888   Y88b                                                   Y8P                   888
        // 888    888                                                                         888
        // 888   d88P 888  888 88888b.       88888b.   8888b.  888  888 888  .d88b.   8888b.  888888 .d88b.  888d888
        // 8888888P"  888  888 888 "88b      888 "88b     "88b 888  888 888 d88P"88b     "88b 888   d88""88b 888P"
        // 888 T88b   888  888 888  888      888  888 .d888888 Y88  88P 888 888  888 .d888888 888   888  888 888
        // 888  T88b  Y88b 888 888  888      888  888 888  888  Y8bd8P  888 Y88b 888 888  888 Y88b. Y88..88P 888
        // 888   T88b  "Y88888 888  888      888  888 "Y888888   Y88P   888  "Y88888 "Y888888  "Y888 "Y88P"  888
        //                                                                       888
        //                                                                  Y8b d88P
        //                                                                   "Y88P"
        start_time = clock();
        run_navigator(pubkey, prvkey, num_sensors, track_filepath, output_filepath, &encoding_params, &serialisation_params);
        end_time = clock();

        // Output only the run time of the navigator in seconds
        printf("%lf\n", (end_time-start_time)/1000000.0);

        // fprintf(stderr, "Navigator finished.\n");
        // fprintf(stderr, "Navigation ended in %lf seconds\n", (end_time-start_time)/1000000.0);

        // Ensure all aggregation keys were sent and free array
        for (int s=0; s<num_sensors; s++){
            MPI_Wait(agg_requests+s, MPI_STATUS_IGNORE);
        }
        free(agg_requests);

        // Free Paillier keys
        free_pubkey(pubkey);
        free_prvkey(prvkey);

    // Remaining processes are the sensors
    } else {

        // If no commandline argument initialise paillier and encoding with defaults
        if (argc == 1){
            sensor_filepath_base = SENSOR_FILEPATH_BASE;
            encoding_params.frac_bits = ENCODING_FRAC_BITS_DEFAULT;
            serialisation_params.paillier_bitsize = PAILLIER_BITSIZE_DEFAULT;
            serialisation_params.paillier_max_key_serialisation_chars = PAILLIER_MAX_KEY_SERIALISATION_CHARS_DEFAULT;
            serialisation_params.paillier_max_enc_serialisation_chars = PAILLIER_MAX_ENC_SERIALISATION_CHARS_DEFAULT;
        
        // If arguments provided, use them instead
        } else if (argc == 6){
            serialisation_params.paillier_bitsize = atoi(argv[1]);
            encoding_params.frac_bits = atoi(argv[2]);
            sensor_filepath_base = argv[5];
            
            serialisation_params.paillier_max_key_serialisation_chars = (serialisation_params.paillier_bitsize/2) + 1;
            serialisation_params.paillier_max_enc_serialisation_chars = (serialisation_params.paillier_bitsize/2) + 1;
        
        // Either all or no arguments must be provided, error otherwise
        } else {
            fprintf(stderr, "Incorrect commandline usage!");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // 8888888b.
        // 888   Y88b
        // 888    888
        // 888   d88P 888  888 88888b.       .d8888b   .d88b.  88888b.  .d8888b   .d88b.  888d888 .d8888b
        // 8888888P"  888  888 888 "88b      88K      d8P  Y8b 888 "88b 88K      d88""88b 888P"   88K
        // 888 T88b   888  888 888  888      "Y8888b. 88888888 888  888 "Y8888b. 888  888 888     "Y8888b.
        // 888  T88b  Y88b 888 888  888           X88 Y8b.     888  888      X88 Y88..88P 888          X88
        // 888   T88b  "Y88888 888  888       88888P'  "Y8888  888  888  88888P'  "Y88P"  888      88888P'



        run_sensor(proc_id, sensor_filepath_base, &encoding_params, &serialisation_params);
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