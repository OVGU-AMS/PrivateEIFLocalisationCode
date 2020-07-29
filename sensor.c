/*
 * 
 */

#include "sensor.h"


// Local functions
void get_bcast_state_var(ciphertext_t *state_var, char *enc_str);
void sensor_input_err_check(int val, int expected, char *msg, int id);


void run_sensor(int id){
    // Key vars
    char key_str[MAX_KEY_SERIALISATION_CHARS];
    pubkey_t *pubkey;
    aggkey_t agg_key;

    // Measurements file var
    FILE *measurements_fp;

    // Get Paillier public key
    MPI_Bcast(key_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, MPI_COMM_WORLD);
    pubkey = deserialise_pubkey(key_str);

    printf("id: %d, key: %s\n", id, key_str);

    // Get sensor's private aggregation key
    MPI_Recv(key_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    deserialise_aggkey(agg_key, key_str);

    printf("id: %d, agg key: %s\n", id, key_str);

    // Open sensor measurements file - TODO currently setup for debugging input only!
    char f_name[100];
    sprintf(f_name, "input/debug_sensor%d.txt", id);
    measurements_fp = fopen(f_name, "r");
    if (measurements_fp == NULL){
        fprintf(stderr, "%d Could not open measurement file!\n", id);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Number of timesteps
    int time_steps;
    sensor_input_err_check(fscanf(measurements_fp, "%d", &time_steps), 1, "Could not read timesteps!", id);
    //printf("%d steps=%d\n", id, time_steps);

    // State dimension
    int state_dim;
    sensor_input_err_check(fscanf(measurements_fp, "%d", &state_dim), 1, "Could not read state dimension!", id);
    //printf("%d state_dim=%d\n", id, state_dim);

    // Sensor location
    double loc_x, loc_y;
    sensor_input_err_check(fscanf(measurements_fp, "%lf %lf", &loc_x, &loc_y), 2, "Could not read sensor location!", id);
    printf("%d loc=(%lf, %lf)\n", id, loc_x, loc_y);

    // State info vars
    char enc_str[MAX_ENC_SERIALISATION_CHARS];
    ciphertext_t *x = NULL;
    ciphertext_t *x2 = NULL;
    ciphertext_t *x3 = NULL;
    ciphertext_t *y = NULL;
    ciphertext_t *xy = NULL;
    ciphertext_t *x2y = NULL;
    ciphertext_t *y2 = NULL;
    ciphertext_t *xy2 = NULL;
    ciphertext_t *y3 = NULL;

    // Result vars
    //c_mtrx_t *hrh;
    //c_mtrx_t *hrz;
    ciphertext_t *partial_sum;
    double weight;

    // Measreument vars
    double measurement;
    double z;
    double inv_R_adj;

    for (int t=0; t<time_steps; t++){
        // Get state variable encryption broadcasts x,x2,x3,y,xy,x2y,y2,xy2,y3 in that order
        get_bcast_state_var(x,   enc_str);
        get_bcast_state_var(x2,  enc_str);
        get_bcast_state_var(x3,  enc_str);
        get_bcast_state_var(y,   enc_str);
        get_bcast_state_var(xy,  enc_str);
        get_bcast_state_var(x2y, enc_str);
        get_bcast_state_var(y2,  enc_str);
        get_bcast_state_var(xy2, enc_str);
        get_bcast_state_var(y3,  enc_str);

        // Get next measurement
        sensor_input_err_check(fscanf(measurements_fp, "%lf", &measurement), 1, "Could not read measurement!", id);
        //printf("%d m=%lf\n", id, measurement);

        // Modified measurement for filter, adjusted to be zero-mean
        z = pow(measurement, 2) - SENSOR_VARIANCE;

        // TODO noise approx may be doable better (overestimate for better consistency?)
        // Measreument noise of the adjusted measreument, approximated by using measurement for distance
        inv_R_adj = 1.0/(2*(2*(pow(measurement, 2) - SENSOR_VARIANCE)*SENSOR_VARIANCE + pow(SENSOR_VARIANCE, 2)));

        // Compute HRH weighted sums

    }


    // Done with measurements, close file and return
    fclose(measurements_fp);
}


void get_bcast_state_var(ciphertext_t *state_var, char *enc_str){
    MPI_Bcast(enc_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, MPI_COMM_WORLD);
    state_var = deserialise_encryption(enc_str);
}

void sensor_input_err_check(int val, int expected, char *msg, int id){
    if (val != expected){
        fprintf(stderr, "%d : %s\n", id, msg);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}