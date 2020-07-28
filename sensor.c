/*
 * 
 */

#include "sensor.h"


// Local functions
void get_bcast_state_var(paillier_ciphertext_t *state_var, char *state_var_enc_str);


void run_sensor(int id){
    // Key vars
    char key_str[MAX_KEY_SERIALISATION_CHARS];
    paillier_pubkey_t *pubkey;
    mpz_t agg_key;

    // Measurements file var
    FILE *measurements_fp;

    // Get Paillier public key
    MPI_Bcast(key_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, MPI_COMM_WORLD);
    pubkey = paillier_pubkey_from_hex(key_str);

    //printf("id: %d, key: %s\n", id, paillier_pubkey_to_hex(pubkey));

    // Get sensor's private aggregation key
    MPI_Recv(key_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    mpz_init_set_str(agg_key, key_str, SERIALISATION_BASE);

    //printf("id: %d, agg key: %s\n", id, key_str);

    // Open sensor measurements file - TODO currently setup for debugging input only!
    char f_name[100];
    sprintf(f_name, "input/debug_sensor%d.txt", id);
    measurements_fp = fopen(f_name, "r");
    if (measurements_fp == NULL){
        printf("%d Could not open measurement file!\n", id);
        exit(0);
    }

    // Number of timesteps
    int time_steps;
    fscanf(measurements_fp, "%d", &time_steps);
    //printf("%d steps=%d\n", id, time_steps);

    // State dimension
    int state_dim;
    fscanf(measurements_fp, "%d", &state_dim);
    //printf("%d state_dim=%d\n", id, state_dim);

    // Sensor location
    double loc_x, loc_y;
    fscanf(measurements_fp, "%lf %lf", &loc_x, &loc_y);
    //printf("%d loc=(%lf, %lf)\n", id, loc_x, loc_y);

    // State info vars
    char state_var_enc_str[MAX_ENC_SERIALISATION_CHARS];
    paillier_ciphertext_t *x = paillier_create_enc_zero();
    paillier_ciphertext_t *x2 = paillier_create_enc_zero();
    paillier_ciphertext_t *x3 = paillier_create_enc_zero();
    paillier_ciphertext_t *y = paillier_create_enc_zero();
    paillier_ciphertext_t *xy = paillier_create_enc_zero();
    paillier_ciphertext_t *x2y = paillier_create_enc_zero();
    paillier_ciphertext_t *y2 = paillier_create_enc_zero();
    paillier_ciphertext_t *xy2 = paillier_create_enc_zero();
    paillier_ciphertext_t *y3 = paillier_create_enc_zero();

    // Result vars
    c_mtrx_t *hrh;
    c_mtrx_t *hrz;
    paillier_ciphertext_t *partial_sum;
    paillier_plaintext_t *weight;

    double measurement;
    for (int t=0; t<time_steps; t++){
        // Get state variable encryption broadcasts x,x2,x3,y,xy,x2y,y2,xy2,y3 in that order
        get_bcast_state_var(x,   state_var_enc_str);
        get_bcast_state_var(x2,  state_var_enc_str);
        get_bcast_state_var(x3,  state_var_enc_str);
        get_bcast_state_var(y,   state_var_enc_str);
        get_bcast_state_var(xy,  state_var_enc_str);
        get_bcast_state_var(x2y, state_var_enc_str);
        get_bcast_state_var(y2,  state_var_enc_str);
        get_bcast_state_var(xy2, state_var_enc_str);
        get_bcast_state_var(y3,  state_var_enc_str);

        // Get next measurement
        fscanf(measurements_fp, "%lf", &measurement);
        //printf("%d m=%lf\n", id, measurement);

    }


    // Done with measurements, close file and return
    fclose(measurements_fp);
}


void get_bcast_state_var(paillier_ciphertext_t *state_var, char *state_var_enc_str){
    MPI_Bcast(state_var_enc_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, MPI_COMM_WORLD);
    mpz_set_str(state_var->c, state_var_enc_str, SERIALISATION_BASE);
    gmp_printf("state var: %Zd\n", state_var->c);
}