/*
 * 
 */

#include "navigator.h"


// Private functions
void broadcast_enc_state_var(pubkey_t *pubkey, ciphertext_t *ct, char *enc_str, double val);
void print_gsl_vector(gsl_vector *v, int d);
void print_gsl_matrix(gsl_matrix *m, int r, int c);
void nav_input_err_check(int val, int expected, char *msg);


void run_navigator(pubkey_t *pubkey, prvkey_t *prvkey){
    // track file var
    FILE *track_fp;

    // Filter state and covariance vars
    gsl_vector *state;
    gsl_matrix *covariance;

    // Open track file for filter init - TODO currently setup for debugging input only!
    char f_name[] = "input/debug_track1.txt";
    track_fp = fopen(f_name, "r");
    if (track_fp == NULL){
        fprintf(stderr, "Could not open track file!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Get number of time steps and the state dimension
    int time_steps;
    int dimension;
    nav_input_err_check(fscanf(track_fp, "%d\n%d", &time_steps, &dimension), 2, "Could not read timesteps and dimension from track file!");
    //printf("0 steps=%d, dimension=%d\n", time_steps, dimension);

    // Initialise state and covariance from track file
    state = gsl_vector_alloc(dimension);
    covariance = gsl_matrix_alloc(dimension, dimension);
    for(int i=0; i<dimension; i++){
        nav_input_err_check(fscanf(track_fp, "%lf", gsl_vector_ptr(state, i)), 1, "Could not read inital state value from track file!");
    }
    for(int i=0; i<dimension; i++){
        for(int j=0; j<dimension; j++){
            nav_input_err_check(fscanf(track_fp, "%lf", gsl_matrix_ptr(covariance, i, j)), 1, "Could not read inital covariance value from track file!");
        }
    }
    print_gsl_vector(state, dimension);
    print_gsl_matrix(covariance, dimension, dimension);

    // Done with init data, close file
    fclose(track_fp);

    // Create state sending variables
    char enc_str[MAX_KEY_SERIALISATION_CHARS];
    ciphertext_t *state_enc = init_ciphertext();

    // Begin tracking
    for(int t=0; t<time_steps; t++){
        // Encrypt and broadcast the state variables x,x2,x3,y,xy,x2y,y2,xy2,y3 in that order
        broadcast_enc_state_var(pubkey, state_enc, enc_str, gsl_vector_get(state, 0));
        broadcast_enc_state_var(pubkey, state_enc, enc_str, pow(gsl_vector_get(state, 0), 2));
        broadcast_enc_state_var(pubkey, state_enc, enc_str, pow(gsl_vector_get(state, 0), 3));
        broadcast_enc_state_var(pubkey, state_enc, enc_str, gsl_vector_get(state, 2));
        broadcast_enc_state_var(pubkey, state_enc, enc_str, gsl_vector_get(state, 0)*gsl_vector_get(state, 2));
        broadcast_enc_state_var(pubkey, state_enc, enc_str, pow(gsl_vector_get(state, 0), 2)*gsl_vector_get(state, 2));
        broadcast_enc_state_var(pubkey, state_enc, enc_str, pow(gsl_vector_get(state, 2), 2));
        broadcast_enc_state_var(pubkey, state_enc, enc_str, gsl_vector_get(state, 0)*pow(gsl_vector_get(state, 2), 2));
        broadcast_enc_state_var(pubkey, state_enc, enc_str, pow(gsl_vector_get(state, 2), 3));
        

    }

    // Free state sending variable
    free_ciphertext(state_enc);
}

// Encrypt serialise and broadcast a state variable
void broadcast_enc_state_var(pubkey_t *pubkey, ciphertext_t *ct, char *enc_str, double val){
    encode_and_enc(pubkey, ct, val, 0);
    serialise_encryption(ct, enc_str);
    MPI_Bcast(enc_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, MPI_COMM_WORLD);
}


void print_gsl_vector(gsl_vector *v, int d){
    fprintf(stderr, "Vector:\n[ ");
    for(int i=0; i<d; i++){
        fprintf(stderr, "%+10.10lf ", gsl_vector_get(v, i));
    }
    fprintf(stderr, "]\n");
}

void print_gsl_matrix(gsl_matrix *m, int r, int c){
    fprintf(stderr, "Matrix:\n");
    for(int i=0; i<r; i++){
        fprintf(stderr, "[ ");
        for(int j=0; j<c; j++){
            fprintf(stderr, "%+10.10lf ", gsl_matrix_get(m, i, j));
        }
        fprintf(stderr, "]\n");
    }
}

void nav_input_err_check(int val, int expected, char *msg){
    if (val != expected){
        fprintf(stderr, "%s\n", msg);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}