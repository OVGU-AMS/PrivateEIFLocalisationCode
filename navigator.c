/*
 * 
 */

#include "navigator.h"


// Private functions
void broadcast_enc_state_var(paillier_pubkey_t *pubkey, 
                                paillier_plaintext_t *state_var, 
                                paillier_ciphertext_t *state_var_enc, 
                                char *state_var_enc_str,
                                double a);
void print_gsl_vector(gsl_vector *v, int d);
void print_gsl_matrix(gsl_matrix *m, int r, int c);


void run_navigator(paillier_pubkey_t *pubkey, paillier_prvkey_t *prvkey){
    // track file var
    FILE *track_fp;

    // Filter state and covariance vars
    gsl_vector *state;
    gsl_matrix *covariance;

    // Open track file for filter init - TODO currently setup for debugging input only!
    char f_name[] = "input/debug_track1.txt";
    track_fp = fopen(f_name, "r");
    if (track_fp == NULL){
        printf("Could not open track file!\n");
        exit(0);
    }

    // Get number of time steps and the state dimension
    int time_steps;
    int dimension;
    fscanf(track_fp, "%d\n%d", &time_steps, &dimension);
    printf("0 steps=%d, dimension=%d\n", time_steps, dimension);

    // Initialise state and covariance from track file
    state = gsl_vector_alloc(dimension);
    covariance = gsl_matrix_alloc(dimension, dimension);
    for(int i=0; i<dimension; i++){
        fscanf(track_fp, "%lf", gsl_vector_ptr(state, i));
    }
    for(int i=0; i<dimension; i++){
        for(int j=0; j<dimension; j++){
            fscanf(track_fp, "%lf", gsl_matrix_ptr(covariance, i, j));
        }
    }
    print_gsl_vector(state, dimension);
    print_gsl_matrix(covariance, dimension, dimension);

    // Done with init data, close file
    fclose(track_fp);

    // Create state sending variables
    char state_var_enc_str[MAX_ENC_SERIALISATION_CHARS];
    paillier_plaintext_t *state_var;
    paillier_ciphertext_t *state_var_enc;
    state_var = paillier_plaintext_from_ui(0);
    state_var_enc = paillier_create_enc_zero();

    for(int t=0; t<time_steps; t++){
        // Encrypt and broadcast the state variables x,x2,x3,y,xy,x2y,y2,xy2,y3 in that order
        broadcast_enc_state_var(pubkey, state_var, state_var_enc, state_var_enc_str, gsl_vector_get(state, 0));
        broadcast_enc_state_var(pubkey, state_var, state_var_enc, state_var_enc_str, pow(gsl_vector_get(state, 0), 2));
        broadcast_enc_state_var(pubkey, state_var, state_var_enc, state_var_enc_str, pow(gsl_vector_get(state, 0), 3));
        broadcast_enc_state_var(pubkey, state_var, state_var_enc, state_var_enc_str, gsl_vector_get(state, 2));
        broadcast_enc_state_var(pubkey, state_var, state_var_enc, state_var_enc_str, gsl_vector_get(state, 0)*gsl_vector_get(state, 2));
        broadcast_enc_state_var(pubkey, state_var, state_var_enc, state_var_enc_str, pow(gsl_vector_get(state, 0), 2)*gsl_vector_get(state, 2));
        broadcast_enc_state_var(pubkey, state_var, state_var_enc, state_var_enc_str, pow(gsl_vector_get(state, 2), 2));
        broadcast_enc_state_var(pubkey, state_var, state_var_enc, state_var_enc_str, gsl_vector_get(state, 0)*pow(gsl_vector_get(state, 2), 2));
        broadcast_enc_state_var(pubkey, state_var, state_var_enc, state_var_enc_str, pow(gsl_vector_get(state, 2), 3));


    }

    // Free state sending variables
    paillier_freeplaintext(state_var);
    paillier_freeciphertext(state_var_enc);
}


void broadcast_enc_state_var(paillier_pubkey_t *pubkey, 
                                paillier_plaintext_t *state_var, 
                                paillier_ciphertext_t *state_var_enc, 
                                char *state_var_enc_str,
                                double a){

    encode_from_dbl(state_var->m, a, MOD_BITS, FRAC_BITS);
    paillier_enc(state_var_enc, pubkey, state_var, paillier_get_rand_devurandom);
    mpz_get_str(state_var_enc_str, SERIALISATION_BASE, state_var_enc->c);
    MPI_Bcast(state_var_enc_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, MPI_COMM_WORLD);
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