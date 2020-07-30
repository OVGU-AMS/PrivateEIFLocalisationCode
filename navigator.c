/*
 * 
 */

#include "navigator.h"


// Private functions
void broadcast_enc_state_var(pubkey_t *pubkey, ciphertext_t *ct, char *enc_str, double val);
void get_c_mtrx_str_from_sensor(int src_sensor, int rows, int cols, char *enc_strs, int send_tag, MPI_Request *r);
void get_c_mtrx_from_enc_str(c_mtrx_t *m, int rows, int cols, char *enc_strs);
void print_gsl_vector(gsl_vector *v, int d);
void print_gsl_matrix(gsl_matrix *m, int r, int c);
void nav_input_err_check(int val, int expected, char *msg);


void run_navigator(pubkey_t *pubkey, prvkey_t *prvkey, int num_sensors){
    // track file var
    FILE *track_fp;

    MPI_Request *hrh_requests = (MPI_Request *)malloc(num_sensors*sizeof(MPI_Request));
    MPI_Request *hrz_requests = (MPI_Request *)malloc(num_sensors*sizeof(MPI_Request));

    c_mtrx_t **hrhs = (c_mtrx_t **)malloc(num_sensors*sizeof(c_mtrx_t *));
    c_mtrx_t **hrzs = (c_mtrx_t **)malloc(num_sensors*sizeof(c_mtrx_t *));
    char **hrh_enc_strs = (char **)malloc(num_sensors*sizeof(char *));
    char **hrz_enc_strs = (char **)malloc(num_sensors*sizeof(char *));

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
    printf("Initial state and covariance:\n");
    print_gsl_vector(state, dimension);
    print_gsl_matrix(covariance, dimension, dimension);
    printf("\n");

    // Done with init data, close file
    fclose(track_fp);

    for (int s=0; s<num_sensors; s++){
        hrhs[s] = c_mtrx_alloc(dimension, dimension);
        hrzs[s] = c_mtrx_alloc(1, dimension);
        hrh_enc_strs[s] = (char *)malloc(dimension*dimension*MAX_ENC_SERIALISATION_CHARS*sizeof(char));
        hrz_enc_strs[s] = (char *)malloc(dimension*MAX_ENC_SERIALISATION_CHARS*sizeof(char));
    }
    

    // Create state sending variables
    char enc_str[MAX_KEY_SERIALISATION_CHARS];
    ciphertext_t *state_enc = init_ciphertext();

    // Begin tracking
    time_steps = 2; // TODO TEMP
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
        
        for (int s=0; s<num_sensors; s++){
            get_c_mtrx_str_from_sensor(s+1, dimension, dimension, hrh_enc_strs[s], 0, hrh_requests+s);
            get_c_mtrx_str_from_sensor(s+1, 1, dimension, hrz_enc_strs[s], 1, hrz_requests+s);
        }

        for (int s=0; s<num_sensors; s++){
            MPI_Wait(hrh_requests+s, MPI_STATUS_IGNORE);
            MPI_Wait(hrz_requests+s, MPI_STATUS_IGNORE);
        }


        // we have error here
        for (int s=0; s<num_sensors; s++){
            get_c_mtrx_from_enc_str(hrhs[s], dimension, dimension, hrh_enc_strs[s]);
            get_c_mtrx_from_enc_str(hrzs[s], 1, dimension, hrz_enc_strs[s]);
            decrypt_mtrx(pubkey, prvkey, hrhs[s], covariance, 1);
            print_gsl_matrix(covariance, dimension, dimension);
        }
        
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

void get_c_mtrx_str_from_sensor(int src_sensor, int rows, int cols, char *enc_strs, int send_tag, MPI_Request *r){
    MPI_Irecv(enc_strs, rows*cols*MAX_ENC_SERIALISATION_CHARS, MPI_CHAR, src_sensor, send_tag, MPI_COMM_WORLD, r);
}

void get_c_mtrx_from_enc_str(c_mtrx_t *m, int rows, int cols, char *enc_strs){
    char *ind;
    for (int r=0; r<rows; r++){
        for (int c=0; c<cols; c++){
            ind = enc_strs+(r*cols*MAX_ENC_SERIALISATION_CHARS+c*MAX_ENC_SERIALISATION_CHARS);
            set_c_mtrx(m, r, c, deserialise_encryption(ind));
        }
    }
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