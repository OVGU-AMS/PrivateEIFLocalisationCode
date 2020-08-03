/*
 * 
 */

#include "navigator.h"


// Private functions
void init_filter_model(gsl_matrix *F, gsl_matrix *Q);
void broadcast_all_enc_state_vars(pubkey_t *pubkey, ciphertext_t *ct, char *enc_str, gsl_vector *state);
void broadcast_enc_state_var(pubkey_t *pubkey, ciphertext_t *ct, char *enc_str, double val);
void get_enc_str_from_sensor(int src_sensor, int rows, int cols, char *enc_strs, int send_tag, MPI_Request *r);
void reset_recieved_flags(int *received_flags, int num_sensors);
int are_all_received(int *received_flags, int num_sensors);
void get_c_mtrx_from_enc_str(c_mtrx_t *m, int rows, int cols, char *enc_strs);
void copy_matrix(c_mtrx_t *to_copy_to, c_mtrx_t *to_copy_from, int rows, int cols);
void print_gsl_vector(gsl_vector *v, int d);
void print_gsl_matrix(gsl_matrix *m, int r, int c);
void nav_input_err_check(int val, int expected, char *msg);


void run_navigator(pubkey_t *pubkey, prvkey_t *prvkey, int num_sensors){
    // track file var
    FILE *track_fp;
    char f_name[100];

    // Sending encrypted state vars
    char enc_str[MAX_KEY_SERIALISATION_CHARS];
    ciphertext_t *state_enc;

    // Vars for sensor encryptions
    c_mtrx_t **hrhs;
    c_mtrx_t **hrzs;
    c_mtrx_t *hrh_sum;
    c_mtrx_t *hrz_sum;

    // Receiving vars
    char **hrh_enc_strs;
    char **hrz_enc_strs;
    MPI_Request *hrh_requests;
    MPI_Request *hrz_requests;
    int *received_hrh;
    int *received_hrz;
    int first_received_hrh;
    int first_received_hrz;
    int test_request;

    // Sim vars
    int time_steps;
    int dimension;

    // Filter vars
    gsl_matrix *F;
    gsl_matrix *Q;
    gsl_vector *state;
    gsl_matrix *covariance;
    gsl_vector *info_vec;
    gsl_matrix *info_mat;
    gsl_vector *sen_hrz_sum;
    gsl_matrix *sen_hrh_sum;
    gsl_permutation *perm;
    int sig_num;

    // 8888888                            888
    //   888                              888
    //   888                              888
    //   888   88888b.  88888b.  888  888 888888
    //   888   888 "88b 888 "88b 888  888 888
    //   888   888  888 888  888 888  888 888
    //   888   888  888 888 d88P Y88b 888 Y88b.
    // 8888888 888  888 88888P"   "Y88888  "Y888
    //                  888
    //                  888
    //                  888

    // Open track file for filter init - TODO currently setup for debugging input only!
    sprintf(f_name, "input/debug_track1.txt");
    track_fp = fopen(f_name, "r");
    if (track_fp == NULL){
        fprintf(stderr, "Could not open track file!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Read number of time steps and the state dimension
    nav_input_err_check(fscanf(track_fp, "%d\n%d", &time_steps, &dimension), 2, "Could not read timesteps and dimension from track file!");
    //printf("0 steps=%d, dimension=%d\n", time_steps, dimension);

    // Allocate state and covariance first as they must be populated from file
    state = gsl_vector_alloc(dimension);
    covariance = gsl_matrix_alloc(dimension, dimension);

    // Read initial state and covariance from track file
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

    // 888b     d888                                        d8888 888 888
    // 8888b   d8888                                       d88888 888 888
    // 88888b.d88888                                      d88P888 888 888
    // 888Y88888P888  .d88b.  88888b.d88b.               d88P 888 888 888  .d88b.   .d8888b
    // 888 Y888P 888 d8P  Y8b 888 "888 "88b             d88P  888 888 888 d88""88b d88P"
    // 888  Y8P  888 88888888 888  888  888            d88P   888 888 888 888  888 888
    // 888   "   888 Y8b.     888  888  888 d8b       d8888888888 888 888 Y88..88P Y88b.
    // 888       888  "Y8888  888  888  888 Y8P      d88P     888 888 888  "Y88P"   "Y8888P




    // Allocate memory for sending the encrypted state
    state_enc = init_ciphertext();

    // Allocate memory for recieved encrypted matrices
    hrhs = (c_mtrx_t **)malloc(num_sensors*sizeof(c_mtrx_t *));
    hrzs = (c_mtrx_t **)malloc(num_sensors*sizeof(c_mtrx_t *));
    for (int s=0; s<num_sensors; s++){
        hrhs[s] = c_mtrx_alloc(dimension, dimension);
        hrzs[s] = c_mtrx_alloc(1, dimension);
    }
    hrh_sum = c_mtrx_alloc(dimension, dimension);
    hrz_sum = c_mtrx_alloc(1, dimension);
    init_c_mtrx(hrh_sum);
    init_c_mtrx(hrz_sum);

    // Allocate memory for recieving serialised encrypted matrices
    hrh_requests = (MPI_Request *)malloc(num_sensors*sizeof(MPI_Request));
    hrz_requests = (MPI_Request *)malloc(num_sensors*sizeof(MPI_Request));
    hrh_enc_strs = (char **)malloc(num_sensors*sizeof(char *));
    hrz_enc_strs = (char **)malloc(num_sensors*sizeof(char *));
    for (int s=0; s<num_sensors; s++){
        hrh_enc_strs[s] = (char *)malloc(dimension*dimension*MAX_ENC_SERIALISATION_CHARS*sizeof(char));
        hrz_enc_strs[s] = (char *)malloc(dimension*MAX_ENC_SERIALISATION_CHARS*sizeof(char));
    }
    received_hrh = (int *)malloc(num_sensors*sizeof(int));
    received_hrz = (int *)malloc(num_sensors*sizeof(int));

    // Allocate for state and covariance which will be read from file
    // state and covariance allocated above for file input
    info_vec = gsl_vector_alloc(dimension);
    info_mat = gsl_matrix_alloc(dimension, dimension);
    F = gsl_matrix_alloc(dimension, dimension);
    Q = gsl_matrix_alloc(dimension, dimension);
    perm = gsl_permutation_alloc(dimension);
    sen_hrh_sum = gsl_matrix_alloc(dimension, dimension);
    sen_hrz_sum = gsl_vector_alloc(dimension);
    init_filter_model(F, Q);

    //  .d8888b.  d8b
    // d88P  Y88b Y8P
    // Y88b.
    //  "Y888b.   888 88888b.d88b.
    //     "Y88b. 888 888 "888 "88b
    //       "888 888 888  888  888
    // Y88b  d88P 888 888  888  888
    //  "Y8888P"  888 888  888  888




    // Sync with all processes before simulation begins
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Begin tracking
    //time_steps = 2; // TODO TEMP
    for(int t=0; t<time_steps; t++){

        // Filter prediction step
        gsl_blas_dgemv(CblasNoTrans, 1.0, F, state, 0.0, state);
        gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, covariance, F, 0.0, covariance);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, covariance, F, 0.0, covariance);
        gsl_matrix_add(covariance, Q);

        // Encrypt and broadcast the state variables
        broadcast_all_enc_state_vars(pubkey, state_enc, enc_str, state);

        // Async receive of all hrh and hrz matrices from all sensors
        for (int s=0; s<num_sensors; s++){
            get_enc_str_from_sensor(s+1, dimension, dimension, hrh_enc_strs[s], 0, hrh_requests+s);
            get_enc_str_from_sensor(s+1, 1, dimension, hrz_enc_strs[s], 1, hrz_requests+s);
        }

        // Construct the hrh and hrz sum matrices as the sensor matrices arrive in order of arrival/poll
        reset_recieved_flags(received_hrh, num_sensors);
        reset_recieved_flags(received_hrz, num_sensors);
        first_received_hrh = 0;
        first_received_hrz = 0;
        // Keep looping sensors until all received
        while (!are_all_received(received_hrh, num_sensors) || !are_all_received(received_hrz, num_sensors)){
            for (int s=0; s<num_sensors; s++){

                // For each sensor where hrh not yet received from, poll the receive request
                if (!received_hrh[s]){
                    MPI_Test(hrh_requests+s, &test_request, MPI_STATUS_IGNORE);
                    if (test_request){
                        received_hrh[s] = 1;
                        get_c_mtrx_from_enc_str(hrhs[s], dimension, dimension, hrh_enc_strs[s]);
                        // If it's the first sensor to send initialise the sum with it's matrix
                        if (!first_received_hrh){
                            first_received_hrh = 1;
                            copy_matrix(hrh_sum, hrhs[s], dimension, dimension);
                        // Otherwise add to the sum
                        } else {
                            add_c_mtrx_c_mtrx(pubkey, hrhs[s], hrh_sum, hrh_sum);
                        }
                    }
                }

                // For each sensor where hrz not yet received from, poll the receive request
                if (!received_hrz[s]){
                    MPI_Test(hrz_requests+s, &test_request, MPI_STATUS_IGNORE);
                    if (test_request){
                        received_hrz[s] = 1;
                        get_c_mtrx_from_enc_str(hrzs[s], 1, dimension, hrz_enc_strs[s]);
                        // If it's the first sensor to send initialise the sum with it's matrix
                        if (!first_received_hrz){
                            first_received_hrz = 1;
                            copy_matrix(hrz_sum, hrzs[s], 1, dimension);
                        // Otherwise add to the sum
                        } else {
                            add_c_mtrx_c_mtrx(pubkey, hrzs[s], hrz_sum, hrz_sum);
                        }
                    }
                }
            }
        }


        // TEMP
        // for (int s=0; s<num_sensors; s++){
        //     decrypt_mtrx(pubkey, prvkey, hrhs[s], covariance, 1);
        //     printf("Matrix from %d:\n", s+1);
        //     print_gsl_matrix(covariance, dimension, dimension);

        //     decrypt_vctr(pubkey, prvkey, hrzs[s], state, 1);
        //     printf("Vector from %d:\n", s+1);
        //     print_gsl_vector(state, dimension);
        // }

        // Decrypt the summed sensor matrices and vectors
        decrypt_mtrx(pubkey, prvkey, hrh_sum, sen_hrh_sum, 1);
        decrypt_vctr(pubkey, prvkey, hrz_sum, sen_hrz_sum, 1);
        // printf("Matrix sum:\n");
        // print_gsl_matrix(sen_hrh_sum, dimension, dimension);
        // printf("Vector sum:\n");
        // print_gsl_vector(sen_hrz_sum, dimension);

        // Filter update step, convert to information filter form then back
        gsl_linalg_LU_decomp(covariance, perm, &sig_num);
        gsl_linalg_LU_invert(covariance, perm, info_mat);
        gsl_blas_dsymv(CblasUpper, 1.0, info_mat, state, 0.0, info_vec);

        gsl_vector_add(info_vec, sen_hrz_sum);
        gsl_matrix_add(info_mat, sen_hrh_sum);

        gsl_linalg_LU_decomp(info_mat, perm, &sig_num);
        gsl_linalg_LU_invert(info_mat, perm, covariance);
        gsl_blas_dsymv(CblasUpper, 1.0, covariance, info_vec, 0.0, state);

        // Output estimated location
        printf("State:\n");
        print_gsl_vector(state, dimension);
        printf("Covariance:\n");
        print_gsl_matrix(covariance, dimension, dimension);

    }

    // 888b     d888                                 8888888888
    // 8888b   d8888                                 888
    // 88888b.d88888                                 888
    // 888Y88888P888  .d88b.  88888b.d88b.           8888888 888d888 .d88b.   .d88b.
    // 888 Y888P 888 d8P  Y8b 888 "888 "88b          888     888P"  d8P  Y8b d8P  Y8b
    // 888  Y8P  888 88888888 888  888  888          888     888    88888888 88888888
    // 888   "   888 Y8b.     888  888  888 d8b      888     888    Y8b.     Y8b.
    // 888       888  "Y8888  888  888  888 Y8P      888     888     "Y8888   "Y8888




    // Free filter vars
    gsl_permutation_free(perm);
    gsl_matrix_free(sen_hrh_sum);
    gsl_vector_free(sen_hrz_sum);
    gsl_matrix_free(Q);
    gsl_matrix_free(F);
    gsl_matrix_free(info_mat);
    gsl_vector_free(info_vec);
    gsl_matrix_free(covariance);
    gsl_vector_free(state);

    // Free receiving vars
    free(received_hrh);
    free(received_hrz);
    for (int s=0; s<num_sensors; s++){
        free(hrh_enc_strs[s]);
        free(hrz_enc_strs[s]);
    }
    free(hrh_enc_strs);
    free(hrz_enc_strs);
    free(hrh_requests);
    free(hrz_requests);

    // Free encrypted storage vars
    c_mtrx_free(hrh_sum);
    c_mtrx_free(hrz_sum);
    for (int s=0; s<num_sensors; s++){
        c_mtrx_free(hrhs[s]);
        c_mtrx_free(hrzs[s]);
    }
    free(hrhs);
    free(hrzs);
    free_ciphertext(state_enc);
    
}

// Initialise values for process model matrices F and Q
void init_filter_model(gsl_matrix *F, gsl_matrix *Q){
    // Noise strength
    double q = 0.01;

    // Time step
    double t = 0.5;

    // F matrix
    gsl_matrix_set(F, 0, 0, 1);
    gsl_matrix_set(F, 0, 1, t);
    gsl_matrix_set(F, 0, 2, 0);
    gsl_matrix_set(F, 0, 3, 0);

    gsl_matrix_set(F, 1, 0, 0);
    gsl_matrix_set(F, 1, 1, 1);
    gsl_matrix_set(F, 1, 2, 0);
    gsl_matrix_set(F, 1, 3, 0);

    gsl_matrix_set(F, 2, 0, 0);
    gsl_matrix_set(F, 2, 1, 0);
    gsl_matrix_set(F, 2, 2, 1);
    gsl_matrix_set(F, 2, 3, t);

    gsl_matrix_set(F, 3, 0, 0);
    gsl_matrix_set(F, 3, 1, 0);
    gsl_matrix_set(F, 3, 2, 0);
    gsl_matrix_set(F, 3, 3, 1);

    // Q matrix
    gsl_matrix_set(Q, 0, 0, pow(t, 3)/3);
    gsl_matrix_set(Q, 0, 1, pow(t, 2)/2);
    gsl_matrix_set(Q, 0, 2, 0);
    gsl_matrix_set(Q, 0, 3, 0);

    gsl_matrix_set(Q, 1, 0, pow(t, 2)/2);
    gsl_matrix_set(Q, 1, 1, t);
    gsl_matrix_set(Q, 1, 2, 0);
    gsl_matrix_set(Q, 1, 3, 0);

    gsl_matrix_set(Q, 2, 0, 0);
    gsl_matrix_set(Q, 2, 1, 0);
    gsl_matrix_set(Q, 2, 2, pow(t, 3)/3);
    gsl_matrix_set(Q, 2, 3, pow(t, 2)/2);

    gsl_matrix_set(Q, 3, 0, 0);
    gsl_matrix_set(Q, 3, 1, 0);
    gsl_matrix_set(Q, 3, 2, pow(t, 2)/2);
    gsl_matrix_set(Q, 3, 3, t);

    gsl_matrix_scale(Q, q);
}

// Encrypt and broadcast the state variables x,x2,x3,y,xy,x2y,y2,xy2,y3 in that order
void broadcast_all_enc_state_vars(pubkey_t *pubkey, ciphertext_t *ct, char *enc_str, gsl_vector *state){
    broadcast_enc_state_var(pubkey, ct, enc_str, gsl_vector_get(state, 0));
    broadcast_enc_state_var(pubkey, ct, enc_str, pow(gsl_vector_get(state, 0), 2));
    broadcast_enc_state_var(pubkey, ct, enc_str, pow(gsl_vector_get(state, 0), 3));
    broadcast_enc_state_var(pubkey, ct, enc_str, gsl_vector_get(state, 2));
    broadcast_enc_state_var(pubkey, ct, enc_str, gsl_vector_get(state, 0)*gsl_vector_get(state, 2));
    broadcast_enc_state_var(pubkey, ct, enc_str, pow(gsl_vector_get(state, 0), 2)*gsl_vector_get(state, 2));
    broadcast_enc_state_var(pubkey, ct, enc_str, pow(gsl_vector_get(state, 2), 2));
    broadcast_enc_state_var(pubkey, ct, enc_str, gsl_vector_get(state, 0)*pow(gsl_vector_get(state, 2), 2));
    broadcast_enc_state_var(pubkey, ct, enc_str, pow(gsl_vector_get(state, 2), 3));
}

// Encrypt serialise and broadcast a state variable
void broadcast_enc_state_var(pubkey_t *pubkey, ciphertext_t *ct, char *enc_str, double val){
    encode_and_enc(pubkey, ct, val, 0);
    serialise_encryption(ct, enc_str);
    MPI_Bcast(enc_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, MPI_COMM_WORLD);
}

// Async recieval of serialised encrypted matrix (or vector) from a sensor
void get_enc_str_from_sensor(int src_sensor, int rows, int cols, char *enc_strs, int send_tag, MPI_Request *r){
    MPI_Irecv(enc_strs, rows*cols*MAX_ENC_SERIALISATION_CHARS, MPI_CHAR, src_sensor, send_tag, MPI_COMM_WORLD, r);
}

// Set all received flags to 0
void reset_recieved_flags(int *received_flags, int num_sensors){
    for (int s=0; s<num_sensors; s++){
        received_flags[s] = 0;
    }
}

// Check if all received flags are true
int are_all_received(int *received_flags, int num_sensors){
    int sum=0;
    for (int s=0; s<num_sensors; s++){
        sum += received_flags[s];
    }
    if (sum == num_sensors){
        return 1;
    } else {
        return 0;
    }
}

// Convert serialised encrypted matrix to encrypted matrix type
void get_c_mtrx_from_enc_str(c_mtrx_t *m, int rows, int cols, char *enc_strs){
    char *ind;
    for (int r=0; r<rows; r++){
        for (int c=0; c<cols; c++){
            ind = enc_strs+(r*cols*MAX_ENC_SERIALISATION_CHARS+c*MAX_ENC_SERIALISATION_CHARS);
            set_c_mtrx(m, r, c, deserialise_encryption(ind));
        }
    }
}

// Initialises matrix and copies values of given matrix to it
void copy_matrix(c_mtrx_t *to_copy_to, c_mtrx_t *to_copy_from, int rows, int cols){
    for (int r=0; r<rows; r++){
        for (int c=0; c<cols; c++){
            copy_encryption(get_c_mtrx(to_copy_to, r, c), get_c_mtrx(to_copy_from, r, c));
        }
    }
}

// Debugging method for printing a gsl vector nicely
void print_gsl_vector(gsl_vector *v, int d){
    fprintf(stderr, "Vector:\n[ ");
    for(int i=0; i<d; i++){
        fprintf(stderr, "%+10.10lf ", gsl_vector_get(v, i));
    }
    fprintf(stderr, "]\n");
}

// Debugging method for printing a gsl matrix nicely
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

// Used to check scanf expected output
void nav_input_err_check(int val, int expected, char *msg){
    if (val != expected){
        fprintf(stderr, "%s\n", msg);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}