/*
 * 
 */

#include "navigator.h"


// Private functions
void init_filter_model(gsl_matrix *F, gsl_matrix *Q);
void broadcast_all_enc_state_vars(pubkey_t *pubkey, ciphertext_t *ct, char *enc_str, gsl_vector *state, encoding_params_t *encoding_params, paillier_serialisation_params_t *serialisation_params);
void broadcast_enc_state_var(pubkey_t *pubkey, ciphertext_t *ct, char *enc_str, double val, encoding_params_t *encoding_params, paillier_serialisation_params_t *serialisation_params);
void print_all_enc_state_vars(gsl_vector *state);
void get_enc_str_from_sensor(int src_sensor, int rows, int cols, char *enc_strs, int send_tag, MPI_Request *r, paillier_serialisation_params_t *serialisation_params);
void reset_recieved_flags(int *received_flags, int num_sensors, int *first_received);
int are_all_received(int *received_flags, int num_sensors);
void poll_sensor_and_aggregate(pubkey_t *pubkey, 
                                int *received_flag, 
                                int *first_received_flag,
                                MPI_Request *request,
                                int dim1, 
                                int dim2, 
                                char *enc_mat_str,
                                c_mtrx_t *enc_mat,
                                c_mtrx_t *enc_mat_sum, 
                                paillier_serialisation_params_t *serialisation_params);
void get_c_mtrx_from_enc_str(c_mtrx_t *m, int rows, int cols, char *enc_strs, paillier_serialisation_params_t *serialisation_params);
void copy_matrix(c_mtrx_t *to_copy_to, c_mtrx_t *to_copy_from, int rows, int cols);
void print_gsl_vector(gsl_vector *v, int d);
void print_gsl_matrix(gsl_matrix *m, int r, int c);
void nav_input_err_check(int val, int expected, char *msg);


void run_navigator(pubkey_t *pubkey, prvkey_t *prvkey, int num_sensors, char *track_filename, char *output_filepath, 
                    encoding_params_t *encoding_params, paillier_serialisation_params_t *serialisation_params){

    // track file var
    FILE *track_fp, *output_fp;

    // Sending encrypted state vars
    char *enc_str;
    ciphertext_t *state_enc;

    // Vars for sensor encryptions
    c_mtrx_t **enc_hrhs;
    c_mtrx_t **enc_hrzs;
    c_mtrx_t *enc_hrh_sum;
    c_mtrx_t *enc_hrz_sum;

    // Receiving vars
    char **hrh_enc_strs;
    char **hrz_enc_strs;
    MPI_Request *hrh_requests;
    MPI_Request *hrz_requests;
    int *hrh_received_flags;
    int *hrz_received_flags;
    int first_received_hrh;
    int first_received_hrz;

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
    gsl_vector *hrz_sum;
    gsl_matrix *hrh_sum;
    gsl_vector *tmp_pred_vec;
    gsl_matrix *tmp_pred_mat;
    gsl_permutation *perm;
    int sig_num;

    // 8888888 .d88888b.
    //   888  d88P" "Y88b
    //   888  888     888
    //   888  888     888
    //   888  888     888
    //   888  888     888
    //   888  Y88b. .d88P
    // 8888888 "Y88888P"




    // Open track file for filter init
    track_fp = fopen(track_filename, "r");
    if (track_fp == NULL){
        fprintf(stderr, "Could not open track file \"%s\"!\n", track_filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // output
    output_fp = fopen(output_filepath, "w");
    if (output_fp == NULL){
        fprintf(stderr, "Could not open output file \"%s\"!\n", output_filepath);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Read number of time steps and the state dimension
    nav_input_err_check(fscanf(track_fp, "%d\n%d", &time_steps, &dimension), 2, "Could not read timesteps and dimension from track file!");
    //fprintf(stderr, "0 steps=%d, dimension=%d\n", time_steps, dimension);

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
    fprintf(stderr, "Initial state and covariance\n");
    print_gsl_vector(state, dimension);
    print_gsl_matrix(covariance, dimension, dimension);
    fprintf(stderr, "\n");

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
    enc_str = (char *)malloc(serialisation_params->paillier_max_enc_serialisation_chars * sizeof(char));
    enc_hrhs = (c_mtrx_t **)malloc(num_sensors*sizeof(c_mtrx_t *));
    enc_hrzs = (c_mtrx_t **)malloc(num_sensors*sizeof(c_mtrx_t *));
    for (int s=0; s<num_sensors; s++){
        enc_hrhs[s] = c_mtrx_alloc(dimension, dimension);
        enc_hrzs[s] = c_mtrx_alloc(1, dimension);
    }
    enc_hrh_sum = c_mtrx_alloc(dimension, dimension);
    enc_hrz_sum = c_mtrx_alloc(1, dimension);
    init_c_mtrx(enc_hrh_sum);
    init_c_mtrx(enc_hrz_sum);

    // Allocate memory for recieving serialised encrypted matrices
    hrh_requests = (MPI_Request *)malloc(num_sensors*sizeof(MPI_Request));
    hrz_requests = (MPI_Request *)malloc(num_sensors*sizeof(MPI_Request));
    hrh_enc_strs = (char **)malloc(num_sensors*sizeof(char *));
    hrz_enc_strs = (char **)malloc(num_sensors*sizeof(char *));
    for (int s=0; s<num_sensors; s++){
        hrh_enc_strs[s] = (char *)malloc(dimension*dimension*(serialisation_params->paillier_max_enc_serialisation_chars)*sizeof(char));
        hrz_enc_strs[s] = (char *)malloc(dimension*(serialisation_params->paillier_max_enc_serialisation_chars)*sizeof(char));
    }
    hrh_received_flags = (int *)malloc(num_sensors*sizeof(int));
    hrz_received_flags = (int *)malloc(num_sensors*sizeof(int));

    // Allocate for state and covariance which will be read from file
    // state and covariance allocated above for file input
    info_vec = gsl_vector_alloc(dimension);
    info_mat = gsl_matrix_alloc(dimension, dimension);
    F = gsl_matrix_alloc(dimension, dimension);
    Q = gsl_matrix_alloc(dimension, dimension);
    perm = gsl_permutation_alloc(dimension);
    hrh_sum = gsl_matrix_alloc(dimension, dimension);
    hrz_sum = gsl_vector_alloc(dimension);
    tmp_pred_vec = gsl_vector_alloc(dimension);
    tmp_pred_mat = gsl_matrix_alloc(dimension, dimension);
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
    for(int t=0; t<time_steps; t++){

        // Filter prediction step
        gsl_blas_dgemv(CblasNoTrans, 1.0, F, state, 0.0, tmp_pred_vec);
        gsl_blas_dcopy(tmp_pred_vec, state);
        gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, covariance, F, 0.0, tmp_pred_mat);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp_pred_mat, F, 0.0, covariance);
        gsl_matrix_add(covariance, Q);

        // Encrypt and broadcast the state variables
        broadcast_all_enc_state_vars(pubkey, state_enc, enc_str, state, encoding_params, serialisation_params);

        // Debugging
        print_all_enc_state_vars(state);

        // Async receive of all hrh and hrz matrices from all sensors
        for (int s=0; s<num_sensors; s++){
            get_enc_str_from_sensor(s+1, dimension, dimension, hrh_enc_strs[s], 0, hrh_requests+s, serialisation_params);
            get_enc_str_from_sensor(s+1, 1, dimension, hrz_enc_strs[s], 1, hrz_requests+s, serialisation_params);
        }

        // Construct the hrh and hrz sum matrices as the sensor matrices arrive, in order of arrival/poll
        reset_recieved_flags(hrh_received_flags, num_sensors, &first_received_hrh);
        reset_recieved_flags(hrz_received_flags, num_sensors, &first_received_hrz);

        while (!are_all_received(hrh_received_flags, num_sensors) || !are_all_received(hrz_received_flags, num_sensors)){
            for (int s=0; s<num_sensors; s++){

                // Keep polling and aggregating until all received
                poll_sensor_and_aggregate(pubkey, 
                                hrh_received_flags+s, 
                                &first_received_hrh, 
                                hrh_requests+s,
                                dimension, 
                                dimension, 
                                hrh_enc_strs[s],
                                enc_hrhs[s],
                                enc_hrh_sum, 
                                serialisation_params);

                poll_sensor_and_aggregate(pubkey, 
                                hrz_received_flags+s, 
                                &first_received_hrz, 
                                hrz_requests+s,
                                1, 
                                dimension,
                                hrz_enc_strs[s],
                                enc_hrzs[s],
                                enc_hrz_sum, 
                                serialisation_params);
            }
        }

        // Decrypt the summed sensor matrices and vectors
        decrypt_mtrx(pubkey, prvkey, enc_hrh_sum, hrh_sum, 1, encoding_params);
        decrypt_vctr(pubkey, prvkey, enc_hrz_sum, hrz_sum, 1, encoding_params);
        fprintf(stderr, "Matrix sum:\n");
        print_gsl_matrix(hrh_sum, dimension, dimension);
        fprintf(stderr, "Vector sum:\n");
        print_gsl_vector(hrz_sum, dimension);

        // Filter update step, convert to information filter form then back
        gsl_linalg_LU_decomp(covariance, perm, &sig_num);
        gsl_linalg_LU_invert(covariance, perm, info_mat);
        gsl_blas_dsymv(CblasUpper, 1.0, info_mat, state, 0.0, info_vec);

        gsl_vector_add(info_vec, hrz_sum);
        gsl_matrix_add(info_mat, hrh_sum);

        gsl_linalg_LU_decomp(info_mat, perm, &sig_num);
        gsl_linalg_LU_invert(info_mat, perm, covariance);
        gsl_blas_dsymv(CblasUpper, 1.0, covariance, info_vec, 0.0, state);

        // Output estimated location
        fprintf(stderr, "time: %d\n", t);
        fprintf(stderr, "State\n");
        print_gsl_vector(state, dimension);
        fprintf(stderr, "Covariance\n");
        print_gsl_matrix(covariance, dimension, dimension);
        fprintf(output_fp, "%lf %lf %lf %lf\n", gsl_vector_get(state, 0), gsl_vector_get(state, 1), gsl_vector_get(state, 2), gsl_vector_get(state, 3));
        fprintf(output_fp, "%lf %lf %lf %lf\n", gsl_matrix_get(covariance, 0, 0), gsl_matrix_get(covariance, 0, 1), gsl_matrix_get(covariance, 0, 2), gsl_matrix_get(covariance, 0, 3));
        fprintf(output_fp, "%lf %lf %lf %lf\n", gsl_matrix_get(covariance, 1, 0), gsl_matrix_get(covariance, 1, 1), gsl_matrix_get(covariance, 1, 2), gsl_matrix_get(covariance, 1, 3));
        fprintf(output_fp, "%lf %lf %lf %lf\n", gsl_matrix_get(covariance, 2, 0), gsl_matrix_get(covariance, 2, 1), gsl_matrix_get(covariance, 2, 2), gsl_matrix_get(covariance, 2, 3));
        fprintf(output_fp, "%lf %lf %lf %lf\n", gsl_matrix_get(covariance, 3, 0), gsl_matrix_get(covariance, 3, 1), gsl_matrix_get(covariance, 3, 2), gsl_matrix_get(covariance, 3, 3));

    }

    // 888b     d888                                 8888888888
    // 8888b   d8888                                 888
    // 88888b.d88888                                 888
    // 888Y88888P888  .d88b.  88888b.d88b.           8888888 888d888 .d88b.   .d88b.
    // 888 Y888P 888 d8P  Y8b 888 "888 "88b          888     888P"  d8P  Y8b d8P  Y8b
    // 888  Y8P  888 88888888 888  888  888          888     888    88888888 88888888
    // 888   "   888 Y8b.     888  888  888 d8b      888     888    Y8b.     Y8b.
    // 888       888  "Y8888  888  888  888 Y8P      888     888     "Y8888   "Y8888


    // Close output file
    fclose(output_fp);

    // Free filter vars
    gsl_permutation_free(perm);
    gsl_matrix_free(tmp_pred_mat);
    gsl_vector_free(tmp_pred_vec);
    gsl_matrix_free(hrh_sum);
    gsl_vector_free(hrz_sum);
    gsl_matrix_free(Q);
    gsl_matrix_free(F);
    gsl_matrix_free(info_mat);
    gsl_vector_free(info_vec);
    gsl_matrix_free(covariance);
    gsl_vector_free(state);

    // Free receiving vars
    free(hrh_received_flags);
    free(hrz_received_flags);
    for (int s=0; s<num_sensors; s++){
        free(hrh_enc_strs[s]);
        free(hrz_enc_strs[s]);
    }
    free(hrh_enc_strs);
    free(hrz_enc_strs);
    free(hrh_requests);
    free(hrz_requests);

    // Free encrypted storage vars
    c_mtrx_free(enc_hrh_sum);
    c_mtrx_free(enc_hrz_sum);
    for (int s=0; s<num_sensors; s++){
        c_mtrx_free(enc_hrhs[s]);
        c_mtrx_free(enc_hrzs[s]);
    }
    free(enc_hrhs);
    free(enc_hrzs);
    free(enc_str);
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
void broadcast_all_enc_state_vars(pubkey_t *pubkey, ciphertext_t *ct, char *enc_str, gsl_vector *state, encoding_params_t *encoding_params, paillier_serialisation_params_t *serialisation_params){
    broadcast_enc_state_var(pubkey, ct, enc_str, gsl_vector_get(state, 0), encoding_params, serialisation_params);
    broadcast_enc_state_var(pubkey, ct, enc_str, pow(gsl_vector_get(state, 0), 2), encoding_params, serialisation_params);
    broadcast_enc_state_var(pubkey, ct, enc_str, pow(gsl_vector_get(state, 0), 3), encoding_params, serialisation_params);
    broadcast_enc_state_var(pubkey, ct, enc_str, gsl_vector_get(state, 2), encoding_params, serialisation_params);
    broadcast_enc_state_var(pubkey, ct, enc_str, gsl_vector_get(state, 0)*gsl_vector_get(state, 2), encoding_params, serialisation_params);
    broadcast_enc_state_var(pubkey, ct, enc_str, pow(gsl_vector_get(state, 0), 2)*gsl_vector_get(state, 2), encoding_params, serialisation_params);
    broadcast_enc_state_var(pubkey, ct, enc_str, pow(gsl_vector_get(state, 2), 2), encoding_params, serialisation_params);
    broadcast_enc_state_var(pubkey, ct, enc_str, gsl_vector_get(state, 0)*pow(gsl_vector_get(state, 2), 2), encoding_params, serialisation_params);
    broadcast_enc_state_var(pubkey, ct, enc_str, pow(gsl_vector_get(state, 2), 3), encoding_params, serialisation_params);
}

// Debugging - Print the state variables x,x2,x3,y,xy,x2y,y2,xy2,y3 in that order
void print_all_enc_state_vars(gsl_vector *state){
    fprintf(stderr, "id: %d sending x = %lf\n", 0, gsl_vector_get(state, 0));
    fprintf(stderr, "id: %d sending x2 = %lf\n", 0, pow(gsl_vector_get(state, 0), 2));
    fprintf(stderr, "id: %d sending x3 = %lf\n", 0, pow(gsl_vector_get(state, 0), 3));
    fprintf(stderr, "id: %d sending y = %lf\n", 0, gsl_vector_get(state, 2));
    fprintf(stderr, "id: %d sending xy = %lf\n", 0, gsl_vector_get(state, 0)*gsl_vector_get(state, 2));
    fprintf(stderr, "id: %d sending x2y = %lf\n", 0, pow(gsl_vector_get(state, 0), 2)*gsl_vector_get(state, 2));
    fprintf(stderr, "id: %d sending y2 = %lf\n", 0, pow(gsl_vector_get(state, 2), 2));
    fprintf(stderr, "id: %d sending xy2 = %lf\n", 0, gsl_vector_get(state, 0)*pow(gsl_vector_get(state, 2), 2));
    fprintf(stderr, "id: %d sending y3 = %lf\n", 0, pow(gsl_vector_get(state, 2), 3));
}

// Encrypt serialise and broadcast a state variable
void broadcast_enc_state_var(pubkey_t *pubkey, ciphertext_t *ct, char *enc_str, double val, encoding_params_t *encoding_params, paillier_serialisation_params_t *serialisation_params){
    encode_and_enc(pubkey, ct, val, 0, encoding_params);
    serialise_encryption(ct, enc_str);
    MPI_Bcast(enc_str, serialisation_params->paillier_max_enc_serialisation_chars, MPI_CHAR, 0, MPI_COMM_WORLD);
}

// Async recieval of serialised encrypted matrix (or vector) from a sensor
void get_enc_str_from_sensor(int src_sensor, int rows, int cols, char *enc_strs, int send_tag, MPI_Request *r, paillier_serialisation_params_t *serialisation_params){
    MPI_Irecv(enc_strs, rows*cols*(serialisation_params->paillier_max_enc_serialisation_chars), MPI_CHAR, src_sensor, send_tag, MPI_COMM_WORLD, r);
}

// Set all received flags to 0
void reset_recieved_flags(int *received_flags, int num_sensors, int *first_received){
    for (int s=0; s<num_sensors; s++){
        received_flags[s] = 0;
    }
    *first_received = 0;
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

// Checking and adding encryptions to information filter sums
void poll_sensor_and_aggregate(pubkey_t *pubkey, 
                                int *received_flag, 
                                int *first_received_flag,
                                MPI_Request *request,
                                int dim1, 
                                int dim2, 
                                char *enc_mat_str,
                                c_mtrx_t *enc_mat,
                                c_mtrx_t *enc_mat_sum,
                                paillier_serialisation_params_t *serialisation_params){
    
    // If matrix not yet received from request, poll the receive request again
    int test_request;
    if (!*received_flag){
        MPI_Test(request, &test_request, MPI_STATUS_IGNORE);
        if (test_request){
            *received_flag = 1;
            get_c_mtrx_from_enc_str(enc_mat, dim1, dim2, enc_mat_str, serialisation_params);
            fprintf(stderr, "Received enc\n");

            // If it's the first sensor to send initialise the sum with it's matrix
            if (!*first_received_flag){
                *first_received_flag = 1;
                copy_matrix(enc_mat_sum, enc_mat, dim1, dim2);
            
            // Otherwise add to the sum
            } else {
                add_c_mtrx_c_mtrx(pubkey, enc_mat, enc_mat_sum, enc_mat_sum);
            }
        }
    }
}

// Convert serialised encrypted matrix to encrypted matrix type
void get_c_mtrx_from_enc_str(c_mtrx_t *m, int rows, int cols, char *enc_strs, paillier_serialisation_params_t *serialisation_params){
    char *ind;
    for (int r=0; r<rows; r++){
        for (int c=0; c<cols; c++){
            ind = enc_strs+(r*cols*(serialisation_params->paillier_max_enc_serialisation_chars)+c*(serialisation_params->paillier_max_enc_serialisation_chars));
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