/*
 * 
 */

#include "sensor.h"

// Container class for all the state information sent from the navigator at each time step
typedef struct EncStateInfoTag {
    ciphertext_t *x;
    ciphertext_t *x2;
    ciphertext_t *x3;
    ciphertext_t *y;
    ciphertext_t *xy;
    ciphertext_t *x2y;
    ciphertext_t *y2;
    ciphertext_t *xy2;
    ciphertext_t *y3;
} enc_state_info_t;

// Debugging
typedef struct StateInfoTag {
    double x;
    double x2;
    double x3;
    double y;
    double xy;
    double x2y;
    double y2;
    double xy2;
    double y3;
} state_info_t;


// Local functions
void get_all_bcast_state_vars(enc_state_info_t *s, char *enc_str, paillier_serialisation_params_t *serialisation_params);
void get_bcast_state_var(ciphertext_t **state_var, char *enc_str, paillier_serialisation_params_t *serialisation_params);
void get_and_print_all_bcast_state_vars(int id, pubkey_t *pubkey, prvkey_t *prvkey, enc_state_info_t *enc_s, state_info_t *s, encoding_params_t *encoding_params);
void encode_mult_then_add(pubkey_t *pubkey, ciphertext_t *sum, ciphertext_t *weighted, ciphertext_t *to_weight, double weight, encoding_params_t *encoding_params);
void encode_then_add(pubkey_t *pubkey, ciphertext_t *sum, ciphertext_t *encrypted, double val, encoding_params_t *encoding_params);
void send_enc_mat(c_mtrx_t *mat, int rows, int cols, char *enc_strs, MPI_Request *req, paillier_serialisation_params_t *serialisation_params);
void sensor_input_err_check(int val, int expected, char *msg, int id);


// The sensor
void run_sensor(int id, char *sensor_filepath_base, encoding_params_t *encoding_params, paillier_serialisation_params_t *serialisation_params){
    // Key vars
    char *key_str;
    pubkey_t *pubkey;
    // prvkey_t *prvkey; // Debugging
    aggkey_t agg_key;

    // Measurements file var
    FILE *measurements_fp;
    char f_name[100];

    // Measreument vars
    double measurement;
    double z;
    double inv_R_adj;

    // Sim vars
    int time_steps;
    int state_dim;
    double loc_x, loc_y;

    // Encrypted state info vars
    char *enc_str;
    enc_state_info_t enc_state;
    // state_info_t state; // Debugging

    // Result vars
    c_mtrx_t *hrh;
    c_mtrx_t *hrz;
    ciphertext_t *mat_elem;
    ciphertext_t *partial_sum;
    // gsl_vector *plain_hrz; // Debugging
    // gsl_matrix *plain_hrh; // Debugging

    // Sending vars
    char *hrh_enc_strs;
    char *hrz_enc_strs;
    MPI_Request hrh_request;
    MPI_Request hrz_request;

    // 888    d8P
    // 888   d8P
    // 888  d8P
    // 888d88K      .d88b.  888  888 .d8888b
    // 8888888b    d8P  Y8b 888  888 88K
    // 888  Y88b   88888888 888  888 "Y8888b.
    // 888   Y88b  Y8b.     Y88b 888      X88
    // 888    Y88b  "Y8888   "Y88888  88888P'
    //                           888
    //                      Y8b d88P
    //                       "Y88P"

    // Serialisation memory alloc, needs to be done earlier to receive keys
    key_str = (char *)malloc(serialisation_params->paillier_max_key_serialisation_chars * sizeof(char));

    // Get Paillier public key
    MPI_Bcast(key_str, serialisation_params->paillier_max_key_serialisation_chars, MPI_CHAR, 0, MPI_COMM_WORLD);
    pubkey = deserialise_pubkey(key_str);
    // fprintf(stderr, "id: %d, pubkey: %s\n", id, key_str); // Debugging

    // Debugging - Get Paillier private key
    // MPI_Bcast(key_str, serialisation_params->paillier_max_key_serialisation_chars, MPI_CHAR, 0, MPI_COMM_WORLD);
    // prvkey = deserialise_prvkey(pubkey, key_str);
    // fprintf(stderr, "id: %d, prvkey: %s\n", id, key_str); // Debugging

    // Get sensor's private aggregation key (blocking)
    MPI_Recv(key_str, serialisation_params->paillier_max_key_serialisation_chars, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // fprintf(stderr, "id: %d, agg key str: %s\n", id, key_str); // Debugging
    deserialise_aggkey(agg_key, key_str);

    // Free the key mpi buffer as it's no longer needed
    free(key_str);
    key_str = NULL;

    // Open sensor measurements file
    sprintf(f_name, sensor_filepath_base, id);
    measurements_fp = fopen(f_name, "r");
    if (measurements_fp == NULL){
        fprintf(stderr, "%d Could not open measurement file \"%s\"!\n", id, f_name);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Read number of timesteps
    sensor_input_err_check(fscanf(measurements_fp, "%d", &time_steps), 1, "Could not read timesteps!", id);
    //fprintf(stderr, "%d steps=%d\n", id, time_steps);

    // Read navigator state dimension
    sensor_input_err_check(fscanf(measurements_fp, "%d", &state_dim), 1, "Could not read state dimension!", id);
    //fprintf(stderr, "%d state_dim=%d\n", id, state_dim);

    // Read sensor location
    sensor_input_err_check(fscanf(measurements_fp, "%lf %lf", &loc_x, &loc_y), 2, "Could not read sensor location!", id);
    //fprintf(stderr, "%d loc=(%lf, %lf)\n", id, loc_x, loc_y);

    // 888b     d888                                        d8888 888 888
    // 8888b   d8888                                       d88888 888 888
    // 88888b.d88888                                      d88P888 888 888
    // 888Y88888P888  .d88b.  88888b.d88b.               d88P 888 888 888  .d88b.   .d8888b
    // 888 Y888P 888 d8P  Y8b 888 "888 "88b             d88P  888 888 888 d88""88b d88P"
    // 888  Y8P  888 88888888 888  888  888            d88P   888 888 888 888  888 888
    // 888   "   888 Y8b.     888  888  888 d8b       d8888888888 888 888 Y88..88P Y88b.
    // 888       888  "Y8888  888  888  888 Y8P      d88P     888 888 888  "Y88P"   "Y8888P




    // Alloc for mpi encryption buffer, note the one for keys has already been allocated and freed above
    enc_str = (char *)malloc((serialisation_params->paillier_max_enc_serialisation_chars) * sizeof(char));

    // Allocate space for result vars now that state dimension is known
    hrh = c_mtrx_alloc(state_dim, state_dim);
    hrz = c_mtrx_alloc(1, state_dim);
    init_c_mtrx(hrh);
    init_c_mtrx(hrz);
    partial_sum = init_ciphertext();
    hrh_enc_strs = (char *)malloc(state_dim*state_dim*(serialisation_params->paillier_max_enc_serialisation_chars)*sizeof(char));
    hrz_enc_strs = (char *)malloc(state_dim*(serialisation_params->paillier_max_enc_serialisation_chars)*sizeof(char));

    // Debugging - Alloc for plaintext matrix and vector
    // plain_hrz = gsl_vector_alloc(state_dim);
    // plain_hrh = gsl_matrix_alloc(state_dim, state_dim);

    // Initialise repeated sends
    MPI_Send_init(hrh_enc_strs, state_dim*state_dim*(serialisation_params->paillier_max_enc_serialisation_chars), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &hrh_request);
    MPI_Send_init(hrz_enc_strs, 1*state_dim*(serialisation_params->paillier_max_enc_serialisation_chars), MPI_CHAR, 0, 1, MPI_COMM_WORLD, &hrz_request);
        
    //  .d8888b.  d8b
    // d88P  Y88b Y8P
    // Y88b.
    //  "Y888b.   888 88888b.d88b.
    //     "Y88b. 888 888 "888 "88b
    //       "888 888 888  888  888
    // Y88b  d88P 888 888  888  888
    //  "Y8888P"  888 888  888  888




    // Sync with all processes before beginning simulation
    MPI_Barrier(MPI_COMM_WORLD);

    // Tracking simulation begins
    for (int t=0; t<time_steps; t++){
        // Get next measurement
        sensor_input_err_check(fscanf(measurements_fp, "%lf", &measurement), 1, "Could not read measurement!", id);
        //printf("%d m=%lf\n", id, measurement);

        // Modified measurement for filter, adjusted to be zero-mean
        z = pow(measurement, 2) - SENSOR_VARIANCE;

        // TODO noise approx may be doable better (overestimate for better consistency?)
        // Measreument noise of the adjusted measreument, approximated by using measurement for distance
        inv_R_adj = 1.0/(2*(2*(pow(measurement + 2*sqrt(SENSOR_VARIANCE), 2))*SENSOR_VARIANCE + pow(SENSOR_VARIANCE, 2)));

        // Get all state variable encryption broadcasts x,x2,x3,y,xy,x2y,y2,xy2,y3
        get_all_bcast_state_vars(&enc_state, enc_str, serialisation_params);

        // get_and_print_all_bcast_state_vars(id, pubkey, prvkey, &enc_state, &state, encoding_params); // Debugging

        // 888    888 8888888b.  888    888
        // 888    888 888   Y88b 888    888
        // 888    888 888    888 888    888
        // 8888888888 888   d88P 8888888888
        // 888    888 8888888P"  888    888
        // 888    888 888 T88b   888    888
        // 888    888 888  T88b  888    888
        // 888    888 888   T88b 888    888




        // hrh[0][0] = (4*invRadj)*x2 + (-8*invRadj*sx)*x + (4*invRadj*sx**2)
        // fprintf(stderr, "id: %d (float) h[0][0] = %lf\n", id, (4*inv_R_adj)*state.x2 + (-8*inv_R_adj*loc_x)*state.x + (4*inv_R_adj*pow(loc_x,2))); // Debugging
        mat_elem = get_c_mtrx(hrh, 0, 0);
        encode_and_mult_enc(pubkey, mat_elem, enc_state.x2, 4*inv_R_adj, 0, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x, -8*inv_R_adj*loc_x, encoding_params);
        encode_then_add(pubkey, mat_elem, partial_sum, 4*inv_R_adj*pow(loc_x, 2), encoding_params);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrh[0][2] = (4*invRadj)*xy + (-4*invRadj*sy)*x + (-4*invRadj*sx)*y + (4*invRadj*sx*sy)
        // fprintf(stderr, "id: %d (float) h[0][2] = %lf\n", id, (4*inv_R_adj)*state.xy + (-4*inv_R_adj*loc_y)*state.x + (-4*inv_R_adj*loc_x)*state.y + (4*inv_R_adj*loc_x*loc_y)); // Debugging
        mat_elem = get_c_mtrx(hrh, 0, 2);
        encode_and_mult_enc(pubkey, mat_elem, enc_state.xy, 4*inv_R_adj, 0, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x, -4*inv_R_adj*loc_y, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y, -4*inv_R_adj*loc_x, encoding_params);
        encode_then_add(pubkey, mat_elem, partial_sum, 4*inv_R_adj*loc_x*loc_y, encoding_params);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrh[2][0] = (4*invRadj)*xy + (-4*invRadj*sy)*x + (-4*invRadj*sx)*y + (4*invRadj*sx*sy)
        // fprintf(stderr, "id: %d (float) h[2][0] = %lf\n", id, (4*inv_R_adj)*state.xy + (-4*inv_R_adj*loc_y)*state.x + (-4*inv_R_adj*loc_x)*state.y + (4*inv_R_adj*loc_x*loc_y)); // Debugging
        mat_elem = get_c_mtrx(hrh, 2, 0);
        copy_encryption(mat_elem, get_c_mtrx(hrh, 0, 2));
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrh[2][2] = (4*invRadj)*y2 + (-8*invRadj*sy)*y + (4*invRadj*sy**2)
        // fprintf(stderr, "id: %d (float) h[2][2] = %lf\n", id, (4*inv_R_adj)*state.y2 + (-8*inv_R_adj*loc_y)*state.y + (4*inv_R_adj*pow(loc_y,2))); // Debugging
        mat_elem = get_c_mtrx(hrh, 2, 2);
        encode_and_mult_enc(pubkey, mat_elem, enc_state.y2, 4*inv_R_adj, 0, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y, -8*inv_R_adj*loc_y, encoding_params);
        encode_then_add(pubkey, mat_elem, partial_sum, 4*inv_R_adj*pow(loc_y, 2), encoding_params);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrh remaining elements are 0
        mat_elem = get_c_mtrx(hrh, 0, 1);
        encrypt_zero(pubkey, mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrh, 0, 3), mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrh, 2, 1), mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrh, 2, 3), mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrh, 1, 0), mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrh, 1, 1), mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrh, 1, 2), mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrh, 3, 0), mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrh, 1, 3), mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrh, 3, 1), mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrh, 3, 2), mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrh, 3, 3), mat_elem);

        // 888    888 8888888b.    .d88                888        .d88          88b.           888    888          88b.
        // 888    888 888   Y88b  d88P"                888       d88P"          "Y88b          888    888          "Y88b
        // 888    888 888    888 d88P                  888      d88P              Y88b         888    888            Y88b
        // 8888888888 888   d88P 888    88888888       88888b.  888    888  888    888   888   8888888888 888  888    888
        // 888    888 8888888P"  888       d88P        888 "88b 888    `Y8bd8P'    888 8888888 888    888 `Y8bd8P'    888
        // 888    888 888 T88b   Y88b     d88P  888888 888  888 Y88b     X88K     d88P   888   888    888   X88K     d88P
        // 888    888 888  T88b   Y88b.  d88P          888  888  Y88b. .d8""8b. .d88P          888    888 .d8""8b. .d88P
        // 888    888 888   T88b   "Y88 88888888       888  888   "Y88 888  888 88P"           888    888 888  888 88P"




        // hrz[0] = (2*invRadj*z)*x + (-2*invRadj*sx*z) + (-2*invRadj*sx**2)*x + (2*invRadj*sx**3) + 
        //          (-2*invRadj*sy**2)*x + (2*invRadj*sx*sy**2) + (2*invRadj)*x3 + (-2*invRadj*sx)*x2 + 
        //          (2*invRadj)*xy2 + (-2*invRadj*sx)*y2
        // fprintf(stderr, "id: %d (float) h[0] = %lf\n", id, (2*inv_R_adj*z)*state.x + (-2*inv_R_adj*loc_x*z) + (-2*inv_R_adj*pow(loc_x,2))*state.x + (2*inv_R_adj*pow(loc_x,3)) 
        //                                             + (-2*inv_R_adj*pow(loc_y,2))*state.x + (2*inv_R_adj*loc_x*pow(loc_y,2)) + (2*inv_R_adj)*state.x3 + (-2*inv_R_adj*loc_x)*state.x2
        //                                             + (2*inv_R_adj)*state.xy2 + (-2*inv_R_adj*loc_x)*state.y2); // Debugging
        mat_elem = get_c_mtrx(hrz, 0, 0);
        encode_and_mult_enc(pubkey, mat_elem, enc_state.x, 2*inv_R_adj*z, 0, encoding_params);
        encode_then_add(pubkey, mat_elem, partial_sum, -2*inv_R_adj*loc_x*z, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x, -2*inv_R_adj*pow(loc_x, 2), encoding_params);
        encode_then_add(pubkey, mat_elem, partial_sum, 2*inv_R_adj*pow(loc_x, 3), encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x, -2*inv_R_adj*pow(loc_y, 2), encoding_params);
        encode_then_add(pubkey, mat_elem, partial_sum, 2*inv_R_adj*loc_x*pow(loc_y, 2), encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x3, 2*inv_R_adj, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x2, -2*inv_R_adj*loc_x, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.xy2, 2*inv_R_adj, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y2, -2*inv_R_adj*loc_x, encoding_params);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrz[2] = (2*invRadj*z)*y + (-2*invRadj*sy*z) + (-2*invRadj*sx**2)*y + (2*invRadj*sy*sx**2) + 
        //          (-2*invRadj*sy**2)*y + (2*invRadj*sy**3) + (2*invRadj)*x2y + (-2*invRadj*sy)*x2 + 
        //          (2*invRadj)*y3 + (-2*invRadj*sy)*y2
        // fprintf(stderr, "id: %d (float) h[2] = %lf\n", id, (2*inv_R_adj*z)*state.y + (-2*inv_R_adj*loc_y*z) + (-2*inv_R_adj*pow(loc_x,2))*state.y + (2*inv_R_adj*loc_y*pow(loc_x,2)) 
        //                                             + (-2*inv_R_adj*pow(loc_y,2))*state.y + (2*inv_R_adj*pow(loc_y,3)) + (2*inv_R_adj)*state.x2y + (-2*inv_R_adj*loc_y)*state.x2
        //                                             + (2*inv_R_adj)*state.y3 + (-2*inv_R_adj*loc_y)*state.y2); // Debugging
        mat_elem = get_c_mtrx(hrz, 0, 2);
        encode_and_mult_enc(pubkey, mat_elem, enc_state.y, 2*inv_R_adj*z, 0, encoding_params);
        encode_then_add(pubkey, mat_elem, partial_sum, -2*inv_R_adj*loc_y*z, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y, -2*inv_R_adj*pow(loc_x, 2), encoding_params);
        encode_then_add(pubkey, mat_elem, partial_sum, 2*inv_R_adj*loc_y*pow(loc_x, 2), encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y, -2*inv_R_adj*pow(loc_y, 2), encoding_params);
        encode_then_add(pubkey, mat_elem, partial_sum, 2*inv_R_adj*pow(loc_y, 3), encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x2y, 2*inv_R_adj, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x2, -2*inv_R_adj*loc_y, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y3, 2*inv_R_adj, encoding_params);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y2, -2*inv_R_adj*loc_y, encoding_params);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrz remaining elements are 0
        mat_elem = get_c_mtrx(hrz, 0, 1);
        encrypt_zero(pubkey, mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrz, 0, 3), mat_elem);

        // Debugging
        // decrypt_mtrx(pubkey, prvkey, hrh, plain_hrh, 1, encoding_params);
        // decrypt_vctr(pubkey, prvkey, hrz, plain_hrz, 1, encoding_params);
        // for(int i=0; i<state_dim; i++){
        //     fprintf(stderr, "id: %d (encrp) h[%d] = %+10.10lf\n", id, i, gsl_vector_get(plain_hrz, i));
        // }
        // for(int i=0; i<state_dim; i++){
        //     for(int j=0; j<state_dim; j++){
        //         fprintf(stderr, "id: %d (encrp) h[%d][%d] = %+10.10lf\n", id, i, j, gsl_matrix_get(plain_hrh, i, j));
        //     }
        // }

        // Add aggregation noise
        add_agg_noise_c_mtrx(pubkey, agg_key, hrh, t, 0);
        add_agg_noise_c_mtrx(pubkey, agg_key, hrz, t, 1);

        // Send encrypted matrix and vector
        send_enc_mat(hrh, state_dim, state_dim, hrh_enc_strs, &hrh_request, serialisation_params);
        send_enc_mat(hrz, 1, state_dim, hrz_enc_strs, &hrz_request, serialisation_params);
        // fprintf(stderr, "\nSensor %d sending : %s\n\n", id, hrh_enc_strs); // Debugging
    }

    // 888b     d888                                 8888888888
    // 8888b   d8888                                 888
    // 88888b.d88888                                 888
    // 888Y88888P888  .d88b.  88888b.d88b.           8888888 888d888 .d88b.   .d88b.
    // 888 Y888P 888 d8P  Y8b 888 "888 "88b          888     888P"  d8P  Y8b d8P  Y8b
    // 888  Y8P  888 88888888 888  888  888          888     888    88888888 88888888
    // 888   "   888 Y8b.     888  888  888 d8b      888     888    Y8b.     Y8b.
    // 888       888  "Y8888  888  888  888 Y8P      888     888     "Y8888   "Y8888




    // Done with repeated sends, wait for last ones to end then free request structs
    MPI_Wait(&hrh_request, MPI_STATUS_IGNORE);
    MPI_Wait(&hrz_request, MPI_STATUS_IGNORE);
    MPI_Request_free(&hrh_request);
    MPI_Request_free(&hrz_request);

    // Free sending vars
    free(hrh_enc_strs);
    free(hrz_enc_strs);

    // Free encrypted matrix and vector vars
    free_ciphertext(partial_sum);
    c_mtrx_free(hrh);
    c_mtrx_free(hrz);

    // Free mpi encryption buffer
    free(enc_str);

    // Done with measurements, close file and finish
    fclose(measurements_fp);
}

// Get state variable encryption broadcasts x,x2,x3,y,xy,x2y,y2,xy2,y3 in that order
void get_all_bcast_state_vars(enc_state_info_t *s, char *enc_str, paillier_serialisation_params_t *serialisation_params){
    get_bcast_state_var(&(s->x),   enc_str, serialisation_params);
    get_bcast_state_var(&(s->x2),  enc_str, serialisation_params);
    get_bcast_state_var(&(s->x3),  enc_str, serialisation_params);
    get_bcast_state_var(&(s->y),   enc_str, serialisation_params);
    get_bcast_state_var(&(s->xy),  enc_str, serialisation_params);
    get_bcast_state_var(&(s->x2y), enc_str, serialisation_params);
    get_bcast_state_var(&(s->y2),  enc_str, serialisation_params);
    get_bcast_state_var(&(s->xy2), enc_str, serialisation_params);
    get_bcast_state_var(&(s->y3),  enc_str, serialisation_params);
}

// Debugging - Decrypt and show broadcasts x,x2,x3,y,xy,x2y,y2,xy2,y3 in that order
void get_and_print_all_bcast_state_vars(int id, pubkey_t *pubkey, prvkey_t *prvkey, enc_state_info_t *enc_s, state_info_t *s, encoding_params_t *encoding_params){
    s->x = dec_and_decode(pubkey, prvkey, enc_s->x, 0, encoding_params);
    s->x2 = dec_and_decode(pubkey, prvkey, enc_s->x2, 0, encoding_params);
    s->x3 = dec_and_decode(pubkey, prvkey, enc_s->x3, 0, encoding_params);
    s->y = dec_and_decode(pubkey, prvkey, enc_s->y, 0, encoding_params);
    s->xy = dec_and_decode(pubkey, prvkey, enc_s->xy, 0, encoding_params);
    s->x2y = dec_and_decode(pubkey, prvkey, enc_s->x2y, 0, encoding_params);
    s->y2 = dec_and_decode(pubkey, prvkey, enc_s->y2, 0, encoding_params);
    s->xy2 = dec_and_decode(pubkey, prvkey, enc_s->xy2, 0, encoding_params);
    s->y3 = dec_and_decode(pubkey, prvkey, enc_s->y3, 0, encoding_params);
    
    fprintf(stderr, "id: %d received x = %lf\n", id, s->x);
    fprintf(stderr, "id: %d received x2 = %lf\n", id, s->x2);
    fprintf(stderr, "id: %d received x3 = %lf\n", id, s->x3);
    fprintf(stderr, "id: %d received y = %lf\n", id, s->y);
    fprintf(stderr, "id: %d received xy = %lf\n", id, s->xy);
    fprintf(stderr, "id: %d received x2y = %lf\n", id, s->x2y);
    fprintf(stderr, "id: %d received y2 = %lf\n", id, s->y2);
    fprintf(stderr, "id: %d received xy2 = %lf\n", id, s->xy2);
    fprintf(stderr, "id: %d received y3 = %lf\n", id, s->y3);
}

// Get individual broadcast state variable encryption
void get_bcast_state_var(ciphertext_t **state_var, char *enc_str, paillier_serialisation_params_t *serialisation_params){
    MPI_Bcast(enc_str, serialisation_params->paillier_max_enc_serialisation_chars, MPI_CHAR, 0, MPI_COMM_WORLD);
    *state_var = deserialise_encryption(enc_str);
    // fprintf(stderr, "got %s\n", enc_str); // Debugging
}

// Convenience method reducing weighted sum step to a shorter call
void encode_mult_then_add(pubkey_t *pubkey, ciphertext_t *sum, ciphertext_t *weighted, ciphertext_t *to_weight, double weight, encoding_params_t *encoding_params){
    encode_and_mult_enc(pubkey, weighted, to_weight, weight, 0, encoding_params);
    add_encs(pubkey, sum, sum, weighted);
}

// Convenience method reducing sum step to a shorter call
void encode_then_add(pubkey_t *pubkey, ciphertext_t *sum, ciphertext_t *encrypted, double val, encoding_params_t *encoding_params){
    encode_and_enc_no_noise(pubkey, encrypted, val, 1, encoding_params);
    add_encs(pubkey, sum, sum, encrypted);
}

// Serialise and send local encrypted matrix to navigator
void send_enc_mat(c_mtrx_t *mat, int rows, int cols, char *enc_strs, MPI_Request *req, paillier_serialisation_params_t *serialisation_params){
    char *ind;
    for (int r=0; r<rows; r++){
        for (int c=0; c<cols; c++){
            ind = enc_strs+(r*cols*(serialisation_params->paillier_max_enc_serialisation_chars)+c*(serialisation_params->paillier_max_enc_serialisation_chars));
            serialise_encryption(get_c_mtrx(mat, r, c), ind);
            // fprintf(stderr, "len enc (%d %d) of (%d %d) %ld\n", r, c, rows, cols, strlen(ind)); // Debugging
        }
    }
    MPI_Wait(req, MPI_STATUS_IGNORE);
    MPI_Start(req);
}

// Used to check scanf's return what they should
void sensor_input_err_check(int val, int expected, char *msg, int id){
    if (val != expected){
        fprintf(stderr, "%d : %s\n", id, msg);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}
