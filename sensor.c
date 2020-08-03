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


// Local functions
void get_all_bcast_state_vars(enc_state_info_t *s, char *enc_str);
void get_bcast_state_var(ciphertext_t **state_var, char *enc_str);
void encode_mult_then_add(pubkey_t *pubkey, ciphertext_t *sum, ciphertext_t *weighted, ciphertext_t *to_weight, double weight);
void encode_then_add(pubkey_t *pubkey, ciphertext_t *sum, ciphertext_t *encrypted, double val);
void send_enc_mat(c_mtrx_t *mat, int rows, int cols, char *enc_strs, MPI_Request *req);
void sensor_input_err_check(int val, int expected, char *msg, int id);


// The sensor
void run_sensor(int id){
    // Key vars
    char key_str[MAX_KEY_SERIALISATION_CHARS];
    pubkey_t *pubkey;
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
    char enc_str[MAX_ENC_SERIALISATION_CHARS];
    enc_state_info_t enc_state;

    // Result vars
    c_mtrx_t *hrh;
    c_mtrx_t *hrz;
    ciphertext_t *mat_elem;
    ciphertext_t *partial_sum;

    // Sending vars
    char *hrh_enc_strs;
    char *hrz_enc_strs;
    MPI_Request hrh_request;
    MPI_Request hrz_request;
    MPI_Request agg_request;

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

    // Get Paillier public key
    MPI_Bcast(key_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, MPI_COMM_WORLD);
    pubkey = deserialise_pubkey(key_str);
    //printf("id: %d, key: %s\n", id, key_str);

    // Get sensor's private aggregation key
    MPI_Irecv(key_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &agg_request);
    deserialise_aggkey(agg_key, key_str);
    //gmp_printf("id: %d, agg key: %Zd\n", id, agg_key);

    // Open sensor measurements file - TODO currently setup for debugging input only!
    sprintf(f_name, "input/debug_sensor%d.txt", id);
    measurements_fp = fopen(f_name, "r");
    if (measurements_fp == NULL){
        fprintf(stderr, "%d Could not open measurement file!\n", id);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Read number of timesteps
    sensor_input_err_check(fscanf(measurements_fp, "%d", &time_steps), 1, "Could not read timesteps!", id);
    //printf("%d steps=%d\n", id, time_steps);

    // Read navigator state dimension
    sensor_input_err_check(fscanf(measurements_fp, "%d", &state_dim), 1, "Could not read state dimension!", id);
    //printf("%d state_dim=%d\n", id, state_dim);

    // Read sensor location
    sensor_input_err_check(fscanf(measurements_fp, "%lf %lf", &loc_x, &loc_y), 2, "Could not read sensor location!", id);
    //printf("%d loc=(%lf, %lf)\n", id, loc_x, loc_y);

    // 888b     d888                                        d8888 888 888
    // 8888b   d8888                                       d88888 888 888
    // 88888b.d88888                                      d88P888 888 888
    // 888Y88888P888  .d88b.  88888b.d88b.               d88P 888 888 888  .d88b.   .d8888b
    // 888 Y888P 888 d8P  Y8b 888 "888 "88b             d88P  888 888 888 d88""88b d88P"
    // 888  Y8P  888 88888888 888  888  888            d88P   888 888 888 888  888 888
    // 888   "   888 Y8b.     888  888  888 d8b       d8888888888 888 888 Y88..88P Y88b.
    // 888       888  "Y8888  888  888  888 Y8P      d88P     888 888 888  "Y88P"   "Y8888P




    // Allocate space for result vars now that state dimension is known
    hrh = c_mtrx_alloc(state_dim, state_dim);
    hrz = c_mtrx_alloc(1, state_dim);
    init_c_mtrx(hrh);
    init_c_mtrx(hrz);
    partial_sum = init_ciphertext();
    hrh_enc_strs = (char *)malloc(state_dim*state_dim*MAX_ENC_SERIALISATION_CHARS*sizeof(char));
    hrz_enc_strs = (char *)malloc(state_dim*MAX_ENC_SERIALISATION_CHARS*sizeof(char));

    // Initialise repeated sends
    MPI_Send_init(hrh_enc_strs, state_dim*state_dim*MAX_ENC_SERIALISATION_CHARS, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &hrh_request);
    MPI_Send_init(hrz_enc_strs, 1*state_dim*MAX_ENC_SERIALISATION_CHARS, MPI_CHAR, 0, 1, MPI_COMM_WORLD, &hrz_request);
        
    //  .d8888b.  d8b
    // d88P  Y88b Y8P
    // Y88b.
    //  "Y888b.   888 88888b.d88b.
    //     "Y88b. 888 888 "888 "88b
    //       "888 888 888  888  888
    // Y88b  d88P 888 888  888  888
    //  "Y8888P"  888 888  888  888




    // Wait for Aggregation key to arrive if not already, and sync with all processes before beginning simulation
    MPI_Wait(&agg_request, MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);

    // Tracking simulation begins
    for (int t=0; t<time_steps; t++){
        // Get next measurement
        sensor_input_err_check(fscanf(measurements_fp, "%lf", &measurement), 1, "Could not read measurement!", id);

        // Modified measurement for filter, adjusted to be zero-mean
        z = pow(measurement, 2) - SENSOR_VARIANCE;

        // TODO noise approx may be doable better (overestimate for better consistency?)
        // Measreument noise of the adjusted measreument, approximated by using measurement for distance
        inv_R_adj = 1.0/(2*(2*(pow(measurement, 2) - SENSOR_VARIANCE)*SENSOR_VARIANCE + pow(SENSOR_VARIANCE, 2)));

        // Get all state variable encryption broadcasts x,x2,x3,y,xy,x2y,y2,xy2,y3
        get_all_bcast_state_vars(&enc_state, enc_str);

        // 888    888 8888888b.  888    888
        // 888    888 888   Y88b 888    888
        // 888    888 888    888 888    888
        // 8888888888 888   d88P 8888888888
        // 888    888 8888888P"  888    888
        // 888    888 888 T88b   888    888
        // 888    888 888  T88b  888    888
        // 888    888 888   T88b 888    888




        // hrh[0][0] = (4*invRadj)*x2 + (-8*invRadj*sx)*x + (4*invRadj*sx**2)
        //printf("sensor %d h[0][0] = %lf\n", id, (4*inv_R_adj)*pow(0.1289492000, 2) + (-8*inv_R_adj*loc_x)*0.1289492000 + (4*inv_R_adj*pow(loc_x,2)));
        mat_elem = get_c_mtrx(hrh, 0, 0);
        encode_and_mult_enc(pubkey, mat_elem, enc_state.x2, 4*inv_R_adj, 0);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x, -8*inv_R_adj*loc_x);
        encode_then_add(pubkey, mat_elem, partial_sum, 4*inv_R_adj*pow(loc_x, 2));
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrh[0][2] = (4*invRadj)*xy + (-4*invRadj*sy)*x + (-4*invRadj*sx)*y + (4*invRadj*sx*sy)
        mat_elem = get_c_mtrx(hrh, 0, 2);
        encode_and_mult_enc(pubkey, mat_elem, enc_state.xy, 4*inv_R_adj, 0);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x, -4*inv_R_adj*loc_y);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y, -4*inv_R_adj*loc_x);
        encode_then_add(pubkey, mat_elem, partial_sum, 4*inv_R_adj*loc_x*loc_y);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrh[2][0] = (4*invRadj)*xy + (-4*invRadj*sy)*x + (-4*invRadj*sx)*y + (4*invRadj*sx*sy)
        mat_elem = get_c_mtrx(hrh, 2, 0);
        copy_encryption(mat_elem, get_c_mtrx(hrh, 0, 2));
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrh[2][2] = (4*invRadj)*y2 + (-8*invRadj*sy)*y + (4*invRadj*sy**2)
        mat_elem = get_c_mtrx(hrh, 2, 2);
        encode_and_mult_enc(pubkey, mat_elem, enc_state.y2, 4*inv_R_adj, 0);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y, -8*inv_R_adj*loc_y);
        encode_then_add(pubkey, mat_elem, partial_sum, 4*inv_R_adj*pow(loc_y, 2));
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
        mat_elem = get_c_mtrx(hrz, 0, 0);
        encode_and_mult_enc(pubkey, mat_elem, enc_state.x, 2*inv_R_adj*z, 0);
        encode_then_add(pubkey, mat_elem, partial_sum, -2*inv_R_adj*loc_x*z);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x, -2*inv_R_adj*pow(loc_x, 2));
        encode_then_add(pubkey, mat_elem, partial_sum, 2*inv_R_adj*pow(loc_x, 3));
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x, -2*inv_R_adj*pow(loc_y, 2));
        encode_then_add(pubkey, mat_elem, partial_sum, 2*inv_R_adj*loc_x*pow(loc_y, 2));
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x3, 2*inv_R_adj);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x2, -2*inv_R_adj*loc_x);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.xy2, 2*inv_R_adj);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y2, -2*inv_R_adj*loc_x);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrz[2] = (2*invRadj*z)*y + (-2*invRadj*sy*z) + (-2*invRadj*sx**2)*y + (2*invRadj*sy*sx**2) + 
        //          (-2*invRadj*sy**2)*y + (2*invRadj*sy**3) + (2*invRadj)*x2y + (-2*invRadj*sy)*x2 + 
        //          (2*invRadj)*y3 + (-2*invRadj*sy)*y2
        mat_elem = get_c_mtrx(hrz, 0, 2);
        encode_and_mult_enc(pubkey, mat_elem, enc_state.y, 2*inv_R_adj*z, 0);
        encode_then_add(pubkey, mat_elem, partial_sum, -2*inv_R_adj*loc_y*z);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y, -2*inv_R_adj*pow(loc_x, 2));
        encode_then_add(pubkey, mat_elem, partial_sum, 2*inv_R_adj*loc_y*pow(loc_x, 2));
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y, -2*inv_R_adj*pow(loc_y, 2));
        encode_then_add(pubkey, mat_elem, partial_sum, 2*inv_R_adj*pow(loc_y, 3));
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x2y, 2*inv_R_adj);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.x2, -2*inv_R_adj*loc_y);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y3, 2*inv_R_adj);
        encode_mult_then_add(pubkey, mat_elem, partial_sum, enc_state.y2, -2*inv_R_adj*loc_y);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrz remaining elements are 0
        mat_elem = get_c_mtrx(hrz, 0, 1);
        encrypt_zero(pubkey, mat_elem);
        refresh_encryption(pubkey, get_c_mtrx(hrz, 0, 3), mat_elem);

        // Add aggregation noise
        add_agg_noise_c_mtrx(pubkey, agg_key, hrh, t, 0);
        add_agg_noise_c_mtrx(pubkey, agg_key, hrz, t, 1);

        // TODO TEMP checking
        // ciphertext_t *x = init_ciphertext();
        // encode_and_enc(pubkey, x, 5.0, 1);
        // set_c_mtrx(hrh, 0, 0, x);

        // Send encrypted matrix and vector
        send_enc_mat(hrh, state_dim, state_dim, hrh_enc_strs, &hrh_request);
        send_enc_mat(hrz, 1, state_dim, hrz_enc_strs, &hrz_request);
        //printf("\nSensor %d sending : %s\n\n", id, hrh_enc_strs);
    }

    // 888b     d888                                 8888888888
    // 8888b   d8888                                 888
    // 88888b.d88888                                 888
    // 888Y88888P888  .d88b.  88888b.d88b.           8888888 888d888 .d88b.   .d88b.
    // 888 Y888P 888 d8P  Y8b 888 "888 "88b          888     888P"  d8P  Y8b d8P  Y8b
    // 888  Y8P  888 88888888 888  888  888          888     888    88888888 88888888
    // 888   "   888 Y8b.     888  888  888 d8b      888     888    Y8b.     Y8b.
    // 888       888  "Y8888  888  888  888 Y8P      888     888     "Y8888   "Y8888




    // Done with repeated sends, free the reused request structs
    MPI_Request_free(&hrh_request);
    MPI_Request_free(&hrz_request);

    // Free sending vars
    free(hrh_enc_strs);
    free(hrz_enc_strs);

    // Free encrypted matrix and vector vars
    free_ciphertext(partial_sum);
    c_mtrx_free(hrh);
    c_mtrx_free(hrz);

    // Done with measurements, close file and finish
    fclose(measurements_fp);
}

// Get state variable encryption broadcasts x,x2,x3,y,xy,x2y,y2,xy2,y3 in that order
void get_all_bcast_state_vars(enc_state_info_t *s, char *enc_str){
    get_bcast_state_var(&(s->x),   enc_str);
    get_bcast_state_var(&(s->x2),  enc_str);
    get_bcast_state_var(&(s->x3),  enc_str);
    get_bcast_state_var(&(s->y),   enc_str);
    get_bcast_state_var(&(s->xy),  enc_str);
    get_bcast_state_var(&(s->x2y), enc_str);
    get_bcast_state_var(&(s->y2),  enc_str);
    get_bcast_state_var(&(s->xy2), enc_str);
    get_bcast_state_var(&(s->y3),  enc_str);
}

// Get individual broadcast state variable encryption
void get_bcast_state_var(ciphertext_t **state_var, char *enc_str){
    MPI_Bcast(enc_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, MPI_COMM_WORLD);
    *state_var = deserialise_encryption(enc_str);
    //printf("got %s\n", enc_str);
}

// Convenience method reducing weighted sum step to a shorter call
void encode_mult_then_add(pubkey_t *pubkey, ciphertext_t *sum, ciphertext_t *weighted, ciphertext_t *to_weight, double weight){
    encode_and_mult_enc(pubkey, weighted, to_weight, weight, 0);
    add_encs(pubkey, sum, sum, weighted);
}

// Convenience method reducing sum step to a shorter call
void encode_then_add(pubkey_t *pubkey, ciphertext_t *sum, ciphertext_t *encrypted, double val){
    encode_and_enc_no_noise(pubkey, encrypted, val, 1);
    add_encs(pubkey, sum, sum, encrypted);
}

// Serialise and send local encrypted matrix to navigator
void send_enc_mat(c_mtrx_t *mat, int rows, int cols, char *enc_strs, MPI_Request *req){
    char *ind;
    for (int r=0; r<rows; r++){
        for (int c=0; c<cols; c++){
            ind = enc_strs+(r*cols*MAX_ENC_SERIALISATION_CHARS+c*MAX_ENC_SERIALISATION_CHARS);
            serialise_encryption(get_c_mtrx(mat, r, c), ind);
            //printf("len enc (%d %d) of (%d %d) %ld\n", r, c, rows, cols, strlen(ind));
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