/*
 * 
 */

#include "sensor.h"


// Local functions
void get_bcast_state_var(ciphertext_t **state_var, char *enc_str);
void send_enc_mat(c_mtrx_t *mat, int rows, int cols, char *enc_strs, int send_tag);
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
    ciphertext_t *x;
    ciphertext_t *x2;
    ciphertext_t *x3;
    ciphertext_t *y;
    ciphertext_t *xy;
    ciphertext_t *x2y;
    ciphertext_t *y2;
    ciphertext_t *xy2;
    ciphertext_t *y3;

    // Result vars
    c_mtrx_t *hrh = c_mtrx_alloc(state_dim, state_dim);
    c_mtrx_t *hrz = c_mtrx_alloc(1, state_dim);
    ciphertext_t *mat_elem;
    ciphertext_t *partial_sum;
    char *hrh_enc_strs = (char *)malloc(state_dim*state_dim*MAX_ENC_SERIALISATION_CHARS*sizeof(char));
    char *hrz_enc_strs = (char *)malloc(state_dim*MAX_ENC_SERIALISATION_CHARS*sizeof(char));

    // Create ciphertext objects in allocated space
    init_c_mtrx(hrh);
    init_c_mtrx(hrz);
    partial_sum = init_ciphertext();

    // Measreument vars
    double measurement;
    double z;
    double inv_R_adj;

    time_steps = 2; // TODO TEMP
    for (int t=0; t<time_steps; t++){
        // Get state variable encryption broadcasts x,x2,x3,y,xy,x2y,y2,xy2,y3 in that order
        get_bcast_state_var(&x,   enc_str);
        get_bcast_state_var(&x2,  enc_str);
        get_bcast_state_var(&x3,  enc_str);
        get_bcast_state_var(&y,   enc_str);
        get_bcast_state_var(&xy,  enc_str);
        get_bcast_state_var(&x2y, enc_str);
        get_bcast_state_var(&y2,  enc_str);
        get_bcast_state_var(&xy2, enc_str);
        get_bcast_state_var(&y3,  enc_str);

        // Get next measurement
        sensor_input_err_check(fscanf(measurements_fp, "%lf", &measurement), 1, "Could not read measurement!", id);
        //printf("%d m=%lf\n", id, measurement);

        // Modified measurement for filter, adjusted to be zero-mean
        z = pow(measurement, 2) - SENSOR_VARIANCE;

        // TODO noise approx may be doable better (overestimate for better consistency?)
        // Measreument noise of the adjusted measreument, approximated by using measurement for distance
        inv_R_adj = 1.0/(2*(2*(pow(measurement, 2) - SENSOR_VARIANCE)*SENSOR_VARIANCE + pow(SENSOR_VARIANCE, 2)));

        // Compute HRH weighted sums

        // hrh[0][0] = (4*invRadj)*x2 + (-8*invRadj*sx)*x + (4*invRadj*sx**2)
        mat_elem = get_c_mtrx(hrh, 0, 0);
        encode_and_mult_enc(pubkey, mat_elem, x2, 4*inv_R_adj, 0);
        encode_and_mult_enc(pubkey, partial_sum, x, -8*inv_R_adj*loc_x, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_enc(pubkey, partial_sum, 4*inv_R_adj*pow(loc_x, 2), 1);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrh[0][2] = (4*invRadj)*xy + (-4*invRadj*sy)*x + (-4*invRadj*sx)*y + (4*invRadj*sx*sy)
        mat_elem = get_c_mtrx(hrh, 0, 2);
        encode_and_mult_enc(pubkey, mat_elem, xy, 4*inv_R_adj, 0);
        encode_and_mult_enc(pubkey, partial_sum, x, -4*inv_R_adj*loc_y, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, partial_sum, y, -4*inv_R_adj*loc_x, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_enc(pubkey, partial_sum, 4*inv_R_adj*loc_x*loc_y, 1);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrh[2][0] = (4*invRadj)*xy + (-4*invRadj*sy)*x + (-4*invRadj*sx)*y + (4*invRadj*sx*sy)
        mat_elem = get_c_mtrx(hrh, 2, 0);
        copy_encryption(mat_elem, get_c_mtrx(hrh, 0, 2));
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrh[2][2] = (4*invRadj)*y2 + (-8*invRadj*sy)*y + (4*invRadj*sy**2)
        mat_elem = get_c_mtrx(hrh, 2, 2);
        encode_and_mult_enc(pubkey, mat_elem, y2, 4*inv_R_adj, 0);
        encode_and_mult_enc(pubkey, partial_sum, y, -8*inv_R_adj*loc_y, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_enc(pubkey, partial_sum, 4*inv_R_adj*pow(loc_y, 2), 1);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrh remaining elements are 0
        encrypt_zero(pubkey, get_c_mtrx(hrh, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrh, 0, 3), get_c_mtrx(hrh, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrh, 2, 1), get_c_mtrx(hrh, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrh, 2, 3), get_c_mtrx(hrh, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrh, 1, 0), get_c_mtrx(hrh, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrh, 1, 1), get_c_mtrx(hrh, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrh, 1, 2), get_c_mtrx(hrh, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrh, 1, 3), get_c_mtrx(hrh, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrh, 3, 0), get_c_mtrx(hrh, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrh, 3, 1), get_c_mtrx(hrh, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrh, 3, 2), get_c_mtrx(hrh, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrh, 3, 3), get_c_mtrx(hrh, 0, 1));

        // Compute HR(z-h(x)+Hx) weighted sums

        // hrz[0] = (2*invRadj*z)*x + (-2*invRadj*sx*z) + (-2*invRadj*sx**2)*x + (2*invRadj*sx**3) + 
        //          (-2*invRadj*sy**2)*x + (2*invRadj*sx*sy**2) + (2*invRadj)*x3 + (-2*invRadj*sx)*x2 + 
        //          (2*invRadj)*xy2 + (-2*invRadj*sx)*y2
        mat_elem = get_c_mtrx(hrz, 0, 0);
        encode_and_mult_enc(pubkey, mat_elem, x, 2*inv_R_adj*z, 0);
        encode_and_enc(pubkey, partial_sum, -2*inv_R_adj*loc_x*z, 1);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, partial_sum, x, -2*inv_R_adj*pow(loc_x, 2), 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_enc(pubkey, partial_sum, 2*inv_R_adj*pow(loc_x, 3), 1);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, partial_sum, x, -2*inv_R_adj*pow(loc_y, 2), 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_enc(pubkey, partial_sum, 2*inv_R_adj*loc_x*pow(loc_y, 2), 1);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, partial_sum, x3, 2*inv_R_adj, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, partial_sum, x2, -2*inv_R_adj*loc_x, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, partial_sum, xy2, 2*inv_R_adj, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, partial_sum, y2, -2*inv_R_adj*loc_x, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrz[2] = (2*invRadj*z)*y + (-2*invRadj*sy*z) + (-2*invRadj*sx**2)*y + (2*invRadj*sy*sx**2) + 
        //          (-2*invRadj*sy**2)*y + (2*invRadj*sy**3) + (2*invRadj)*x2y + (-2*invRadj*sy)*x2 + 
        //          (2*invRadj)*y3 + (-2*invRadj*sy)*y2
        mat_elem = get_c_mtrx(hrz, 0, 2);
        encode_and_mult_enc(pubkey, mat_elem, y, 2*inv_R_adj*z, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_enc(pubkey, partial_sum, -2*inv_R_adj*loc_y*z, 1);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, mat_elem, y, -2*inv_R_adj*pow(loc_x, 2), 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_enc(pubkey, partial_sum, 2*inv_R_adj*loc_y*pow(loc_x, 2), 1);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, mat_elem, y, -2*inv_R_adj*pow(loc_y, 2), 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_enc(pubkey, partial_sum, 2*inv_R_adj*pow(loc_y, 3), 1);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, mat_elem, x2y, 2*inv_R_adj, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, mat_elem, x2, -2*inv_R_adj*loc_y, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, mat_elem, y3, 2*inv_R_adj, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        encode_and_mult_enc(pubkey, mat_elem, y2, -2*inv_R_adj*loc_y, 0);
        add_encs(pubkey, mat_elem, mat_elem, partial_sum);
        refresh_encryption(pubkey, mat_elem, mat_elem);

        // hrz remaining elements are 0
        encrypt_zero(pubkey, get_c_mtrx(hrz, 0, 1));
        refresh_encryption(pubkey, get_c_mtrx(hrz, 0, 3), get_c_mtrx(hrz, 0, 1));

        // Add aggregation noise

        // TODO TEMP checking
        ciphertext_t *x = init_ciphertext();
        encode_and_enc(pubkey, x, 5.0, 1);
        set_c_mtrx(hrh, 0, 0, x);

        // Send encrypted matrix and vector
        send_enc_mat(hrh, state_dim, state_dim, hrh_enc_strs, 0);
        send_enc_mat(hrz, 1, state_dim, hrz_enc_strs, 1);
    }


    // Done with measurements, close file and return
    fclose(measurements_fp);
}


void get_bcast_state_var(ciphertext_t **state_var, char *enc_str){
    MPI_Bcast(enc_str, MAX_KEY_SERIALISATION_CHARS, MPI_CHAR, 0, MPI_COMM_WORLD);
    *state_var = deserialise_encryption(enc_str);
}

void send_enc_mat(c_mtrx_t *mat, int rows, int cols, char *enc_strs, int send_tag){
    char *ind;
    for (int r=0; r<rows; r++){
        for (int c=0; c<cols; c++){
            ind = enc_strs+(r*cols*MAX_ENC_SERIALISATION_CHARS+c*MAX_ENC_SERIALISATION_CHARS);
            serialise_encryption(get_c_mtrx(mat, r, c), ind);
            printf("len enc (%d %d) of (%d %d) %ld\n", r, c, rows, cols, strlen(ind));
        }
    }
    MPI_Send(enc_strs, rows*cols*MAX_ENC_SERIALISATION_CHARS, MPI_CHAR, 0, send_tag, MPI_COMM_WORLD);
}

void sensor_input_err_check(int val, int expected, char *msg, int id){
    if (val != expected){
        fprintf(stderr, "%d : %s\n", id, msg);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}