/*
 *
 */

#include "key_setup.h"

// Private declarations
void init_rand_agg(gmp_randstate_t rand, paillier_get_rand_t get_rand, int bytes);



// Wrapped for the Paillier public and private key generation
void gen_phe_keys(int paillier_bitsize, paillier_pubkey_t **pubkey, paillier_prvkey_t **prvkey){
    paillier_keygen(paillier_bitsize, pubkey, prvkey, paillier_get_rand_devrandom);
}


// Broadcasting of the Paillier public key to all sensors (starting from process index 1)
void dist_phe_key(int num_sensors, paillier_pubkey_t *pubkey){
    // Serialize key using the Paillier library function provided
    // (note character buffer len +1 for null termination)
    char *key = paillier_pubkey_to_hex(pubkey);
    int key_len = strlen(key)+1;

    // Send key to all
    MPI_Bcast(key, key_len, MPI_CHAR, 0, MPI_COMM_WORLD);

    // Free character array allocated in paillier_pubkey_to_hex
    free(key);
    key = NULL;
}


// Generating and distributing aggregation private keys (starting from process index 1)
void gen_dist_agg_keys(int num_sensors, paillier_pubkey_t *pubkey){
    // Key variables
    char *key_str;
    mpz_t k;
    mpz_t k_last;

    // Init and allocate memeory for large integers
    mpz_init(k);
    mpz_init(k_last);

    // Required for large random number generation. Copied from paillier library
    gmp_randstate_t rnd;
    init_rand_agg(rnd, paillier_get_rand_devrandom, (pubkey->bits)*2 / 8 + 1);

    // Generate all but one key as random values less than the Paillier modulus
    for (int s=1; s<num_sensors; s++){
        do
            mpz_urandomb(k, rnd, (pubkey->bits)*2);
        while(mpz_cmp(k, pubkey->n_squared) >= 0);

        gmp_printf("\n%Zd\n", k);

        // Keep track of sum of keys
        mpz_add(k_last, k_last, k);

        // Convert key to string and send it to respective sensor process
        // (note character buffer len +1 for null termination)
        key_str = mpz_get_str(NULL, SERIALISATION_BASE, k);
        MPI_Send(key_str, strlen(key_str)+1, MPI_CHAR, s, 0, MPI_COMM_WORLD);

        // Free string allocated in mpz_get_str
        free(key_str);
        key_str = NULL;
    }

    // Last key computed and negation of sum of all other keys mod the Paillier modulus
    // i.e. such that the sum of all aggregation keys is (0 % N^2)
    mpz_neg(k_last, k_last);
    mpz_mod(k_last, k_last, pubkey->n_squared);

    gmp_printf("\n%Zd\n", k_last);

    // Convert key to string and send it to the final sensor process
    // (note character buffer len +1 for null termination)
    key_str = mpz_get_str(NULL, SERIALISATION_BASE, k_last);
    MPI_Send(key_str, strlen(key_str)+1, MPI_CHAR, num_sensors, 0, MPI_COMM_WORLD);

    // Free string allocated in mpz_get_str
    free(key_str);
    key_str = NULL;
}



// Copied from libpaillier-X.X library for getting random mpz values for aggregation
void init_rand_agg(gmp_randstate_t rnd, paillier_get_rand_t get_rand, int bytes){
	void* buf;
	mpz_t s;

	buf = malloc(bytes);
	get_rand(buf, bytes);

	gmp_randinit_default(rnd);
	mpz_init(s);
	mpz_import(s, bytes, 1, 1, 0, 0, buf);
	gmp_randseed(rnd, s);
	mpz_clear(s);

	free(buf);
}
