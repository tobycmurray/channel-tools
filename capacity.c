/* capacity.c

   Calculate Shannon capacity of the given channel matrix.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "channel_algorithms.h"
#include "sparse.h"
#include "log.h"

int
main(int argc, char *argv[]) {
    csc_mat_t *Q;
    csc_errno_t e;
    float c, epsilon, e_obs;
    FILE *in;
    int quiet= 0;
#ifdef CAP_BENCH
    struct timespec start, end;
#endif

    if(argc < 3) {
        fprintf(stderr, "Usage: %s <channel_matrix> <precision> [-q]\n",
                argv[0]);
        return 1;
    }

    epsilon= strtof(argv[2], NULL);

    if(argc > 3 && !strcmp(argv[3], "-q")) quiet= 1;

    if(!quiet) printf("Loading channel matrix...");
    fflush(stdout);
    in= fopen(argv[1], "rb");
    if(!in) { perror("fopen"); return 1; }

    Q= csc_load_binary(in, &e);
    if(!Q) { csc_perror(e, "csc_load_binary"); return 1; }

    fclose(in);
    if(!quiet) printf(" done.\n");

    if(!quiet) printf("Writing log table...");
    fflush(stdout);
    write_log_table();
    if(!quiet) printf(" done.\n");

    csc_prune_cols(Q);

    if(!quiet) printf("Finding capacity with target precision %.3e...", epsilon);
    fflush(stdout);
#ifdef CAP_BENCH
    if(clock_gettime(CLOCK_REALTIME, &start)) {
        perror("clock_gettime"); abort();
    }
#endif
    c= ba_phased(Q, epsilon, &e_obs);
#ifdef CAP_BENCH
    if(clock_gettime(CLOCK_REALTIME, &end)) {
        perror("clock_gettime"); abort();
    }
#endif
    if(!quiet) printf(" done.\n");

    if(!quiet)
        printf("Channel capacity is %e(+%e,-0) bits per usage.\n", c, e_obs);
    else
        printf("%.12e %.12e\n", c, c + e_obs);

#ifdef CAP_BENCH
    {
        double iv= end.tv_sec   + end.tv_nsec*1e-9
                 - start.tv_sec - start.tv_nsec*1e-9;

        printf("Converged in %.2fs and %d iterations: %.2e/s\n",
               iv, iterations, iterations / iv);
    }
#endif

    return 0;
}
