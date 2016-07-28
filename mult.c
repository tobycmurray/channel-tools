/* mult.c

   Measure matrix multiplication speed.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "sparse.h"

#define REPS 10

int
main(int argc, char *argv[]) {
    FILE *in;
    csc_mat_t *M;
    csc_errno_t e;
    dv_t *x, *y;
    struct timespec start, end;
    double iv;
    int i;

    if(argc < 2) {
        fprintf(stderr, "Usage: %s <matrix_filename>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    in= fopen(argv[1], "rb");
    if(!in) { perror("fopen"); exit(EXIT_FAILURE); }

    M= csc_load_binary(in, &e);
    if(!M) { csc_perror(e, "csc_load_binary"); exit(EXIT_FAILURE); }
    fclose(in);

    x= dv_new(M->nrow);
    if(!x) { perror("dv_new"); exit(EXIT_FAILURE); }

    y= dv_new(M->ncol);
    if(!y) { perror("dv_new"); exit(EXIT_FAILURE); }

    if(!csc_check(M, 1)) abort();
    csc_stats(M);

    dv_uniform(x, 1.0);

    if(STRIDE_OF(M) > 1) {
        if(M->flags & CSC_F_CFREE) {
            clock_gettime(CLOCK_REALTIME, &start);
            for(i= 0; i < REPS; i++) csc_mult_cf(y, x, M);
            clock_gettime(CLOCK_REALTIME, &end);
        }
        else {
            clock_gettime(CLOCK_REALTIME, &start);
            for(i= 0; i < REPS; i++) csc_str_mult_nv(y, x, M);
            clock_gettime(CLOCK_REALTIME, &end);
        }
    }
    else {
        clock_gettime(CLOCK_REALTIME, &start);
        for(i= 0; i < REPS; i++) mult_csc_dv(y, x, M);
        clock_gettime(CLOCK_REALTIME, &end);
    }

    iv= end.tv_sec   + end.tv_nsec*1e-9
      - start.tv_sec - start.tv_nsec*1e-9;
    printf("%d reps %lld entries %.2esec %.2fGFLOPs\n", REPS,
           (long long int)M->nnz, iv, (REPS * 2 * M->nnz) / iv / 1e9);

    dv_destroy(y);
    dv_destroy(x);
    csc_mat_destroy(M);

    return 0;
}
