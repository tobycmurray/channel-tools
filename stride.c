/* stride.c

   Stride the given matrix.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <stdio.h>
#include <stdlib.h>

#include "sparse.h"

int
main(int argc, char *argv[]) {
    FILE *in, *out;
    csc_mat_t *M;
    csc_errno_t e;
    int stride;
    int span= -1;

    if(argc < 4) {
        fprintf(stderr, "Usage: %s <stride> <input> <output> [<span>]\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }

    stride= atoi(argv[1]);
    if(stride < 0 || stride > CSC_MAX_STRIDE) {
        fprintf(stderr, "Valid stride range: 0-%d\n", CSC_MAX_STRIDE);
        exit(EXIT_FAILURE);
    }

    in= fopen(argv[2], "rb");
    if(!in) { perror("fopen"); exit(EXIT_FAILURE); }

    out= fopen(argv[3], "wb");
    if(!out) { perror("fopen"); exit(EXIT_FAILURE); }

    if(argc > 4) {
        span= atoi(argv[4]);

        if(span < 0 || span > CSC_MAX_SPAN) {
            fprintf(stderr, "Valid span range: 0-%d\n", CSC_MAX_SPAN);
            exit(EXIT_FAILURE);
        }
    }

    printf("Reading %s...", argv[2]);
    fflush(stdout);
    M= csc_load_binary(in, &e);
    if(!M) { csc_perror(e, "csc_load_binary"); exit(EXIT_FAILURE); }
    fclose(in);
    printf(" done.\n");

    if(!csc_check(M, 1)) abort();
    csc_stats(M);

    if(STRIDE_OF(M) > 1) {
        fprintf(stderr, "Input already strided @ %d\n", STRIDE_OF(M));
        exit(EXIT_FAILURE);
    }

    printf("Striding to %d...", 1 << stride);
    fflush(stdout);
    csc_stride(M, stride);
    if(!csc_check(M, 1)) abort();
    printf(" done.\n");

    if(span != -1) {
        printf("Making collision-free with row span %d...", 1 << span);
        fflush(stdout);
        csc_make_cfree(M, span);
        if(!csc_check(M, 1)) abort();
        printf(" done.\n");
    }

    csc_stats(M);

    printf("Writing %s...", argv[3]);
    e= csc_store_binary(M, out);
    if(e != E_CSC_SUCCESS) {
        fflush(stdout);
        csc_perror(e, "csc_store_binary");
        exit(EXIT_FAILURE);
    }
    printf(" done.\n");

    csc_mat_destroy(M);

    return 0;
}
