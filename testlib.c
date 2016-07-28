/* testlib.c

   Helpers for both test_hist and test_sparse.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "sparse.h"
#include "testlib.h"

bsc_hist_t *
bsc_random(int nrows, int ncols, int nent, int quiet) {
    bsc_hist_t *H;
    int i;

    if(!quiet) {
        printf("Generating a %dx%dx%d bsc histogram...", nrows, ncols, nent);
        fflush(stdout);
    }

    H= bsc_hist_new();

    for(i= 0; i < nent; i++) {
        int c, r;

        c= random() % ncols;
        r= random() % nrows;
        bsc_hist_count(H, c, r, 1);
    }
    if(!quiet) printf(" done.\n");

    return H;
}

dv_t *
dv_random(int length, float mag, int quiet) {
    dv_t *v;
    int i;

    assert(0 <= length);

    if(!quiet) {
        printf("Generating random vector of length %d...", length);
        fflush(stdout);
    }

    v= dv_new(length);
    for(i= 0; i < length; i++)
        v->entries[i]= (mag * random()) / RAND_MAX;
    if(!quiet) printf(" done.\n");

    return v;
}
