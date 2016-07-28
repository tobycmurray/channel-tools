/* speed_sparse.c

   Measure speed of sparse matrix operations.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "sparse.h"
#include "testlib.h"

#define ST_DENSE_ROWS 100
#define ST_SPARSE_ROWS 10000
#define ST_COLS 10000
#define ST_VSPARSE_COLS 100000
#define ST_ENTRIES 10000000
#define ST_BLOCKSIZE 4

void
random_build(void) {
    bsc_hist_t *H;
    csc_mat_t *M;
    int *rows_d, *rows_s, *cols, *cols_vs;
    int i;
    struct timespec start, end;
    double interval;

    printf("*** Testing build speed for uniform random data\n\n");

    printf("Generating samples...");
    fflush(stdout);
    rows_d= malloc(ST_ENTRIES * sizeof(int));
    if(!rows_d) abort();
    rows_s= malloc(ST_ENTRIES * sizeof(int));
    if(!rows_s) abort();
    cols= malloc(ST_ENTRIES * sizeof(int));
    if(!cols) abort();
    cols_vs= malloc(ST_ENTRIES * sizeof(int));
    if(!cols_vs) abort();
    for(i= 0; i < ST_ENTRIES; i++) {
        rows_d[i]=  random() % ST_DENSE_ROWS;
        rows_s[i]=  random() % ST_SPARSE_ROWS;
        cols[i]=    random() % ST_COLS;
        cols_vs[i]= random() % ST_VSPARSE_COLS;
    }
    printf(" done.\n");

    printf("Building dense bsc histogram...");
    fflush(stdout);
    if(clock_gettime(CLOCK_REALTIME, &start)) abort();
    H= bsc_hist_new();
    for(i= 0; i < ST_ENTRIES; i++)
        bsc_hist_count(H, cols[i], rows_d[i], 1);
    if(clock_gettime(CLOCK_REALTIME, &end)) abort();
    printf(" done.\n");

    interval= (end.tv_sec + end.tv_nsec * 1e-9) -
              (start.tv_sec + start.tv_nsec * 1e-9);
    printf("Took %.3es, %.3ecounts/s\n", interval, ST_ENTRIES/interval);

    bsc_stats(H);

    printf("Building channel matrix...");
    fflush(stdout);
    if(clock_gettime(CLOCK_REALTIME, &start)) abort();
    M= bsc_normalise(H);
    if(clock_gettime(CLOCK_REALTIME, &end)) abort();
    printf(" done.\n");

    interval= (end.tv_sec + end.tv_nsec * 1e-9) -
              (start.tv_sec + start.tv_nsec * 1e-9);
    printf("Took %.3es, %.3eentries/s\n", interval, M->nnz/interval);

    csc_stats(M);
    printf("\n");

    bsc_hist_destroy(H);
    csc_mat_destroy(M);

    printf("Building sparse bsc histogram...");
    fflush(stdout);
    if(clock_gettime(CLOCK_REALTIME, &start)) abort();
    H= bsc_hist_new();
    for(i= 0; i < ST_ENTRIES; i++)
        bsc_hist_count(H, cols[i], rows_s[i], 1);
    if(clock_gettime(CLOCK_REALTIME, &end)) abort();
    printf(" done.\n");

    interval= (end.tv_sec + end.tv_nsec * 1e-9) -
              (start.tv_sec + start.tv_nsec * 1e-9);
    printf("Took %.3es, %.3ecounts/s\n", interval, ST_ENTRIES/interval);

    bsc_stats(H);

    printf("Building channel matrix...");
    fflush(stdout);
    if(clock_gettime(CLOCK_REALTIME, &start)) abort();
    M= bsc_normalise(H);
    if(clock_gettime(CLOCK_REALTIME, &end)) abort();
    printf(" done.\n");

    interval= (end.tv_sec + end.tv_nsec * 1e-9) -
              (start.tv_sec + start.tv_nsec * 1e-9);
    printf("Took %.3es, %.3eentries/s\n", interval, M->nnz/interval);

    csc_stats(M);
    printf("\n");

    bsc_hist_destroy(H);
    csc_mat_destroy(M);
}

static const int nr_sizes[]= { 1000,     10000, 0};
static const int nc_sizes[]= { 1000,     10000, 0};
static const int ne_sizes[]= {10000, 100000000, 0};

#define RM_REPS 10
#define RM_STRIDE 7

void
random_mult(void) {
    int nr, nc, ne, i, j;
    bsc_hist_t *H;
    csc_mat_t *M;
    dv_t *x, *y;
    struct timespec start, end;
    double interval;

    printf("*** Testing multiplication speed.\n\n");

    for(i= 0; nr_sizes[i] > 0; i++) {
        nr= nr_sizes[i];
        nc= nc_sizes[i];
        ne= ne_sizes[i];

        printf("Using %d rows, %d cols, %d samples.\n", nr, nc, ne);

        H= bsc_random(nr, nc, ne, 0);

        bsc_stats(H);

        M= bsc_normalise(H);
        bsc_hist_destroy(H);

        csc_stats(M);


        x= dv_random(M->nrow, 1.0, 0);
        y= dv_new(M->ncol);
        printf("Measuring non-strided speed... "); fflush(stdout);
        if(clock_gettime(CLOCK_REALTIME, &start)) abort();
        for(j= 0; j < RM_REPS; j++) mult_csc_dv(y, x, M);
        if(clock_gettime(CLOCK_REALTIME, &end)) abort();

        interval= (end.tv_sec + end.tv_nsec * 1e-9) -
                  (start.tv_sec + start.tv_nsec * 1e-9);
        printf("%d iterations took %.3es, %.3fGFLOPs\n\n", RM_REPS,
               interval, RM_REPS*M->nnz/interval*2*1e-9);

        dv_destroy(x); dv_destroy(y);

        printf("Striding matrix at %d...", 1<<RM_STRIDE); fflush(stdout);
        csc_stride(M, RM_STRIDE);
        printf(" done.\n");

        x= dv_random(M->nrow, 1.0, 0);
        y= dv_new(M->ncol);

        printf("Measuring strided speed... "); fflush(stdout);
        if(clock_gettime(CLOCK_REALTIME, &start)) abort();
        for(j= 0; j < RM_REPS; j++) csc_str_mult_nv(y, x, M);
        if(clock_gettime(CLOCK_REALTIME, &end)) abort();

        interval= (end.tv_sec + end.tv_nsec * 1e-9) -
                  (start.tv_sec + start.tv_nsec * 1e-9);
        printf("%d iterations took %.3es, %.3fGFLOPs\n\n", RM_REPS,
               interval, RM_REPS*M->nnz/interval*2*1e-9);

        csc_mat_destroy(M);
    }
}

int
main(int argc, char *argv[]) {
    struct timespec t;
    unsigned int seed;

    if(argc > 1) {
        seed= strtoul(argv[1], NULL, 10);
    }
    else {
        if(clock_gettime(CLOCK_REALTIME, &t)) {
            perror("clock_gettime");
            abort();
        }
        seed= t.tv_sec + t.tv_nsec;
    }
    printf("Seed: %u\n", seed);
    srandom(seed);

    //random_build();

    random_mult();

    return 0;
}
