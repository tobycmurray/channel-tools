/* test_sparse.c

   Test sparse matrix operations.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <assert.h>
#include <malloc.h>
#include <mcheck.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "sparse.h"
#include "testlib.h"

#ifdef NDEBUG
#error Building tests without asserts is meaningless.
#endif

#define RC_NROWS  10000
#define RC_NCOLS  10000
#define RC_NENT   100000
#define RC_STRIDE 2
#define RC_CFSPAN 2

#define PROB_DELTA 1e-6

void
bsc_check_total(bsc_hist_t *H) {
    int c, r, t;
    int *row_totals;

    row_totals= calloc(H->end_row, sizeof(int));
    if(!row_totals) { perror("calloc"); abort(); }

    t= 0;
    for(c= 0; c < H->end_col; c++) {
        for(r= H->start_rows[c]; r < H->end_rows[c]; r++) {
            row_totals[r]+= H->entries[c][r - H->start_rows[c]];
            t+= H->entries[c][r - H->start_rows[c]];
        }
    }

    for(r= 0; r < H->end_row; r++) {
        assert(row_totals[r] == H->row_total[r]);
    }

    assert(t == H->total);
}

void
check_alloc(bsc_hist_t *H) {
    int c, n;

    n= 0;
    for(c= 0; c < H->end_col; c++) {
        if(H->end_rows[c] > H->start_rows[c])
            n+= H->end_rows[c] - H->start_rows[c];
    }

    assert(n == H->nalloc);
}

void
leak_test(void) {
    bsc_hist_t *H;
    int i;
    struct mallinfo mi_start, mi_end;

    mtrace();

    mi_start= mallinfo();
    printf("Testing for leaks in bsc_hist_destroy()...");
    fflush(stdout);
    for(i= 0; i < 10; i++) {
        H= bsc_random(RC_NROWS, RC_NCOLS, RC_NENT, 1);
        bsc_hist_destroy(H);
    }
    printf(" done.\n");
    mi_end= mallinfo();
    assert(mi_start.uordblks == mi_end.uordblks);

    muntrace();
}

void
check_row_prob(csc_mat_t *M) {
    int c, i, r;
    double *rp;
    int *row_nz;
    
    rp= calloc(M->nrow, sizeof(double));
    if(!rp) { perror("calloc"); abort(); }

    row_nz= calloc(M->nrow, sizeof(int));
    if(!row_nz) { perror("calloc"); abort(); }

    if(STRIDE_OF(M) > 1) {
        for(c= 0; c < M->ncol / STRIDE_OF(M); c++) {
            for(i= M->si[c]; i < M->si[c+1]; i++) {
                int row;

                if(M->flags & CSC_F_CFREE)
                    row= M->rows[i/STRIDE_OF(M)] + M->row_offsets[i];
                else
                    row= M->rows[i];

                assert(row >= 0 && row < M->nrow);
                rp[row]+= M->entries[i];
                if(M->entries[i] != 0) row_nz[row]++;
            }
        }
    }
    else {
        for(c= 0; c < M->ncol; c++) {
            for(i= M->ci[c]; i < M->ci[c+1]; i++) {
                rp[M->rows[i]]+= M->entries[i];
                row_nz[M->rows[i]]++;
            }

        }
    }

    for(r= 0; r < M->nrow; r++) {
        assert(row_nz[r] == 0 || abs(1.0 - rp[r]) <= PROB_DELTA);
    }
}

#define RIGGED_COLS 1000

void
rigged_create_test(void) {
    bsc_hist_t *H;
    int c, r;

    printf("Generating a histogram with known counts...");

    H= bsc_hist_new();

    for(c= 0; c < RIGGED_COLS; c++) {
        int i;

        for(i= 0; i < 10; i++) {
            for(r= RIGGED_COLS - c; r < RIGGED_COLS; r++) {
                bsc_hist_count(H, c, r, 1);
            }
        }
    }

    assert(H->end_col= RIGGED_COLS);
    for(c= 0; c < H->end_col; c++) {
        for(r= H->start_rows[c]; r < H->end_rows[c]; r++) {
            int ri= r - H->start_rows[c];
            if(r < RIGGED_COLS - c || r >= RIGGED_COLS) {
                assert(H->entries[c][ri] == 0);
            }
            else {
                assert(H->entries[c][ri] == 10);
            }
        }
    }

    bsc_hist_destroy(H);

    printf(" done.\n");
}

void
random_create_test(void) {
    bsc_hist_t *H;
    csc_mat_t *M;

    H= bsc_random(RC_NROWS, RC_NCOLS, RC_NENT, 0);

    printf("\tChecking totals...");
    fflush(stdout);
    assert(H->total= RC_NENT);
    bsc_check_total(H);
    printf(" done.\n");

    printf("\tChecking alloc counts...");
    fflush(stdout);
    check_alloc(H);
    printf(" done.\n");

    printf("\tGenerating probability matrix...");
    fflush(stdout);
    M= bsc_normalise(H);
    assert(M);
    printf(" done.\n");

    printf("\tChecking probability matrix...");
    fflush(stdout);
    if(!csc_check(M, 1)) abort();
    printf(" done.\n");

    printf("\tChecking row probabilities...");
    fflush(stdout);
    check_row_prob(M);
    printf(" done.\n");

    printf("\tStriding (%d)...", 1<<RC_STRIDE);
    fflush(stdout);
    csc_stride(M, RC_STRIDE);
    printf(" done.\n");

    printf("\tRe-checking probability matrix...");
    fflush(stdout);
    if(!csc_check(M, 1)) abort();
    printf(" done.\n");

    printf("\tRe-checking row probabilities...");
    fflush(stdout);
    check_row_prob(M);
    printf(" done.\n");

    printf("\tMaking matrix collision-free (%d)...", 1<<RC_CFSPAN);
    fflush(stdout);
    csc_make_cfree(M, RC_CFSPAN);
    printf(" done.\n");

    printf("\tRe-checking probability matrix...");
    fflush(stdout);
    if(!csc_check(M, 1)) abort();
    printf(" done.\n");

    printf("\tRe-checking row probabilities...");
    fflush(stdout);
    check_row_prob(M);
    printf(" done.\n");

    csc_mat_destroy(M);
    bsc_hist_destroy(H);
}

void
save_load_test(void) {
    bsc_hist_t *H;
    csc_mat_t *M, *N;
    csc_errno_t e;
    FILE *tmp;
    int i;

    printf("Testing binary save and load of csc matrices.\n");

    H= bsc_random(RC_NROWS, RC_NCOLS, RC_NENT, 0);
    assert(H);

    printf("Generating probability matrix...");
    fflush(stdout);
    M= bsc_normalise(H);
    assert(M);
    printf(" done.\n");

    printf("\tChecking probability matrix...");
    fflush(stdout);
    if(!csc_check(M, 1)) abort();
    printf(" done.\n");

    tmp= tmpfile();
    if(!tmp) { perror("tmpfile"); abort(); }

    printf("\tWriting to disk...");
    fflush(stdout);
    e= csc_store_binary(M, tmp);
    if(e != E_CSC_SUCCESS) {
        csc_perror(e, "csc_store_binary");
        abort();
    }
    printf(" done.\n");

    printf("\tRe-reading...");
    rewind(tmp);
    N= csc_load_binary(tmp, &e);
    if(!N) {
        csc_perror(e, "csc_load_binary");
        abort();
    }
    printf(" done.\n");

    fclose(tmp);

    printf("\tVerifying...");
    assert(M->nrow == N->nrow);
    assert(M->ncol == N->ncol);
    assert(M->nnz == N->nnz);
    for(i= 0; i < M->ncol + 1; i++)
        assert(M->ci[i] == N->ci[i]);

    for(i= 0; i < M->nnz; i++) {
        assert(M->rows[i] == N->rows[i]);
        assert(M->entries[i] == N->entries[i]);
    }
    printf(" done.\n");

    csc_mat_destroy(M);
    bsc_hist_destroy(H);
}

void
test_mult(void) {
    bsc_hist_t *H;
    csc_mat_t *M;
    dv_t *x, *x2, *y;
    int i;

    printf("Testing matrix-vector multiplication...");
    fflush(stdout);
    x= dv_new(RC_NCOLS);
    x2= dv_new(RC_NCOLS);
    y= dv_new(RC_NROWS);
    if(!x || !x2 || !y) { perror("dv_new"); abort(); }
    dv_uniform(y, 1.0);

    H= bsc_random(RC_NROWS, RC_NCOLS, RC_NENT, 1);
    M= bsc_normalise(H);
    bsc_hist_destroy(H);

    if(!csc_check(M, 1)) abort();

    mult_csc_dv(x, y, M);

    for(i= 0; i < RC_NCOLS; i++) {
        int j;
        float s= 0.0;

        for(j= M->ci[i]; j < M->ci[i+1]; j++)
            s+= M->entries[j];

        assert(abs(s - x->entries[i]) - PROB_DELTA);
    }
    printf(" done.\n");

#if 0
    printf("Testing strided (%d) matrix-vector multiplication...",
           RC_STRIDE);
    fflush(stdout);
    csc_stride(M, RC_STRIDE);
    if(!csc_check(M, 1)) abort();
    csc_str_mult_nv(x2, y, M);

    for(i= 0; i < x->length; i++)
        assert(x->entries[i] == x2->entries[i]);
    printf(" done.\n");

    printf("Testing strided (%d) collision-free multiplication...",
           RC_STRIDE);
    fflush(stdout);
    csc_make_cfree(M, RC_CFSPAN);
    if(!csc_check(M, 1)) abort();
    csc_mult_cf(x2, y, M);

    for(i= 0; i < x->length; i++)
        assert(x->entries[i] == x2->entries[i]);
    printf(" done.\n");
#endif

    csc_mat_destroy(M);
    dv_destroy(x);
    dv_destroy(x2);
    dv_destroy(y);
}

void
test_prune(void) {
    bsc_hist_t *H;
    csc_mat_t *M;

    printf("Testing column pruning on random matrix:\n");
    fflush(stdout);
    /* Setting nent=ncol almost guarantees an empty column. */
    H= bsc_random(RC_NROWS, RC_NCOLS, RC_NCOLS, 1);
    M= bsc_normalise(H);
    bsc_hist_destroy(H);
    printf("\tChecking probability matrix...");
    fflush(stdout);
    if(!csc_check(M, 1)) abort();
    printf(" done.\n");
    printf("\tHad %d columns\n", M->ncol);
    csc_prune_cols(M);
    printf("\tNow have %d columns\n", M->ncol);
    check_row_prob(M);
    printf("\tRow sums still match\n");
    csc_mat_destroy(M);
    printf("done.\n");
}

int
main(int argc, char *argv[]) {
    struct timespec t;
    unsigned int seed;
    int i;

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

    leak_test();

    rigged_create_test();

    for(i= 0; i < 4; i++) random_create_test();

    save_load_test();

    test_mult();

    test_prune();

    return 0;
}
