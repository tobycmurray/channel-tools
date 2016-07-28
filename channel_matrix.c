/* channel_matrix.c

   Build a channel matrix from the given samples.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sparse.h"

#define MAXLINE 1024

int
main(int argc, char *argv[]) {
    bsc_hist_t *H;
    csc_mat_t *M;
    csc_errno_t e;
    int c, r;
    FILE *out;
    int rmin= INT_MIN, rmax= INT_MAX;
    int climit= -1;
    int *counts= NULL;
    int discard= 0;
    size_t malformed= 0, out_of_range= 0;
    char buf[MAXLINE];

    if(argc != 2 && argc != 4 && argc != 5 && argc != 6) {
        printf("Usage: %s <output_filename> [<row. min> <row. max> "
               "[<count limit> [<discard>]]]\n",
                argv[0]);
        return 1;
    }

    out= fopen(argv[1], "wb");
    if(!out) {
        perror("fopen");
        return 1;
    }

    if(argc >= 4) {
        rmin= atoi(argv[2]);
        rmax= atoi(argv[3]);
    }

    if(argc >= 5) {
        climit= atoi(argv[4]);
        counts= calloc(rmax - rmin + 1, sizeof(int));
        if(!counts) {
            perror("calloc");
            exit(EXIT_FAILURE);
        }
    }

    if(argc >= 6) {
        discard= atoi(argv[5]);
    }

    printf("Building histogram...");
    fflush(stdout);
    H= bsc_hist_new();
    while(fgets(buf, MAXLINE, stdin)) {
        int n;

        n= sscanf(buf, "%d %d\n", &r, &c);

        if(n != 2) {
            malformed++;
            continue;
        }

        if(r < rmin || rmax < r) {
            out_of_range++;
            continue;
        }

        if(climit >= 0) {
            int i= r - rmin;

            if(counts[i] >= climit + discard)
                continue;

            counts[i]++;

            if(counts[i] <= discard)
                continue;
        }

        bsc_hist_count(H, c, r, 1);
    }
    printf(" done.\n");
    bsc_stats(H);

    if(!bsc_check(H, 1)) abort();

    printf("Building matrix\n");
    M= bsc_normalise(H);
    csc_stats(M);
    bsc_hist_destroy(H);

    if(!csc_check(M, 1)) abort();

    /* Check the row totals. */
    {
        float *row_tot= malloc(M->nrow * sizeof(float));
        if(!row_tot) { perror("malloc"); abort(); }
        bzero(row_tot, M->nrow * sizeof(float));

        for(c= 0; c < M->ncol; c++) {
            int i;

            for(i= M->ci[c]; i < M->ci[c+1]; i++) {
                row_tot[M->rows[i]]+= M->entries[i];
            }
        }

        for(r= 0; r < M->nrow; r++) {
            if(fabs(row_tot[r] - 1.0) > 1e-3 && row_tot[r] != 0.0) {
                fprintf(stderr, "Row %d sums to %.12e\n", r, row_tot[r]);
                abort();
            }
        }
        free(row_tot);
    }

    printf("Writing matrix\n");
    e= csc_store_binary(M, out);
    if(e != E_CSC_SUCCESS) {
        csc_perror(e, "csc_write_binary");
        return 1;
    }

    fclose(out);

    printf("%lu malformed entries, %lu columns out of range\n",
            malformed, out_of_range);

    return 0;
}
