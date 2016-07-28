/* row_average.c

   Average each row of the channel matrix, and print it.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "sparse.h"

int
main(int argc, char *argv[]) {
    csc_mat_t *Q;
    FILE *in;
    int quiet;
    csc_errno_t e;
    double *row_avg;
    int c, r;

    if(argc < 2) {
        fprintf(stderr, "Usage: %s <channel_matrix> [-q]\n", argv[0]);
        return 1;
    }

    if(argc > 2 && !strcmp(argv[2], "-q")) quiet= 1;

    if(!quiet) printf("Loading channel matrix...");
    fflush(stdout);
    in= fopen(argv[1], "rb");
    if(!in) { perror("fopen"); return 1; }

    Q= csc_load_binary(in, &e);
    if(!Q) { csc_perror(e, "csc_load_binary"); return 1; }

    fclose(in);
    if(!quiet) printf(" done.\n");

    row_avg= calloc(Q->nrow, sizeof(double));
    if(!row_avg) {
        perror("calloc");
        return 1;
    }

    for(c= 0; c < Q->ncol; c++) {
        int64_t i;

        for(i= Q->ci[c]; i < Q->ci[c+1]; i++) {
            r= Q->rows[i];
            double p= Q->entries[i];

            assert(0 <= r);
            assert(r < Q->nrow);

            row_avg[r]+= p * c;
        }
    }

    for(r= 0; r < Q->nrow; r++)
        printf("%d %.12e\n", r, row_avg[r]);

    free(row_avg);
    csc_mat_destroy(Q);

    return 0;
}
