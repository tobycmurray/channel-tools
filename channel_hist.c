/* channel_hist.c

   Calculate and display the histogram of the input samples.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "sparse.h"

#define MAXLINE 1024

int
main(int argc, char *argv[]) {
    bsc_hist_t *H;
    int c, r;
    FILE *out;
    int ne_cols= 0;
    int min_row= INT_MAX;
    int max_row= 0;
    int cmin= INT_MIN, cmax= INT_MAX;
    int climit= -1;
    int *counts;
    size_t malformed= 0, out_of_range= 0;
    int discard= 0;
    char buf[MAXLINE];

    if(argc != 2 && argc != 4 && argc != 5 && argc != 6) {
        printf("Usage: %s <output_filename> [<col. min> <col. max> "
               "[<count limit> [<discard>]]]\n",
                argv[0]);
        return 1;
    }

    out= fopen(argv[1], "w");
    if(!out) {
        perror("fopen");
        return 1;
    }

    if(argc >= 4) {
        cmin= atoi(argv[2]);
        cmax= atoi(argv[3]);
    }

    if(argc >= 5) {
        climit= atoi(argv[4]);
        counts= calloc(cmax - cmin + 1, sizeof(int));
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

        n= sscanf(buf, "%d %d\n", &c, &r);

        if(n != 2) {
            malformed++;
            continue;
        }

        if(c < cmin || cmax < c) {
            out_of_range++;
            continue;
        }

        if(climit >= 0) {
            int i= c - cmin;

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

    for(c= 0; c < H->end_col; c++) {
        if(H->start_rows[c] < H->end_rows[c]) ne_cols++;
        if(H->start_rows[c] < H->end_rows[c]) {
            if(H->start_rows[c] < min_row)
                min_row= H->start_rows[c];
            if(max_row < H->end_rows[c])
                max_row= H->end_rows[c];
        }
    }
    assert(min_row <= max_row);

    fprintf(out, "%d %d %d\n", ne_cols, min_row, max_row);

    for(c= 0; c < H->end_col; c++) {
        for(r= H->start_rows[c]; r < H->end_rows[c]; r++) {
            int ri= r - H->start_rows[c];
            if(H->entries[c][ri] > 0) {
                fprintf(out, "%d %d %d\n", c, r, H->entries[c][ri]);
            }
        }
    }

    fclose(out);

    fprintf(stderr, "%lu malformed entries, %lu columns out of range\n",
            malformed, out_of_range);

    return 0;
}

