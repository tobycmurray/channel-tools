/* analyse.c

   Read samples and display histogram statistics.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "sparse.h"

int
main(void) {
    bsc_hist_t *H;
    int c, r;
    int min_start, max_start, min_end, max_end, min_range, max_range;
    int nempty;
    uint64_t sum_ranges;

    H= bsc_hist_new();
    if(!H) { perror("bsc_hist_new"); abort(); }
    while(scanf("%d %d\n", &c, &r) == 2) bsc_hist_count(H, c, r, 1);
    fprintf(stderr, "  %d columns, %llu entries, %.2fMb\n",
                    H->end_col, (long long unsigned int)H->nnz,
                    bsc_size(H) / (1024.0 * 1024));

    nempty= 0;
    min_start= INT_MAX; max_start= -1;
    min_end=   INT_MAX; max_end=   -1;
    min_range= INT_MAX; max_range= -1;
    sum_ranges= 0;
    for(c= 0; c < H->end_col; c++) {
        int range;

        if(H->end_rows[c] <= H->start_rows[c]) {
            nempty++;
            continue;
        }

        if(min_start > H->start_rows[c]) min_start= H->start_rows[c];
        if(H->start_rows[c] > max_start) max_start= H->start_rows[c];
        if(H->end_rows[c] > max_end) max_end= H->end_rows[c];
        if(min_end > H->end_rows[c]) min_end= H->end_rows[c];

        range= H->end_rows[c] - H->start_rows[c];
        if(range > max_range) max_range= range;
        if(min_range > range) min_range= range;

        sum_ranges+= range;
    }
    printf("Range of start rows is (%d-%d)\n", min_start, max_start);
    printf("Range of end rows is (%d-%d)\n", min_end, max_end);
    printf("Range of row spans is (%d-%d)\n", min_range, max_range);
    printf("%d empty columns (%.2f%% filled)\n", nempty,
           (100.0 * (H->end_col - nempty)) / H->end_col);
    {
        uint64_t ncol= H->end_col, nfull= H->end_col - nempty;
        uint64_t est_sum= sum_ranges * ncol / nfull;
        printf("Sum of ranges is %llu (est. %llu or %.2fMb if complete)\n",
               (unsigned long long)sum_ranges,
               (unsigned long long)est_sum,
               (double)est_sum * sizeof(int) / (1024 * 1024));
    }

    return 0;
}
