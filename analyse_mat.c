/* analyse_mat.c

   Calculate the effect of different strides on the matrix.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "sparse.h"

int
stride_init(csc_mat_t *M, int stride, int c_base, int *i) {
    int s;
    int k= -1, r_min= INT_MAX; /* Sentinels. */

    /* Find the least starting row. */
    for(s= 0; s < stride; s++) {
        int c= c_base + s;

        /* Start the column off. */
        i[s]= M->ci[c];

        /* Only considering non-empty columns. */
        if(M->ci[c] < M->ci[c+1]) {
            /* Assume row values don't interfere with the sentinel. */
            assert(M->rows[M->ci[c]] < INT_MAX);

            /* Find the smallest starting row. */
            if(M->rows[M->ci[c]] < r_min) {
                k= s;
                r_min= M->rows[i[s]+1];
            }
        }
    }

    if(k != -1) return i[k];
    else        return -1; /* It's possible that all columns were empty. */
}

/* Iterate along 'stride' merged columns.  Returns -1 on completion. */
int
stride_step(csc_mat_t *M, int stride, int c_base, int *i, int j) {
    int s;
    int k= -1, r_min= INT_MAX; /* Sentinels. */

    for(s= 0; s < stride; s++) {
        int c= c_base + s;
        /* Only consider columns not yet exhausted. */
        if(i[s] < M->ci[c+1]) {
            /* j should have been the previous minimum. */
            assert(M->rows[i[s]+1] >= M->rows[j]);

            /* Assume row values don't interfere with the sentinel. */
            assert(M->rows[i[s]+1] < INT_MAX);

            /* Find the smallest row we could advance to. */
            if(M->rows[i[s]+1] < r_min) {
                k= s;
                r_min= M->rows[i[s]+1];
            }
        }
    }

    /* If we found a column to advance, do it. */
    if(k != -1) {
        assert(i[k] < M->ci[c_base + k + 1]);
        i[k]++;
        return i[k];
    }
    else {
        return -1; /* Stop iteration. */
    }
}

void
stride_sep(csc_mat_t *M, int stride) {
    int c_base;
    long long int span= 0;

    assert(stride > 0);

    for(c_base= 0; c_base < M->ncol; c_base+= stride) {
        int s;
        int start= INT_MAX, end= -1;

        /* Find the lowest start row and highest end row. */
        for(s= 0; s < stride && c_base + s < M->ncol; s++) {
            int c= c_base + s;

            assert(M->rows[M->ci[c]] < INT_MAX);

            if(M->rows[M->ci[c]] < start) start= M->rows[M->ci[c]];
            if(end < M->rows[M->ci[c+1]-1]) end= M->rows[M->ci[c+1]-1];
        }

        /* Ignore completely empty blocks. */
        if(end < start) continue;

        span+= end - start;
    }

    printf("Average element separation with stride %d: %.2f rows, %.2fB\n",
           stride, ((double)span) / M->nnz,
           sizeof(float) * ((double)span) / M->nnz);
}

void
stride_stats(csc_mat_t *M) {
    int s;
    int min_rows= INT_MAX, max_rows= 0;
    int n_under_4=  0;
    int n_under_8=  0;
    int n_under_16= 0;
    int n_under_32= 0;
    int blocks= 0;
    int partial= 0;
    int collisions= 0;
    int cf_blocks;

    for(s= 0; s < M->ncol / STRIDE_OF(M); s++) {
        int i;

        for(i= M->si[s]; i < M->si[s+1] - ((M->si[s+1]-M->si[s]) % STRIDE_OF(M));
            i+= STRIDE_OF(M)) {
            int rows= M->rows[i+STRIDE_OF(M)-1] - M->rows[i] + 1;
            char cols[STRIDE_OF(M)];
            char collision= 0;
            int j;

            if(!(M->flags & CSC_F_CFREE)) {
                for(j= 0; j < STRIDE_OF(M); j++)
                    cols[j]= 0;

                for(j= 0; j < STRIDE_OF(M); j++) {
                    if(cols[M->sc[i+j]]) collision= 1;
                    else                 cols[M->sc[i+j]]= 1;
                }
                if(collision) collisions++;
            }

            if(rows < min_rows) min_rows= rows;
            if(rows > max_rows) max_rows= rows;
            if(rows <= 4)  n_under_4++;
            if(rows <= 8)  n_under_8++;
            if(rows <= 16) n_under_16++;
            if(rows <= 32) n_under_32++;

            blocks++;
        }

        if(i < M->si[s+1]) partial++;
    }

    printf("%d blocks\n", blocks);
    printf("Min/max rows in a full block: %d %d\n", min_rows, max_rows);
    printf("% 8d blocks with row span <=  4 (%.2f%%)\n", n_under_4,
           (double)n_under_4 / blocks * 100.0);
    printf("% 8d blocks with row span <=  8 (%.2f%%)\n", n_under_8,
           (double)n_under_8 / blocks * 100.0);
    printf("% 8d blocks with row span <= 16 (%.2f%%)\n", n_under_16,
           (double)n_under_16 / blocks * 100.0);
    printf("% 8d blocks with row span <= 32 (%.2f%%)\n", n_under_32,
           (double)n_under_32 / blocks * 100.0);
    printf("% 8d partial blocks (%.2f%%)\n", partial,
           (double)partial / blocks * 100.0);
    if(!(M->flags & CSC_F_CFREE)) {
        printf("% 8d column collisions (%.2f%%)\n", collisions,
               (double)collisions / blocks * 100.0);
    }

    if(!(M->flags & CSC_F_CFREE)) {
        cf_blocks= estimate_cfree(M, 1);
        printf("Would need %.2f%% more (%d total) blocks to be collision-free"
               " with span 1\n",
               (double)(cf_blocks - blocks) / blocks * 100, cf_blocks);
        cf_blocks= estimate_cfree(M, 4);
        printf("Would need %.2f%% more (%d total) blocks to be collision-free"
               " with span 4\n",
               (double)(cf_blocks - blocks) / blocks * 100, cf_blocks);
        cf_blocks= estimate_cfree(M, 8);
        printf("Would need %.2f%% more (%d total) blocks to be collision-free"
               " with span 8\n",
               (double)(cf_blocks - blocks) / blocks * 100, cf_blocks);
        cf_blocks= estimate_cfree(M, 16);
        printf("Would need %.2f%% more (%d total) blocks to be collision-free"
               " with span 16\n",
               (double)(cf_blocks - blocks) / blocks * 100, cf_blocks);
        cf_blocks= estimate_cfree(M, 32);
        printf("Would need %.2f%% more (%d total) blocks to be collision-free"
               " with span 32\n",
               (double)(cf_blocks - blocks) / blocks * 100, cf_blocks);
        cf_blocks= estimate_cfree(M, 64);
        printf("Would need %.2f%% more (%d total) blocks to be collision-free"
               " with span 64\n",
               (double)(cf_blocks - blocks) / blocks * 100, cf_blocks);
    }
}

void
run_length(csc_mat_t *M) {
    int c;
    int64_t total_run= 0;
    int shortest_longest= INT_MAX;

    for(c= 0; c < M->ncol; c++) {
        int i, run_start= INT_MIN, run_end= INT_MIN;
        int longest_start= INT_MIN, longest_end= INT_MIN;

        for(i= M->ci[c]; i < M->ci[c+1]; i++) {
            if(M->rows[i] == run_end+1) {
                /* Extend run. */
                run_end++;
            }
            else {
                /* Start new run. */
                run_start= M->rows[i];
                run_end= run_start;
            }
            if(run_end - run_start > longest_end - longest_start) {
                longest_start= run_start;
                longest_end= run_end;
            }
        }

        total_run+= longest_end - longest_start;

        if(longest_end - longest_start < shortest_longest) {
            shortest_longest= longest_end - longest_start;
        }
    }

    printf("Average contiguous per column: %.2f\n", 
           (float)total_run / M->ncol);
    printf("Fraction contiguous: %.2f%%\n",
           (float)total_run / M->nnz * 100.0);
    printf("Shortest run: %d\n", shortest_longest);
}

int
main(int argc, char *argv[]) {
    FILE *in;
    csc_mat_t *M;
    csc_errno_t e;
    int r;

    if(argc < 2) {
        fprintf(stderr, "Usage: %s <matrix_filename>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    in= fopen(argv[1], "rb");
    if(!in) { perror("fopen"); exit(EXIT_FAILURE); }

    M= csc_load_binary(in, &e);
    if(!M) { csc_perror(e, "csc_load_binary"); exit(EXIT_FAILURE); }
    fclose(in);

    if(!csc_check(M, 1)) abort();
    csc_stats(M);

    printf("Average column length: %.2f\n", (double)M->nnz / M->ncol);

    if(STRIDE_OF(M) == 1) {
        for(r= 0; r <= 10; r++) stride_sep(M, 1<<r);
        run_length(M);
    }
    else {
        printf("Average stride length: %.2f\n",
               (double)M->nnz / (M->ncol/STRIDE_OF(M)));
        stride_stats(M);
    }

    csc_mat_destroy(M);

    return 0;
}
