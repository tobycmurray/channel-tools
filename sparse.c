/* sparse.c

   A library for large sparse histograms, and compressed sparse matrices,
   tailored to support the generation of channel matrices for side-channel
   analysis.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <assert.h>
#include <limits.h>
#include <malloc.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sparse.h"

#undef DEBUG_SPARSE

#ifdef DEBUG_SPARSE
    #define D(f,...) printf("%s: " f, __func__, __VA_ARGS__)
#else
    #define D(f,...)
#endif

/* Allocate and return an empty histogram. */
bsc_hist_t*
bsc_hist_new(void) {
    bsc_hist_t *H;

    H= (bsc_hist_t *)calloc(1, sizeof(bsc_hist_t));
    if(!H) { perror("malloc"); abort(); }

    return H;
}

/* Deallocate a histogram. */
void
bsc_hist_destroy(bsc_hist_t *H) {
    int c;
    assert(H);
    if(H->row_total) free(H->row_total);
    if(H->entries) {
        for(c= 0; c < H->end_col; c++) {
            if(H->entries[c]) free(H->entries[c]);
        }
        free(H->entries);
    }
    if(H->start_rows) free(H->start_rows);
    if(H->end_rows) free(H->end_rows);
    free(H);
}

/* Extend histogram to include (at least) column c. */
void
bsc_extend(bsc_hist_t *H, int c) {
    int new_end;
    int i;

    D("%p, %d\n", H, c);

    assert(H);
    assert(0 <= c);
    assert(c < INT_MAX);

    /* Nothing to do */
    if(c < H->end_col) return;

    /* Round from c up the next block size. */
    new_end= c+1;

    assert(H->end_col < new_end);

    /* Extend arrays.  New columns start at INT_MAX... */
    H->start_rows= (int *)realloc(H->start_rows, new_end * sizeof(int));
    if(!H->start_rows) { perror("realloc"); abort(); }
    for(i= H->end_col; i < new_end; i++) H->start_rows[i]= INT_MAX;

    /* and end before 0.  Sentinels. */
    H->end_rows= (int *)realloc(H->end_rows, new_end * sizeof(int));
    if(!H->end_rows) { perror("realloc"); abort(); }
    memset(H->end_rows + H->end_col, 0, (new_end - H->end_col) * sizeof(int));

    H->entries= (int **)realloc(H->entries, new_end * sizeof(int *));
    if(!H->entries) { perror("realloc"); abort(); }
    memset(H->entries + H->end_col, 0, (new_end - H->end_col) * sizeof(int *));

    H->end_col= new_end;
}

/* Extend column c to include all existing rows, plus (at least) row r. */
void
bsc_extend_col(bsc_hist_t *H, int c, int r) {
    int new_start, new_end;
    int *new_col;

    D("%p, %d, %d\n", H, c, r);

    assert(H);
    assert(0 <= c);
    assert(c < H->end_col);
    assert(0 <= r);

    /* Extend the row total array. */
    if(r >= H->end_row) {
        H->row_total= (int *)realloc(H->row_total, (r+1) * sizeof(int));
        if(!H->row_total) { perror("realloc"); abort(); }
        memset(H->row_total + H->end_row, 0,
                (r+1 - H->end_row) * sizeof(int));

        H->end_row= r + 1;
    }

    /* Nothing more to do */
    if(H->start_rows[c] <= r && r < H->end_rows[c]) return;

    new_start= H->start_rows[c];
    new_end=   H->end_rows[c];

    if(r < H->start_rows[c]) {
        /* Round r down to block size */
        new_start= r - (r % ROW_BLOCK);
    }
    if(H->end_rows[c] <= r) {
        /* Round r up to block size */
        new_end=   r + ROW_BLOCK;
        new_end-=  new_end % ROW_BLOCK;
    }
    assert(new_start <= H->start_rows[c]);
    assert(H->end_rows[c] <= new_end);
    assert(new_start <= new_end);

    /* Reallocate column */
    new_col= (int *)calloc(new_end - new_start, sizeof(int));
    if(!new_col) { perror("calloc"); abort(); }

    /* Copy existing entries */
    if(H->start_rows[c] < H->end_rows[c]) { /* Guard against sentinels. */
        assert(H->entries[c] != (void *)0);
        memcpy(new_col + (H->start_rows[c] - new_start),
               H->entries[c],
               (H->end_rows[c] - H->start_rows[c]) * sizeof(int));
        H->nalloc+= (H->start_rows[c] - new_start) +
                    (new_end - H->end_rows[c]);
    }
    else {
        H->nalloc+= ROW_BLOCK;
    }

    /* Free the old column */
    free(H->entries[c]);

    /* Update bookkeeping */
    H->start_rows[c]= new_start;
    H->end_rows[c]=   new_end;
    H->entries[c]=    new_col;
}

/* Increment the count at (r,c), extending as required. */
void
bsc_hist_count(bsc_hist_t *H, int c, int r, int n) {
    assert(H);
    assert(n > 0);

    D("%p, %d, %d, %d\n", H, c, r, n);

    /* Extend as required */
    bsc_extend(H, c);
    bsc_extend_col(H, c, r);
    assert(c < H->end_col);
    assert(H->start_rows[c] <= r);
    assert(r < H->end_rows[c]);

    /* Increment the appropriate entry */
    if(H->entries[c][r - H->start_rows[c]] == 0) H->nnz++;
    H->entries[c][r - H->start_rows[c]]+= n;

    /* Increment the totals */
    H->total+= n;
    H->row_total[r]+= n;
}

/* Normalize a histogram, to give a conditional probability matrix. */
csc_mat_t *
bsc_normalise(bsc_hist_t *H) {
    int c, r, i;
    csc_mat_t *M;

    assert(H);

    M= (csc_mat_t *)calloc(1, sizeof(csc_mat_t));
    if(!M) { perror("calloc"); abort(); }

    /* Allocate and initialise */
    M->ver= CSC_VERSION;
    M->flags= 0; /* No stride. */
    M->nrow= H->end_row;
    M->ncol= H->end_col;
    M->nnz=  H->nnz;
    M->ci=   memalign(64, (M->ncol + 1) * sizeof(int));
    if(!M->ci) { perror("malloc"); abort(); }
    M->rows= memalign(64, M->nnz * sizeof(int));
    if(!M->rows) { perror("malloc"); abort(); }
    M->entries= memalign(64, M->nnz * sizeof(float));
    if(!M->entries) { perror("malloc"); abort(); }

    i= 0;
    for(c= 0; c < H->end_col; c++) {
        /* Record the start of column. */
        M->ci[c]= i;
        for(r= H->start_rows[c]; r < H->end_rows[c]; r++) {
            int n;

            n= H->entries[c][r - H->start_rows[c]];
            if(n > 0) { /* Drop zeros */
                assert(i < M->nnz);
                assert(0 < H->row_total[r]);

                /* Record row number. */
                M->rows[i]= r;
                /* Scale and store count. */
                M->entries[i]= (float)n / H->row_total[r];
                /* Increment in destination array. */
                i++;
            }
        }
    }
    M->ci[c]= i;

    return M;
}

/* Size (in bytes) of the allocated histogram. */
uint64_t
bsc_size(bsc_hist_t *H) {
    int cols= H->end_col;
    int rows= H->end_row;
    assert(H);
           /* The root struct. */
    return sizeof(bsc_hist_t) +
           /* start_rows and end_rows */
           2 * cols * sizeof(int) +
           /* row_total */
           1 * rows * sizeof(int) +
           /* entries */
           cols * sizeof(int *) +
           /* The entries themselves. */
           H->nalloc * sizeof(int);
}

#define FAIL(s,...) { if(verbose) { fprintf(stderr, s, ##__VA_ARGS__); } \
                      return 0; }
#define SUCCEED() { return 1; }

/* Sanity check the histogram structure. */
int
bsc_check(bsc_hist_t *H, int verbose) {
    int r, c;
    int64_t nnz;
    int total;
    int *row_total;

    if(!H) FAIL("H is null\n");
    if(H->nalloc < 0) FAIL("nalloc (%ld) is negative\n", H->nalloc);
    if(H->nnz < 0) FAIL("nnz (%ld) is negative\n", H->nnz);
    if(H->nnz > H->nalloc) FAIL("nnz (%ld) > nalloc (%ld)\n",
            H->nnz, H->nalloc);
    if(H->end_col < 0) FAIL("end_col (%d) is negative\n", H->end_col);
    if(H->end_row < 0) FAIL("end_row (%d) is negative\n", H->end_row);
    if(!H->start_rows) FAIL("start_rows is NULL\n");
    for(c= 0; c < H->end_col; c++) {
        if(H->start_rows[c] < 0)
            FAIL("start_rows[%d] (%d) is negative\n", c, H->start_rows[c]);
    }
    row_total= calloc(H->end_row, sizeof(int));
    if(!row_total) {
        perror("calloc");
        abort();
    }
    nnz= 0; total= 0;
    for(c= 0; c < H->end_col; c++) {
        for(r= H->start_rows[c]; r < H->end_rows[c]; r++) {
            int e= H->entries[c][r - H->start_rows[c]];
            if(e < 0)
                FAIL("entries[%d][%d] (%d) < 0", c, r - H->start_rows[c], e);
            if(e > 0) {
                if(r > H->end_row)
                    FAIL("%d,%d (%d) > end_row (%d)\n", c, r, e, H->end_row);
                nnz++;
                total+= e;
                row_total[r]+= e;
            }
        }
    }
    if(H->nnz != nnz) FAIL("nnz (%ld) != calculated nnz (%ld)\n",
            H->nnz, nnz);
    if(H->total != total) FAIL("total (%d) != calculated total (%d)\n",
            H->total, total);
    for(r= 0; r < H->end_row; r++) {
        if(H->row_total[r] != row_total[r])
            FAIL("row_total[%d] (%d) != calculated row_total[%d] (%d)\n",
                    r, H->row_total[r], r, row_total[r]);
    }

    SUCCEED();
}

/* Sanity check the matrix structure. */
int
csc_check(csc_mat_t *M, int verbose) {
    int64_t i;

    if(!M) FAIL("M is null\n");
    if(M->ver > CSC_VERSION) FAIL("Unknown version %d\n", M->ver);
    if(M->flags & ~CSC_M_ALLFLAGS) FAIL("Unknown flags\n");
    if(M->nrow < 0) FAIL("nrow is negative\n");
    if(M->ncol < 0) FAIL("ncol is negative\n");
    if(M->nnz < 0)  FAIL("nnz is negative\n");
    if(STRIDE_OF(M) > 1) {
        int32_t last= 0;
        if(M->ci) FAIL("Stride is %d but ci is not NULL.\n", STRIDE_OF(M));
        if(!M->si) FAIL("Stride is %d but si is NULL.\n", STRIDE_OF(M));
        if(M->flags & CSC_F_CFREE) {
            if(!M->row_offsets)
                FAIL("CFREE is true but row_offsets is NULL\n");
            if(M->sc) FAIL("CFREE is true d but sc is not NULL.\n");
            for(i= 0; i < M->nnz; i++) {
                if(M->row_offsets[i] >= SPAN_OF(M))
                    FAIL("Row offset %d too large at element %lld\n",
                         M->row_offsets[i], (long long int)i);
            }
        }
        else {
            if(M->row_offsets)
                FAIL("CFREE is false but row_offsets is not NULL\n");
            if(!M->sc) FAIL("Stride is %d but sc is NULL.\n", STRIDE_OF(M));
            for(i= 0; i < M->nnz; i++) {
                if(M->sc[i] > STRIDE_OF(M))
                    FAIL("Stride column %d is greater than stride width"
                          " %d at element %lld\n", M->sc[i], STRIDE_OF(M),
                          (long long int)i);
            }
        }
        for(i= 0; i < M->ncol / STRIDE_OF(M) + 1; i++) {
            int j;
            int32_t last_row= 0;
            if(M->si[i] < 0 || M->si[i] > M->nnz)
                FAIL("Stride index %d out of range at element %lld\n",
                      M->si[i], (long long int)i);
            if(M->si[i] < last)
                FAIL("Stride index %d out of order at element %lld\n",
                        M->si[i], (long long int)i);
            last= M->si[i];

            if(M->flags & CSC_F_CFREE) {
                if(M->si[i] % STRIDE_OF(M))
                    FAIL("Index for stride %lld (%d) not a multiple of"
                         " stride\n", (long long int)i, M->si[i]);
                for(j= M->si[i]; j < M->si[i+1]; j+= STRIDE_OF(M)) {
                    if(M->rows[j/STRIDE_OF(M)] < last_row)
                        FAIL("Row %d out of order at element %lld\n",
                             M->rows[j/STRIDE_OF(M)], (long long int)j);

                    last_row= M->rows[j/STRIDE_OF(M)];
                }
            }
            else {
                for(j= M->si[i]; j < M->si[i+1]; j++) {
                    if(M->rows[j] < last_row)
                        FAIL("Row %d out of order at element %lld\n",
                             M->rows[j], (long long int)j);

                    last_row= M->rows[j];
                }
            }
        }
    }
    else {
        int32_t last = 0;
        if(M->flags & CSC_F_CFREE)
            FAIL("Stride is %d but flag CFREE set\n", STRIDE_OF(M));
        if(!M->ci) FAIL("Stride is %d but ci is NULL.\n", STRIDE_OF(M));
        for(i= 0; i < M->ncol+1; i++) {
            if(M->ci[i] < 0 || M->ci[i] > M->nnz)
                FAIL("Column index %d out of range at element %lld\n",
                      M->ci[i], (long long int)i);
            if(M->ci[i] < last)
                FAIL("Column index %d out of order at element %lld\n",
                        M->ci[i], (long long int)i);
            last= M->ci[i];
        }
        if(M->si) FAIL("Stride is %d but si is not NULL.\n", STRIDE_OF(M));
        if(M->sc) FAIL("Stride is %d but sc is not NULL.\n", STRIDE_OF(M));
    }
    if(!M->rows) FAIL("rows is NULL\n");
    if(M->flags & CSC_F_CFREE) {
        for(i= 0; i < M->nnz/STRIDE_OF(M); i++) {
            if(M->rows[i] < 0 || M->rows[i] >= M->nrow)
                FAIL("Invalid row number %d in entry %lld\n",
                     M->rows[i], (long long int)i);
        }
    }
    else {
        for(i= 0; i < M->nnz; i++) {
            if(M->rows[i] < 0 || M->rows[i] >= M->nrow)
                FAIL("Invalid row number %d in entry %lld\n",
                     M->rows[i], (long long int)i);
        }
    }
    if(!M->entries) FAIL("entries is NULL\n");

    SUCCEED();
}


/* Size (in bytes) of the allocated matrix. */
uint64_t
csc_size(csc_mat_t *M) {
    /* Strided matrix. */
    if(STRIDE_OF(M) > 1) {
               /* root struct */
        return sizeof(csc_mat_t) +
               /* si */
               (M->ncol/STRIDE_OF(M)+1) * sizeof(int32_t) +
               /* sc */
               (M->nnz) * sizeof(uint8_t) +
               /* rows */
               (M->nnz) * sizeof(int32_t) +
               /* entries */
               (M->nnz) * sizeof(float);
    }
    /* Unstrided matrix. */
    else {
               /* root struct */
        return sizeof(csc_mat_t) +
               /* ci */
               (M->ncol+1) * sizeof(int32_t) +
               /* rows */
               (M->nnz) * sizeof(int32_t) +
               /* entries */
               (M->nnz) * sizeof(float);
    }
}

/* Deallocate matrix. */
void
csc_mat_destroy(csc_mat_t *M) {
    assert(M);
    if(STRIDE_OF(M) > 1) {
        assert(M->si);
        free(M->si);
        if(M->flags & CSC_F_CFREE) {
            assert(M->row_offsets);
            free(M->row_offsets);
        }
        else {
            assert(M->sc);
            free(M->sc);
        }
    }
    else {
        assert(M->ci);
        free(M->ci);
    }
    assert(M->rows);
    free(M->rows);
    assert(M->entries);
    free(M->entries);
    free(M);
}

static const char *csc_errstr[] = {
    "Success",
    "File truncated",
    "Bad magic bytes",
    "System error - check errno",
};

/* Readable errors. */
void
csc_perror(csc_errno_t e, const char *s) {
    if(e < 0 || e >= E_CSC_COUNT) return;

    if(e == E_CSC_ERRNO) perror(s);
    else fprintf(stderr, "%s: %s\n", s, csc_errstr[e]);
}

/* Serialise into FD f */
csc_errno_t
csc_store_binary(csc_mat_t *M, FILE *f) {
    if(fputs(CSC_MAGIC, f) == EOF) return E_CSC_ERRNO;
    if(fwrite(&M->ver,    sizeof(int32_t), 1, f) != 1)
        return E_CSC_ERRNO;
    if(fwrite(&M->flags,  sizeof(int32_t), 1, f) != 1)
        return E_CSC_ERRNO;
    if(fwrite(&M->nrow,   sizeof(int32_t), 1, f) != 1)
        return E_CSC_ERRNO;
    if(fwrite(&M->ncol,   sizeof(int32_t), 1, f) != 1)
        return E_CSC_ERRNO;
    if(fwrite(&M->nnz,    sizeof(int64_t), 1, f) != 1)
        return E_CSC_ERRNO;
    if(STRIDE_OF(M) > 1) {
        if(fwrite(M->si, sizeof(int32_t),
                  M->ncol/STRIDE_OF(M)+1, f) != M->ncol/STRIDE_OF(M)+1)
            return E_CSC_ERRNO;
        if(M->flags & CSC_F_CFREE) {
            if(fwrite(M->row_offsets, sizeof(uint8_t), M->nnz, f) != M->nnz)
                return E_CSC_ERRNO;
        }
        else {
            if(fwrite(M->sc, sizeof(uint8_t), M->nnz, f) != M->nnz)
                return E_CSC_ERRNO;
        }
    }
    else {
        if(fwrite(M->ci, sizeof(int32_t), M->ncol+1, f) != M->ncol+1)
            return E_CSC_ERRNO;
    }
    if(fwrite(M->rows,    sizeof(int32_t), M->nnz, f) != M->nnz)
        return E_CSC_ERRNO;
    if(fwrite(M->entries, sizeof(float), M->nnz, f)   != M->nnz)
        return E_CSC_ERRNO;

    return E_CSC_SUCCESS;
}

/* Deserialise from FD f */
csc_mat_t *
csc_load_binary(FILE *f, csc_errno_t *e) {
    csc_mat_t *M;
    size_t n;
    char magic[12];

    /* Check for magic */
    n= fread(magic, 1, 12, f);
    if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
    if(n < 12) { *e= E_CSC_TRUNC; return NULL; }
    if(strncmp(magic, CSC_MAGIC, 12)) { *e= E_CSC_BADMAGIC; return NULL; }

    /* Allocate root struct */
    M= (csc_mat_t *)calloc(1,sizeof(csc_mat_t));
    if(!M) { *e= E_CSC_ERRNO; return NULL; }

    /* Read header */
    n= 0;
    n+= fread(&M->ver,   sizeof(int32_t), 1, f);
    if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
    n+= fread(&M->flags, sizeof(int32_t), 1, f);
    if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
    n+= fread(&M->nrow,  sizeof(int32_t), 1, f);
    if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
    n+= fread(&M->ncol,  sizeof(int32_t), 1, f);
    if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
    n+= fread(&M->nnz,   sizeof(int64_t), 1, f);
    if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
    if(n != 5) { free(M); *e= E_CSC_TRUNC; return NULL; }

    /* Allocate tables */
    M->rows= memalign(64, M->nnz * sizeof(int32_t));
    if(!M->rows) { free(M); *e= E_CSC_ERRNO; return NULL; }
    M->entries= memalign(64, M->nnz * sizeof(int32_t));
    if(!M->entries) { free(M->rows); free(M); *e= E_CSC_ERRNO; return NULL; }
    if(STRIDE_OF(M) > 1) {
        M->si= memalign(64, (M->ncol/STRIDE_OF(M) + 1) * sizeof(int32_t));
        if(!M->si) { free(M->entries); free(M->rows); free(M);
                     *e= E_CSC_ERRNO; return NULL; }
        if(M->flags & CSC_F_CFREE) {
            M->row_offsets= memalign(64, M->nnz * sizeof(uint8_t));
            if(!M->row_offsets) { free(M->si); free(M->entries);
                                  *e= E_CSC_ERRNO; free(M->rows); free(M);
                                  return NULL; }
        }
        else {
            M->sc= memalign(64, M->nnz * sizeof(uint8_t));
            if(!M->sc) { free(M->si); free(M->entries); *e= E_CSC_ERRNO;
                         free(M->rows); free(M); return NULL; }
        }
    }
    else {
        M->ci= memalign(64, (M->ncol + 1) * sizeof(int32_t));
        if(!M->ci) { free(M->entries); free(M->rows); free(M);
                     *e= E_CSC_ERRNO; return NULL; }
    }

    /* Read tables */
    if(STRIDE_OF(M) > 1) {
        n= 0;
        n+= fread(M->si, sizeof(int32_t), M->ncol/STRIDE_OF(M) + 1, f);
        if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
        if(M->flags & CSC_F_CFREE) {
            n+= fread(M->row_offsets, sizeof(uint8_t), M->nnz, f);
        }
        else {
            n+= fread(M->sc, sizeof(uint8_t), M->nnz, f);
        }
        if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
        n+= fread(M->rows, sizeof(int32_t), M->nnz, f);
        if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
        n+= fread(M->entries, sizeof(float),   M->nnz, f);
        if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
        if(n != M->ncol/STRIDE_OF(M) + 1 + 3 * M->nnz) {
            free(M->si); free(M->sc); free(M->entries);
            free(M->rows); free(M);
            *e= E_CSC_TRUNC; return NULL;
        }
    }
    else {
        n= 0;
        n+= fread(M->ci, sizeof(int32_t), M->ncol + 1, f);
        if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
        n+= fread(M->rows,    sizeof(int32_t), M->nnz,      f);
        if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
        n+= fread(M->entries, sizeof(float),   M->nnz,      f);
        if(ferror(f)) { *e= E_CSC_ERRNO; return NULL; }
        if(n != M->ncol+1 + 2 * M->nnz) {
            free(M->ci); free(M->entries);
            free(M->rows); free(M); *e= E_CSC_TRUNC; return NULL;
        }
    }

    return M;
}

/* Prune (remove) empty columns. */
void
csc_prune_cols(csc_mat_t *M) {
    int c_old, c_new;

    assert(STRIDE_OF(M) == 1);

    c_new= 0;
    for(c_old= 0; c_old < M->ncol; c_old++) {
        /* Only preserve nonempty columns. */
        if(M->ci[c_old] != M->ci[c_old+1]) {
            assert(c_new <= c_old);
            /* Re-index the column.  We never move anything to the right. */
            M->ci[c_new]= M->ci[c_old];
            c_new++;
        }
    }
    /* Move the sentinel. */
    M->ci[c_new]= M->ci[c_old];

    /* Shrink the array. */
    M->ncol= c_new;
    M->ci= realloc(M->ci, (M->ncol+1) * sizeof(int32_t));
    if(!M->ci) { perror("realloc"); abort(); }
}

void
mult_csc_dv_range(dv_t *y, dv_t *x, csc_mat_t *A, int start, int end) {
    int c, i;

    assert(x); assert(A); assert(y);
    assert(STRIDE_OF(A) == 1);
    assert(A->nrow == x->length);
    assert(A->ncol == y->length);
    assert(0 <= start && start <= end && end <= A->ncol);

#ifdef __AVX__
#if 0
    mult_csc_dv_avx(y, x, A, 0, A->ncol);
    return;
#endif
#endif

    /* Step through the columns */
    for(c= start; c < end; c++) {
        float tmp= 0.0;
        /* Step through the non-zero entries of this column. */
        for(i= A->ci[c]; i < A->ci[c+1]; i++) {
            /* Multiply by the appropriate entry of y. */
            tmp+= A->entries[i] * x->entries[A->rows[i]];
        }
        y->entries[c]= tmp;
    }
}

/* Compute y= x*A, where A is a Compressed-Sparse-Column sparse matrix,
   and y and x are dense vectors. */
void
mult_csc_dv(dv_t *y, dv_t *x, csc_mat_t *A) {
    mult_csc_dv_range(y, x, A, 0, A->ncol);
}

/* Threaded multiplication. */
struct mult_args {
    dv_t *y, *x;
    csc_mat_t *A;
    int start, end;
};

void *
mult_thread(void *_args) {
    struct mult_args *args= (struct mult_args *)_args;

    mult_csc_dv_range(args->y, args->x, args->A, args->start, args->end);

    return NULL;
}

void
mult_csc_dv_p(dv_t *y, dv_t *x, csc_mat_t *A, int n) {
    pthread_t threads[n];
    pthread_attr_t attr;
    struct mult_args args[n];
    int i;
    int block= A->ncol / n;

    assert(0 < n);

    if(pthread_attr_init(&attr))
        { perror("pthread_attr_init"); abort(); }

    for(i= 0; i < n - 1; i++) {
        args[i].y= y;
        args[i].x= x;
        args[i].A= A;
        args[i].start= i * block;
        args[i].end= (i+1) * block;

        if(pthread_create(&threads[i], &attr, mult_thread, &args[i]))
            { perror("pthread_create"); abort(); }
    }
    args[n-1].y= y;
    args[n-1].x= x;
    args[n-1].A= A;
    args[n-1].start= (n-1) * block;
    args[n-1].end= A->ncol;

    if(pthread_create(&threads[n-1], &attr, mult_thread, &args[n-1]))
        { perror("pthread_create"); abort(); }

    for(i= 0; i < n; i++) {
        if(pthread_join(threads[i], NULL))
            { perror("pthread_join"); abort(); }
    }

    if(pthread_attr_destroy(&attr))
        { perror("pthread_attr_destroy"); abort(); }
}

/* Fallback non-vectorised multiply for strided matrix. */
void
csc_str_mult_nv(dv_t *y, dv_t *x, csc_mat_t *A) {
    int s;
    int stride;

    assert(x); assert(A); assert(y);
    assert(A->nrow == x->length);
    assert(A->ncol == y->length);
    assert(A->ncol % STRIDE_OF(A) == 0);

#ifdef __AVX__
    if(STRIDE_OF(A) == 4) {
        csc_str_mult_nv_4(y, x, A);
        return;
    }
    else if(STRIDE_OF(A) == 8) {
        csc_str_mult_nv_8(y, x, A);
        return;
    }
#endif

    stride= STRIDE_OF(A);

    for(s= 0; s < A->ncol/stride; s++) {
        int i;
        float tmp[stride];

        memset(tmp, 0, stride * sizeof(float));

        for(i= A->si[s]; i < A->si[s+1]; i++) {
            assert(A->sc[i] < stride);

            tmp[A->sc[i]]+= A->entries[i] * x->entries[A->rows[i]];
        }

        for(i= 0; i < stride; i++) {
            y->entries[s * stride + i]= tmp[i];
        }
    }
}

/* Collision-free multiplication. */
void
csc_mult_cf(dv_t *y, dv_t *x, csc_mat_t *A) {
    int s;
    int stride;

    assert(x); assert(A); assert(y);
    assert(A->nrow == x->length);
    assert(A->ncol == y->length);
    assert(A->ncol % STRIDE_OF(A) == 0);
    assert(A->nnz % STRIDE_OF(A) == 0);
    assert(A->flags & CSC_F_CFREE);

#ifdef __AVX__
    if(STRIDE_OF(A) == 4 && SPAN_OF(A) == 4) {
        csc_mult_cf_4_4(y, x, A);
        return ;
    }
#endif

    stride= STRIDE_OF(A);

    dv_zero(y);

    for(s= 0; s < A->ncol/stride; s++) {
        int i;
        float tmp[stride];

        memset(tmp, 0, stride * sizeof(float));

        assert(A->si[s] % stride == 0);
        assert(A->si[s+1] % stride == 0);

        for(i= A->si[s]; i < A->si[s+1]; i++) {
            int row= A->rows[i/stride] + A->row_offsets[i];

            tmp[i%stride]+= A->entries[i] * x->entries[row];
        }

        for(i= 0; i < stride; i++) {
            y->entries[s * stride + i]= tmp[i];
        }
    }
}

int
find_row(int32_t *rows, int s, int e, int32_t r) {
    assert(s <= e);
    assert(s < INT_MAX/2);
    assert(e < INT_MAX/2);

    /* Binary search. */
    while(s < e) {
        int i= (s + e) / 2;

        if(rows[i] <= r) s= i+1;
        else             e= i;
    }

    assert(s == e);
    return e;
}

void
sort_stride(csc_mat_t *M, int s) {
    int stride_start= M->si[s];
    int stride_end=   M->si[s+1];
    int i;

    printf("Stride %d, length %d\n", s, stride_end - stride_start);

    /* An insertion sort for the moment. */
    for(i= stride_start + 1; i < stride_end; i++) {
        /* INV: M->rows[stride_start,i) is sorted. */

        /* Find a spot for element i in [stride_start,i]. */
        int j= find_row(M->rows, stride_start, i, M->rows[i]);

        /* This guard is technically unnecessary. */
        if(j != i) {
            /* Save the values to move. */
            uint8_t sc_i=      M->sc[i];
            int32_t rows_i=    M->rows[i];
            float   entries_i= M->entries[i];

            /* Shift the array tails along. */
            memmove(&M->sc[j+1],      &M->sc[j],      (i-j) * sizeof(uint8_t));
            memmove(&M->rows[j+1],    &M->rows[j],    (i-j) * sizeof(int32_t));
            memmove(&M->entries[j+1], &M->entries[j], (i-j) * sizeof(float));

            /* Put the saved values at j. */
            M->sc[j]=      sc_i;
            M->rows[j]=    rows_i;
            M->entries[j]= entries_i;
        }
    }
}

void
merge_stride(csc_mat_t *M, int s, int width) {
    int stride_start= M->si[s];
    int stride_end=   M->si[s+1];
    int len= stride_end - stride_start;
    uint8_t *sc_tmp=      malloc(len * sizeof(uint8_t));
    int32_t *rows_tmp=    malloc(len * sizeof(int32_t));;
    float   *entries_tmp= malloc(len * sizeof(float));;
    int     *sc_i=        malloc(width * sizeof(int));
    int i;

    if(width == 0) return;

    if(!sc_tmp || !rows_tmp || !entries_tmp)
        { perror("malloc"); exit(EXIT_FAILURE); }

    /* printf("Merging stride %d, length %d, width %d\n", s, len, width); */

    /* Initialise the column indices. */
    for(i= 0; i < width; i++) sc_i[i]= M->ci[s * STRIDE_OF(M) + i];

    /* Merge columns element-by-element into temporary space. */ 
    for(i= 0; i < len; i++) {
        int j;
        int min_j= -1, min_row= INT_MAX;

        /* Find the smallest remaining head. */
        for(j= 0; j < width; j++) {
            /* Skip exhausted columns. */
            if(sc_i[j] == M->ci[s * STRIDE_OF(M) + j + 1]) continue;

            if(M->rows[sc_i[j]] < min_row) {
                min_j= j;
                min_row= M->rows[sc_i[j]];
            }
        }
        assert(min_j != -1); /* We must have found something. */

        /* Copy the elements. */
        sc_tmp[i]=      M->sc[sc_i[min_j]];
        rows_tmp[i]=    M->rows[sc_i[min_j]];
        entries_tmp[i]= M->entries[sc_i[min_j]];

        /* Cut the head off the column we just used. */
        sc_i[min_j]++;
    }

    /* We should have exhausted all columns. */
    for(i= 0; i < width; i++)
        assert(sc_i[i] == M->ci[s * STRIDE_OF(M) + i + 1]);

    free(sc_i);

    /* Write back the temporary workspace. */
    memcpy(&M->sc[stride_start],      sc_tmp,     len * sizeof(uint8_t));
    memcpy(&M->rows[stride_start],    rows_tmp,    len * sizeof(int32_t));
    memcpy(&M->entries[stride_start], entries_tmp, len * sizeof(float));

    free(entries_tmp); free(rows_tmp); free(sc_tmp);
}

/* Convert a matrix to strided format, in-place. */
void
csc_stride(csc_mat_t *M, int r_stride) {
    int s, c;
    int ncol_new;
    int stride= 1<<r_stride;

    assert(M);
    assert(0 <= r_stride && r_stride <= CSC_MAX_STRIDE);
    assert(STRIDE_OF(M) == 1);

    /* Set the stride. */
    M->flags &= ~CSC_M_STRIDERADIX;
    M->flags |= r_stride & CSC_M_STRIDERADIX;

    /* Round up the number of columns. */
    ncol_new= ROUND_UP(M->ncol, stride);

    /* Allocate the stride structures. */
    M->si= malloc(sizeof(int32_t) * (ncol_new / stride + 1));
    M->sc= malloc(sizeof(uint8_t) * M->nnz);
    if(!M->si || !M->sc) { perror("malloc"); exit(EXIT_FAILURE); }

    /* Strides start at the beginning of the first column contained. */
    for(s= 0; s < ncol_new / stride; s++) {
        M->si[s]= M->ci[s * stride];
    }
    /* The sentinel must be moved explicitly. */
    M->si[ncol_new / stride]= M->ci[M->ncol];

    /* The entries all start out grouped by column. */
    for(c= 0; c < M->ncol; c++) {
        int i;

        for(i= M->ci[c]; i < M->ci[c+1]; i++) M->sc[i]= c % stride;
    }

    /* Now sort the buggers */
    for(s= 0; s < M->ncol / STRIDE_OF(M); s++)
        merge_stride(M, s, STRIDE_OF(M));
    /* The last stride may be incomplete in the original index. */
    merge_stride(M, s, M->ncol % STRIDE_OF(M));

    /* We're now done with the column indices (they're now
       invalid anyway), so free them. */
    free(M->ci); M->ci= NULL;

    M->ncol= ncol_new;
}

/* Estimate the number of blocks required to make the matrix
   collision-free. */
int
estimate_cfree(csc_mat_t *M, int row_span) {
    int s;
    int cf_blocks= 0;

    for(s= 0; s < M->ncol / STRIDE_OF(M); s++) {
        int i;
        /* Force non-empty strides to start a new block. */
        int block_start_row = -row_span;
        int cols[STRIDE_OF(M)];

        for(i= M->si[s]; i < M->si[s+1]; i++) {
            if(M->rows[i] - block_start_row >= row_span || cols[M->sc[i]]) {
                int j;

                for(j= 0; j < STRIDE_OF(M); j++) cols[j]= 0;
                block_start_row= M->rows[i];
                cf_blocks++;
            }
            cols[M->sc[i]]= 1;
        }
    }

    return cf_blocks;
}

/* Re-lay-out the matrix to be collision free for columns within stride blocks,
   and so that stride blocks have row spans of less than 2^'row_span'.  We
   insert zeroes to pad out the blocks. */
void
csc_make_cfree(csc_mat_t *M, int row_span) {
    float *new_entries;
    int32_t *new_rows;
    int cf_blocks;
    int s, i, block_index;

    assert(M);

    /* Nothing to do. */
    if(STRIDE_OF(M) == 1 || M->flags & CSC_F_CFREE) return;

    M->flags&= ~CSC_M_ROWSPAN;
    M->flags |= (row_span << CSC_I_ROWSPAN) & CSC_M_ROWSPAN;

    /* How many blocks will we need? */
    cf_blocks= estimate_cfree(M, SPAN_OF(M));

    new_entries= calloc(cf_blocks * STRIDE_OF(M), sizeof(float));
    new_rows= calloc(cf_blocks, sizeof(int32_t));
    M->row_offsets= calloc(cf_blocks * STRIDE_OF(M), sizeof(uint8_t));
    if(!new_entries || !new_rows || !M->row_offsets)
        { perror("calloc"); exit(EXIT_FAILURE); }

    block_index= -1;
    for(s= 0; s < M->ncol / STRIDE_OF(M); s++) {
        /* Force non-empty strides to start a new block. */
        int block_start_row= -SPAN_OF(M);
        int cols[STRIDE_OF(M)];
        int j;
        int start_of_stride= 1;
        int start_block;

        /* Empty strides take the index of the last allocated block. */
        if(block_index == -1) start_block= 0;
        else start_block= block_index;

        for(i= M->si[s]; i < M->si[s+1]; i++) {
            int new_index;

            /* If we've exceeded the row span or had a column collision,
               start a new block. */
            if(M->rows[i] - block_start_row >= SPAN_OF(M) || cols[M->sc[i]]) {
                for(j= 0; j < STRIDE_OF(M); j++) cols[j]= 0;
                block_start_row= M->rows[i];
                assert(block_start_row < M->nrow);
                block_index++;

                assert(block_index >= 0 && block_index < cf_blocks);

                /* Save the new start row. */
                new_rows[block_index]= block_start_row;

                /* If this is the first new block for the stride,
                   record its location. */
                if(start_of_stride) {
                    start_block= block_index;
                    start_of_stride= 0;
                }
            }
            cols[M->sc[i]]= 1;

            assert(block_start_row != -1);
            assert(block_index != -1);
            new_index= block_index * STRIDE_OF(M) + M->sc[i];

            /* Copy the entry. */
            new_entries[new_index]= M->entries[i];

            /* Save the row offset. */
            assert(M->rows[i] < M->nrow);
            assert(M->rows[i] - block_start_row < SPAN_OF(M));
            M->row_offsets[new_index]= M->rows[i] - block_start_row;
        }

        /* Record the new stride start. */
        M->si[s]= start_block * STRIDE_OF(M);
    }
    /* Insert a new sentinel. */
    M->si[s]= (block_index+1) * STRIDE_OF(M);

    /* Column number is now implicit. */
    free(M->sc);
    M->sc= NULL;

    /* Save the new arrays. */
    free(M->rows); M->rows= new_rows;
    free(M->entries); M->entries= new_entries;

    M->nnz= cf_blocks * STRIDE_OF(M);
    M->flags |= CSC_F_CFREE;
}

/* Cache-align structures. */
void
csc_align(csc_mat_t *M, int n) {
    uint64_t nnz_new;
    int c;
    int32_t *new_rows;
    float *new_entries;
    uint64_t j;

    assert(STRIDE_OF(M) == 1);
    assert(0 < n);

    nnz_new= 0;
    for(c= 0; c < M->ncol; c++) {
        nnz_new+= ROUND_UP(M->ci[c+1] - M->ci[c], n);
    }

    new_rows= memalign(64, nnz_new * sizeof(int32_t));
    new_entries= memalign(64, nnz_new * sizeof(float));
    if(!new_rows || !new_entries) {
        perror("memalign");
        exit(EXIT_FAILURE);
    }

    memset(new_rows, 0xff, nnz_new * sizeof(int32_t));
    memset(new_entries, 0, nnz_new * sizeof(int32_t));

    j= 0;
    for(c= 0; c < M->ncol; c++) {
        int i;
        int new_ci= j;

        for(i= M->ci[c]; i < M->ci[c+1]; i++) {
            new_rows[j]= M->rows[i];
            new_entries[j]= M->entries[j];
            j++;
        }
        j= ROUND_UP(j, n);
        M->ci[c]= new_ci;
    }
    M->ci[c]= nnz_new;

    free(M->rows); M->rows= new_rows;
    free(M->entries); M->entries= new_entries;
    M->nnz= nnz_new;
}

/* Allocate a vector. */
dv_t *
dv_new(int length) {
    dv_t *x;

    assert(0 <= length);

    x= (dv_t *)malloc(sizeof(dv_t));
    if(!x) { perror("malloc"); abort(); }
    x->length= length;
    x->entries= memalign(64, length * sizeof(float));
    if(!x->entries) { perror("malloc"); abort(); }

    return x;
}

/* Deallocate a vector. */
void
dv_destroy(dv_t *v) {
    assert(v);
    assert(v->entries);
    free(v->entries);
    free(v);
}

/* Zero entries. */
void
dv_zero(dv_t *v) {
    memset(v->entries, 0, v->length * sizeof(float));
}

/* Initialise uniformly. */
void
dv_uniform(dv_t *v, float x) {
    int i;

    for(i= 0; i < v->length; i++) v->entries[i]= x;
}

/* Dot product. */
float
dv_dot(dv_t *u, dv_t *v) {
    float x= 0.0;
    int i;

    assert(u); assert(v);
    assert(u->length == v->length);

    for(i= 0; i < u->length; i++) x+= u->entries[i] * v->entries[i];

    return x;
}

/* Maximum entry. */
float
dv_max(dv_t *u) {
    float x= -INFINITY;
    int i;

    for(i= 0; i < u->length; i++) {
        if(u->entries[i] > x) x = u->entries[i];
    }

    return x;
}

/* Normalise such that entries sum to 1. */
void
dv_normalise(dv_t *v) {
    float Ps;
    int i;

    Ps= 0;
    for(i= 0; i < v->length; i++) Ps+= v->entries[i];
    for(i= 0; i < v->length; i++) v->entries[i]/= Ps;
}

/* Print statistics. */
void
bsc_stats(bsc_hist_t *H) {
    printf("bsc_hist: %dc x %dr, %llunz, %.2fMb, %.2fb/e\n",
           H->end_col, H->end_row,
           (long long unsigned int)H->nnz, bsc_size(H) / (1024.0 * 1024),
           (double)bsc_size(H) / H->nnz);
}

/* Print statistics. */
void
csc_stats(csc_mat_t *M) {
    printf("csc_mat: %dc x %dr, %llunz, %.2f%% dense, "
           "%.2fMb, %.2fb/e, stride %d",
           M->ncol, M->nrow,
           (long long unsigned int)M->nnz,
           100.0 * (double)M->nnz /
           ((long long unsigned int)M->nrow * M->ncol),
           csc_size(M) / (1024.0 * 1024),
           (double)csc_size(M) / M->nnz, STRIDE_OF(M));
    if(M->flags & CSC_F_CFREE) printf(", cfree span %d", SPAN_OF(M));
    printf("\n");
}
