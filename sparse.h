/* sparse.h

   A library for large sparse histograms, and compressed sparse matrices,
   tailored to support the generation of channel matrices for side-channel
   analysis.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#ifndef __SPARSE_H
#define __SPARSE_H

#include <stdint.h>
#include <stdio.h>

#define ROUND_UP(x,n) (((x)+(n)-1)-(((x)+(n)-1)%n))
#define ROUND_DOWN(x,n) ((x)-((x)%(n)))

/*** Block-Sparse-Column (BSC) Histograms. ***/

/* Quantum for allocating column space. */
#define ROW_BLOCK 16

/* 2-D BSC histogram.  For each nonempty column, stores all entries between
 * the least and greatest rows observed. */
typedef struct bsc_hist {
    int end_col;     /* First unallocated column. */
    int end_row;     /* First row unallocated in all columns. */
    int *start_rows; /* First non-empty row by column block. */
    int *end_rows;   /* First unallocated row by column block. */
    int **entries;   /* Entries indexed by block and index, then row. */
    int *row_total;  /* Totals by (allocated) row. */
    int total;       /* Total of all counts. */
    int64_t nalloc;  /* Number of allocated entries. */
    int64_t nnz;     /* Number of non-zero entries. */
} bsc_hist_t;

/*** Compressed-Sparse-Column Matrices. ***/

#define CSC_MAX_STRIDE 7
#define CSC_MAX_SPAN 7
#define CSC_VERSION 1
#define CSC_MAGIC "CSC_MATRIXLE"

#define CSC_M_STRIDERADIX 0x7
#define CSC_F_CFREE       0x8
#define CSC_I_ROWSPAN     4
#define CSC_M_ROWSPAN     (0x7 << CSC_I_ROWSPAN)
#define CSC_M_ALLFLAGS (CSC_M_ROWSPAN | CSC_F_CFREE | CSC_M_STRIDERADIX)

#define STRIDE_OF(M) (1 << ((M)->flags & CSC_M_STRIDERADIX))
#define SPAN_OF(M) (1 << (((M)->flags & CSC_M_ROWSPAN) >> CSC_I_ROWSPAN))

#define CSC_ROW_INVALID 0xff

/* CSC matrix.  For each nonempty column, stores a list of row,value pairs (as
 * two separate arrays).  Often acheives a compression of greater than 10^4
 * relative to BSC.  Supports a number of optional data formatting options to
 * improve cache locality, or vectorisation. */
typedef struct csc_matrix {
    /* We use explicit sizes to allow serialisation. */
    int32_t ver;    /* Version. */
    int32_t flags;  /* Bits 2:0 -> log2(stride length), 0 means no stride. */
    int32_t nrow;   /* Number of Rows */
    int32_t ncol;   /* Number of Columns */
    int64_t nnz;    /* Number of Non-Zero entries */

    /* Only one of {ci,si} is used, depending on whether there is a stride.
     * The other is NULL. */
    int32_t *ci;    /* start-of-Column Indicies in entries
                       (size ncol+1) */
    int32_t *si;    /* start-of-Stride Indicies in entries
                       (size (ncol/STRIDE)+1) */

    /* Only used with a stride, otherwise NULL. */
    uint8_t *sc;    /* Column number within the stride (size nnz) */

    int32_t *rows;  /* Row numbers of entries (size nnz, or
                       nnz/stride if collision-free) */

    uint8_t *row_offsets;
                    /* Offset from the row number of the containing block.
                       (size nnz), NULL if not collision-free. */

    float *entries; /* Entry values (size nnz) */
} csc_mat_t;

/* A vector, stored densely. */
typedef struct dense_vec {
    int length;
    float *entries;
} dv_t;

typedef enum {
    E_CSC_SUCCESS  = 0,
    E_CSC_TRUNC    = 1,
    E_CSC_BADMAGIC = 2,
    /* Pass through libc errno. */
    E_CSC_ERRNO    = 3,
    E_CSC_COUNT
} csc_errno_t;

#ifdef __AVX__
/* Optional vectorised routines for x86 CPUs. */
#include "sparse_avx.h"
#endif

/*** Histograms ***/

/* Create. */
bsc_hist_t *bsc_hist_new(void);
/* Destroy. */
void bsc_hist_destroy(bsc_hist_t *M);
/* Allocate space for (at least) all columns up to c. */
void bsc_extend(bsc_hist_t *M, int c);
/* Allocate space for (at least) all rows up to r in column c. */
void bsc_extend_col(bsc_hist_t *M, int c, int r);
/* Count samples. Calls bsc_extend and bsc_extend_col as required. */
void bsc_hist_count(bsc_hist_t *M, int c, int r, int n);
/* Generate a CSC matrix by normalising such that each row sums to 1. */
csc_mat_t *bsc_normalise(bsc_hist_t *H);
/* Return size (in bytes). */
uint64_t bsc_size(bsc_hist_t *H);
/* Print statistics. */
void bsc_stats(bsc_hist_t *H);
/* Sanity checking. */
int bsc_check(bsc_hist_t *H, int verbose);

/*** Matrices ***/

/* Destroy.  Matrices are created by normalising a histogram. */
void csc_mat_destroy(csc_mat_t *M);
/* Serialise to fd f */
csc_errno_t csc_store_binary(csc_mat_t *M, FILE *f);
/* De-serialise from fd f */
csc_mat_t *csc_load_binary(FILE *f, csc_errno_t *e);
/* Stride the matrix, by interleaving r_stride adjacent columns. */
void csc_stride(csc_mat_t *M, int r_stride);
/* Left-multiply A by dense row vector x, storing in y. y = x * A.  y must be
 * preallocated, with y->length == A->ncol and x->length == A->nrow. */
void mult_csc_dv(dv_t *y, dv_t *x, csc_mat_t *A);
/* Ditto, for a strided matrix. */
void csc_str_mult_nv(dv_t *y, dv_t *x, csc_mat_t *A);
/* Prune any empty columns. */
void csc_prune_cols(csc_mat_t *M);
/* Make the matrix collision free.  This guarantees that within a stride
 * (across n columns), a column appears at most once.  This allows easier
 * vectorisation.  Can also ensure that blocks have limited row span. */
void csc_make_cfree(csc_mat_t *M, int row_span);
/* Return the number of blocks required to make this matrix collision-free.
 * Called internally by csc_make_cfree. */
int estimate_cfree(csc_mat_t *M, int row_span);
/* Multiplication optimised for collision-free matrices. */
void csc_mult_cf(dv_t *y, dv_t *x, csc_mat_t *A);
/* Multithreaded multiplication. */
void mult_csc_dv_p(dv_t *y, dv_t *x, csc_mat_t *A, int n);
/* Ensure that internal structures are cache-line aligned. */
void csc_align(csc_mat_t *M, int n);
/* Return size (in bytes). */
uint64_t csc_size(csc_mat_t *M);
/* Print statistics. */
void csc_stats(csc_mat_t *M);
/* Sanity checking. */
int csc_check(csc_mat_t *M, int verbose);
/* Verbose error description. */
void csc_perror(csc_errno_t e, const char *s);

/*** Vectors ***/

/* Create. */
dv_t *dv_new(int length);
/* Destroy. */
void dv_destroy(dv_t *v);
/* Initialise with a uniform value. */
void dv_uniform(dv_t *v, float x);
/* Initialise with zero. */
void dv_zero(dv_t *v);
/* Dot product. */
float dv_dot(dv_t *u, dv_t *v);
/* Point-wise maximum. */
float dv_max(dv_t *u);
/* Normalise. */
void dv_normalise(dv_t *v);

#endif /* __SPARSE_H */
