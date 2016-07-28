/* sparse_avx.h

   (not very well) AVX-optimised csc matrix operations.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#ifndef __SPARSE_AVX_H
#define __SPARSE_AVX_H

void csc_str_mult_nv_4(dv_t *y, dv_t *x, csc_mat_t *A);
void csc_str_mult_nv_8(dv_t *y, dv_t *x, csc_mat_t *A);
void csc_mult_cf_4_4(dv_t *y, dv_t *x, csc_mat_t *A);

#endif /* __SPARSE_AVX_H */
