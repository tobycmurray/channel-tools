/* testlib.c

   Helpers for both test_hist and test_sparse.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#ifndef __TESTLIB_H
#define __TESTLIB_H

bsc_hist_t *bsc_random(int nrows, int ncols, int nent, int quiet);
dv_t       *dv_random(int length, float mag, int quiet);

#endif /* __TESTLIB_H */
