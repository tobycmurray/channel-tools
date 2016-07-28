/* channel_algorithms.h

   Implementations of the Arimoto-Blahut algorithm (and improvements), with
   varying optimisations, using sparse matrices.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#ifndef __CHANNEL_ALGORITHMS_H
#define __CHANNEL_ALGORITHMS_H

#include "sparse.h"

#ifdef CAP_BENCH
extern int iterations;
#endif

/* The basic algorithm, in single-precision arithmetic. */
float blahut_arimoto(csc_mat_t *Q, float epsilon, float *e_obs);

/* The basic algorithm, using double-precision arithmetic internally. */
float blahut_arimoto_precise(csc_mat_t *Q, float epsilon, float *e_obs);

/* The squeezed algorithm of Yu, in double precision -
 *   Yaming Yu, "Squeezing the Arimoto–Blahut Algorithm for Faster
 *   Convergence," Information Theory, IEEE Transactions on , vol.56, no.7,
 *   pp.3149,3157, July 2010
 *   doi: 10.1109/TIT.2010.2048452
 */
float blahut_arimoto_precise_squeezed(csc_mat_t *Q, float epsilon, float *e_obs);

/* The squeezed algorithm in long double. */
float blahut_arimoto_ld_squeezed(csc_mat_t *Q, float epsilon, float *e_obs);

/* The squeezed algorithm, with adaptively increased precision.  This is
 * generally the best choice. */
float ba_phased(csc_mat_t *Q, float epsilon, float *e_obs);

#endif /* __CHANNEL_ALGORITHMS_H */
