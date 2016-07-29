/* bayes_inference.h

  Implementation of simple inference that updates the belief about a
  channel's input, given the channel matrix and an observation on that channel.

  This code is experimental, and error-handling is primitive.
*/

/* Written by Toby Murray */

#ifndef __BAYES_INFERENCE_H
#define __BAYES_INFERENCE_H

#include "sparse.h"

/* For a channel characterised by the matrix M, and an initial
 * probability distribution on inputs, belief, 
 * given a new observation on the channel c (denoting a column of M), 
 * update belief by applying Bayes' rule via the channel matrix M.
 * The vector scratch is used as temporary working space for this function
 * to hold the calculated output probability distribution from the initial
 * belief. Its initial contents will be overwritten. 
 * Requires: belief->length == M->nrow,  c < M->ncol, 
 *           scratch->length == M->ncol */   
void bayes_update(csc_mat_t *M, dv_t *belief, unsigned int c,
                  dv_t *scratch);

#endif
