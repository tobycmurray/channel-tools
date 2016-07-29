/* bayes_inference.c

  Implementation of simple inference that updates the belief about a
  channel's input, given the channel matrix and an observation on that channel.

  This code is experimental, and error-handling is primitive.
*/

/* Written by Toby Murray */

#include <assert.h>

#include "bayes_inference.h"

void bayes_update(csc_mat_t *M, dv_t *belief, unsigned int c,
                  dv_t *scratch){
  assert(belief->length == M->nrow);
  assert(c < M->ncol);
  assert(scratch->length == M->ncol);
  assert(STRIDE_OF(M) == 1); /* TODO: handle strided matrices */

  /* calculate the probability distribution on outputs (columns) */
  mult_csc_dv(scratch, belief, M);

  int32_t i;
  for(i = M->ci[c]; i < M->ci[c+1]; i++) {
    int32_t r = M->rows[i];
    float prob = M->entries[i];
    
    /* prob is p(c|r). We update belief[r] to contain the probability of
       r given the observation c, i.e. p(r|c) = p(c|r)*p(r)/p(c) */
    belief->entries[r] *= (prob / scratch->entries[c]);
  }

  /* normalise belief -- TODO: necessary? */
  dv_normalise(belief);
}
