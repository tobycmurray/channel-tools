/* guess.c

   Given a series of outputs on a channel, output most likely guess
   at the input as well as the probability of each input.

   This code is experimental, and error-handling is primitive.
*/

/* Written by Toby Murray */


#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "sparse.h"
#include "bayes_inference.h"

int
main(int argc, char *argv[]) {
    FILE *in;
    csc_mat_t *M;
    csc_errno_t e;
    dv_t *rows; dv_t*cols;

    if(argc < 2) {
        fprintf(stderr, "Usage: %s <matrix_filename>\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }

    in= fopen(argv[1], "rb");
    if(!in) { perror("fopen"); exit(EXIT_FAILURE); }

    M= csc_load_binary(in, &e);
    if(!M) { csc_perror(e, "csc_load_binary"); exit(EXIT_FAILURE); }
    fclose(in);

    if(!csc_check(M, 1)) abort();
    printf("# ");
    csc_stats(M);

    rows = dv_new(M->nrow);
    cols = dv_new(M->ncol);
    if (!rows || !cols){
      perror("couldn't allocate vectors"); exit(EXIT_FAILURE);
    }

    /* begin with uniform probability distribution */
    dv_uniform(rows, (1.0 / M->nrow));

    unsigned int output = 0;
    unsigned int used_outputs = 0;
    while (scanf("%d",&output) == 1){
      if (output >= M->ncol){
	perror("Output out of range. Ignoring.");
	continue;
      }

      bayes_update(M, rows, output, cols);
      used_outputs++;
    }

    /* print out the final distribution */
    int32_t i;
    float max = 0.0;
    int32_t maxi = 0;
    for (i=0; i<rows->length; i++){
      printf("%d %f\n",i,rows->entries[i]);
      if (rows->entries[i] > max){
	max = rows->entries[i];
	maxi = i;
      }
    }

    printf("# Guess after %u outputs: %d, with probability: %f\n",used_outputs,maxi,max);

    dv_destroy(cols);
    dv_destroy(rows);
    csc_mat_destroy(M);

    return 0;
}
