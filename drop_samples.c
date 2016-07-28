/* drop_samples.c

   Drop n samples from every modulation.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <stdio.h>
#include <stdlib.h>

#define MAX_LINE 1024

int
main(int argc, char *argv[]) {
    int in_old= -1, count= 0, discard;
    int in, out;
    char buf[MAX_LINE];

    if(argc < 2) {
        fprintf(stderr, "Usage: %s <samples to discard per modulation>\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }

    discard= atoi(argv[1]);

    while(fgets(buf, MAX_LINE, stdin)) {
        if(sscanf(buf, "%d %d", &in, &out) == 2) {
            if(in != in_old) {
                in_old= in;
                count= 0;
            }
            count++;
            if(count > discard) printf("%d %d\n", in, out);
        }
    }

    return EXIT_SUCCESS;
}
