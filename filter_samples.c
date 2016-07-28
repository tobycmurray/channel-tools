/* filter_samples.c

   Filter out-of-range samples.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <stdio.h>
#include <stdlib.h>

#define LINE_MAX 1024

int
main(int argc, char *argv[]) {
    int cmin, cmax, rmin, rmax;
    int c, r;
    int discarded= 0;
    char buf[LINE_MAX];

    if(argc < 5) {
        fprintf(stderr, "Usage: %s <cmin> <cmax> <rmin> <rmax>\n", argv[0]);
        return 1;
    }

    cmin= atoi(argv[1]);
    cmax= atoi(argv[2]);
    rmin= atoi(argv[3]);
    rmax= atoi(argv[4]);

    while(fgets(buf, LINE_MAX, stdin)) {
        if(scanf("%d %d\n", &c, &r) == 2) {
            if(cmin <= c && c <= cmax &&
               rmin <= r && r <= rmax) {
                printf("%d %d\n", c, r);
            }
            else discarded++;
        }
    }

    fprintf(stderr, "Discarded %d samples\n", discarded);

    return 0;
}
