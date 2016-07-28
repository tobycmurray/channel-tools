/* smooth.c

   For each input, output Sum([n-r,n+r]).

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <stdio.h>
#include <stdlib.h>

#define LINE_MAX 1024

int
main(int argc, char *argv[]) {
    int n, count, i;
    int *xs;
    double *ys;
    char buf[LINE_MAX];

    if(argc < 2) {
        fprintf(stderr, "Usage: %s <averaging radius>\n", argv[0]);
        return 1;
    }

    n= atoi(argv[1]);

    xs= malloc((2*n + 1) * sizeof(int));
    ys= malloc((2*n + 1) * sizeof(double));
    if(!xs || !ys) {
        perror("malloc");
        return 1;
    }

    i= 0; count= 0;
    while(fgets(buf, LINE_MAX, stdin)) {
        int x;
        double y;

        if(sscanf(buf, "%d %lf", &x, &y) != 2) continue;

        xs[i]= x;
        ys[i]= y;
        i= (i+1) % (2*n + 1);
        count++;

        if(count >= 2*n + 1) {
            int j;
            double z= 0.0;

            for(j= 0; j < 2*n + 1; j++)
                z+= ys[(i + j) % (2*n + 1)];

            printf("%d %.12e\n", xs[(i+n) % (2*n + 1)], z / (2*n + 1));
        }
    }

    free(ys);
    free(xs);

    return 0;
}
