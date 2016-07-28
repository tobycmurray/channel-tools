/* log_bench.c

   Benchmark logarithm implementations.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "log.h"

#define REPS 100000000

#define TEST_STEPS 1000000

int
main(int argc, char *argv[]) {
    struct timespec start, end;
    int i;
    double iv;
    float x, y;
    float emax= 0.0;
    __v4sf w, z;
    __v8sf w8, z8;

    write_log_table();

    x= 0;
    for(i= 0; i < 2 * TEST_STEPS; i++) {
        float e;
        x+= 1.0 / TEST_STEPS;
        e= fabsf(log2f_table(x) - log2f(x));
        if(e > emax) emax= e;
    }
    printf("Max error between 0 and 2: %e\n", emax);

    x= 1.0; y= 0.0;
    clock_gettime(CLOCK_REALTIME, &start);
    for(i= 0; i < REPS; i++) {
        y+= log2f(x);
        x+= 1.0;
    }
    clock_gettime(CLOCK_REALTIME, &end);
    printf("GLIBC log2f: %.2f\n", y);

    iv= end.tv_sec   + end.tv_nsec*1e-9
      - start.tv_sec - start.tv_nsec*1e-9;
    printf("%d reps %.2esec %.2e/s\n", REPS, iv, REPS / iv);

    x= 1.0; y= 0.0;
    clock_gettime(CLOCK_REALTIME, &start);
    for(i= 0; i < REPS; i++) {
        y+= log2f_table(x);
        x+= 1.0;
    }
    clock_gettime(CLOCK_REALTIME, &end);
    printf("log2f_table: %.2f\n", y);

    iv= end.tv_sec   + end.tv_nsec*1e-9
      - start.tv_sec - start.tv_nsec*1e-9;
    printf("%d reps %.2esec %.2e/s\n", REPS, iv, REPS / iv);

    w[0]= 1.0; w[1]= 1.1; w[2]= 1.2; w[3]= 1.3;
    z[0]= 0.0; z[1]= 0.0; z[2]= 0.0; z[3]= 0.0;
    clock_gettime(CLOCK_REALTIME, &start);
    for(i= 0; i < REPS; i++) {
        z+= log2f_table_v(w);
        w+= 1.0;
    }
    clock_gettime(CLOCK_REALTIME, &end);
    printf("log2f_table_v: %.2f %.2f %.2f %.2f\n",
           z[0], z[1], z[2], z[3]);

    iv= end.tv_sec   + end.tv_nsec*1e-9
      - start.tv_sec - start.tv_nsec*1e-9;
    printf("%d reps %.2esec %.2e/s\n", REPS, iv, 4 * REPS / iv);

    w8[0]= 1.0; w8[1]= 1.1; w8[2]= 1.2; w8[3]= 1.3;
    z8[0]= 0.0; z8[1]= 0.0; z8[2]= 0.0; z8[3]= 0.0;
    clock_gettime(CLOCK_REALTIME, &start);
    for(i= 0; i < REPS; i++) {
        z8+= log2f_table_v8(w8);
        w8+= 1.0;
    }
    clock_gettime(CLOCK_REALTIME, &end);
    printf("log2f_table_v8: %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",
           z8[0], z8[1], z8[2], z8[3], z8[4], z8[5], z8[6], z8[7]);

    iv= end.tv_sec   + end.tv_nsec*1e-9
      - start.tv_sec - start.tv_nsec*1e-9;
    printf("%d reps %.2esec %.2e/s\n", REPS, iv, 8 * REPS / iv);

    return 0;
}
