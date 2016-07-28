/* mvec.c

   Measure raw arithmetic speed.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef __AVX__
#error "Need AVX"
#endif

#include <immintrin.h>

#define REPS 1000000000ULL

__v8sf
test_reg(void) {
    struct timespec start, end;
    double iv;
    long long int i;
    __v8sf x= {1,1,1,1,1,1,1,1},
           y= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z= {0,0,0,0,0,0,0,0};
    __v8sf x2= {1,1,1,1,1,1,1,1},
           y2= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z2= {0,0,0,0,0,0,0,0};
    __v8sf x3= {1,1,1,1,1,1,1,1},
           y3= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z3= {0,0,0,0,0,0,0,0};
    __v8sf x4= {1,1,1,1,1,1,1,1},
           y4= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z4= {0,0,0,0,0,0,0,0};
    __v8sf x5= {1,1,1,1,1,1,1,1},
           y5= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z5= {0,0,0,0,0,0,0,0};
    __v8sf x6= {1,1,1,1,1,1,1,1},
           y6= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z6= {0,0,0,0,0,0,0,0};
    __v8sf x7= {1,1,1,1,1,1,1,1},
           y7= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z7= {0,0,0,0,0,0,0,0};
    __v8sf x8= {1,1,1,1,1,1,1,1},
           y8= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z8= {0,0,0,0,0,0,0,0};
    __v8sf x9= {1,1,1,1,1,1,1,1},
           y9= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z9= {0,0,0,0,0,0,0,0};
    __v8sf x10= {1,1,1,1,1,1,1,1},
           y10= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z10= {0,0,0,0,0,0,0,0};
    __v8sf x11= {1,1,1,1,1,1,1,1},
           y11= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z11= {0,0,0,0,0,0,0,0};
    __v8sf x12= {1,1,1,1,1,1,1,1},
           y12= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z12= {0,0,0,0,0,0,0,0};

    printf("Testing register-register arithmetic speed...");
    fflush(stdout);

    if(clock_gettime(CLOCK_REALTIME, &start)) {
        perror("clock_gettime"); exit(EXIT_FAILURE);
    }
    asm volatile("foo:nop");
    for(i= 0; i < REPS; i++) {
        z= __builtin_ia32_addps256(z, __builtin_ia32_mulps256(y,x));
        z2= __builtin_ia32_addps256(z2, __builtin_ia32_mulps256(y2,x2));
        z3= __builtin_ia32_addps256(z3, __builtin_ia32_mulps256(y3,x3));
        z4= __builtin_ia32_addps256(z4, __builtin_ia32_mulps256(y4,x4));
        z5= __builtin_ia32_addps256(z5, __builtin_ia32_mulps256(y5,x5));
        z6= __builtin_ia32_addps256(z6, __builtin_ia32_mulps256(y6,x6));
        z7= __builtin_ia32_addps256(z7, __builtin_ia32_mulps256(y7,x7));
        z8= __builtin_ia32_addps256(z8, __builtin_ia32_mulps256(y8,x8));
        z9= __builtin_ia32_addps256(z9, __builtin_ia32_mulps256(y9,x9));
        z10= __builtin_ia32_addps256(z10, __builtin_ia32_mulps256(y10,x10));
        z11= __builtin_ia32_addps256(z11, __builtin_ia32_mulps256(y11,x11));
        z12= __builtin_ia32_addps256(z12, __builtin_ia32_mulps256(y12,x12));
    }
    asm volatile("bar:nop");
    if(clock_gettime(CLOCK_REALTIME, &end)) {
        perror("clock_gettime"); exit(EXIT_FAILURE);
    }

    iv= end.tv_sec + end.tv_nsec*1e-9
      - start.tv_sec - start.tv_nsec*1e-9;
    printf("%.3e REPS %.3es %.3fGFLOPS\n", (double)REPS, iv,
           (12 * 2 * 8 * REPS) / iv / 1e9);

    return z + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11 + z12;
}

__v8sf
test_rm(void) {
    __v8sf x[256]; /* 8kB */
    __v8sf y= {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5},
           z= {0,0,0,0,0,0,0,0},
           z2= {0,0,0,0,0,0,0,0},
           z3= {0,0,0,0,0,0,0,0},
           z4= {0,0,0,0,0,0,0,0},
           z5= {0,0,0,0,0,0,0,0},
           z6= {0,0,0,0,0,0,0,0},
           z7= {0,0,0,0,0,0,0,0},
           z8= {0,0,0,0,0,0,0,0};
    struct timespec start, end;
    double iv;
    long long int i;

    printf("Testing register-memory arithmetic speed...");
    fflush(stdout);

    for(i= 0; i < 256*8; i++) ((float *)x)[i]= 1.0;

    if(clock_gettime(CLOCK_REALTIME, &start)) {
        perror("clock_gettime"); exit(EXIT_FAILURE);
    }
    asm volatile("foo_rm:");
    for(i= 0; i < REPS; i++) {
        z= __builtin_ia32_addps256(z, __builtin_ia32_mulps256(y,x[(8*i)%256]));
        z2= __builtin_ia32_addps256(z2, __builtin_ia32_mulps256(y,x[(8*i+1)%256]));
        z3= __builtin_ia32_addps256(z3, __builtin_ia32_mulps256(y,x[(8*i+2)%256]));
        z4= __builtin_ia32_addps256(z4, __builtin_ia32_mulps256(y,x[(8*i+3)%256]));
        z5= __builtin_ia32_addps256(z5, __builtin_ia32_mulps256(y,x[(8*i+4)%256]));
        z6= __builtin_ia32_addps256(z6, __builtin_ia32_mulps256(y,x[(8*i+5)%256]));
        z7= __builtin_ia32_addps256(z7, __builtin_ia32_mulps256(y,x[(8*i+6)%256]));
        z8= __builtin_ia32_addps256(z8, __builtin_ia32_mulps256(y,x[(8*i+7)%256]));
    }
    asm volatile("bar_rm:");
    if(clock_gettime(CLOCK_REALTIME, &end)) {
        perror("clock_gettime"); exit(EXIT_FAILURE);
    }

    iv= end.tv_sec + end.tv_nsec*1e-9
      - start.tv_sec - start.tv_nsec*1e-9;
    printf("%.3e REPS %.3es %.3fGFLOPS\n", (double)REPS, iv,
           (8 * 2 * 8 * REPS) / iv / 1e9);

    return z + z2 + z3 + z4 + z5 + z6 + z7 + z8;
}

__v8sf
test_mm(void) {
    __v8sf x[128]; /* 4kB */
    __v8sf y[128]; /* 4kB */
    __v8sf z= {0,0,0,0,0,0,0,0},
           z2= {0,0,0,0,0,0,0,0},
           z3= {0,0,0,0,0,0,0,0},
           z4= {0,0,0,0,0,0,0,0},
           z5= {0,0,0,0,0,0,0,0},
           z6= {0,0,0,0,0,0,0,0},
           z7= {0,0,0,0,0,0,0,0},
           z8= {0,0,0,0,0,0,0,0};
    struct timespec start, end;
    double iv;
    long long int i;

    printf("Testing memory-memory arithmetic speed...");
    fflush(stdout);

    for(i= 0; i < 128*8; i++) ((float *)x)[i]= 1.0;
    for(i= 0; i < 128*8; i++) ((float *)y)[i]= 1.0;

    if(clock_gettime(CLOCK_REALTIME, &start)) {
        perror("clock_gettime"); exit(EXIT_FAILURE);
    }
    asm volatile("foo_mm:");
    for(i= 0; i < REPS; i++) {
        z= __builtin_ia32_addps256(z,
                __builtin_ia32_mulps256(y[(8*i)%128],x[(8*i)%128]));
        z2= __builtin_ia32_addps256(z2,
                __builtin_ia32_mulps256(y[(8*i+1)%128],x[(8*i+1)%128]));
        z3= __builtin_ia32_addps256(z3,
                __builtin_ia32_mulps256(y[(8*i+2)%128],x[(8*i+2)%128]));
        z4= __builtin_ia32_addps256(z4,
                __builtin_ia32_mulps256(y[(8*i+3)%128],x[(8*i+3)%128]));
        z5= __builtin_ia32_addps256(z5,
                __builtin_ia32_mulps256(y[(8*i+4)%128],x[(8*i+4)%128]));
        z6= __builtin_ia32_addps256(z6,
                __builtin_ia32_mulps256(y[(8*i+5)%128],x[(8*i+5)%128]));
        z7= __builtin_ia32_addps256(z7,
                __builtin_ia32_mulps256(y[(8*i+6)%128],x[(8*i+6)%128]));
        z8= __builtin_ia32_addps256(z8,
                __builtin_ia32_mulps256(y[(8*i+7)%128],x[(8*i+7)%128]));
    }
    asm volatile("bar_mm:");
    if(clock_gettime(CLOCK_REALTIME, &end)) {
        perror("clock_gettime"); exit(EXIT_FAILURE);
    }

    iv= end.tv_sec + end.tv_nsec*1e-9
      - start.tv_sec - start.tv_nsec*1e-9;
    printf("%.3e REPS %.3es %.3fGFLOPS\n", (double)REPS, iv,
           (8 * 2 * 8 * REPS) / iv / 1e9);

    return z + z2 + z3 + z4 + z5 + z6 + z7 + z8;
}

#define STRIDE 1
#define LENGTH (1<<25)
#define BIGREPS 10

__v8sf
test_big(long long int ymask) {
    long long int i, j;
    struct timespec start, end;
    __v8sf *x, *y;
    __v8sf z[STRIDE], zf= {0,0,0,0,0,0,0,0};
    double iv;
    int r;

    printf("Testing memory-memory arithmetic speed with large vectors...\n");
    printf("\tone working set restricted to %lldkB\n", ((ymask+1)*32)/1024);

    r=  posix_memalign((void **)&x, sizeof(__v8sf), LENGTH * sizeof(__v8sf));
    r|= posix_memalign((void **)&y, sizeof(__v8sf), LENGTH * sizeof(__v8sf));
    if(r) { fprintf(stderr, "posix_memalign failed"); abort(); }

    for(j= 0; j < LENGTH * 8ULL; j++) {
        ((float *)x)[j]= 1.0;
    }
    for(j= 0; j < LENGTH * 8ULL; j++) {
        ((float *)y)[j]= 0.5;
    }

    for(j= 0; j < STRIDE * 8; j++) {
        ((float *)z)[j]= 0.0;
    }

    clock_gettime(CLOCK_REALTIME, &start);
    asm volatile("foo_big:");
    for(r= 0; r < BIGREPS; r++) {
        for(i= 0; i < LENGTH; i+= STRIDE) {
            for(j=0; j < STRIDE; j++) {
                z[j]= __builtin_ia32_addps256(z[j],
                        __builtin_ia32_mulps256(y[(i+j)&ymask],x[i+j]));
            }
        }
    }
    asm volatile("end_big:");
    clock_gettime(CLOCK_REALTIME, &end);

    iv= end.tv_sec + end.tv_nsec*1e-9
      - start.tv_sec - start.tv_nsec*1e-9;
    printf("\t%.3e REPS %.3es %.3fGFLOPS\n", BIGREPS * (double)LENGTH, iv,
           (BIGREPS * 2 * 8ULL * LENGTH) / iv / 1e9);

    free(x);
    free(y);

    for(j= 0; j < STRIDE; j++) zf+= z[j];
    return zf;
}

int
main(void) {
    test_reg();
    test_rm();
    test_mm();
    test_big(0xffffffff);
    test_big(0x1fff);
    test_big(0xff);

    return 0;
}
