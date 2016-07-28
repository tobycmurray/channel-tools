/* channel_algorithms.c

   Implementations of the Arimoto-Blahut algorithm (and improvements), with
   varying optimisations, using sparse matrices.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "log.h"
#include "sparse.h"

#undef DEBUG_CAPACITY

#ifdef DEBUG_CAPACITY
#define D(fmt,...) fprintf(stderr,fmt,## __VA_ARGS__)
#else
#define D(fmt,...)
#endif

#ifdef CAP_BENCH
int iterations;
#endif

/*** The Basic Algorithm ***/

/* Arimoto, S., "An algorithm for computing the capacity of arbitrary discrete
 * memoryless channels," Information Theory, IEEE Transactions on , vol.18,
 * no.1, pp.14,20, Jan 1972
 * doi: 10.1109/TIT.1972.1054753
 */

/* Blahut, R.E., "Computation of channel capacity and rate-distortion
 * functions," Information Theory, IEEE Transactions on , vol.18, no.4,
 * pp.460,473, Jul 1972
 * doi: 10.1109/TIT.1972.1054855 */

/* Precalculates the logs of table entries, and keeps results as logs wherever
 * possible to avoid precision loss. */
float
blahut_arimoto(csc_mat_t *Q, float epsilon, float *e_obs) {
    dv_t *p, *q, *c;
    float Il, Iu, e;
    float Il_last= 0, Iu_last= INFINITY;
    float Il_best= 0, e_best= INFINITY;
    int i, col, row;
    float *logQ, *logq;

    D("\n");

    p= dv_new(Q->nrow);
    q= dv_new(Q->ncol);
    c= dv_new(Q->nrow);
    if(!p || !q || !c) { perror("dv_new"); abort(); }

    D("Precalculating logs...");
    logQ= malloc(Q->nnz * sizeof(float));
    if(!logQ) { perror("malloc"); abort(); }
    /* We precalculate single-precision logs.  It's slightly faster on modern
     * hardware with large caches. */
    for(i= 0; i < Q->nnz; i++) logQ[i]= log2f_table(Q->entries[i]);
    D(" done.\n");

    logq= malloc(q->length * sizeof(float));
    if(!logq) { perror("malloc"); abort(); }

    /* Start with a uniform input distribution. */
    dv_uniform(p, 1.0 / Q->nrow);

#ifdef CAP_BENCH
    iterations= 0;
#endif

    /* Invariant: Il <= I <= Iu */
    while(1) {
#ifdef CAP_BENCH
        iterations++;
#endif

        /* Calculate q - output distribution. */
        mult_csc_dv(q, p, Q);

        /* Arithmetic imprecision causes the sum to wander. */
        dv_normalise(q);

        /* Precalculate log(q). */
        for(i= 0; i < q->length; i++) logq[i]= log2f_table(q->entries[i]);

        /* Calculate log(c). */
        dv_zero(c);
        for(col= 0; col < Q->ncol; col++) {
            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                row= Q->rows[i];

                assert(Q->entries[i] > 0);

                /* Use the log tables. */
                c->entries[row]+= Q->entries[i] * (logQ[i] - logq[col]);

            }
        }

        /* Find the largest element of c.  This is Iu, as log is monotonic. */
        Iu= dv_max(c);

        {
            float tmp= 0.0;

            /* Update Il.  By taking out the largest element of c, the
             * rest are all less than 1. */
            for(row= 0; row < p->length; row++) {
                tmp+= p->entries[row] * exp2f(c->entries[row] - Iu);
            }
            Il= Iu + log2f_table(tmp);
        }

        e= Iu - Il;

        /* Once we reach the precision limit, we will lose monotonicity, so
         * keep track of the best result we've seen. */
        if(e < e_best) {
            Il_best= Il; e_best= e;
        }

        D("Il= %.6e Iu= %.6e e= %.6e\n", Il, Iu, e);
        if(e < epsilon || (Il==Il_last && Iu==Iu_last)) break;

        Il_last= Il; Iu_last= Iu;

        /* Update input distribution.  Here again, we avoid overflow. */
        for(i= 0; i < p->length; i++) {
            p->entries[i]*= exp2f(c->entries[i] - Il);
        }

        /* Again, prevent the sum from wandering. */
        dv_normalise(p);
    } 

    free(logq);
    free(logQ);
    free(p);
    free(q);
    free(c);

    if(e_obs) *e_obs= e_best;
    return Il_best;
}

/*** The Double-Precision Algorithm ***/

float
blahut_arimoto_precise(csc_mat_t *Q, float epsilon, float *e_obs) {
    double *p, *q, *c;
    double Il, Iu, e;
    double Il_best, e_best= INFINITY;
    int i, col, row;
    double *logQ, *logq;

    D("\n");

    p= malloc(Q->nrow * sizeof(double));
    q= malloc(Q->ncol * sizeof(double));
    c= malloc(Q->nrow * sizeof(double));
    if(!p || !q || !c) { perror("malloc"); abort(); }

    D("Precalculating logs...");
    logQ= malloc(Q->nnz * sizeof(double));
    if(!logQ) { perror("malloc"); abort(); }
    /* We can't precalculate all DP logs. */
    for(i= 0; i < Q->nnz; i++) logQ[i]= log2((double)Q->entries[i]);
    D(" done.\n");

    logq= malloc(Q->ncol * sizeof(double));
    if(!logq) { perror("malloc"); abort(); }

    /* Start with a uniform input distribution. */
    for(i= 0; i < Q->nrow; i++) p[i]= 1.0 / Q->nrow;

#ifdef CAP_BENCH
    iterations= 0;
#endif

    while(1) {
#ifdef CAP_BENCH
        iterations++;
#endif

        /* Calculate q - output distribution. */
        {
            int c;

            for(c= 0; c < Q->ncol; c++) {
                int i;

                q[c]= 0.0;
                for(i= Q->ci[c]; i < Q->ci[c+1]; i++)
                    q[c]+= (double)Q->entries[i] * p[Q->rows[i]];
            }
        }

        /* Precalculate log(q). */
        for(i= 0; i < Q->ncol; i++) logq[i]= log2(q[i]);

        /* Calculate log(c). */
        bzero(c, Q->nrow * sizeof(double));
        for(col= 0; col < Q->ncol; col++) {
            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                row= Q->rows[i];

                assert(Q->entries[i] > 0);

                /* Use the log tables. */
                c[row]+= (double)Q->entries[i] * (logQ[i] - logq[col]);

            }
        }

        /* Find the largest element of c.  This is Iu, as log is monotonic. */
        Iu= -INFINITY;
        for(i= 0; i < Q->nrow; i++) {
            if(c[i] > Iu) Iu= c[i];
        }

        {
            double tmp= 0.0;

            /* Update Il.  By taking out the largest element of c, the
             * rest are all less than 1. */
            for(row= 0; row < Q->nrow; row++) {
                tmp+= p[row] * exp2(c[row] - Iu);
            }
            Il= Iu + log2(tmp);
        }

        e= Iu - Il;

        if(e < e_best) {
            Il_best= Il; e_best= e;
        }

        D("Il= %.6e Iu= %.6e e= %.6e\n", Il, Iu, e);
        if(e < epsilon) break;

        /* Update input distribution.  Here again, we avoid overflow. */
        for(i= 0; i < Q->nrow; i++) {
            p[i]*= exp2(c[i] - Il);
        }
    } 

    free(logq);
    free(logQ);
    free(p);
    free(q);
    free(c);

    if(e_obs) *e_obs= (float)e_best;
    return (float)Il_best;
}

/*** The Squeezed Algorithm ***/

/* Yaming Yu, "Squeezing the Arimoto–Blahut Algorithm for Faster Convergence,"
 * Information Theory, IEEE Transactions on , vol.56, no.7, pp.3149,3157, July
 * 2010
 * doi: 10.1109/TIT.2010.2048452
*/

/* In double precision. */
float
blahut_arimoto_precise_squeezed(csc_mat_t *Q, float epsilon, float *e_obs) {
    double *p, *q, *z;
    double Il, Iu, e;
    double Il_best, e_best= INFINITY;
    int i, col, row;
    double *logQ, *logq;
    double lambda;
    double scale_den;

    D("\n");

    p= malloc(Q->nrow * sizeof(double));
    q= malloc(Q->ncol * sizeof(double));
    z= malloc(Q->nrow * sizeof(double));
    if(!p || !q || !z) { perror("malloc"); abort(); }

    D("Precalculating logs...");
    logQ= malloc(Q->nnz * sizeof(double));
    if(!logQ) { perror("malloc"); abort(); }
    for(i= 0; i < Q->nnz; i++) logQ[i]= log2((double)Q->entries[i]);
    D(" done.\n");

    logq= malloc(Q->ncol * sizeof(double));
    if(!logq) { perror("malloc"); abort(); }

    /* Start with a uniform input distribution. */
    for(i= 0; i < Q->nrow; i++) p[i]= 1.0 / Q->nrow;

    /* Find the greatest accelerating/squeezing factor guaranteeing monotone
     * convergence. */
    {
        double som= 0.0;

        for(col= 0; col < Q->ncol; col++) {
            double col_min= INFINITY;

            if(Q->ci[col] >= Q->ci[col+1]) continue;

            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                if(Q->entries[i] < col_min)
                    col_min= Q->entries[i];
            }

            som+= col_min;
        }

        D("som= %.e\n", som);
        lambda= 1.0 / (1.0 - som);
        D("lambda= %.e\n", lambda);
    }

#ifdef CAP_BENCH
    iterations= 0;
#endif

    while(1) {
#ifdef CAP_BENCH
        iterations++;
#endif

        /* Calculate q - output distribution. */
        {
            int c;

            for(c= 0; c < Q->ncol; c++) {
                int i;

                q[c]= 0.0;
                for(i= Q->ci[c]; i < Q->ci[c+1]; i++)
                    q[c]+= (double)Q->entries[i] * p[Q->rows[i]];
            }
        }

        /* Precalculate log(q). */
        for(i= 0; i < Q->ncol; i++) logq[i]= log2(q[i]);

        /* Calculate KL divergence. */
        bzero(z, Q->nrow * sizeof(double));
        for(col= 0; col < Q->ncol; col++) {
            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                row= Q->rows[i];

                assert(Q->entries[i] > 0);

                /* Use the log tables. */
                z[row]+= (double)Q->entries[i] * (logQ[i] - logq[col]);

            }
        }

        /* Find the largest element of c.  This is Iu, as log is monotonic. */
        Iu= -INFINITY;
        for(i= 0; i < Q->nrow; i++) {
            if(z[i] > Iu) Iu= z[i];
        }

        {
            double tmp= 0.0;

            /* Update Il.  By taking out the largest element of c, the
             * rest are all less than 1. */
            for(row= 0; row < Q->nrow; row++) {
                tmp+= p[row] * exp2(z[row] - Iu);
            }
            Il= Iu + log2(tmp);
        }

        e= Iu - Il;

        if(e < e_best) {
            Il_best= Il; e_best= e;
        }

        D("Il= %.6e Iu= %.6e e= %.6e\n", Il, Iu, e);
        if(e < epsilon) break;

        /* Find the denominator of the probability scale equation. */
        scale_den= 0.0;
        for(i= 0; i < Q->nrow; i++)
            scale_den+= p[i] * exp2(lambda * z[i]);

        /* Scale the input probabilities. */
        for(i= 0; i < Q->nrow; i++)
            /* The squeezing factor is applied here. */
            p[i]= p[i] * exp2(lambda * z[i]) / scale_den;
    } 

    free(logq);
    free(logQ);
    free(p);
    free(q);
    free(z);

    if(e_obs) *e_obs= (float)e_best;
    return (float)Il_best;
}

/* In long-double precision. */
float
blahut_arimoto_ld_squeezed(csc_mat_t *Q, float epsilon, float *e_obs) {
    long double *p, *q, *z;
    long double Il, Iu, e;
    long double Il_best, e_best= INFINITY;
    int i, col, row;
    long double *logQ, *logq;
    long double lambda;
    long double scale_den;

    D("\n");

    p= malloc(Q->nrow * sizeof(long double));
    q= malloc(Q->ncol * sizeof(long double));
    z= malloc(Q->nrow * sizeof(long double));
    if(!p || !q || !z) { perror("malloc"); abort(); }

    D("Precalculating logs...");
    logQ= malloc(Q->nnz * sizeof(long double));
    if(!logQ) { perror("malloc"); abort(); }
    for(i= 0; i < Q->nnz; i++) logQ[i]= log2l((long double)Q->entries[i]);
    D(" done.\n");

    logq= malloc(Q->ncol * sizeof(long double));
    if(!logq) { perror("malloc"); abort(); }

    /* Start with a uniform input distribution. */
    for(i= 0; i < Q->nrow; i++) p[i]= 1.0 / Q->nrow;

    /* Find the greatest accelerating/squeezing factor guaranteeing monotone
     * convergence. */
    {
        long double som= 0.0;

        for(col= 0; col < Q->ncol; col++) {
            long double col_min= INFINITY;

            if(Q->ci[col] >= Q->ci[col+1]) continue;

            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                if(Q->entries[i] < col_min)
                    col_min= Q->entries[i];
            }

            som+= col_min;
        }

        D("som= %.le\n", (double)som);
        lambda= 1.0 / (1.0 - som);
        D("lambda= %.le\n", (double)lambda);
    }

#ifdef CAP_BENCH
    iterations= 0;
#endif

    while(1) {
#ifdef CAP_BENCH
        iterations++;
#endif

        /* Calculate q - output distribution. */
        {
            int c;

            for(c= 0; c < Q->ncol; c++) {
                int i;

                q[c]= 0.0;
                for(i= Q->ci[c]; i < Q->ci[c+1]; i++)
                    q[c]+= (long double)Q->entries[i] * p[Q->rows[i]];
            }
        }

        /* Precalculate log(q). */
        for(i= 0; i < Q->ncol; i++) logq[i]= log2l(q[i]);

        /* Calculate KL divergence. */
        bzero(z, Q->nrow * sizeof(long double));
        for(col= 0; col < Q->ncol; col++) {
            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                row= Q->rows[i];

                assert(Q->entries[i] > 0);

                /* Use the log tables. */
                z[row]+= (long double)Q->entries[i] * (logQ[i] - logq[col]);

            }
        }

        /* Find the largest element of c.  This is Iu, as log is monotonic. */
        Iu= -INFINITY;
        for(i= 0; i < Q->nrow; i++) {
            if(z[i] > Iu) Iu= z[i];
        }

        {
            long double tmp= 0.0;

            /* Update Il.  By taking out the largest element of c, the
             * rest are all less than 1. */
            for(row= 0; row < Q->nrow; row++) {
                tmp+= p[row] * exp2l(z[row] - Iu);
            }
            Il= Iu + log2l(tmp);
        }

        e= Iu - Il;

        if(e < e_best) {
            Il_best= Il; e_best= e;
        }

        D("Il= %.6le Iu= %.6le e= %.6le\n", (double)Il, (double)Iu, (double)e);
        if(e < epsilon) break;

        /* Find the denominator of the probability scale equation. */
        scale_den= 0.0;
        for(i= 0; i < Q->nrow; i++)
            scale_den+= p[i] * exp2l(lambda * z[i]);

        /* Scale the input probabilities. */
        for(i= 0; i < Q->nrow; i++)
            p[i]= p[i] * exp2l(lambda * z[i]) / scale_den;
    } 

    free(logq);
    free(logQ);
    free(p);
    free(q);
    free(z);

    if(e_obs) *e_obs= (float)e_best;
    return (float)Il_best;
}

/*** The Adaptive-Precision or Squeezed Algorithm ***/

/* We allow this many iterations for recovery after a precision bump, before
 * looking for further nonmonotonicity. */
#define GRACE_PERIOD 25

/* After this many nonmonotonic steps, we increase precision. */
#define NONMONOTONIC_TOLERANCE 4

/* We bump up underflowed probabilities by this much. */
#define JIGGER 1e-3

/* The fastest, single-precision phase. */
float
single_phase(csc_mat_t *Q, float epsilon, float *e_obs, float *p,
        float lambda) {
    float *q, *z;
    float Il, Iu, e, e_last= *e_obs;
    float Il_best= 0.0, e_best= INFINITY;
    int col, row;
    float *logQ, *logq;
    float tmp;
    int increases= 0, iterations= 0;

    D("Single-precision phase\n");

    q= malloc(Q->ncol * sizeof(float));
    z= malloc(Q->nrow * sizeof(float));
    if(!q || !z) { perror("malloc"); abort(); }

    D("Precalculating logs...");
    logQ= malloc(Q->nnz * sizeof(float));
    if(!logQ) { perror("malloc"); abort(); }
    {
        int i;
        for(i= 0; i < Q->nnz; i++) logQ[i]= log2f(Q->entries[i]);
    }
    D(" done.\n");

    logq= malloc(Q->ncol * sizeof(float));
    if(!logq) { perror("malloc"); abort(); }

    while(1) {
        int limited= 0;

        /* Calculate q - output distribution. */
        for(col= 0; col < Q->ncol; col++) {
            int i;

            q[col]= 0.0;
            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                if(p[Q->rows[i]] == 0) limited= 1;
                q[col]+= Q->entries[i] * p[Q->rows[i]];
            }
        }
        if(limited) {
            D("Reached precision limit\n");
            break;
        }

        /* Precalculate log(q). */
        for(col= 0; col < Q->ncol; col++) {
            if(q[col] == 0) limited= 1;
            logq[col]= log2f(q[col]);
        }
        if(limited) {
            D("Reached precision limit\n");
            break;
        }

        /* Calculate KL divergence. */
        bzero(z, Q->nrow * sizeof(float));
        for(col= 0; col < Q->ncol; col++) {
            int i;

            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                row= Q->rows[i];

                assert(Q->entries[i] > 0);

                /* Use the log tables. */
                z[row]+= Q->entries[i] * (logQ[i] - logq[col]);
            }
        }

        /* Find the largest element of c.  This is Iu, as log is monotonic. */
        Iu= -INFINITY;
        for(col= 0; col < Q->nrow; col++) {
            if(z[col] > Iu) Iu= z[col];
        }

        /* Update Il.  By taking out the largest element of c, the
         * rest are all less than 1. */
        tmp= 0.0;
        for(row= 0; row < Q->nrow; row++) {
            tmp+= p[row] * exp2f(z[row] - Iu);
        }
        Il= Iu + log2f(tmp);

        e= Iu - Il;

        if(e < e_best) {
            Il_best= Il; e_best= e;
        }

        D("Il= %.6e Iu= %.6e e= %.6e\n", Il, Iu, e);
        if(e < epsilon) {
            D("Terminated in single-precision phase\n");
            break;
        }
        if(e >= e_last) {
            if(iterations > GRACE_PERIOD) {
                increases++;
                if(increases > NONMONOTONIC_TOLERANCE) {
                    D("Reached precision limit\n");
                    break;
                }
            }
        }
        iterations++;

        e_last= e;

        /* Scale the input probabilities. */
        tmp= 0.0;
        for(row= 0; row < Q->nrow; row++) {
            p[row]*= exp2f(lambda * z[row]);
            tmp+= p[row];
        }
        for(row= 0; row < Q->nrow; row++)
            p[row]/= tmp;
    } 

    free(logq);
    free(logQ);
    free(z);
    free(q);

    if(e_obs) *e_obs= e_best;
    return Il_best;
}

/* The slower double-precision phase. */
float
double_phase(csc_mat_t *Q, float epsilon, float *e_obs, double *p,
        float lambda, float Il_best_prev, float e_best_prev) {
    double *q, *z;
    double Il, Iu, e, e_last= INFINITY;
    double Il_best= Il_best_prev, e_best= e_best_prev;
    int col, row;
    double *logQ, *logq;
    double tmp;
    int increases= 0, iterations= 0;

    D("Double-precision phase\n");

    q= malloc(Q->ncol * sizeof(double));
    z= malloc(Q->nrow * sizeof(double));
    if(!q || !z) { perror("malloc"); abort(); }

    D("Precalculating logs...");
    logQ= malloc(Q->nnz * sizeof(double));
    if(!logQ) { perror("malloc"); abort(); }
    {
        int i;
        for(i= 0; i < Q->nnz; i++) logQ[i]= log2(Q->entries[i]);
    }
    D(" done.\n");

    logq= malloc(Q->ncol * sizeof(double));
    if(!logq) { perror("malloc"); abort(); }

    while(1) {
        int limited= 0;

        /* Calculate q - output distribution. */
        for(col= 0; col < Q->ncol; col++) {
            int i;

            q[col]= 0.0;
            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                if(p[Q->rows[i]] == 1) limited= 1;
                q[col]+= Q->entries[i] * p[Q->rows[i]];
            }
        }
        if(limited) {
            D("1 Reached precision limit\n");
            break;
        }

        /* Precalculate log(q). */
        for(col= 0; col < Q->ncol; col++) {
            if(q[col] == 0) limited= 1;
            logq[col]= log2(q[col]);
        }
        if(limited) {
            D("2 Reached precision limit\n");
            break;
        }

        /* Calculate KL divergence. */
        bzero(z, Q->nrow * sizeof(double));
        for(col= 0; col < Q->ncol; col++) {
            int i;

            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                row= Q->rows[i];

                assert(Q->entries[i] > 0);

                /* Use the log tables. */
                z[row]+= (double)Q->entries[i] * (logQ[i] - logq[col]);

            }
        }

        /* Find the largest element of c.  This is Iu, as log is monotonic. */
        Iu= -INFINITY;
        for(col= 0; col < Q->nrow; col++) {
            if(z[col] > Iu) Iu= z[col];
        }

        /* Update Il.  By taking out the largest element of c, the
         * rest are all less than 1. */
        tmp= 0.0;
        for(row= 0; row < Q->nrow; row++) {
            tmp+= p[row] * exp2(z[row] - Iu);
        }
        Il= Iu + log2(tmp);

        e= Iu - Il;

        if(e < e_best) {
            Il_best= Il; e_best= e;
        }

        D("Il= %.6le Iu= %.6le e= %.6e\n", Il, Iu, e);
        if(e < epsilon) {
            D("Terminated in double-precision phase\n");
            break;
        }
        if(e >= e_last) {
            if(iterations > GRACE_PERIOD) {
                increases++;
                if(increases > NONMONOTONIC_TOLERANCE) {
                    D("Reached precision limit\n");
                    break;
                }
            }
        }
        iterations++;

        e_last= e;

        /* Find the denominator of the probability scale equation. */
        tmp= 0.0;
        for(row= 0; row < Q->nrow; row++)
            tmp+= p[row] * exp2(lambda * z[row]);

        /* Scale the input probabilities. */
        for(row= 0; row < Q->nrow; row++)
            p[row]= p[row] * exp2(lambda * z[row]) / tmp;
    } 

    free(logq);
    free(logQ);
    free(z);
    free(q);

    if(e_obs) *e_obs= e_best;
    return Il_best;
}

/* The slowest, long-double, phase. */
float
long_double_phase(csc_mat_t *Q, float epsilon, float *e_obs,
        long double *p, float lambda, float Il_best_prev, float e_best_prev) {
    long double *q, *z;
    long double Il, Iu, e, e_last= INFINITY;
    long double Il_best= Il_best_prev, e_best= e_best_prev;
    int col, row;
    long double *logQ, *logq;
    long double tmp;
    int increases= 0, iterations= 0;

    D("Extended-precision phase\n");

    q= malloc(Q->ncol * sizeof(long double));
    z= malloc(Q->nrow * sizeof(long double));
    if(!q || !z) { perror("malloc"); abort(); }

    D("Precalculating logs...");
    logQ= malloc(Q->nnz * sizeof(long double));
    if(!logQ) { perror("malloc"); abort(); }
    {
        int i;
        for(i= 0; i < Q->nnz; i++) logQ[i]= log2l(Q->entries[i]);
    }
    D(" done.\n");

    logq= malloc(Q->ncol * sizeof(long double));
    if(!logq) { perror("malloc"); abort(); }

    while(1) {
        int limited= 0;

        /* Calculate q - output distribution. */
        for(col= 0; col < Q->ncol; col++) {
            int i;

            q[col]= 0.0;
            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                if(p[Q->rows[i]] == 0) limited= 1;
                q[col]+= Q->entries[i] * p[Q->rows[i]];
            }
        }
        if(limited) {
            D("1 Reached precision limit\n");
            break;
        }

        /* Precalculate log(q). */
        for(col= 0; col < Q->ncol; col++) {
            if(q[col] == 0) limited= 1;
            logq[col]= log2l(q[col]);
        }
        if(limited) {
            D("2 Reached precision limit\n");
            break;
        }

        /* Calculate KL divergence. */
        bzero(z, Q->nrow * sizeof(long double));
        for(col= 0; col < Q->ncol; col++) {
            int i;

            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                row= Q->rows[i];

                assert(Q->entries[i] > 0);

                /* Use the log tables. */
                z[row]+= (long double)Q->entries[i] * (logQ[i] - logq[col]);

            }
        }

        /* Find the largest element of c.  This is Iu, as log is monotonic. */
        Iu= -INFINITY;
        for(col= 0; col < Q->nrow; col++) {
            if(z[col] > Iu) Iu= z[col];
        }

        /* Update Il.  By taking out the largest element of c, the
         * rest are all less than 1. */
        tmp= 0.0;
        for(row= 0; row < Q->nrow; row++) {
            tmp+= p[row] * exp2l(z[row] - Iu);
        }
        Il= Iu + log2l(tmp);

        e= Iu - Il;

        if(e < e_best) {
            Il_best= Il; e_best= e;
        }

        D("Il= %.6le Iu= %.6le e= %.6e\n", (double)Il, (double)Iu, (double)e);
        if(e < epsilon) {
            D("Terminated in Extended-precision phase\n");
            break;
        }
        if(e >= e_last) {
            if(iterations > GRACE_PERIOD) {
                increases++;
                if(increases > NONMONOTONIC_TOLERANCE) {
                    D("Reached precision limit\n");
                    break;
                }
            }
        }
        iterations++;

        e_last= e;

        /* Find the denominator of the probability scale equation. */
        tmp= 0.0;
        for(row= 0; row < Q->nrow; row++)
            tmp+= p[row] * exp2l(lambda * z[row]);

        /* Scale the input probabilities. */
        for(row= 0; row < Q->nrow; row++)
            p[row]= p[row] * exp2l(lambda * z[row]) / tmp;
    } 

    free(logq);
    free(logQ);
    free(z);
    free(q);

    if(e_obs) *e_obs= e_best;
    return Il_best;
}

/* The phased implementation. */
float
ba_phased(csc_mat_t *Q, float epsilon, float *e_obs) {
    float *ps;
    double *pd;
    long double *pld;
    float som= 0.0;
    float lambda;
    float e_phase= INFINITY, Il;
    int col, i;

    ps= malloc(Q->nrow * sizeof(float));
    if(!ps) { perror("malloc"); abort(); }

    /* Find the greatest accelerating/squeezing factor guaranteeing monotone
     * convergence. */
    for(col= 0; col < Q->ncol; col++) {
        float col_min= INFINITY;

        /* If there's a zero in this column, the minimum is zero. */
        if(Q->ci[col+1] - Q->ci[col+1] < Q->nrow) continue;

        for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
            if(Q->entries[i] < col_min)
                col_min= Q->entries[i];
        }

        som+= col_min;
    }

    if(som == INFINITY || som == 1.0)
        lambda= 1.0;
    else
        lambda= 1.0 / (1.0 - som);
    D("som= %.e\n", som);
    D("lambda= %.e\n", lambda);

    /* Start with a uniform input distribution. */
    for(i= 0; i < Q->nrow; i++) ps[i]= 1.0 / Q->nrow;

    /* Start with single precision. */
    Il= single_phase(Q, epsilon, &e_phase, ps, lambda);

    /* Terminate if we've already hit the target. */
    if(e_phase < epsilon) {
        *e_obs= e_phase;
        free(ps);
        return Il;
    }

    /* Otherwise, increase the precision. */
    pd= malloc(Q->nrow * sizeof(double));
    if(!pd) { perror("malloc"); abort(); }
    {
        double sum= 0.0;
        for(i= 0; i < Q->nrow; i++) {
            /* If we had an underflow, bump it up a little. */
            if(ps[i] == 0.0)
                pd[i]= ((double)JIGGER) / Q->nrow;
            else
                pd[i]= (double)ps[i];
            sum+= pd[i];
        }
        for(i= 0; i < Q->nrow; i++) pd[i]/= sum;
    }
    free(ps);

    /* Continue in double precision. */
    Il= double_phase(Q, epsilon, &e_phase, pd, lambda, Il, e_phase);

    /* Again, terminate if good enough. */
    if(e_phase < epsilon) {
        *e_obs= e_phase;
        free(pd);
        return Il;
    }

    /* Otherwise increase precision again. */
    pld= malloc(Q->nrow * sizeof(long double));
    if(!pld) { perror("malloc"); abort(); }
    {
        long double sum= 0.0;
        for(i= 0; i < Q->nrow; i++) {
            if(pd[i] == 0.0)
                pld[i]= ((long double)JIGGER) / Q->nrow;
            else
                pld[i]= (long double)pd[i];
            sum+= pld[i];
        }
        for(i= 0; i < Q->nrow; i++) pld[i]/= sum;
    }
    free(pd);

    /* Our last shot is long double precision, which may be no better than
     * double, depending on platform, but is usually slow. */
    Il= long_double_phase(Q, epsilon, &e_phase, pld, lambda, Il, e_phase);

    free(pld);

    *e_obs= e_phase;
    return Il;
}
