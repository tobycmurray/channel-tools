/* sample_error.c

   Run a monte-carlo simulation to approximate the distribution of capacities
   expected if the given matrix had zero bandwidth.  Used to generate
   confidence intervals.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "log.h"
#include "sparse.h"
#include "dSFMT-src-2.2.1/dSFMT.h"
#include "channel_algorithms.h"

/* Sample from the given distribution. */
int
sample(dv_t *cp, dsfmt_t *rng) {
    /* Choose x in [0,1). */
    float x= dsfmt_genrand_close_open(rng);
    int s= 0, e= cp->length;

    /* Binary search for last point with cumulative probability <= x. */
    while(s+1 < e) {
        int m= (s+e)/2;

        if(x < cp->entries[m]) e= m;
        else                   s= m;
    }

    return s;
}

/* Generate a matrix by taking 'samples' samples per column from the given
 * distributions, using lcp for the left-hand side and rcp for the right. */
csc_mat_t *
sampled_matrix(dv_t *lcp, dv_t *rcp, int rows, int samples, dsfmt_t *rng) {
    bsc_hist_t *H;
    csc_mat_t *M;
    int r;

    H= bsc_hist_new();
    for(r= 0; r < rows; r++) {
        int i;

        for(i= 0; i < samples; i++) {
            int c;
            if(r < rows/2)
                c= sample(lcp, rng);
            else
                c= sample(rcp, rng);
            bsc_hist_count(H, c, r, 1);
        }
    }
    M= bsc_normalise(H);
    bsc_hist_destroy(H);
    return M;
}

/* Average along the columns of the matrix. */
dv_t *
average_cols(csc_mat_t *M) {
    dv_t *v;
    int c;
    float Ps= 0.0;

    v= dv_new(M->ncol);
    if(!v) {
        perror("dv_new");
        exit(EXIT_FAILURE);
    }
    dv_zero(v);

    /* Sum along every column. */
    for(c= 0; c < M->ncol; c++) {
        int64_t i;

        for(i= M->ci[c]; i < M->ci[c+1]; i++) {
            if(M->entries[i] == 0) continue;
            v->entries[c]+= M->entries[i];
        }
    }

    /* Rescale to ensure total prob = 1 */
    for(c= 0; c < v->length; c++) Ps+= v->entries[c];
    for(c= 0; c < v->length; c++) v->entries[c]/= Ps;

    return v;
}

/* Calculate the cumulative probability. */
dv_t *
accumulate_prob(dv_t *p) {
    dv_t *cp;
    float Pc= 0.0;
    int i;

    cp= dv_new(p->length);
    if(!cp) {
        perror("dv_new");
        exit(EXIT_FAILURE);
    }

    for(i= 0; i < p->length; i++) {
        Pc+= p->entries[i];
        cp->entries[i]= Pc;
    }

    return cp;
}

/* Don't let the workers trample each other's output. */
pthread_mutex_t output_lock= PTHREAD_MUTEX_INITIALIZER;

/* Generate n noisy matrices, and print their capacity. */
void
noisy_matrices(float cap, dv_t *lcp, dv_t *rcp, float epsilon, int n,
        int rows, int samples, dsfmt_t *rng) {
    int i;

    for(i= 0; i < n; i++) {
        /* Draw the samples. */
        csc_mat_t *M= sampled_matrix(lcp, rcp, rows, samples, rng);
        float I, e;
        /* We may have empty columns, which the ABA doesn't like. */
        csc_prune_cols(M);
        /* Find the capacity (mutual information). */
        I= ba_phased(M, epsilon, &e);
        /* Deallocate the matrix. */
        csc_mat_destroy(M);
        /* Print the capacity. */
        pthread_mutex_lock(&output_lock);
        printf("%.12e %.12e %.12e\n", cap, I, e);
        /* Very important, as the simulation is long-running and may be
         * interrupted.  You don't want to lose an hour's calculation in the
         * buffer. */
        fflush(stdout);
        pthread_mutex_unlock(&output_lock);
    }
}

/* A worker thread's context. */
struct job {
    float cap;
    dv_t *lcp;
    dv_t *rcp;
    float epsilon;
    int n;
    int rows;
    int samples;
    dsfmt_t rng;
    pthread_t thread;
};

/* Worker threads simply call noisy_matrices. */
static void *
worker(void *arg) {
    struct job *j= (struct job *)arg;
    noisy_matrices(j->cap, j->lcp, j->rcp, j->epsilon, j->n, j->rows,
            j->samples, &j->rng);
    return (void *)0;
}

/* Take a distribution and cut it in half. */
void
split_prob(dv_t *prob, dv_t **top_half, dv_t **bot_half) {
    int i;

    assert(top_half);
    assert(bot_half);

    *top_half= dv_new(prob->length);
    *bot_half= dv_new(prob->length);
    if(!top_half || !bot_half) {
        perror("dv_new");
        exit(EXIT_FAILURE);
    }

    /* Put the first n/2 values into bot_half. */
    dv_zero(*bot_half);
    for(i= 0; i < prob->length/2; i++) {
        (*bot_half)->entries[i]= prob->entries[i];
    }

    /* Put the second n/2 values into top_half. */
    dv_zero(*top_half);
    for(i= prob->length/2; i < prob->length; i++) {
        (*top_half)->entries[i]= prob->entries[i];
    }

    /* We now have two completely disjoint sub-distributions.  Added together,
     * they'll give the original. */
}

/* Blend two distributions. */
dv_t *
join_prob(dv_t *left, dv_t *right, float alpha) {
    int i;
    float Ps= 0.0;
    dv_t *prob;

    assert(left->length == left->length);
    assert(0 <= alpha);
    assert(alpha <= 1);

    prob= dv_new(left->length);
    if(!prob) {
        perror("dv_new");
        exit(EXIT_FAILURE);
    }

    /* Put the weighted sum of the sub-distributions into prob. */
    for(i= 0; i < prob->length; i++) {
        prob->entries[i]=
            alpha * left->entries[i] + (1-alpha) * right->entries[i];
        Ps+= prob->entries[i];
    }

    /* Rescale */
    for(i= 0; i < prob->length; i++) prob->entries[i]/= Ps;

    return prob;
}

float
plogp(float p) {
    assert(0.0 <= p);
    assert(p <= 1.0);
    if(p == 0.0) return 0.0;
    return p * log2f(p);
}

float
h(float p) {
    return -(plogp(p) + plogp(1 - p));
}

/* Find the mutual information of the two-column matrix [Pa:Pb], given input
 * distribution (pa, 1-pa). */
float
I(dv_t *Pa, dv_t *Pb, float pa) {
    /* Initial entropy */
    float Hi= h(pa);
    /* Expected final entropy. */
    float EHf= 0.0;
    int o;

    /* For all outputs. */
    for(o= 0; o < Pa->length; o++) {
        /* The probability of the output. */
        float po= pa * Pa->entries[o] + (1 - pa) * Pb->entries[o];
        /* The posterior distribution (Bayes' rule). */
        float pa_o= (Pa->entries[o] / po) * pa;
        /* The weighted entropy of the posterior. */
        EHf+= po * h(pa_o);
    }

    /* The mutual information. */
    return Hi - EHf > 0 ? Hi - EHf : 0;
}

/* Find the capacity by maximising over pa. */
float
find_binary_cap(dv_t *Pa, dv_t *Pb, float epsilon) {
    float Il, Iu, pa;
    int n= Pa->length;

    assert(Pa->length == Pb->length);

    pa= 0.5;
    Il= 0.0;
    Iu= 1.0;

    while(Iu - Il > epsilon) {
        float ca, cb;
        float lg_ca, lg_cb;
        int i;

        //printf(">> %e %e %e\n", pa, Il, Iu); fflush(stdout);

        lg_ca= 0.0;
        lg_cb= 0.0;
        for(i= 0; i < n; i++) {
            float po= pa * Pa->entries[i] + (1-pa) * Pb->entries[i];
            lg_ca+= Pa->entries[i] * log2f(Pa->entries[i] / po);
            lg_cb+= Pb->entries[i] * log2f(Pb->entries[i] / po);
        }
        ca= exp2f(lg_ca);
        cb= exp2f(lg_cb);

        Il= log2f(pa * ca + (1-pa) * cb);
        if(ca > cb) Iu= log2f(ca);
        else        Iu= log2f(cb);

        pa*= ca / (pa * ca + (1-pa) * cb);
    }

    //printf(">> %e %e %e\n", pa, Il, Iu); fflush(stdout);

    return Il;
}

int
main(int argc, char *argv[]) {
    FILE *in;
    csc_mat_t *M;
    csc_errno_t e;
    dv_t *avg_prob;
    dv_t *top_half, *bot_half;
    int samples, runs, runs_per_job;
    float epsilon;
    int cap_steps, start_step, end_step;
    int nthreads;
    struct job *jobs;
    int seed;
    int quiet= 0;
    int i, s;

    srandom(time(NULL));

    if(argc < 8) {
        fprintf(stderr,
                "Usage: %s <channel matrix> <samples per column> "
                "<runs> <max error> <capacity steps> <start step> "
                "<end step> [-q]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    in= fopen(argv[1], "rb");
    if(!in) { perror("fopen"); exit(EXIT_FAILURE); }

    samples= atoi(argv[2]);
    runs= atoi(argv[3]);
    epsilon= atof(argv[4]);
    cap_steps= atoi(argv[5]);
    start_step= atoi(argv[6]);
    end_step= atoi(argv[7]);

    if(argc > 8 && !strcmp(argv[8], "-q"))
        quiet= 1;

    M= csc_load_binary(in, &e);
    if(!M) { csc_perror(e, "csc_load_binary"); exit(EXIT_FAILURE); }
    fclose(in);

    if(!csc_check(M, 1)) abort();
    csc_prune_cols(M);
    if(!quiet) csc_stats(M);

    if(STRIDE_OF(M) > 1) {
        fprintf(stderr, "Can't handle a strided matrix\n");
        exit(EXIT_FAILURE);
    }

    if(!quiet) fprintf(stderr, "Averaging matrix columns...");
    if(!quiet) fflush(stderr);
    avg_prob= average_cols(M);
    split_prob(avg_prob, &top_half, &bot_half);
    if(!quiet) fprintf(stderr, " done.\n");

    nthreads= sysconf(_SC_NPROCESSORS_ONLN);
    if(nthreads < 1) nthreads= 1;

    write_log_table();

    jobs= malloc(nthreads * sizeof(struct job));
    if(!jobs) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    if(!quiet) fprintf(stderr, "Generating noisy matrices...");
    if(!quiet) fflush(stderr);

    for(s= start_step; s < end_step; s++) {
        dv_t *left_prob, *right_prob,
             *left_cprob, *right_cprob;
        float alpha= 0.5 + s * (0.5 / cap_steps);
        float cap;
        int runs_step= runs;

        left_prob= join_prob(top_half, bot_half, alpha);
        right_prob= join_prob(top_half, bot_half, 1.0 - alpha);
        left_cprob= accumulate_prob(left_prob);
        right_cprob= accumulate_prob(right_prob);

        /* Find the capacity of the step matrix. */
        cap= find_binary_cap(left_prob, right_prob, epsilon);

        runs_per_job= runs_step / nthreads;
        for(i= 0; i < nthreads; i++) {
            jobs[i].cap= cap;
            jobs[i].lcp= left_cprob;
            jobs[i].rcp= right_cprob;
            jobs[i].epsilon= epsilon;
            jobs[i].n= runs_per_job;
            jobs[i].rows= M->nrow;
            jobs[i].samples= samples;
            seed= random();
            dsfmt_init_gen_rand(&jobs[i].rng, seed);
            runs_step-= runs_per_job;
        }
        jobs[nthreads-1].n+= runs_step;
        for(i= 0; i < nthreads; i++) {
            if(pthread_create(&jobs[i].thread, NULL, worker, &jobs[i])) {
                perror("pthread_create");
                exit(EXIT_FAILURE);
            }
        }
        for(i= 0; i < nthreads; i++) {
            if(pthread_join(jobs[i].thread, NULL)) {
                perror("pthread_join");
                exit(EXIT_FAILURE);
            }
        }

        dv_destroy(left_prob);
        dv_destroy(right_prob);
        dv_destroy(left_cprob);
        dv_destroy(right_cprob);
    }
    if(!quiet) fprintf(stderr, " done.\n");

    dv_destroy(avg_prob);
    dv_destroy(top_half);
    dv_destroy(bot_half);
    csc_mat_destroy(M);

    return 0;
}
