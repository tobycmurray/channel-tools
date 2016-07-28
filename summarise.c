/* summarise.c

   Calculate statistics for the input samples.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXLINE 1024

int
main(int argc, char *argv[]) {
    int r, c;
    int min, max;
    int res_min, res_max;
    int rmin, rmax;
    int cmin, cmax;
    size_t *counts, total_count, out_of_range, malformed;
    FILE *oor_log, *mal_log;
    char buf[MAXLINE];

    if(argc < 7) {
        fprintf(stderr, "Usage: %s <input min> <input max> "
                "<output min> <output max> <out of range "
                "logfile> <malformed logfile>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    min= atoi(argv[1]);
    max= atoi(argv[2]);
    res_min= atoi(argv[3]);
    res_max= atoi(argv[4]);
    oor_log= fopen(argv[5], "w");
    mal_log= fopen(argv[6], "w");
    if(!oor_log || !mal_log) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    if(min < 0 || max < 0 || max < min) {
        fprintf(stderr, "Invalid range: %d - %d\n", min, max);
        exit(EXIT_FAILURE);
    }

    counts= calloc(max - min + 1, sizeof(size_t));
    if(!counts) {
        perror("calloc");
        exit(EXIT_FAILURE);
    }
    total_count= 0;
    out_of_range= 0;
    malformed= 0;
    rmin= INT_MAX;
    rmax= INT_MIN;
    cmin= INT_MAX;
    cmax= INT_MIN;

    while(fgets(buf, MAXLINE, stdin)) {
        int n;

        n= sscanf(buf, "%d %d\n", &c, &r);
        if(n != 2) {
            malformed++;
            fprintf(mal_log, "%s\n", buf);
            continue;
        }

        total_count++;

        if(c < cmin) cmin= c;
        if(cmax < c) cmax= c;

        if(c < min || max < c) {
            out_of_range++;
            fprintf(oor_log, "%s", buf);
            continue;
        }

        if(r < rmin) rmin= r;
        if(rmax < r) rmax= r;

        if(r < res_min || res_max < r) {
            out_of_range++;
            fprintf(oor_log, "%s", buf);
            continue;
        }

        counts[c - min]++;
    }

    printf("%lu %lu %lu %d %d %d %d", total_count, out_of_range, malformed,
            rmin, rmax, cmin, cmax);
    for(c= 0; c <= max - min; c++)
        printf(" %lu", counts[c]);
    printf("\n");

    fclose(mal_log);
    fclose(oor_log);

    return 0;
}
