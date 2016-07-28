/* log.c

   Table-based single-precision logarithm.

   This code is experimental, and error-handling is primitive.
*/

/* Copyright 2013, NICTA.  See COPYRIGHT for license details. */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include <immintrin.h>

#define TABLE_BITS 23
#define TABLE_ENTRIES (1<<(TABLE_BITS))
#define TABLE_BASE ((float)1.0)
#define TABLE_INCREMENT ((float)1.0/TABLE_ENTRIES)

float log_table[TABLE_ENTRIES];

void
write_log_table(void) {
    int i;

    for(i= 0; i < TABLE_ENTRIES; i++)
        log_table[i]= log2f(TABLE_BASE + i*TABLE_INCREMENT);
}

float
log2f_table(float x) {
    __v4sf y;
    uint32_t raw;
    int exp, mantissa;
    float v;

    /* Get the float as a uint32. */
    y[0]= x;
    raw= ((__v4si) y)[0];

    /* Extract mantissa and exponent. */
    exp= (raw >> 23) & 0xff;
    mantissa= raw & 0x7fffff;

    /* Lookup the base value using the most-significant bits. */
    v= log_table[mantissa >> (23 - TABLE_BITS)];

    /* Add in the (de-biased) exponent. */
    return v + (exp - 127);
}

#ifdef __SSE2__
__v4sf
log2f_table_v(__v4sf x) {
    __v4si raw;
    __v4si exp, mantissa;
    __v4sf v;

    /* Get the float as a uint32. */
    raw= (__v4si)x;

    /* Extract mantissa and exponent. */
    exp= ((raw >> 23) & 0xff) - 127;
    mantissa= raw & 0x7fffff;

    /* Lookup the base value using the most-significant bits. */
    v[0]= log_table[mantissa[0] >> (23 - TABLE_BITS)];
    v[1]= log_table[mantissa[1] >> (23 - TABLE_BITS)];
    v[2]= log_table[mantissa[2] >> (23 - TABLE_BITS)];
    v[3]= log_table[mantissa[3] >> (23 - TABLE_BITS)];

    /* Add in the (de-biased) exponent. */
    return v + __builtin_ia32_cvtdq2ps(exp);
}
#endif /* __SSE2__ */

#ifdef __AVX__
/* This is awful on AVX1.  They seem to have left all the useful
   instructions for AVX2! */
__v8sf
log2f_table_v8(__v8sf x) {
    __v8si raw;
    __v8si exp, mantissa;
    __v8sf v;

    /* Get the float as a uint32. */
    raw= (__v8si)x;

    /* Extract mantissa and exponent. */
    exp= ((raw >> 23) & 0xff) - 127;
    mantissa= raw & 0x7fffff;

    /* Lookup the base value using the most-significant bits. */
    v[0]= log_table[mantissa[0] >> (23 - TABLE_BITS)];
    v[1]= log_table[mantissa[1] >> (23 - TABLE_BITS)];
    v[2]= log_table[mantissa[2] >> (23 - TABLE_BITS)];
    v[3]= log_table[mantissa[3] >> (23 - TABLE_BITS)];
    v[4]= log_table[mantissa[4] >> (23 - TABLE_BITS)];
    v[5]= log_table[mantissa[5] >> (23 - TABLE_BITS)];
    v[6]= log_table[mantissa[6] >> (23 - TABLE_BITS)];
    v[7]= log_table[mantissa[7] >> (23 - TABLE_BITS)];

    /* Add in the (de-biased) exponent. */
    return v + __builtin_ia32_cvtdq2ps256(exp);
}
#endif /* __AVX__ */
