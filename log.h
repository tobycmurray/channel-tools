#ifndef __LOG_H
#define __LOG_H

#include <x86intrin.h>

void write_log_table(void);
float log2f_table(float x);

#ifdef __SSE2__
__v4sf log2f_table_v(__v4sf x);
#endif

#ifdef __AVX__
__v8sf log2f_table_v8(__v8sf x);
#endif

#endif
