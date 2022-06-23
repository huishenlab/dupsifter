/* Utility functions for dupsifter
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2022 Jacob.Morrison@vai.org
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _DS_UTILS_H
#define _DS_UTILS_H

#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "sam.h"

// Is a primary alignment (i.e., not secondary or supplementary)
static inline uint8_t is_primary(bam1_t *b) {
    bam1_core_t *c = &b->core;
    if (!((c->flag & BAM_FSECONDARY) || (c->flag & BAM_FSUPPLEMENTARY))) {
        return 1;
    }

    return 0;
}

// Read is flagged as unmapped
static inline uint8_t is_unmapped(bam1_t *b) {
    bam1_core_t *c = &b->core;
    if (c->flag & BAM_FUNMAP) {
        return 1;
    }

    return 0;
}

// Read is flagged as PCR duplicate
static inline uint8_t is_duplicate(bam1_t *b) {
    bam1_core_t *c = &b->core;
    if (c->flag & BAM_FDUP) {
        return 1;
    }

    return 0;
}

// Read is flagged as read 1
static inline uint8_t is_first_read(bam1_t *b) {
    bam1_core_t *c = &b->core;
    if (c->flag & BAM_FREAD1) {
        return 1;
    }

    return 0;
}

// Read is flagged as read 2
static inline uint8_t is_second_read(bam1_t *b) {
    bam1_core_t *c = &b->core;
    if (c->flag & BAM_FREAD2) {
        return 1;
    }

    return 0;
}

// Find sum of base qualities for read
static inline uint32_t total_qual(const bam1_t *b) {
    if (!b) { return 0; }

    uint32_t total = 0;
    uint8_t *qual = bam_get_qual(b);

    int i;
    for (i=0; i<b->core.l_qseq; ++i) { total += qual[i]; }

    return total;
}

// Exit with error message and EXIT_FAILURE return
// Pulled from wzfatal() in zwdzwd/utils/wzmisc.h on GitHub
static inline void fatal_error(const char *err, ...) {
    va_list args;
    va_start(args, err);
    vfprintf(stderr, err, args);
    va_end(args);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

// Timing functions
// CPU time usage
static inline double get_cpu_runtime() {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);

    double secs = r.ru_utime.tv_sec + r.ru_stime.tv_sec;
    double msec = 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);

    return secs + msec;
}

// Walltime
static inline double get_current_time() {
    struct timeval  tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);

    return tp.tv_sec + 1e-6*tp.tv_usec;
}

#endif /* _DS_UTILS_H */
