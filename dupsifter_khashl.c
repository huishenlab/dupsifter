/* mark duplicates in bisulfite-converted reads
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2016 Wanding.Zhou@vai.org
 *               2022 Jacob.Morrison@vai.org
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
 *
 *
 * This utility is largely based on samblaster's methodology for marking
 * duplicates, as it is fast and memory efficient. However, there are some
 * cytosine-conversion-aware things that must be taken into account that
 * samblaster does not take into account.
 *
 * Original Header from Samblaster samblaster.cpp
 *
 * -*- mode: C++ ; indent-tabs-mode: nil ; c-file-style: "stroustrup" -*-
 *
 *  Project: samblaster
 *           Fast mark duplicates in read-ID grouped SAM file.
 *           Also, optionally pull discordants, splitters, and/or unmappend/clipped reads.
 *  Author:  Greg Faust (gf4ea@virginia.edu)
 *  Date:    October 2013
 *  Cite:    SAMBLASTER: fast duplicate marking and structural variant read extraction
 *           GG Faust, IM Hall
 *           Bioinformatics 30 (17), 2503-2505
 *           https://academic.oup.com/bioinformatics/article/30/17/2503/2748175
 *
 *  File:    samblaster.cpp  code file for the main routine and most of the other code.
 *
 *  License Information:
 *
 *  Copyright 2013-2020 Gregory G. Faust
 *
 *  Licensed under the MIT license (the "License");
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at http://opensource.org/licenses/MIT
 *
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include "sam.h"
// TODO: Test whether khash or khashl is better to use
#include "khashl.h"
#include "util.h"
#include "refcache.h"

//---------------------------------------------------------------------------------------------------------------------
// Config struct and initialization function
typedef struct {
    char     *reffn;               /* reference file name */
    char     *infn;                /* name of input file */
    char     *outfn;               /* name of output file */
    char      out_mode[6];         /* write mode of output file */
    uint32_t  min_base_qual;       /* threshold for high quality bases */
    uint8_t   rm_dup;              /* flag to remove duplicates */
    uint32_t  cnt_pe_dup;          /* number of PE reads marked as dup */
    uint32_t  cnt_pe;              /* number of PE reads */
    uint32_t  cnt_se_dup;          /* number of SE reads marked as dup */
    uint32_t  cnt_se;              /* number of SE reads */
    uint32_t  cnt_dangle;          /* number of dangling PE reads */
    uint8_t   verbose;             /* level of messages to print */
    uint32_t  max_length;          /* max read length allowed */
    uint8_t   single_end;          /* process all reads as single-end */
} mkdup_conf_t;

void mkdup_conf_init(mkdup_conf_t *conf) {
    strcpy(conf->out_mode, "w");

    conf->outfn = (char *)"-";
    conf->min_base_qual = 0;
    conf->rm_dup = 0;
    conf->cnt_pe_dup = 0;
    conf->cnt_pe = 0;
    conf->cnt_se_dup = 0;
    conf->cnt_se = 0;
    conf->cnt_dangle = 0;
    conf->verbose = 0;
    conf->max_length = 10000;
    conf->single_end = 0;
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
// Structs for parsing CIGAR string and the signature, plus any init functions that are needed
typedef struct {
    uint32_t rlen; /* reference length */
    uint32_t qlen; /* query (read) length */
    uint32_t sclp; /* soft/hard clips at start of read */
    uint32_t eclp; /* soft/hard clips at end of read */
} parsed_cigar_t;

typedef struct {
    uint16_t bin_start; /* bin number for starting position, taken from forward strand read */
    uint16_t bin_end;   /* bin number for ending position, taken from reverse strand read */
    uint32_t pos_start; /* bin position for starting position, taken from forward strand read */
    uint32_t pos_end;   /* bin position for ending position, taken from reverse strand read */
} signature_t;

signature_t signature_init() {
    // Chances of the bin or position being the max value are very low, so set default value to max value
    signature_t sig = {0};
    sig.bin_start = UINT16_MAX;
    sig.bin_end   = UINT16_MAX;
    sig.pos_start = UINT32_MAX;
    sig.pos_end   = UINT32_MAX;

    return sig;
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
// Initialize hash maps and define functions needed for initialization
khint_t hash_sig(signature_t s) {
    if (s.pos_start == UINT32_MAX) {
        return kh_hash_uint32(s.pos_end);
    } else {
        return kh_hash_uint32(s.pos_start);
    }
}

int sig_equal(signature_t s1, signature_t s2) {
    if ((s1.bin_start == s2.bin_start) &&
        (s1.pos_start == s2.pos_start) &&
        (s1.bin_end   == s2.bin_end  ) &&
        (s1.pos_end   == s2.pos_end  )) {
        return 1;
    } else {
        return 0;
    }
}

KHASHL_SET_INIT(, SigHash, sh, signature_t, hash_sig, sig_equal)
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
// Genome bins/Hash map struct, init and destroy functions
typedef struct {
    uint32_t  n_contigs;    /* number of contigs */
    uint64_t  total_length; /* total length of genome */
    uint32_t *length;       /* length of contigs */
    uint64_t *offset;       /* offset from start of contigs */
    uint32_t  n_bins;       /* number of bins in genome */
} mkdup_bins_t;

mkdup_bins_t *mkdup_bins_init() {
    return (mkdup_bins_t *)calloc(1, sizeof(mkdup_bins_t));
}

void destroy_mkdup_bins(mkdup_bins_t *b) {
    if (b->length != NULL) { free(b->length); }
    if (b->offset != NULL) { free(b->offset); }
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
// Padding functions
//
// These functions account for any reads that may run off the start or end of a chromosome, which may happen due to
// soft clipped reads
static inline uint32_t pad_chrom_length(uint32_t length, uint32_t max_length) {
    return length + (2 * max_length);
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
// Bisulfite strand functions
uint8_t infer_bsstrand(refcache_t *rs, bam1_t *b, uint32_t min_base_qual) {
    bam1_core_t *c = &b->core;
    uint32_t rpos = c->pos+1, qpos = 0, nC2T = 0, nG2A = 0;
    uint32_t i, op, oplen;
    char rb, qb;
    uint32_t j;
    for (i=0; i<c->n_cigar; ++i) {
        op    = bam_cigar_op(bam_get_cigar(b)[i]);
        oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
        switch(op) {
            case BAM_CMATCH:
                for (j=0; j<oplen; ++j) {
                    rb = refcache_getbase_upcase(rs, rpos+j);
                    qb = bscall(b, qpos+j);
                    if (bam_get_qual(b)[qpos+j] < min_base_qual) { continue; }
                    if (rb == 'C' && qb == 'T') { nC2T++; }
                    if (rb == 'G' && qb == 'A') { nG2A++; }
                }
                rpos += oplen;
                qpos += oplen;
                break;
            case BAM_CINS:
                qpos += oplen;
                break;
            case BAM_CDEL:
                rpos += oplen;
                break;
            case BAM_CSOFT_CLIP:
                qpos += oplen;
                break;
            case BAM_CHARD_CLIP:
                qpos += oplen;
                break;
            default:
                fatal_error("Unknown cigar, %u\n", op);
        }
    }
    if (nC2T >= nG2A) { return 0; }
    else              { return 1; }
}

uint8_t get_bsstrand(refcache_t *rs, bam1_t *b, uint32_t min_base_qual, uint8_t allow_u) {
    uint8_t *s;
  
    /* bwa-meth flag has highest priority */
    s = bam_aux_get(b, "YD");
    if (s) {
        s++;
        if      (*s == 'f')            { return 0; }
        else if (*s == 'r')            { return 1; }
        else if (*s == 'u' && allow_u) { return 2; }
    }

    /* bsmap flag */
    s = bam_aux_get(b, "ZS");
    if (s) {
        s++;
        if      (*s == '+') { return 0; }
        else if (*s == '-') { return 1; }
    }

    /* bismark flag */
    s = bam_aux_get(b, "XG");
    if (s) {
        s++;
        if      (strcmp((char*)s, "CT")==0) { return 0; }
        else if (strcmp((char*)s, "GA")==0) { return 1; }
    }

    /* otherwise, guess the bsstrand from nCT and nGA */
    return infer_bsstrand(rs, b, min_base_qual);
}

uint8_t determine_bsstrand_from_pair(refcache_t *rs, bam_hdr_t *hdr, bam1_t *one, bam1_t *two, mkdup_conf_t *conf, uint8_t allow_u) {
    refcache_fetch(rs, hdr->target_name[one->core.tid], one->core.pos > 100 ? one->core.pos-100 : 1, one->core.pos + one->core.l_qseq + 100);
    int8_t bss1 = one ? get_bsstrand(rs, one, conf->min_base_qual, allow_u) : -1;

    refcache_fetch(rs, hdr->target_name[two->core.tid], two->core.pos > 100 ? two->core.pos-100 : 1, two->core.pos + two->core.l_qseq + 100);
    int8_t bss2 = two ? get_bsstrand(rs, two, conf->min_base_qual, allow_u) : -1;

    if (bss1 == bss2 && bss1 >= 0) {
        return bss1;
    } else {
        // Resolve which bisulfite strand to use
        if      (bss1 < 0 && bss2 >= 0) { return bss2; }
        else if (bss2 < 0 && bss1 >= 0) { return bss1; }
        else if (bss1 < 0 && bss2 <  0) {
            // Rare case, assume OT/CTOT
            if (conf->verbose) {
                fprintf(stderr, "[%s] Warning: No valid reads to determine bisulfite strand info. Assuming OT/CTOT strand.\n", __func__);
                fflush(stderr);
            }
            return 0;
        }
        else {
            if (conf->verbose) {
                fprintf(stderr, "[%s] Warning: Inconsistent bisulfite strand info. ", __func__);
                fprintf(stderr, "Taking strand from read with higher total base quality.");
                fflush(stderr);
            }
            uint8_t out = total_qual(one) > total_qual(two) ? bss1 : bss2;
            return out;
        }
    }
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
// Following samblaster's method for handling contigs (particularly genomes with many, many contigs), concatenate all
// contigs in the genome (using the BAM header to define the contigs) into a single full-genome contig. Then, partition
// the full-genome contig out into roughly evenly sized bins (see https://github.com/GregoryFaust/samblaster/issues/21)
// for rationale for this choice of binning.

// Value that samblaster uses for resizing length and offset arrays
// This value may be able to be optimized, but use it for now
#define BIN_RESIZE 32768

// To recap samblaster, 27 bits was large enough to hold human chromosome 1 offset
// This value could be optimized to balance memory/time trade-offs, but for now
// this can stay as it is
#define BIN_SHIFT ((uint64_t)27)
#define BIN_MASK  ((uint64_t)((1 << BIN_SHIFT) - 1))

mkdup_bins_t *prepare_hash_bins(bam_hdr_t *hdr, mkdup_conf_t *conf) {

    mkdup_bins_t *b = mkdup_bins_init();
    b->n_contigs = hdr->n_targets;

    uint64_t sum_length = 0;

    int32_t i;
    for (i=0; i<hdr->n_targets; i++) {
        // Resize length and offset arrays
        if (i % BIN_RESIZE == 0) {
            b->length = (uint32_t *)realloc(b->length, (i+BIN_RESIZE)*sizeof(uint32_t));
            b->offset = (uint64_t *)realloc(b->offset, (i+BIN_RESIZE)*sizeof(uint64_t));
        }

        // Fill length and offset values for specified contig
        b->length[i] = pad_chrom_length(hdr->target_len[i], conf->max_length);
        b->offset[i] = sum_length;

        sum_length += (uint64_t)(b->length[i]);
    }

    b->total_length = sum_length;

    b->n_bins = (sum_length >> BIN_SHIFT) + 1;
    if (b->n_bins >= BIN_RESIZE) {
        // NOTE: exiting here may cause a memory leak, depending on how clean up is handled
        fprintf(stderr, "[%s:%d] Error: Too many contigs found in BAM header\n", __func__, __LINE__);
        exit(1);
    }

    return b;
}

//---------------------------------------------------------------------------------------------------------------------
parsed_cigar_t parse_cigar(bam1_t *b) {
    bam1_core_t *c = &b->core;
    uint32_t i, op, oplen;

    uint32_t rlen = 0, qlen = 0, sclip = 0, eclip = 0;
    uint8_t first_op = 1;
    for (i=0; i<c->n_cigar; ++i) {
        op    = bam_cigar_op(bam_get_cigar(b)[i]);
        oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
        switch(op) {
            case BAM_CMATCH:
                rlen += oplen;
                qlen += oplen;
                first_op = 0;
                break;
            case BAM_CINS:
                qlen += oplen;
                first_op = 0;
                break;
            case BAM_CDEL:
                rlen += oplen;
                first_op = 0;
                break;
            case BAM_CREF_SKIP:
                rlen += oplen;
                first_op = 0;
                break;
            case BAM_CSOFT_CLIP:
                if (first_op) { sclip += oplen; }
                else          { eclip += oplen; }
                break;
            case BAM_CHARD_CLIP:
                if (first_op) { sclip += oplen; }
                else          { eclip += oplen; }
                break;
            case BAM_CEQUAL:
                rlen += oplen;
                qlen += oplen;
                first_op = 0;
                break;
            case BAM_CDIFF:
                rlen += oplen;
                qlen += oplen;
                first_op = 0;
                break;
            default:
                fatal_error("Unknown cigar operation: %u\n", op);
        }
    }

    parsed_cigar_t out = {0};
    out.rlen = rlen;
    out.qlen = qlen;
    out.sclp = sclip;
    out.eclp = eclip;

    return out;
}

signature_t create_signature(bam1_t *forward, bam1_t *reverse, mkdup_conf_t *conf, mkdup_bins_t *bins) {
    // Variables to store things
    parsed_cigar_t cigar;
    uint32_t total_length;
    uint64_t position;

    // Variables for output
    signature_t sig = signature_init();

    // Handle forward strand read
    if (forward != NULL) {
        cigar = parse_cigar(forward);

        total_length = cigar.sclp + cigar.qlen + cigar.eclp;
        if (total_length > conf->max_length) {
            fatal_error("Error: read with length %u is longer than max read length %u\n",
                    total_length, conf->max_length);
        }

        position = bins->offset[forward->core.tid] + forward->core.pos - cigar.sclp;
        sig.bin_start = position >> BIN_SHIFT;
        sig.pos_start = position &  BIN_MASK;
    }

    // Handle reverse strand read
    if (reverse != NULL) {
        cigar = parse_cigar(reverse);

        total_length = cigar.sclp + cigar.qlen + cigar.eclp;
        if (total_length > conf->max_length) {
            fatal_error("Error: read with length %u is longer than max read length %u\n",
                    total_length, conf->max_length);
        }

        position = bins->offset[reverse->core.tid] + reverse->core.pos + cigar.rlen + cigar.eclp;
        sig.bin_end = position >> BIN_SHIFT;
        sig.pos_end = position &  BIN_MASK;
    }

    return sig;
}

void check_if_dup(bam1_t *curr, bam1_t *last, bam_hdr_t *hdr, mkdup_conf_t *conf, refcache_t *rs, mkdup_bins_t *bins,
        SigHash *ot_map, SigHash *ot_for, SigHash *ot_rev,
        SigHash *ob_map, SigHash *ob_for, SigHash *ob_rev) {
    // If unable to determine OT/CTOT or OB/CTOB, defaults to OT/CTOT
    uint8_t bss = determine_bsstrand_from_pair(rs, hdr, last, curr, conf, 0);

    // Ignore case where both reads are not primary alignments
    if (!is_primary(curr) && !is_primary(last)) {
        return;
    }

    if (is_unmapped(curr) && is_unmapped(last)) {
        return;
    }

    bam1_t *forward = bam_init1();
    bam1_t *reverse = bam_init1();

    uint8_t is_forward = 0, is_reverse = 0;

    if (conf->single_end) {
        if (bam_is_rev(curr)) {
            forward = NULL;
            reverse = curr;
            is_reverse = 1;
        } else {
            forward = curr;
            reverse = NULL;
            is_forward = 1;
        }
    } else {
        if (curr->core.flag & BAM_FPROPER_PAIR) {
            if (bam_is_rev(curr) && bam_is_mrev(last)) {
                forward = last;
                reverse = curr;
            } else {
                forward = curr;
                reverse = last;
            }
            is_forward = 1;
            is_reverse = 1;
        } else {
            if (is_unmapped(curr) && bam_is_rev(last)) {
                forward = NULL;
                reverse = last;
                is_reverse = 1;
            } else if (is_unmapped(curr) && bam_is_mrev(last)) {
                forward = last;
                reverse = NULL;
                is_forward = 1;
            } else if (is_unmapped(last) && bam_is_rev(curr)) {
                forward = NULL;
                reverse = curr;
                is_reverse = 1;
            } else if (is_unmapped(last) && bam_is_mrev(curr)) {
                forward = curr;
                reverse = NULL;
                is_forward = 1;
            }
        }
    }

    signature_t sig = create_signature(forward, reverse, conf, bins);

    int not_a_dup;
    if (!is_forward && !is_reverse) {
        return;
    } else if (is_forward && is_reverse) {
        if (bss) {
            sh_put(ob_map, sig, &not_a_dup);
        } else {
            sh_put(ot_map, sig, &not_a_dup);
        }
        if (!not_a_dup) { // signature found!
            // Both reads need to marked as duplicates in this case
            forward->core.flag |= BAM_FDUP;
            reverse->core.flag |= BAM_FDUP;
        }
    } else if (is_forward && !is_reverse) {
        if (bss) {
            sh_put(ob_for, sig, &not_a_dup);
        } else {
            sh_put(ot_for, sig, &not_a_dup);
        }
        if (!not_a_dup) { // signature found!
            // Only the forward strand read needs to be marked as a duplicate
            // The reverse strand either doesn't exist (SE) or is unmapped (PE)
            forward->core.flag |= BAM_FDUP;
        }
    } else if (!is_forward && is_reverse) {
        if (bss) {
            sh_put(ob_rev, sig, &not_a_dup);
        } else {
            sh_put(ot_rev, sig, &not_a_dup);
        }
        if (!not_a_dup) { // signature found!
            // Only the reverse strand read needs to be marked as a duplicate
            // The forward strand either doesn't exist (SE) or is unmapped (PE)
            reverse->core.flag |= BAM_FDUP;
        }
    }
}

int markdups(mkdup_conf_t *conf) {

    // Read input file and BAM header
    htsFile *infh = hts_open(conf->infn, "r");
    if (infh == NULL) {
        fprintf(stderr, "[%s:%d] Error: Unable to read from %s\n",
                __func__, __LINE__, strcmp(conf->infn, "-") == 0 ? "stdin" : conf->infn);
        fflush(stderr);
        return 1;
    }
    bam_hdr_t *hdr = sam_hdr_read(infh);

    // Prepare hash table bins based on input file
    mkdup_bins_t *bins = prepare_hash_bins(hdr, conf);

    // Open output file
    htsFile *outfh = hts_open(conf->outfn, conf->out_mode);
    if (outfh == NULL) {
        fprintf(stderr, "[%s:%d] Error: Unable to write to %s\n",
                __func__, __LINE__, strcmp(conf->outfn, "-") == 0 ? "stdout" : conf->outfn);
        fflush(stderr);

        hts_close(infh);

        return 1;
    }

    // Write header to output file
    // TODO: Add new PG line to header before writing
    if (sam_hdr_write(outfh, hdr) < 0) {
        fprintf(stderr, "[%s:%d] Error: Header writing failed. Abort.\n", __func__, __LINE__);
        fflush(stderr);

        hts_close(outfh);
        bam_hdr_destroy(hdr);
        hts_close(infh);

        return 1;
    }

    // Initialize variables for reading
    refcache_t *rs = init_refcache(conf->reffn, 1000, 1000);
    int ret;

    uint32_t n_read_ids      = 0; // Total number of read IDs processed
    uint32_t n_matching_ids  = 0; // Number of PE read IDs processed that were next to their mate (i.e., name sorted)
    uint8_t  is_coord_sorted = 0; // File appears to be coordinate sorted (1 if true, 0 if false)
    uint8_t  is_paired       = 0; // Read is paired-end (1 if true, 0 if false) - used for running in single-end mode
    uint8_t  successful_run  = 1; // Run had no errors crop up

    // Prepare hash maps
    SigHash *ot_map = sh_init();
    SigHash *ot_for = sh_init();
    SigHash *ot_rev = sh_init();
    SigHash *ob_map = sh_init();
    SigHash *ob_for = sh_init();
    SigHash *ob_rev = sh_init();

    if (conf->single_end) { // process file as single-end BAM
        bam1_t *curr = bam_init1();

        while ((ret = sam_read1(infh, hdr, curr)) >= 0) {
            // Check for PE read
            if (curr->core.flag & BAM_FPAIRED) {
                successful_run = 0;
                is_paired = 1;
                break;
            }

            // Remove existing duplicate info
            curr->core.flag &= ~BAM_FDUP;

            check_if_dup(curr, NULL, hdr, conf, rs, bins, ot_map, ot_for, ot_rev, ob_map, ob_for, ob_rev);

            n_read_ids++;

            // Print out read if not marked as a duplicate or user requested they be printed
            if (!is_duplicate(curr) || (is_duplicate(curr) && !conf->rm_dup)) {
                sam_write1(outfh, hdr, curr);
            }
        }

        bam_destroy1(curr);
    } else { // process file as paired-end BAM
        // Read in first read - primes loop through remaining reads as well as
        // checking to see if there are any reads in the file
        bam1_t *last = bam_init1();
        if ((ret = sam_read1(infh, hdr, last)) < 0) {
            fprintf(stderr, "[%s:%d] Error: No reads found in %s\n", __func__,
                    __LINE__, strcmp(conf->infn, "-") ? conf->infn : "stdin");
            fflush(stderr);

            bam_destroy1(last);
            hts_close(outfh);
            bam_hdr_destroy(hdr);
            hts_close(infh);

            return 1;
        }
        n_read_ids = 1;

        // Remove existing duplicate info
        last->core.flag &= ~BAM_FDUP;

        // FIXME: If there is only 1 read in the file, then no reads will be written
        //        to output. Need to figure out how to handle this edge case.
        bam1_t *curr = bam_init1();
        while ((ret = sam_read1(infh, hdr, curr)) >= 0) {
            // Remove existing duplicate info
            curr->core.flag &= ~BAM_FDUP;

            if (strcmp(bam_get_qname(last), bam_get_qname(curr)) == 0) {
                n_matching_ids++;
                check_if_dup(curr, last, hdr, conf, rs, bins, ot_map, ot_for, ot_rev, ob_map, ob_for, ob_rev);

                // Print lines if not duplicates or user requests that duplicates be printed
                if (!is_duplicate(last) || (is_duplicate(last) && !conf->rm_dup)) {
                    sam_write1(outfh, hdr, last);
                }
                if (!is_duplicate(curr) || (is_duplicate(curr) && !conf->rm_dup)) {
                    sam_write1(outfh, hdr, curr);
                }

                last = curr;
                curr = bam_init1();
            } else {
                last = curr;
                curr = bam_init1();
                n_read_ids++;

                if (n_read_ids > 100 && n_matching_ids < 50) {
                    successful_run = 0;
                    is_coord_sorted = 1;
                    break;
                }
            }
        }

        bam_destroy1(curr);
        bam_destroy1(last);
    }

    // Verify PE was not coordinate sorted
    if (!conf->single_end && !is_coord_sorted && n_read_ids <= 100 && n_read_ids > n_matching_ids) {
        successful_run = 0;
        is_coord_sorted = 1;
    }

    // Print some processing metrics
    if (successful_run) {
        fprintf(stderr, "[%s] processing mode: %s\n", __func__, (conf->single_end) ? "single-end" : "paired-end");
        fprintf(stderr, "[%s] processed number of reads ids: %u\n", __func__, n_read_ids);
        fflush(stderr);
    } else {
        if (is_coord_sorted) {
            fprintf(stderr, "ERROR: THIS APPEARS TO BE A COORDINATE SORTED BAM. PLEASE NAME SORT (samtools sort -n) AND RETRY.\n");
        }
        if (is_paired) {
            fprintf(stderr, "ERROR: PAIRED-END READ FOUND WHEN RUNNING IN SINGLE-END MODE (-s). MIXED BAMS ARE NOT ALLOWED.\n");
        }
        fprintf(stderr, "Error occurred. Deleting output BAM (if not streaming)\n");
    }

    // Clean up
    sh_destroy(ob_rev);
    sh_destroy(ob_for);
    sh_destroy(ob_map);
    sh_destroy(ot_rev);
    sh_destroy(ot_for);
    sh_destroy(ot_map);

    free_refcache(rs);

    hts_close(outfh);
    destroy_mkdup_bins(bins);
    bam_hdr_destroy(hdr);
    hts_close(infh);

    // Delete output BAM file if there was an error
    if (!successful_run && strcmp(conf->outfn, "-") != 0) {
        if (remove(conf->outfn) == 0) {
            fprintf(stderr, "Successfully deleted %s\n", conf->outfn);
        } else {
            fprintf(stderr, "Unable to delete %s\n", conf->outfn);
        }
    }

    if (successful_run) {
        return 0;
    } else {
        return 1;
    }
}

//---------------------------------------------------------------------------------------------------------------------
static int usage(mkdup_conf_t *conf) {
    fprintf(stderr, "\n");
    fprintf(stderr, "dupsifter [options] <ref.fa> [in.bam]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -o STR    name of output file [stdout]\n");
    fprintf(stderr, "Input options:\n");
    fprintf(stderr, "    -s        run for single-end data\n");
    fprintf(stderr, "    -l INT    maximum read length for paired end duplicate-marking [%d]\n", conf->max_length);
    fprintf(stderr, "    -b INT    minimum base quality [30].\n");
    fprintf(stderr, "    -r        toggle to remove marked duplicate.\n");
    fprintf(stderr, "    -v        print extra messages\n");
    fprintf(stderr, "    -h        this help.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note, [in.bam] must be name sorted. If not provided, assume the input is stdin.\n");
    fprintf(stderr, "Another note, assumes either ALL reads are paired-end (default) or single-end.\n");
    fprintf(stderr, "If a singleton read is found in paired-end mode, the code will break nicely.\n");
    fprintf(stderr, "\n");

    return 1;
}

int main(int argc, char *argv[]) {

    int c;
    mkdup_conf_t conf = {0};
    mkdup_conf_init(&conf);

    // TODO: Add functionality for adding mate tags (MC and MQ) to reads
    if (argc < 1) { return usage(&conf); }
    while ((c = getopt(argc, argv, ":b:l:o:rsqvh")) >= 0) {
        switch (c) {
            case 'b': conf.min_base_qual = atoi(optarg); break;
            case 'l': conf.max_length = (uint32_t)atoi(optarg); break;
            case 'o': conf.outfn = optarg; break;
            case 'r': conf.rm_dup = 1; break;
            case 's': conf.single_end = 1; break;
            case 'v': conf.verbose = 1; break;
            case 'h': return usage(&conf);
            case ':':
                fprintf(stderr, "Option needs an argument: -%c\n", optopt);
                return usage(&conf);
            case '?':
                fprintf(stderr, "Unrecognized option: -%c\n", optopt);
                return usage(&conf);
            default:
                fprintf(stderr, "[%s:%d] Unrecognized command: -%c\n", __func__, __LINE__, c);
                return usage(&conf);
        }
    }

    conf.reffn = optind < argc ? argv[optind++] : NULL;
    conf.infn  = optind < argc ? argv[optind++] : "-";
    if (!conf.reffn) {
        fprintf(stderr, "Please provide a reference\n");
        return usage(&conf);
    }

    // Set up input file
    if (strcmp(conf.infn, "-") == 0) {
        fprintf(stderr, "Reading input from stdin\n");
    } else {
        fprintf(stderr, "Reading input from %s\n", conf.infn);
    }

    // Set up output file
    if (strcmp(conf.outfn, "-") == 0) {
        fprintf(stderr, "Writing output to stdout\n");
    } else {
        fprintf(stderr, "Writing output to %s\n", conf.outfn);
    }

    // Set write mode, auto-detect from file name, assume no extra compression
    // or other extra additions to the write mode
    sam_open_mode(conf.out_mode+1, conf.outfn, NULL);

    markdups(&conf);

    return 0;
}
