/* Mark duplicates for both WGS and WGBS reads
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
#include <getopt.h>
#include "sam.h"
#include "khashl.h"
#include "util.h"
#include "refcache.h"
#include "version.h"

//---------------------------------------------------------------------------------------------------------------------
// Config struct, initialization function, and output function
#define DEF_STAT_NAME "dupsifter.stat"
#define TAG_STAT_NAME ".dupsifter.stat"

typedef struct {
    char     *reffn;               /* reference file name */
    char     *infn;                /* name of input file */
    char     *outfn;               /* name of output file */
    char     *statfn;               /* name of statistics file */
    char      out_mode[6];         /* write mode of output file */
    char     *arg_list;            /* input argument list for PG tag */
    uint32_t  min_base_qual;       /* threshold for high quality bases */
    uint32_t  max_length;          /* max read length allowed */
    uint8_t   rm_dup;              /* flag to remove duplicates */
    uint8_t   is_wgs;              /* process reads as WGS instead of WGBS */
    uint8_t   verbose;             /* level of messages to print */
    uint8_t   single_end;          /* process all reads as single-end */
    uint8_t   add_mate_tags;       /* add MC and MQ tags to mate reads */

    uint32_t  n_reads_read;        /* Number of reads read from file */
    uint32_t  cnt_id_both_map;     /* number of PE reads */
    uint32_t  cnt_id_forward;      /* number of reads where only one read is mapped and it's on the forward strand */
    uint32_t  cnt_id_reverse;      /* number of reads where only one read is mapped and it's on the reverse strand */
    uint32_t  cnt_id_both_dup;     /* number of PE reads marked as dup */
    uint32_t  cnt_id_forward_dup;  /* number of duplicates on the forward strand */
    uint32_t  cnt_id_reverse_dup;  /* number of duplicates on the reverse strand */
    uint32_t  cnt_id_no_map;       /* number of reads with both reads unmapped (PE) or the read is unmapped (SE) */
    uint32_t  cnt_id_no_prim;      /* number of reads with no primary reads */
    uint32_t  cnt_id_n_primary;    /* number of primary alignments */
    uint32_t  cnt_id_n_sec_sup;    /* number of secondary and supplementary alignments */
} ds_conf_t;

void ds_conf_init(ds_conf_t *conf) {
    strcpy(conf->out_mode, "w");

    conf->outfn = (char *)"-";
    conf->arg_list = NULL;
    conf->statfn = NULL;
    conf->min_base_qual = 0;
    conf->max_length = 10000;
    conf->rm_dup = 0;
    conf->is_wgs = 0;
    conf->verbose = 0;
    conf->single_end = 0;
    conf->add_mate_tags = 0;

    conf->n_reads_read = 0;
    conf->cnt_id_both_map = 0;
    conf->cnt_id_forward = 0;
    conf->cnt_id_reverse = 0;
    conf->cnt_id_both_dup = 0;
    conf->cnt_id_forward_dup = 0;
    conf->cnt_id_reverse_dup = 0;
    conf->cnt_id_no_map = 0;
    conf->cnt_id_no_prim = 0;
    conf->cnt_id_n_primary = 0;
    conf->cnt_id_n_sec_sup = 0;
}

void ds_conf_destroy(ds_conf_t *conf) {
    free(conf->arg_list);
    if (conf->statfn != NULL) { free(conf->statfn); }
}

void ds_conf_print(ds_conf_t *conf) {
    FILE *statfh = fopen(conf->statfn, "w");

    fprintf(statfh, "[dupsifter] processing mode: %s\n", (conf->single_end) ? "single-end" : "paired-end");
    fprintf(statfh, "[dupsifter] number of individual reads processed: %u\n", conf->n_reads_read);
    fprintf(statfh, "[dupsifter] number of reads with both reads mapped: %u\n", conf->cnt_id_both_map);
    fprintf(statfh, "[dupsifter] number of reads with only one read mapped to the forward strand: %u\n", conf->cnt_id_forward);
    fprintf(statfh, "[dupsifter] number of reads with only one read mapped to the reverse strand: %u\n", conf->cnt_id_reverse);
    fprintf(statfh, "[dupsifter] number of reads with both reads marked as duplicates: %u\n", conf->cnt_id_both_dup);
    fprintf(statfh, "[dupsifter] number of reads on the forward strand marked as duplicates: %u\n", conf->cnt_id_forward_dup);
    fprintf(statfh, "[dupsifter] number of reads on the reverse strand marked as duplicates: %u\n", conf->cnt_id_reverse_dup);
    fprintf(statfh, "[dupsifter] number of individual primary-alignment reads: %u\n", conf->cnt_id_n_primary);
    fprintf(statfh, "[dupsifter] number of individual secondary- and supplementary-alignment reads: %u\n", conf->cnt_id_n_sec_sup);
    fprintf(statfh, "[dupsifter] number of reads with no reads mapped: %u\n", conf->cnt_id_no_map);
    fprintf(statfh, "[dupsifter] number of reads with no primary reads: %u\n", conf->cnt_id_no_prim);

    fclose(statfh);
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
// Struct for storing all reads matching a given read name
typedef struct bam1_chain bam1_chain_t;
struct bam1_chain {
    bam1_t *read;
    struct bam1_chain *next;
};

bam1_chain_t *bam1_chain_init(bam1_t *b) {
    bam1_chain_t *out = (bam1_chain_t *)malloc(sizeof(bam1_chain_t));
    if (b == NULL) { out->read = bam_init1(); }
    else           { out->read = b; }
    out->next = NULL;

    return out;
}

void destroy_bam1_chain(bam1_chain_t *b) {
    if (b == NULL) { return; }

    struct bam1_chain *temp = b;
    while (temp != NULL) {
        struct bam1_chain *next = temp->next;
        if (temp->read != NULL) { bam_destroy1(temp->read); }
        free(temp);
        temp = next;
    }
    
    if (temp != NULL) { free(temp); }
}

void add_new_read_to_chain(bam1_chain_t **head, bam1_t *b) {
    bam1_chain_t *new  = bam1_chain_init(b);
    bam1_chain_t *last = *head;

    if (*head == NULL) {
        *head = new;
        return;
    }

    while (last->next != NULL) { last = last->next; }
    last->next = new;

    return;
}

uint32_t get_count(bam1_chain_t *bc) {
    uint32_t count = 0;
    bam1_chain_t *curr = bc;
    while (curr != NULL) {
        count++;
        curr = curr->next;
    }

    return count;
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

// Read orientation possibilities (r1-r2: forward-forward, reverse-reverse, forward-reverse, reverse-forward)
#define OO_FF 0
#define OO_RR 1
#define OO_FR 2
#define OO_RF 3

// Is read single (or paired with unmapped mate) or paired?
#define S_N 0
#define S_Y 1

// Is read1 in paired reads the leftmost read (single/paired with unmapped mate always L_N)?
#define L_N 0
#define L_Y 1

typedef struct {
    uint8_t  l_or_s; /* packed variable, (from leftmost bit) bit 5 = is r1 leftmost read, bit 6/7 = orientation, bit 8 = is single */
    uint16_t r1_bin; /* bin number for r1 position */
    uint16_t r2_bin; /* bin number for r2 position */
    uint32_t r1_pos; /* bin position for r1 position */
    uint32_t r2_pos; /* bin position for r2 position */
} signature_t;

signature_t signature_init() {
    // Chances of the bin or position being the max value are very low, so set default value to max value
    signature_t sig = {0};
    sig.l_or_s = 128;
    sig.r1_bin = UINT16_MAX;
    sig.r2_bin = UINT16_MAX;
    sig.r1_pos = UINT32_MAX;
    sig.r2_pos = UINT32_MAX;

    return sig;
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
// Initialize hash maps and define functions needed for initialization
khint_t hash_sig(signature_t s) {
    khint_t to_hash;

    // Create input value for hashing function
    // Pull off last 16 bits of r1_pos
    uint16_t pos = s.r1_pos & 0xffff;
    // Pull off last 12 bits of r1_bin
    uint16_t bin = s.r1_bin & 0x0fff;

    // Bit pattern in to_hash
    // PPPP PPPP PPPP PPPP BBBB BBBB BBBB LOOS
    // |_________________| |____________| |__|
    //          |                |         |
    // |-----------------| |------------| |--|
    //       Position           Bin       r1_leftmost (L)/orientation (OO)/is_single (S)
    to_hash = (pos << 16) + (bin << 4) + (s.l_or_s & 0xf);

    return kh_hash_uint32(to_hash);
}

int sig_equal(signature_t s1, signature_t s2) {
    if (s1.l_or_s != s2.l_or_s) { return 0; }
    if (s1.r1_bin != s2.r1_bin) { return 0; }
    if (s1.r1_pos != s2.r1_pos) { return 0; }
    if (s1.r2_bin != s2.r2_bin) { return 0; }
    if (s1.r2_pos != s2.r2_pos) { return 0; }

    return 1;
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
} ds_bins_t;

ds_bins_t *ds_bins_init() {
    return (ds_bins_t *)calloc(1, sizeof(ds_bins_t));
}

void destroy_ds_bins(ds_bins_t *b) {
    if (b->length != NULL) { free(b->length); }
    if (b->offset != NULL) { free(b->offset); }
    free(b);
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
// Padding functions
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
                fatal_error("[dupsifter] Error: Unknown cigar, %u\n", op);
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

uint8_t determine_bsstrand(refcache_t *rs, bam_hdr_t *hdr, bam1_t *one, bam1_t *two,
        ds_conf_t *conf, uint8_t allow_u) {
    int8_t bss1 = -1, bss2 = -1;
    if (one && !is_unmapped(one)) {
        refcache_fetch(rs, hdr->target_name[one->core.tid],
                one->core.pos > 100 ? one->core.pos-100 : 1,
                one->core.pos + one->core.l_qseq + 100);
        bss1 = get_bsstrand(rs, one, conf->min_base_qual, allow_u);
    }

    if (two && !is_unmapped(two)) {
        refcache_fetch(rs, hdr->target_name[two->core.tid],
                two->core.pos > 100 ? two->core.pos-100 : 1,
                two->core.pos + two->core.l_qseq + 100);
        bss2 = get_bsstrand(rs, two, conf->min_base_qual, allow_u);
    }

    if (bss1 == bss2 && bss1 >= 0) {
        return (uint8_t)bss1;
    } else {
        // Resolve which bisulfite strand to use
        if      (bss1 < 0 && bss2 >= 0) { return (uint8_t)bss2; }
        else if (bss2 < 0 && bss1 >= 0) { return (uint8_t)bss1; }
        else if (bss1 < 0 && bss2 <  0) {
            // Rare case, assume OT/CTOT
            if (conf->verbose) {
                fprintf(stderr,
                        "[%s] Warning: No valid reads to determine bisulfite strand info. Assuming OT/CTOT strand.\n",
                        __func__);
                fflush(stderr);
            }
            return 0;
        }
        else {
            if (conf->verbose) {
                fprintf(stderr, "[%s] Warning: Inconsistent bisulfite strand info. ", __func__);
                fprintf(stderr, "Taking strand from read with higher total base quality.\n");
                fflush(stderr);
            }
            uint8_t out = total_qual(one) > total_qual(two) ? (uint8_t)bss1 : (uint8_t)bss2;
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

ds_bins_t *prepare_hash_bins(bam_hdr_t *hdr, ds_conf_t *conf) {

    ds_bins_t *b = ds_bins_init();
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
                fatal_error("[dupsifter] Error: Unknown cigar operation: %u\n", op);
        }
    }

    parsed_cigar_t out = {0};
    out.rlen = rlen;
    out.qlen = qlen;
    out.sclp = sclip;
    out.eclp = eclip;

    return out;
}

signature_t create_single_signature(bam1_t *read, ds_conf_t *conf, ds_bins_t *bins) {
    if (read == NULL) {
        fatal_error("[dupsifter] Error: Read information for single or mate-unmapped read is missing\n");
    }

    parsed_cigar_t cigar = parse_cigar(read);
    uint32_t total_length = cigar.sclp + cigar.qlen + cigar.eclp;
    if (total_length > conf->max_length) {
        fatal_error("[dupsifter] Error: Read with length %u is longer than max read length %u\n",
                total_length, conf->max_length);
    }

    uint64_t position;
    uint8_t  orientation;

    if (bam_is_rev(read)) {
        position = bins->offset[read->core.tid] + read->core.pos + cigar.rlen + cigar.eclp;
        orientation = OO_RR;
    } else {
        position = bins->offset[read->core.tid] + read->core.pos - cigar.sclp;
        orientation = OO_FF;
    }

    signature_t sig = signature_init();
    sig.r1_bin = position >> BIN_SHIFT;
    sig.r1_pos = position &  BIN_MASK;
    sig.l_or_s = (L_N << 3) + (orientation << 1) + S_Y;

    return sig;
}

signature_t create_paired_signature(bam1_t *read1, bam1_t *read2, ds_conf_t *conf, ds_bins_t *bins) {
    if (read1 == NULL || read2 == NULL) {
        fatal_error("[dupsifter] Error: Read information for paired end read is missing\n");
    }

    parsed_cigar_t cigar1 = parse_cigar(read1);
    parsed_cigar_t cigar2 = parse_cigar(read2);
    
    uint32_t total_length_1 = cigar1.sclp + cigar1.qlen + cigar1.eclp;
    uint32_t total_length_2 = cigar2.sclp + cigar2.qlen + cigar2.eclp;

    if (total_length_1 > conf->max_length) {
        fatal_error("[dupsifter] Error: Read with length %u is longer than max read length %u\n",
                total_length_1, conf->max_length);
    }
    if (total_length_2 > conf->max_length) {
        fatal_error("[dupsifter] Error: Read with length %u is longer than max read length %u\n",
                total_length_2, conf->max_length);
    }

    uint64_t pos1, pos2, beg1, beg2, end1, end2;
    uint32_t ref1, ref2;
    uint8_t  r1_leftmost, orientation;

    beg1 = bins->offset[read1->core.tid] + read1->core.pos - cigar1.sclp;
    beg2 = bins->offset[read2->core.tid] + read2->core.pos - cigar2.sclp;
    end1 = bins->offset[read1->core.tid] + read1->core.pos + cigar1.rlen + cigar1.eclp;
    end2 = bins->offset[read2->core.tid] + read2->core.pos + cigar2.rlen + cigar2.eclp;

    ref1 = read1->core.tid;
    ref2 = read2->core.tid;

    // Work out if read1 is leftmost in pair
    if (ref1 != ref2) {
        r1_leftmost = ref1 < ref2;
    } else {
        if (bam_is_rev(read1) == bam_is_rev(read1)) {
            if (!bam_is_rev(read1)) {
                r1_leftmost = beg1 <= beg2;
            } else {
                r1_leftmost = end1 <= end2;
            }
        } else {
            if (bam_is_rev(read1)) {
                r1_leftmost = end1 <= end2;
            } else {
                r1_leftmost = beg1 <= beg2;
            }
        }
    }
    r1_leftmost = r1_leftmost ? L_Y : L_N;

    // Determine orientation of reads
    if (bam_is_rev(read1) == bam_is_rev(read2)) {
        if (!bam_is_rev(read1)) {
            orientation = OO_FF;
        } else {
            orientation = OO_RR;
        }
    } else {
        if (!bam_is_rev(read1)) {
            orientation = OO_FR;
        } else {
            orientation = OO_RF;
        }
    }

    // Find final positions
    switch (orientation) {
        case OO_FF:
            pos1 = beg1;
            pos2 = beg2;
            break;
        case OO_RR:
            pos1 = end1;
            pos2 = end2;
            break;
        case OO_FR:
            pos1 = beg1;
            pos2 = end2;
            break;
        case OO_RF:
            pos1 = end1;
            pos2 = beg2;
            break;
        default:
            fatal_error("[dupsifter] Error: Unknown orientation found.\n");
    }

    signature_t sig = signature_init();
    sig.r1_bin = pos1 >> BIN_SHIFT;
    sig.r1_pos = pos1 &  BIN_MASK;
    sig.r1_bin = pos2 >> BIN_SHIFT;
    sig.r1_pos = pos2 &  BIN_MASK;
    sig.l_or_s = (r1_leftmost << 3) + (orientation << 1) + S_N;

    return sig;
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
// Functions for adding mate (MC/MQ) tags
static uint8_t create_cigar_string(bam1_t *b, kstring_t *str) {
    // Empty CIGAR gets a "*"
    if (b->core.n_cigar == 0) {
        return (kputc('*', str) == EOF) ? 0 : 1;
    }

    uint32_t *cigar = bam_get_cigar(b);
    uint32_t i;

    for (i=0; i<b->core.n_cigar; i++) {
        if (kputw(bam_cigar_oplen(cigar[i]), str) == EOF) { return 0; }
        if (kputc(bam_cigar_opchr(cigar[i]), str) == EOF) { return 0; }
    }

    return 1;
}

void add_MCMQ(bam1_chain_t *bc, bam1_t *b_read, bam1_t *b_mate) {
    if (is_unmapped(b_read)) { return; }

    int mask = (b_mate->core.flag & (0x40 | 0x80));

    kstring_t mc = { 0, 0, NULL };
    if (!create_cigar_string(b_read, &mc)) { return; }

    uint32_t mq = b_read->core.qual;

    for (bam1_chain_t *next = bc; next != NULL; next = next->next) {
        bam1_t *read = next->read;

        if (read->core.flag & mask) {
            if (bam_aux_get(read, "MC") == NULL) {
                bam_aux_append(read, "MC", 'Z', ks_len(&mc)+1, (uint8_t *)ks_str(&mc));
            }

            if (bam_aux_get(read, "MQ") == NULL) {
                bam_aux_append(read, "MQ", 'i', sizeof(uint32_t), (uint8_t *)&mq);
            }
        }
    }

    free(mc.s);
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
void problem_chain(bam1_chain_t *bc, uint32_t count) {
    fatal_error("[dupsifter] Error: Can't find read 1 and/or read 2 in %u reads with read ID: %s. Are these reads coordinate sorted?\n",
            count, bam_get_qname(bc->read));
}

void mark_dup(bam1_chain_t *bc, bam_hdr_t *hdr, ds_conf_t *conf, refcache_t *rs, ds_bins_t *bins,
        SigHash *ot_map, SigHash *ot_for, SigHash *ot_rev,
        SigHash *ob_map, SigHash *ob_for, SigHash *ob_rev) {
    bam1_t *r1 = NULL;
    bam1_t *r2 = NULL;
    uint32_t count = 0;
    for (bam1_chain_t *curr = bc; curr != NULL; curr = curr->next) {
        count++;

        // Remove any duplicate marks that have already been applied
        curr->read->core.flag &= ~BAM_FDUP;

        // Secondary and supplementary alignments are not used to determine duplicates
        if (!is_primary(curr->read)) { conf->cnt_id_n_sec_sup++; continue; }

        // If unpaired, set as r1 for creating signature
        if      (!is_paired(curr->read))     { r1 = curr->read; }
        // Set reads 1 and 2
        else if (is_first_read(curr->read))  { r1 = curr->read; }
        else if (is_second_read(curr->read)) { r2 = curr->read; }

        conf->cnt_id_n_primary++;
    }

    if (r1 == NULL && r2 == NULL) {
        conf->cnt_id_no_prim++;
        fprintf(stderr, "[dupsifter] Warning: No valid primary alignments found for read ID: %s", bam_get_qname(bc->read));
        return;
    }

    // Check for singleton (either SE read or PE read with unmapped mate) reads
    uint8_t is_single = 0;
    if (r1 == NULL || r2 == NULL) {
        if (r1 == NULL) { swap_bam1_pointers(&r1, &r2); }

        // If there is only one read that says it's paired, but is unmapped or
        // it's mate is mapped, then something went wrong
        if (is_paired(r1) && (is_unmapped(r1) || !is_mate_unmapped(r1))) {
            problem_chain(bc, count);
        }

        // If the only read is unmapped, then it isn't a duplicate
        if (is_unmapped(r1) && !conf->single_end) {
            problem_chain(bc, count);
        } else if (is_unmapped(r1) && conf->single_end) {
            conf->cnt_id_no_map++;
            return;
        }

        is_single = 1;
    } else {
        // Add MC and MQ tags if desired
        if (conf->add_mate_tags) {
            add_MCMQ(bc, r1, r2);
            add_MCMQ(bc, r2, r1);
        }

        // Don't mark reads as duplicates if both are unmapped
        if (is_unmapped(r1) && is_unmapped(r2)) {
            conf->cnt_id_no_map++;
            return;
        }

        // Check for unmapped mate
        is_single = (is_unmapped(r1) || is_unmapped(r2));
        if (is_unmapped(r1) && !is_unmapped(r2)) {
            swap_bam1_pointers(&r1, &r2);
        }
    }

    // If unable to determine OT/CTOT or OB/CTOB, defaults to OT/CTOT
    // For WGS case, only use OT/CTOT
    uint8_t bss = 0;
    if (!conf->is_wgs) {
        bss = determine_bsstrand(rs, hdr, r1, r2, conf, 0);
    }

    uint8_t is_forward = 0; // used only for single reads to determine strand
    signature_t sig = signature_init();

    if (is_single) {
        if (bam_is_rev(r1)) {
            conf->cnt_id_reverse++;
        } else {
            conf->cnt_id_forward++;
            is_forward = 1;
        }

        sig = create_single_signature(r1, conf, bins);
    } else {
        conf->cnt_id_both_map++;

        sig = create_paired_signature(r1, r2, conf, bins);
    }

    // Do duplicate marking, either all reads in chain will marked as dups or none will
    int not_a_dup;
    if (!is_single) {
        if (bss) {
            sh_put(ob_map, sig, &not_a_dup);
        } else {
            sh_put(ot_map, sig, &not_a_dup);
        }
        if (!not_a_dup) { // signature found!
            // Both reads need to marked as duplicates in this case
            conf->cnt_id_both_dup++;
            for (bam1_chain_t *curr = bc; curr != NULL; curr = curr->next) {
                curr->read->core.flag |= BAM_FDUP;
            }
        }
    } else {
        if (is_forward) {
            if (bss) {
                sh_put(ob_for, sig, &not_a_dup);
            } else {
                sh_put(ot_for, sig, &not_a_dup);
            }
            if (!not_a_dup) { // signature found!
                conf->cnt_id_forward_dup++;
                for (bam1_chain_t *curr = bc; curr != NULL; curr = curr->next) {
                    curr->read->core.flag |= BAM_FDUP;
                }
            }
        } else {
            if (bss) {
                sh_put(ob_rev, sig, &not_a_dup);
            } else {
                sh_put(ot_rev, sig, &not_a_dup);
            }
            if (!not_a_dup) { // signature found!
                conf->cnt_id_reverse_dup++;
                for (bam1_chain_t *curr = bc; curr != NULL; curr = curr->next) {
                    curr->read->core.flag |= BAM_FDUP;
                }
            }
        }
    }
}

int dupsifter(ds_conf_t *conf) {
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
    ds_bins_t *bins = prepare_hash_bins(hdr, conf);

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
    if (sam_hdr_add_pg(hdr, "dupsifter", "VN", dupsifter_version(),
                conf->arg_list ? "CL" : NULL, conf->arg_list ? conf->arg_list : NULL, NULL)) {
        fprintf(stderr, "[dupsifter] Error: Failed to add PG tag to output header\n");
        fflush(stderr);

        hts_close(outfh);
        bam_hdr_destroy(hdr);
        hts_close(infh);

        return 1;
    }
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

    // Prepare hash maps
    SigHash *ot_map = sh_init();
    SigHash *ot_for = sh_init();
    SigHash *ot_rev = sh_init();
    SigHash *ob_map = sh_init();
    SigHash *ob_for = sh_init();
    SigHash *ob_rev = sh_init();

    // Read in first read - primes loop through remaining reads as well as
    // checking to see if there are any reads in the file
    bam1_t *curr = bam_init1();
    if ((ret = sam_read1(infh, hdr, curr)) < 0) {
        fprintf(stderr, "[%s:%d] Error: No reads found in %s\n", __func__,
                __LINE__, strcmp(conf->infn, "-") ? conf->infn : "stdin");
        fflush(stderr);

        bam_destroy1(curr);
        hts_close(outfh);
        bam_hdr_destroy(hdr);
        hts_close(infh);

        return 1;
    }
    conf->n_reads_read = 1;
    bam1_chain_t *bc_curr = bam1_chain_init(curr);

    bam1_t *next = bam_init1();
    while ((ret = sam_read1(infh, hdr, next)) >= 0) {
        conf->n_reads_read++;
        if (strcmp(bam_get_qname(next), bam_get_qname(bc_curr->read)) != 0) {
            mark_dup(bc_curr, hdr, conf, rs, bins,
                    ot_map, ot_for, ot_rev, ob_map, ob_for, ob_rev);
            for (bam1_chain_t *curr = bc_curr; curr != NULL; curr = curr->next) {
                // If the first read in the chain is a duplicate then all others will be,
                // so if we want to remove dupliates, then we can break out right away
                if (is_duplicate(curr->read) && conf->rm_dup) { break; }
                sam_write1(outfh, hdr, curr->read);
            }
            destroy_bam1_chain(bc_curr);
            bc_curr = bam1_chain_init(next);
        } else {
            add_new_read_to_chain(&bc_curr, next);
        }
        next = bam_init1();
    }

    if (get_count(bc_curr) > 0) {
        mark_dup(bc_curr, hdr, conf, rs, bins, 
                ot_map, ot_for, ot_rev, ob_map, ob_for, ob_rev);
        for (bam1_chain_t *curr = bc_curr; curr != NULL; curr = curr->next) {
            // If the first read in the chain is a duplicate then all others will be,
            // so if we want to remove dupliates, then we can break out right away
            if (is_duplicate(curr->read) && conf->rm_dup) { break; }
            sam_write1(outfh, hdr, curr->read);
        }
    }

    // Clean up
    bam_destroy1(next);
    destroy_bam1_chain(bc_curr);

    sh_destroy(ob_rev);
    sh_destroy(ob_for);
    sh_destroy(ob_map);
    sh_destroy(ot_rev);
    sh_destroy(ot_for);
    sh_destroy(ot_map);

    free_refcache(rs);

    hts_close(outfh);
    destroy_ds_bins(bins);
    bam_hdr_destroy(hdr);
    hts_close(infh);

    return 0;
}
//---------------------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------------------
static int print_version() {
    fprintf(stderr, "Program: dupsifter\n");
    fprintf(stderr, "Version: %s\n", DUPSIFTER_VERSION);
    fprintf(stderr, "Contact: Jacob Morrison <jacob.morrison@vai.org>\n");

    return 0;
}

static int usage(ds_conf_t *conf) {
    fprintf(stderr, "\n");
    print_version();
    fprintf(stderr, "\n");
    fprintf(stderr, "dupsifter [options] <ref.fa> [in.bam]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -o, --output STR             name of output file [stdout]\n");
    fprintf(stderr, "    -O, --stats-output STR       name of file to write statistics to (see Note 3 for details)\n");
    fprintf(stderr, "Input options:\n");
    fprintf(stderr, "    -s, --single-end             run for single-end data\n");
    fprintf(stderr, "    -m, --add-mate-tags          add MC and MQ mate tags to mate reads\n");
    fprintf(stderr, "    -W, --wgs-only               process WGS reads instead of WGBS\n");
    fprintf(stderr, "    -l, --max-read-length INT    maximum read length for paired end duplicate-marking [%u]\n", conf->max_length);
    fprintf(stderr, "    -b, --min-base-qual INT      minimum base quality [%u]\n", conf->min_base_qual);
    fprintf(stderr, "    -r, --remove-dups            toggle to remove marked duplicate\n");
    fprintf(stderr, "    -v, --verbose                print extra messages\n");
    fprintf(stderr, "    -h, --help                   this help\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Note 1, [in.bam] must be name sorted. If not provided, assume the input is stdin.\n");
    fprintf(stderr, "Note 2, assumes either ALL reads are paired-end (default) or single-end.\n");
    fprintf(stderr, "    If a singleton read is found in paired-end mode, the code will break nicely.\n");
    fprintf(stderr, "Note 3, defaults to dupsifter.stat if streaming or (-o basename).dupsifter.stat \n");
    fprintf(stderr, "    if the -o option is provided. If -o and -O are provided, then -O will be used.\n");
    fprintf(stderr, "\n");

    return 1;
}

int main(int argc, char *argv[]) {

    int c;
    ds_conf_t conf = {0};
    ds_conf_init(&conf);

    static const struct option loptions[] = {
        {"output"         , required_argument, NULL, 'o'},
        {"stats-output"   , required_argument, NULL, 'O'},
        {"single-end"     , no_argument      , NULL, 's'},
        {"add-mate-tags"  , no_argument      , NULL, 'm'},
        {"wgs-only"       , no_argument      , NULL, 'W'},
        {"max-read-length", required_argument, NULL, 'l'},
        {"min-base-qual"  , required_argument, NULL, 'b'},
        {"remove-dups"    , no_argument      , NULL, 'r'},
        {"verbose"        , no_argument      , NULL, 'v'},
        {"help"           , no_argument      , NULL, 'h'},
        {"version"        , no_argument      , NULL, 1},
        {NULL, 0, NULL, 0}
    };

    if (argc < 1) { return usage(&conf); }
    while ((c = getopt_long(argc, argv, "b:l:o:O:Wrsmqvh", loptions, NULL)) >= 0) {
        switch (c) {
            case 'b': conf.min_base_qual = (uint32_t)atoi(optarg); break;
            case 'l': conf.max_length = (uint32_t)atoi(optarg); break;
            case 'o': conf.outfn = optarg; break;
            case 'O': conf.statfn = strdup(optarg); break;
            case 'W': conf.is_wgs = 1; break;
            case 'r': conf.rm_dup = 1; break;
            case 's': conf.single_end = 1; break;
            case 'm': conf.add_mate_tags = 1; break;
            case 'v': conf.verbose = 1; break;
            case 'h': return usage(&conf);
            case 1: print_version(); return 0;
            default:
                return usage(&conf);
        }
    }

    conf.reffn = optind < argc ? argv[optind++] : NULL;
    conf.infn  = optind < argc ? argv[optind++] : "-";
    if (!conf.reffn) {
        fprintf(stderr, "Please provide a reference\n");
        return usage(&conf);
    }

    // Setup statistics filename
    if (conf.statfn == NULL) {
        int len = (strcmp(conf.outfn, "-") != 0) ? strlen(conf.outfn)+strlen(TAG_STAT_NAME)+1 : strlen(DEF_STAT_NAME);

        if (strcmp(conf.outfn, "-") != 0) {
            conf.statfn = (char *)calloc(len, sizeof(char));
            strncpy(conf.statfn, conf.outfn, strlen(conf.outfn)-4);
            strcat(conf.statfn, TAG_STAT_NAME);
        } else {
            conf.statfn = strdup(DEF_STAT_NAME);
        }
    }

    if (conf.verbose) {
        // Input file location
        if (strcmp(conf.infn, "-") == 0) {
            fprintf(stderr, "Reading input from stdin\n");
        } else {
            fprintf(stderr, "Reading input from %s\n", conf.infn);
        }

        // Output file location
        if (strcmp(conf.outfn, "-") == 0) {
            fprintf(stderr, "Writing output to stdout\n");
        } else {
            fprintf(stderr, "Writing output to %s\n", conf.outfn);
        }

        // Metrics file location
        fprintf(stderr, "Writing stats to %s\n", conf.statfn);
    }

    // Create argument list string
    if (!(conf.arg_list = stringify_argv(argc, argv))) {
        fatal_error("[dupsifter] Error: Unable to create argument list for PG string\n");
    }

    // Set write mode, auto-detect from file name, assume no extra compression
    // or other extra additions to the write mode
    sam_open_mode(conf.out_mode+1, conf.outfn, NULL);

    double t1 = get_current_time();
    double c1 = get_cpu_runtime();
    dupsifter(&conf);
    double c2 = get_cpu_runtime();
    double t2 = get_current_time();

    // Print output of run
    ds_conf_print(&conf);
    fprintf(stderr, "[dupsifter:%s] Wall time: %.3f seconds, CPU time: %.3f seconds\n",
            __func__, t2-t1, c2-c1);

    ds_conf_destroy(&conf);

    return 0;
}
