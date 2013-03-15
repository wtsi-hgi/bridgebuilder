/*
 * binnie_process.h - main binnie data processing
 *
 * Copyright (c) 2013 Genome Research Ltd. 
 * Author: Joshua C. Randall <jcrandall@alum.mit.edu>
 *
 * This file is part of BridgeBuilder. 
 *
 * BridgeBuilder is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 3 of the License, or (at your option) any later 
 * version. 
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  
 * 
 * You should have received a copy of the GNU General Public License along with 
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef BINNIE_PROCESS_H
#define BINNIE_PROCESS_H

/* gnulib headers */
#include <stdbool.h>
#include "gl_xlist.h"

/* binnie includes */
#include "binnie.h"

/* htslib for sam/bam processing */
#include <htslib/sam.h>

typedef enum {
  unchanged = BINNIE_UNCHANGED,
  bridged   = BINNIE_BRIDGED,
  remap     = BINNIE_REMAP,
} binnie_bin_t;

typedef struct {
  bool        bam_read_present;
  bam1_t      *bam_read;
} binnie_read_t;

typedef struct binnie_binned_read_t_ {
  binnie_read_t *br;
  binnie_bin_t  bin;
  int32_t       expected_mate_count;
  int32_t       mate_count;
  int32_t       original_refid;
  int32_t       original_pos;
  struct binnie_binned_read_t_ *next_mate;
  struct binnie_binned_read_t_ *prev_mate;
} binnie_binned_read_t;


bool binnie_process(int buffer_size, int max_buffer_bases, samFile *original_in_fp, samFile *bridge_in_fp, samFile *unchanged_out_fp, samFile *bridged_out_fp, samFile *remap_out_fp);

binnie_binned_read_t *binnie_read_bin(binnie_read_t *original_read, binnie_read_t *bridge_read);

void fixup_bridge_from_original (binnie_read_t *bridge_read, binnie_read_t *original_read);

void binnie_read_buffer (binnie_binned_read_t *bbr, gl_list_t output_buffer);

int32_t br_get_refid (const binnie_read_t *br);

int32_t br_get_pos (const binnie_read_t *br);

int32_t br_get_mapq (const binnie_read_t *br);

int32_t br_get_segment_index (const binnie_read_t *br);

int32_t br_get_num_segments (const binnie_read_t *br);

char *br_get_read_group (const binnie_read_t *br);

char *br_get_qname (const binnie_read_t *br);

bool br_equals(const binnie_read_t *br1, const binnie_read_t *br2);

char *br_get_uid_alloc(const binnie_read_t *br);

void br_dispose(const binnie_read_t *br);

binnie_read_t *br_init();

char *bbr_get_bin_name(binnie_binned_read_t *bbr);

static bool bbr_equals(const void *elt1, const void *elt2);

size_t bbr_hashcode(const void *elt);

void bbr_dispose(const void *elt);

binnie_binned_read_t *bbr_init(binnie_read_t *br);


#endif
