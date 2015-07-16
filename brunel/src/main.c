// Copyright (c) 2013 Genome Research Limited.
//
// This file is part of Brunel which is part of BridgeBuilder.
//
// Brunel is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program. If not, see L<http://www.gnu.org/licenses/>.

// Linux fixes
#define _GNU_SOURCE
#ifndef SIZE_T_MAX
#define SIZE_T_MAX      (~((size_t) 0))
#endif

#include "config.h"

#include <htslib/sam.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <unistd.h>

struct parsed_opts {
    char* output_header_name;
    size_t input_count;
    char** input_name;
    char** input_trans_name;
    char* output_name;
};

typedef struct parsed_opts parsed_opts_t;

struct state {
    size_t input_count;
    samFile** input_file;
    bam_hdr_t** input_header;
    int32_t** input_trans;
    samFile* output_file;
    bam_hdr_t* output_header;
};

typedef struct state state_t;

void cleanup_state(state_t* status);
void cleanup_opts(parsed_opts_t* opts);


parsed_opts_t* parse_args(int argc, char** argv) {
    if (argc < 3) {
        dprintf(STDERR_FILENO, "Arguments should be: brunel <newheader.sam> <input1.bam[:trans_tbl.txt]> <input2.bam[:trans_tbl.txt]> [<inputX.bam[:trans_tbl.txt]> ...] <output.bam>\r\n");
        return NULL;
    }
    
    parsed_opts_t* retval = malloc(sizeof(parsed_opts_t));
    if (! retval ) return NULL;
    
    retval->output_header_name = strdup(argv[1]);

    retval->input_count = argc-3;
    retval->input_name = (char**)calloc(retval->input_count,sizeof(char*));
    retval->input_trans_name = (char**)calloc(retval->input_count,sizeof(char*));
    size_t i = 0;
    for (; i < retval->input_count; i++) {
        char* temp = strdup(argv[i+2]);
        char* sep = temp;
        retval->input_name[i] = strsep(&sep, ":");
        retval->input_trans_name[i] = sep;
    }

    retval->output_name = strdup(argv[i+2]);

    return retval;
}

int* build_translation_file(const char* trans_name, bam_hdr_t* file_header, bam_hdr_t* replace_header) {
    FILE* trans_file = fopen(trans_name, "r");
    
    int file_entries = file_header->n_targets;
    int replace_entries = replace_header->n_targets;
    
    int *trans = malloc(sizeof(int)*file_entries);
    
    char* linepointer = NULL;
    size_t read = 0;
    
    int counter = file_entries;
    
    while (!feof(trans_file) && !ferror(trans_file) && counter > 0) {
        getline(&linepointer, &read, trans_file);
        char* sep = linepointer;
        strsep(&sep, "\t");
        if (sep == NULL) break;
        char* two = sep;
        strsep(&two, "\t\n");

        // lookup tid of original and replacement
        int i = 0;
        for ( ; i < file_entries; i++ ) {
            char* item = file_header->target_name[i];
            if (!strcmp(item,linepointer)) { break; }
        }

        if(i < file_entries) {
            // A tid requires replacement
            int j = 0;
            for ( ; j < replace_entries; j++ ) {
                char* item = replace_header->target_name[j];
                if (!strcmp(item,sep)) { break; }
            }
            trans[i] = j;
        }
        counter--;
    }
    free(linepointer);
    
    fclose(trans_file);
    return trans;
}

int* build_translation( bam_hdr_t* file_header, bam_hdr_t* replace_header ) {
    int file_entries = file_header->n_targets;
    int replace_entries = replace_header->n_targets;
    
    int *trans = malloc(sizeof(int)*file_entries);
    bool exact_match = true;
    for (int i = 0; i < file_entries; i++) {
        if (!strcmp(file_header->target_name[i], replace_header->target_name[i])) {
            trans[i] = i;
        } else {
            exact_match = false;
            bool error = true;
            for (int j = 0; j < replace_entries; j++) {
                if (!strcmp(file_header->target_name[i], replace_header->target_name[j])) {
                    trans[i] = j;
                    error = false;
                    break;
                }
            }
            if (error) {
	      dprintf(STDERR_FILENO, "Translation table entry missing for entry %d. file SQ: [%s]. output SQ: [%s]\n", i, file_header->target_name[i], replace_header->target_name[i]);
                exit(-1);
            }
        }
    }

    if (!exact_match) {
        return trans;
    } else {
        free(trans);
        return NULL;
    }
}

state_t* init(parsed_opts_t* opts) {
    state_t* retval = malloc(sizeof(state_t));
    if (!retval) {
        dprintf(STDERR_FILENO, "Out of memory\n");
        return NULL;
    }

    // TODO: add option to merge headers instead of just asking for one
    // TODO: create RG translation table
    samFile* hdr_load = sam_open(opts->output_header_name, "r", 0);
    if (!hdr_load) {
        dprintf(STDERR_FILENO, "Could not open header file\n");
        return NULL;
    }
    retval->output_header = sam_hdr_read(hdr_load);
    sam_close(hdr_load);
    if (!(retval->output_header->n_targets > 0)) {
      dprintf(STDERR_FILENO, "Header has no SQ targets, pointless to proceed!\n");
      return NULL;
    }
    retval->output_file = sam_open(opts->output_name, "wb", 0);
    
    if (retval->output_file == NULL) {
        printf("Could not open output file: %s\r\n", opts->output_name);
        cleanup_state(retval);
        return NULL;
    }

    retval->input_count = opts->input_count;
    
    // Open files
    retval->input_trans = (int**)calloc(opts->input_count, sizeof(int*));
    retval->input_file = (samFile**)calloc(opts->input_count, sizeof(samFile*));
    retval->input_header = (bam_hdr_t**)calloc(opts->input_count, sizeof(bam_hdr_t*));
    if (!retval->input_file || !retval->input_header) {
        dprintf(STDERR_FILENO, "Out of memory");
        free(retval);
        return NULL;
    }
    for (size_t i = 0; i < opts->input_count; i++) {
        retval->input_file[i] = sam_open(opts->input_name[i], "rb", 0);
        if (retval->input_file[i] == NULL) {
            dprintf(STDERR_FILENO, "Could not open input file: %s\r\n", opts->input_name[i]);
            return NULL;
        }
        retval->input_header[i] = sam_hdr_read(retval->input_file[i]);
        if (opts->input_trans_name[i] != NULL)
        {
            retval->input_trans[i] = build_translation_file(opts->input_trans_name[i], retval->input_header[i], retval->output_header);
        }
        else
        {
            retval->input_trans[i] = build_translation(retval->input_header[i], retval->output_header);
        }
    }

    return retval;
}

size_t selectRead( bam1_t **file_read, size_t input_count )
{
    assert(input_count != 0);
    // need to find element with lowest tid and pos
    size_t min = SIZE_T_MAX;
    // horrible hack of the day
    // treat the int32_t tid as a uint32_t so that -1 (aka unmapped) is treated as UINT32_MAX
    uint32_t tid_min;
    int32_t pos_min;
    // load initial value
    size_t i = 0;
    for (; i < input_count; i++) {
        if (file_read[i] != NULL) {
            tid_min = (uint32_t)file_read[i]->core.tid;
            pos_min = file_read[i]->core.pos;
            min = i;
            break;
        }
    }
    assert(min != SIZE_T_MAX); // No valid files?

    // then resume our search
    for (;i < input_count; i++) {
        if (file_read[i] != NULL) {
            // To complicate matters tid == -1 is a special value which should always go last
            if ((tid_min > (uint32_t)file_read[i]->core.tid ) ||
                (tid_min == (uint32_t)file_read[i]->core.tid && pos_min > file_read[i]->core.pos)) {
                tid_min = (uint32_t)file_read[i]->core.tid;
                pos_min = file_read[i]->core.pos;
                min = i;
            }
        }
    }
    
    assert(min != SIZE_T_MAX); // No valid files?
    
    return min;
}

bool merge(state_t* opts) {
    if (sam_hdr_write(opts->output_file, opts->output_header) != 0) {
        dprintf(STDERR_FILENO, "Could not write output file header\n");
        return false;
    }
    
    bam1_t** file_read = calloc(opts->input_count, sizeof(bam1_t*));
    size_t files_to_merge = opts->input_count;
    // initialise the first read for each input file
    for (size_t i = 0; i < opts->input_count; i++) {
        file_read[i] = bam_init1();  
        // Read the first record
        if (sam_read1(opts->input_file[i], opts->input_header[i], file_read[i]) < 0) {
            // Nothing more to read?  Ignore this file
            bam_destroy1(file_read[i]);
            file_read[i] = NULL;
            files_to_merge--;
        } else {
            if (opts->input_trans[i]) {
                // Translate the tid and mate tid but only if they're not null values
                if (file_read[i]->core.tid != -1) {
                    file_read[i]->core.tid = opts->input_trans[i][file_read[i]->core.tid];
                }
                if (file_read[i]->core.mtid != -1) {
                    file_read[i]->core.mtid = opts->input_trans[i][file_read[i]->core.mtid];
                }
            }
        }
    }

    while (files_to_merge > 0) {
        size_t i = selectRead(file_read, opts->input_count);
        // Write the read out and replace it with the next one to process
        sam_write1(opts->output_file, opts->output_header, file_read[i]);
        if (sam_read1(opts->input_file[i], opts->input_header[i], file_read[i]) < 0) {
            // Nothing more to read?  Ignore this file in future
            bam_destroy1(file_read[i]);
            file_read[i] = NULL;
            files_to_merge--;
        } else {
            if (opts->input_trans[i]) {
                // Translate the tid and mate tid but only if they're not null values
                if (file_read[i]->core.tid != -1) {
                    file_read[i]->core.tid = opts->input_trans[i][file_read[i]->core.tid];
                }
                if (file_read[i]->core.mtid != -1) {
                    file_read[i]->core.mtid = opts->input_trans[i][file_read[i]->core.mtid];
                }
            }
        }
    }

    // Clean up
    for (size_t i = 0; i < opts->input_count; i++) {
        if (file_read[i]) { bam_destroy1(file_read[i]); }
    }

    return true;
}

void cleanup_state(state_t* status) {
    sam_close(status->output_file);
    for (size_t i = 0; i < status->input_count; i++) {
        sam_close(status->input_file[i]);
    }
    free(status->input_file);
}

void cleanup_opts(parsed_opts_t* opts) {
    free(opts->output_name);
    for (size_t i = 0; i < opts->input_count; i++) {
        free(opts->input_name[i]);
    }
    free(opts->input_name);
}

int main(int argc, char** argv) {

    parsed_opts_t* opts = parse_args(argc, argv);
    if (!opts ) return -1;
    state_t* status = init(opts);
    if (!status) return -1;
    
    if (!merge(status)) return -1;
    
    cleanup_state(status);
    cleanup_opts(opts);
      
    return 0;
}
