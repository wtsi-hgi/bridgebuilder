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

#include <htslib/sam.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <unistd.h>

struct parsed_opts {
    size_t input_count;
    char** input_name;
    char* output_name;
};

typedef struct parsed_opts parsed_opts_t;

struct state {
    size_t input_count;
    samFile** input_file;
    bam_hdr_t** input_header;
    int** input_trans;
    samFile* output_file;
    bam_hdr_t* output_header;
};

typedef struct state state_t;

void cleanup_state(state_t* status);
void cleanup_opts(parsed_opts_t* opts);


parsed_opts_t* parse_args(int argc, char** argv) {
    if (argc < 3) {
        dprintf(STDERR_FILENO, "Arguments should be: merge <input1.bam> <input2.bam> [<inputX.bam> ...] <output.bam>\r\n");
        return NULL;
    }
    
    parsed_opts_t* retval = malloc(sizeof(parsed_opts_t));
    if (! retval ) return NULL;

    retval->input_count = argc-2;
    retval->input_name = (char**)calloc(retval->input_count,sizeof(char*));
    size_t i = 0;
    for (; i < retval->input_count; i++) {
        retval->input_name[i] = strdup(argv[i+1]);
    }

    retval->output_name = strdup(argv[i+1]);

    return retval;
}

state_t* init(parsed_opts_t* opts) {
    state_t* retval = malloc(sizeof(state_t));
    if (!retval) {
        dprintf(STDERR_FILENO, "Out of memory");
        return NULL;
    }
    
    retval->input_count = opts->input_count;
    
    // Open files
    retval->input_file = (samFile**)calloc(opts->input_count, sizeof(samFile*));
    retval->input_header = (bam_hdr_t**)calloc(opts->input_count, sizeof(bam_hdr_t*));
    for (size_t i = 0; i < opts->input_count; i++) {
        retval->input_file[i] = sam_open(opts->input_name[i], "rb", 0);
        if (retval->input_file[i] == NULL) {
            dprintf(STDERR_FILENO, "Could not open input file: %s\r\n", opts->input_name[i]);
            return NULL;
        }
        retval->input_header[i] = sam_hdr_read(retval->input_file[i]);
    }
    
    // TODO: merge headers instead of just taking first one
    // TODO: create SQ translation table
    // TODO: create RG translation table
    retval->output_header = retval->input_header[0];

    retval->output_file = sam_open(opts->output_name, "wb", 0);
    
    if (retval->output_file == NULL) {
        printf("Could not open output file: %s\r\n", opts->output_name);
        cleanup_state(retval);
        return NULL;
    }
    
    return retval;
}

size_t selectRead( bam1_t **file_read, size_t input_count )
{
    // need to find element with lowest tid and pos
    size_t min = SIZE_T_MAX;
    int32_t tid_min = INT32_MAX, pos_min = INT32_MAX;
    
    for (size_t i = 0; i < input_count; i++) {
        if (file_read[i] !=NULL) {
            if (tid_min > file_read[i]->core.tid ||
                (tid_min == file_read[i]->core.tid && pos_min > file_read[i]->core.pos)) {
                tid_min = file_read[i]->core.tid;
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
        dprintf(STDERR_FILENO, "Could not write output file header");
        return false;
    }
    
    bam1_t** file_read = calloc(opts->input_count, sizeof(bam1_t*));
    size_t files_to_merge = opts->input_count;
    for (size_t i = 0; i < opts->input_count; i++) {
        file_read[i] = bam_init1();
        if (sam_read1(opts->input_file[i], opts->input_header[i], file_read[i]) < 0) {
            bam_destroy1(file_read[i]);
            file_read[i] = NULL;
            files_to_merge--;
        }
    }

    while (files_to_merge > 0) {
        size_t i = selectRead(file_read, opts->input_count);
        sam_write1(opts->output_file, opts->output_header, file_read[i]);
        if (sam_read1(opts->input_file[i], opts->input_header[i], file_read[i]) < 0) {
            bam_destroy1(file_read[i]);
            file_read[i] = NULL;
            files_to_merge--;
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
