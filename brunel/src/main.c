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

struct parsed_opts {
    size_t input_count;
    char** input_name;
    samFile** input_file;
    bam_hdr_t** input_header;
    char* output_name;
    samFile* output_file;
    bam_hdr_t* output_header;
};

typedef struct parsed_opts parsed_opts_t;

parsed_opts_t* parse_args(int argc, char** argv) {
    if (argc < 3) {
        printf("Arguments should be: merge <input1.bam> <input2.bam> [<inputX.bam> ...] <output.bam>\r\n");
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

bool init(parsed_opts_t* opts) {
    // Open files
    opts->input_file = (samFile**)calloc(opts->input_count, sizeof(samFile*));
    opts->input_header = (bam_hdr_t**)calloc(opts->input_count, sizeof(bam_hdr_t*));
    for (size_t i = 0; i < opts->input_count; i++) {
        opts->input_file[i] = sam_open(opts->input_name[i], "rb", 0);
        if (opts->input_file[i] == NULL) {
            printf("Could not open input file: %s\r\n", opts->input_name[i]);
            return false;
        }
        opts->input_header[i] = sam_hdr_read(opts->input_file[i]);
    }
    
    // TODO: merge headers instead of just taking first one
    // TODO: create SQ translation table
    // TODO: create RG translation table
    opts->output_header = opts->input_header[0];

    opts->output_file = sam_open(opts->output_name, "wb", 0);
    
    if (opts->output_file == NULL) {
        printf("Could not open output file: %s\r\n", opts->output_name);
        return false;
    }
    
    return true;
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

bool merge(parsed_opts_t* opts) {
    if (sam_hdr_write(opts->output_file, opts->output_header) != 0) {
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
    for (size_t i = 0; i < opts->input_count; i++) { if (file_read[i]) { bam_destroy1(file_read[i]); } }

    return true;
}

void cleanup(parsed_opts_t* opts) {
    sam_close(opts->output_file);
    free(opts->output_name);
    for (size_t i = 0; i < opts->input_count; i++) {
        sam_close(opts->input_file[i]);
        free(opts->input_name[i]);
    }
}

int main(int argc, char** argv) {

    parsed_opts_t* opts = parse_args(argc, argv);
    if (!opts || !init(opts)) return -1;
    
    if (!merge(opts)) return -1;
    
    cleanup(opts);
      
    return 0;
}
