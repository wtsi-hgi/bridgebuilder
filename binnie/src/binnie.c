/*
 * binnie - processes original and bridge remapped BAM into four "bins"
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <stdint.h>
#include <errno.h>
#include <error.h>
#include <getopt.h>

#include "config.h"

#include "gettext.h"
#include <locale.h> 

#include <htslib/sam.h>

static int verbosity;
static int debug_flag;

char* file_name;
char* replace_name;
char* trans_name;
char* out_name;

samFile* file_iter;
samFile* replace_iter;
samFile* out_file;

bam_hdr_t* file_header;
bam_hdr_t* replace_header;

int* trans;


char* usage() 
{
  return gettext("Usage: binnie [options] <original(bam|sam)> <bridge(bam|sam)>");
}

/* 
 * binnie procedure
 * ----------------
 * 
 * 1. Individual Reads processed into a buffer containing the read data and a result bin: 
 * --------------------------------------
 * Original  Bridge    Bin
 * --------------------------------------
 * Unmapped  Unmapped  Unchanged
 * Unmapped  MQ >= 0   Bridged
 * MQ == 0   Unmapped  Unchanged
 * MQ == 0   MQ == 0   Unchanged*
 * MQ == 0   MQ > 0    Remap
 * MQ > 0    Unmapped  Unchanged
 * MQ > 0    MQ == 0   Remap
 * MQ > 0    MQ > 0    Remap
 * Deleted   (any)     Remap
 * --------------------------------------
 *
 * 2. Check if the read's pair is already in the buffer -- if it is, then update the 
 * other read to note that it's pair is here and change the result bin for both reads: 
 * ------------------------------------------
 * Result_1   Result_2   Bin_1  Bin_2
 * ------------------------------------------
 * Remap      (any)      Remap  Remap
 * (any)      Remap      Remap  Remap
 * Unchanged  Bridged    Remap  Remap
 * Bridged    Unchanged  Remap  Remap
 * ------------------------------------------
 *
 * 3. When the buffer is full, start writing to output bins, but first perform one final check: 
 * if a read's pair has not been added to buffer, then change its bin to Remap.
 *
 */

int main(int argc, char** argv) 
{
  int c;
  
  /* init globals */
  verbosity = 0;
  debug_flag = 0;
  
  /* setup gettext */
  setlocale (LC_ALL, "");
  bindtextdomain (PACKAGE, LOCALEDIR);
  textdomain (PACKAGE);
  
  while (1)
    {
      static struct option long_options[] =
	{
	  {"verbose",		optional_argument,	0,	'v'},
	  {"debug",		no_argument,		0,	'd'},
	  {"unchanged_out",	required_argument,	0,	'u'},
	  {"bridged_out",	required_argument,	0,	'b'},
	  {"remap_out",		required_argument,	0,	'r'},
	  {0, 0, 0, 0}
	};
      int option_index = 0;
      
      c = getopt_long(argc, argv, "vdu:b:r:", long_options, &option_index);

      if (c < 0)
	break;
      
      switch (c)
	{
	case 'v': 
	  break;
	case 'd':
	  break;
	case 'u':
	  unchanged_out_file = optarg;
	  break;
	case 'b':
	  bridged_out_file = optarg;
	  break;
	case 'r':
	  remap_out_file = optarg;
	  break;
	case '?':
	  /* getopt_long will have already printed an error */
	  break;
	default:
	  abort ();
	}
    }

  if (optind + 2 != argc) {
    perror(usage());
    return -1;
  } else {
    original_in_file = argv[optind++];
    bridge_in_file = argv[optind++];
  }
  
  file_name = argv[1];
  replace_name = argv[2];
  trans_name = argv[3];
  out_name = argv[4];

  // Open files
  file_iter = sam_open(file_name, "rb", 0);
  replace_iter = sam_open(replace_name, "r", 0);
  out_file = sam_open(out_name, "wb", 0);
  
  if (file_iter == NULL || replace_iter == NULL || out_file == NULL) {
    //TODO: usage();
    return -1;
  }
  
  // Iterate through BAM
  file_header = sam_hdr_read(file_iter);
  replace_header = sam_hdr_read(replace_iter);
  
  // build translation
  FILE* trans_file = fopen(trans_name, "r");
  
  int file_entries = file_header->n_targets;
  int replace_entries = replace_header->n_targets;
  
  //trans = new int[file_entries];
  trans = calloc(file_entries, sizeof(int));
  if (trans == NULL)
    {
      //todo close files
      perror("doh.\n");
      return -2;
    }
  
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
    int j = 0;
    for ( ; j < replace_entries; j++ ) {
      char* item = replace_header->target_name[j];
      if (!strcmp(item,sep)) { break; }
    }
    
    trans[i] = j;
    counter--;
  }
  free(linepointer);
  
  fclose(trans_file);
  // done building translation

  // parse
  if (sam_hdr_write(out_file, replace_header) != 0) {
    //throw new std::runtime_error("IEEEE!");
    //todo: could not write sam header
    perror("IEEEE!\n");
    return -3;
  }
  
  bam1_t* file_read = bam_init1();
  
  while (sam_read1(file_iter, file_header, file_read) >= 0) {
    file_read->core.tid = trans[file_read->core.tid];
    sam_write1(out_file, file_header, file_read);
  }
  
  // Clean up
  if (file_read) { bam_destroy1(file_read); }
  free(trans);
  // done parsing

  sam_close(out_file);
  sam_close(replace_iter);
  sam_close(file_iter);
  
  return 0;
}

