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

#include "config.h"

#include <err.h>
#include <errno.h>
#include <getopt.h>
#include <locale.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

/* gnulib headers */
#include "error.h"
#include "progname.h"
#include "xalloc.h"

#include "gettext.h"

/* binnie includes */
#include "binnie_log.h"
#include "binnie_files.h"

#include <htslib/sam.h>

/* filename of input (BAM/SAM) consisting of original reads */
 char *original_in_file;

/* filename of input (BAM/SAM) consisting of original reads mapped to the bridge (in original order) */
 char *bridge_in_file;

/* filename of output bin (BAM/SAM) of original reads which have not been changed*/
 char *unchanged_out_file;

/* filename of output bin (BAM/SAM) of reads that have been newly-mapped to bridge */
 char *bridged_out_file;

/* filename of output bin (BAM/SAM) of reads that must be remapped to the full reference */
 char *remap_out_file;

/* pointers to opened samFile-s corresponding to above filenames */
 samFile *original_in_fp;
 samFile *bridge_in_fp;
 samFile *unchanged_out_fp;
 samFile *bridged_out_fp;
 samFile *remap_out_fp;

/* suffix to add to original input file to get an output name if unchanged_out_file was not specified */
 const char *unchanged_out_suffix = "_unchanged.bam";

/* suffix to add to original input file to get an output name if bridged_out_file was not specified */
 const char *bridged_out_suffix = "_bridged.bam";

/* suffix to add to original input file to get an output name if remap_out_file was not specified */
 const char *remap_out_suffix = "_remap.bam";


void print_usage() 
{
  fprintf(stderr, gettext("Usage: %s [options] <original(bam|sam)> <bridge(bam|sam)>\n"), program_name);
}

void print_help() 
{
  print_usage();
  fprintf(stderr, gettext("Options: \n"));
  fprintf(stderr, gettext("  -u, --unchanged_out    Filename of output bin (.bam/.sam) for original reads which did not map to bridge\n"));
  fprintf(stderr, gettext("  -b, --bridged_out      Filename of output bin (.bam/.sam) for reads that have been newly mapped to bridge\n"));
  fprintf(stderr, gettext("  -r, --remap-out        Filaname of output bin (.bam/.sam) for reads that need remapping against the full reference\n"));
  fprintf(stderr, gettext("  -h, --help             Print short help message and exit\n"));
  fprintf(stderr, gettext("  -v, --verbose[=level]  Increase/Set level of verbosity (-vvv sets level 3 as does --verbose=3)\n"));
  fprintf(stderr, gettext("  -d, --debug            Print debugging messages to stderr (also sets -v 3)\n"));
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

  /* setup progname */
  set_program_name (argv[0]);

  /* init globals */
  verbosity = 0;
  debug_flag = 0;
  unchanged_out_file = NULL;
  bridged_out_file = NULL;
  remap_out_file = NULL;
  
  blog(9, "main: started");

  /* setup gettext */
  //setlocale (LC_ALL, "");
  //bindtextdomain (PACKAGE, LOCALEDIR);
  //textdomain (PACKAGE);


  /* get command-line options */
  while (1)
    {
      int c;
      int option_index;
      static struct option binnie_options[] =
	{
 	  {"verbose",	        optional_argument,	0,	 0 },
 	  {"verbose",           no_argument,		0,	'v'},
	  {"debug",		no_argument,		0,	'd'},
	  {"help",		no_argument,		0,	'h'},
	  {"unchanged_out",	required_argument,	0,	'u'},
	  {"bridged_out",	required_argument,	0,	'b'},
	  {"remap_out",		required_argument,	0,	'r'},
	  {0, 0, 0, 0}
	};
      option_index = 0;
      
      c = getopt_long(argc, argv, "vdhu:b:r:", binnie_options, &option_index);

      if (c < 0)
	break;
      
      switch (c)
	{
	case 0:
	  if (!strcmp(binnie_options[option_index].name, "verbose")) {
	    if (optarg)
	      verbosity = atoi(optarg);
	    else 
	      verbosity++;
	  }
	  break;
	case 'v': 
	  verbosity++;
	  break;
	case 'd':
	  debug_flag = 1;
	  break;
	case 'h':
	  print_help();
	  exit(0);
	  break;
	case 'u':
	  unchanged_out_file = xstrdup(optarg);
	  break;
	case 'b':
	  bridged_out_file = xstrdup(optarg);
	  break;
	case 'r':
	  remap_out_file = xstrdup(optarg);
	  break;
	case '?':
	  /* getopt_long will have already printed an error */
	  print_usage();
	  break;
	default:
	  error(0, 0, gettext("unhandled option [-%c]"), c);
	  print_usage();
	}
    }

  if (verbosity > 0) 
    {
      error(0, 0, gettext("verbosity set to %u"), verbosity);
    }

  if (debug_flag != 0)
    {
      error(0, 0, gettext("printing debugging messages"));
    }

  /* get remaining command-line arguments (original and bridge input file names) */
  if (optind + 2 != argc) {
    print_usage();
    err(1, gettext("two filenames should be given as arguments following the options"));
  } else {
    original_in_file = xstrdup(argv[optind++]);
    blog(3, gettext("original_in_file set to %s"), original_in_file);
    bridge_in_file = xstrdup(argv[optind++]);
    blog(3, gettext("bridge_in_file set to %s"), bridge_in_file);
  }
  

  /* open BAM/SAM input files */
  original_in_fp = binnie_open_in(original_in_file);
  bridge_in_fp = binnie_open_in(bridge_in_file);

  /* check if input files are open */
  if (original_in_fp <= 0 || bridge_in_fp <= 0) 
    {
      /* print error and exit */
      err(2, gettext("could not open one or more input files"));
    }
  else 
    {
      blog(1, gettext("input files opened"));
      blog(2, gettext("\toriginal=[%s]"), original_in_fp->fn);
      blog(2, gettext("\tbridge=[%s]"), bridge_in_fp->fn);
    }

  /* name output files after input if they aren't specified */
  if (unchanged_out_file == NULL) 
    {
      blog(9, gettext("overriding unchanged_out_file with original_in_file + unchanged_out_suffix"));
      unchanged_out_file = xmalloc(strlen(original_in_file) + strlen(unchanged_out_suffix) + 1);
      strcat(unchanged_out_file, original_in_file);
      strcat(unchanged_out_file, unchanged_out_suffix);
    }
  blog(3, gettext("unchanged_out_file set to %s"), unchanged_out_file);

  if (bridged_out_file == NULL) 
    {
      blog(9, gettext("overriding unchanged_out_file with original_in_file + bridged_out_suffix"));
      bridged_out_file = xmalloc(strlen(original_in_file) + strlen(bridged_out_suffix) + 1);
      strcat(bridged_out_file, original_in_file);
      strcat(bridged_out_file, bridged_out_suffix);
    }
  blog(3, gettext("bridged_out_file set to %s"), bridged_out_file);

  if (remap_out_file == NULL) 
    {
      blog(9, gettext("overriding unchanged_out_file with original_in_file + remap_out_suffix"));
      remap_out_file = xmalloc(strlen(original_in_file) + strlen(remap_out_suffix) + 1);
      strcat(remap_out_file, original_in_file);
      strcat(remap_out_file, remap_out_suffix);
    }
  blog(3, gettext("remap_out_file set to %s"), remap_out_file);
  
  
  /* open BAM/SAM output files */
  unchanged_out_fp = binnie_open_out(unchanged_out_file);
  bridged_out_fp = binnie_open_out(bridged_out_file);
  remap_out_fp = binnie_open_out(remap_out_file);

  /* check if output files are open */
  if (unchanged_out_fp <= 0 || bridged_out_fp <= 0 || remap_out_fp <= 0)
    {
      /* print error and exit */
      err(3, gettext("could not open one or more output files"));
    }
  else 
    {
      blog(1, gettext("output files opened"));
      blog(2, gettext("\tunchanged=[%s]"), unchanged_out_fp->fn);
      blog(2, gettext("\tbridged=[%s]"), bridged_out_fp->fn);
      blog(2, gettext("\tremap=[%s]"), remap_out_fp->fn);
    }


  /* read BAM/SAM headers */
  blog(3, gettext("main: reading BAM/SAM headers"));
  bam_hdr_t *original_header = sam_hdr_read(original_in_fp);
  blog(9, gettext("main: original has %d targets"), original_header->n_targets);
  bam_hdr_t *bridge_header = sam_hdr_read(bridge_in_fp);
  blog(9, gettext("main: bridge has %d targets"), bridge_header->n_targets);


  /* write BAM/SAM headers */
  blog(3, gettext("main: writing BAM/SAM headers"));
  sam_hdr_write(unchanged_out_fp, original_header);
  sam_hdr_write(bridged_out_fp, original_header);
  sam_hdr_write(remap_out_fp, original_header);
  

  /* initialize BAM/SAM read buffers */
  blog(3, gettext("main: initializing BAM/SAM read buffers"));
  bam1_t* original_read = bam_init1();
  bam1_t* bridge_read = bam_init1();


  /* reading original and bridge */
  blog(3, gettext("main: beginning read processing"));
  int success = 1;
  do 
    {
      if (sam_read1(original_in_fp, original_header, original_read) < 0)
	{
	  error(0, errno, gettext("failed to read from original"));
	  success = 0;
	}
      if (sam_read1(bridge_in_fp, bridge_header, bridge_read) < 0) 
	{
	  error(0, errno, gettext("failed to read from bridge"));
	  success = 0;
	}
      
      if (success)
	{
	  /* for now just write original to all three files */
	  sam_write1(unchanged_out_fp, original_header, original_read);
	  sam_write1(bridged_out_fp, original_header, original_read);
	  sam_write1(remap_out_fp, original_header, original_read);
	}
    } while(success);
  
  
  blog(3, gettext("main: destroying BAM/SAM read buffers"));
  if (original_read) 
    bam_destroy1(original_read); 
  if (bridge_read) 
    bam_destroy1(bridge_read); 

  blog(3, gettext("main: closing open files"));
  binnie_close(original_in_fp);
  binnie_close(bridge_in_fp);
  binnie_close(unchanged_out_fp);
  binnie_close(bridged_out_fp);
  binnie_close(remap_out_fp);

  blog(3, gettext("main: freeing memory"));
  free(original_in_file);
  free(bridge_in_file);
  free(unchanged_out_file);
  free(bridged_out_file);
  free(remap_out_file);

  blog(9, gettext("main: returning"));
  return 0;
}

