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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

/* gnulib headers */
#include "error.h"
#include "progname.h"
#include "xalloc.h"
#include "version-etc.h"

/* internationalisation */
#include "gettext.h"

/* binnie includes */
#include "binnie.h"
#include "binnie_log.h"
#include "binnie_files.h"
#include "binnie_process.h"

/* copyright notice for --version output (%s is symbol and %d is year) */
const char version_etc_copyright[] = "Copyright %s %d Genome Research Limited";

/* max size of output buffer (in number of reads or 0 for unlimited) */
unsigned int buffer_size;

/* max size of output buffer (in number of mapped bases or 0 for unlimited) */
unsigned int max_buffer_bases;

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
  fprintf(stderr, gettext("  -u, --unchanged_out          Filename of output bin (.bam/.sam) for original reads which did not map to bridge\n"));
  fprintf(stderr, gettext("  -b, --bridged_out            Filename of output bin (.bam/.sam) for reads that have been newly mapped to bridge\n"));
  fprintf(stderr, gettext("  -r, --remap-out              Filename of output bin (.bam/.sam) for reads that need remapping against the full reference\n"));
  fprintf(stderr, gettext("  -s, --buffer_size            Size of output buffer (in reads) [default: %d]\n"), BINNIE_DEFAULT_BUFFER_SIZE);
  fprintf(stderr, gettext("  -m, --max_buffer_bases       Size of output buffer (in bases) [default: %d]\n"), BINNIE_DEFAULT_BUFFER_BASES);
  fprintf(stderr, gettext("  -i, --ignore_rg              Ignore read group (RG) when matching reads between original and bridge\n"));
  fprintf(stderr, gettext("  -a, --allow_sorted_unmapped  Allow reads with flag 0x4 set to be sorted according to their refid and pos\n"));
  fprintf(stderr, gettext("  -h, --help                   Print short help message and exit\n"));
  fprintf(stderr, gettext("  -v, --verbose[=level]        Increase/Set level of verbosity (-vvv sets level 3 as does --verbose=3)\n"));
#ifdef DEBUG
  fprintf(stderr, gettext("  -d, --debug                  Print debugging messages to stderr (also sets -v 3)\n"));
#endif
  fprintf(stderr, gettext("  -V, --version                Print version information to stdout and exit\n"));
}

/* 
 * main
 * ----
 *
 * Sets defaults, processes command-line options and arguments, opens 
 * input and output files, calls binnie_process to do the processing, 
 * cleans up, and exits.
 *
 */
int main(int argc, char **argv) 
{

  /* setup progname */
  set_program_name (argv[0]);

  /* init globals */
  verbosity = 0;
#ifdef DEBUG
  debug_flag = false;
#endif
  ignore_rg = false;
  allow_sorted_unmapped = false;
  buffer_size = BINNIE_DEFAULT_BUFFER_SIZE;
  max_buffer_bases = BINNIE_DEFAULT_BUFFER_BASES;
  unchanged_out_file = NULL;
  bridged_out_file = NULL;
  remap_out_file = NULL;
  
  DLOG("main: started");

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
	  {"unchanged_out",		required_argument,	0,	'u'},
	  {"bridged_out",		required_argument,	0,	'b'},
	  {"remap_out",			required_argument,	0,	'r'},
	  {"buffer_size",		required_argument,	0,	's'},
	  {"max_buffer_bases",  	required_argument,	0,	'm'},
	  {"ignore_rg",         	no_argument,            0,      'i'},
	  {"allow_sorted_unmapped",    	no_argument,            0,      'a'},
	  {"help",			no_argument,		0,	'h'},
 	  {"verbose",	        	optional_argument,	0,	 0 },
 	  {"verbose",           	no_argument,		0,	'v'},
#ifdef DEBUG
	  {"debug",			no_argument,		0,	'd'},
#endif
 	  {"version",           	no_argument,		0,	'V'},
	  {0, 0, 0, 0}
	};
      option_index = 0;
      
      c = getopt_long(argc, argv, "u:b:r:s:m:iahvdV", binnie_options, &option_index);

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
	case 's':
	  buffer_size = atoi(optarg);
	  break;
	case 'm':
	  max_buffer_bases = atoi(optarg);
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
	case 'i':
	  ignore_rg = true;
	  break;
	case 'a':
	  allow_sorted_unmapped = true;
	  break;
	case 'h':
	  print_help();
	  exit(BINNIE_EXIT_SUCCESS);
	  break;
	case 'v': 
	  verbosity++;
	  break;
#ifdef DEBUG
	case 'd':
	  debug_flag = true;
	  break;
#endif
	case 'V':
	  version_etc(stdout, NULL, PACKAGE_NAME, PACKAGE_VERSION, "Joshua C. Randall", (char *) 0);
	  exit(BINNIE_EXIT_SUCCESS);
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

#ifdef DEBUG
  if (debug_flag)
    {
      error(0, 0, gettext("printing debugging messages"));
    }
#endif

  if (ignore_rg)
    {
      blog(0, gettext("ignoring read group (RG) when matching original and bridge-mapped reads"));
    }

  if (allow_sorted_unmapped)
    {
      blog(0, gettext("allowing reads with 0x4 flag set to be sorted according to their refid and pos"));
    }

  if (buffer_size > 0)
    {
      blog(0, gettext("buffer size set to %d reads"), buffer_size);
    }

  if (max_buffer_bases > 0)
    {
      blog(0, gettext("max buffer bases set to %d bases"), max_buffer_bases);
    }

  /* get remaining command-line arguments (original and bridge input file names) */
  if (optind + 2 != argc) {
    print_usage();
    errx(BINNIE_EXIT_ERR_ARGS, gettext("two filenames should be given as arguments following the options"));
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
      err(BINNIE_EXIT_ERR_IN_FILES, gettext("could not open one or more input files"));
    }
  else 
    {
      blog(0, gettext("input files opened"));
      blog(1, gettext("\toriginal=[%s]"), original_in_fp->fn);
      blog(1, gettext("\tbridge=[%s]"), bridge_in_fp->fn);
    }

  /* name output files after input if they aren't specified */
  if (unchanged_out_file == NULL) 
    {
      DLOG(gettext("overriding unchanged_out_file with original_in_file + unchanged_out_suffix"));
      unchanged_out_file = xmalloc(strlen(original_in_file) + strlen(unchanged_out_suffix) + 1);
      strcat(unchanged_out_file, original_in_file);
      strcat(unchanged_out_file, unchanged_out_suffix);
    }
  blog(3, gettext("unchanged_out_file set to %s"), unchanged_out_file);

  if (bridged_out_file == NULL) 
    {
      DLOG(gettext("overriding unchanged_out_file with original_in_file + bridged_out_suffix"));
      bridged_out_file = xmalloc(strlen(original_in_file) + strlen(bridged_out_suffix) + 1);
      strcat(bridged_out_file, original_in_file);
      strcat(bridged_out_file, bridged_out_suffix);
    }
  blog(3, gettext("bridged_out_file set to %s"), bridged_out_file);

  if (remap_out_file == NULL) 
    {
      DLOG(gettext("overriding unchanged_out_file with original_in_file + remap_out_suffix"));
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
      err(BINNIE_EXIT_ERR_OUT_FILES, gettext("could not open one or more output files"));
    }
  else 
    {
      blog(0, gettext("output files opened"));
      blog(1, gettext("\tunchanged=[%s]"), unchanged_out_fp->fn);
      blog(1, gettext("\tbridged=[%s]"), bridged_out_fp->fn);
      blog(1, gettext("\tremap=[%s]"), remap_out_fp->fn);
    }


  /* process data */
  blog(1, gettext("beginning binnie processing"));
  binnie_process(buffer_size, max_buffer_bases, original_in_fp, bridge_in_fp, unchanged_out_fp, bridged_out_fp, remap_out_fp);


  /* clean up */
  blog(1, gettext("cleaning up"));

  blog(2, gettext("closing open files"));
  binnie_close(original_in_fp);
  binnie_close(bridge_in_fp);
  binnie_close(unchanged_out_fp);
  binnie_close(bridged_out_fp);
  binnie_close(remap_out_fp);

  blog(2, gettext("freeing memory"));
  free(original_in_file);
  free(bridge_in_file);
  free(unchanged_out_file);
  free(bridged_out_file);
  free(remap_out_file);


  blog(1, gettext("finished!"));
  DLOG(gettext("main: returning"));
  return 0;
}

