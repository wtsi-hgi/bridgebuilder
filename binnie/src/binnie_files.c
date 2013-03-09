/*
 * binnie_files.c - handles BAM/SAM input and output (via htslib)
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

#include <errno.h>
#include <getopt.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h> 

#include "error.h"
#include "xalloc.h"

#include <htslib/sam.h>

#include "binnie_files.h"
#include "binnie_log.h"

/*
 * binnie_open_out
 *
 * Opens the file named FILENAME for writing, with the mode (BAM/SAM)
 * depending on the file extension.
 *
 * Returns: a pointer to the opened samFile, or < 0 on error.
 */
samFile *binnie_open_out(const char *filename) 
{
  samFile *fp;
  int filename_len;

  DLOG("binnie_open_out: filename=[%s]", filename);

  if(!filename) 
    {
      error(0, 0, "binnie_open_out: null filename");
      return fp;
    }

  /*
   * check whether filename is bam or sam
   */
  filename_len = strlen(filename);
  if ( !strcasecmp(".bam", filename + filename_len - 4) )
    {
      fp = sam_open(filename, "wb", 0);
      if (!fp) 
	{
	  error(0, errno, "binnie_open_out: error opening [%s] as bam", filename);
	}
    }
  else if ( !strcasecmp(".sam", filename + filename_len - 4) )
    {
      fp = sam_open(filename, "w", 0);
      if (fp == NULL) 
	{
	  error(0, errno, "binnie_open_out: error opening [%s] as sam", filename);
	}
    }
  else
    {
      error(0, 0, "binnie_open_out: filename [%s] does not end in .bam or .sam", filename);
      return fp;
    }

  blog(3, "binnie_open_out: opened fp->fn=[%s]", fp->fn);

  DLOG("binnie_open_out: returning fp=[%u] for filename=[%s]", fp, filename);
  return fp;
}


/*
 * binnie_open_in
 *
 * Opens the file named FILENAME for reading, with the mode (BAM/SAM)
 * depending on the file extension.
 *
 * Returns: a pointer to the opened samFile, or < 0 on error.
 */
samFile* binnie_open_in(const char *filename) 
{
  samFile *fp;
  int filename_len;

  DLOG("binnie_open_in: filename=[%s]", filename);
  
  if(!filename) 
    {
      error(0, 0, "binnie_open_in: null filename");
      return fp;
    }

  /*
   * check whether filename is bam or sam
   */
  filename_len = strlen(filename);
  if ( !strcasecmp(".bam", filename + filename_len - 4) )
    {
      fp = sam_open(filename, "rb", 0);
      if (fp == NULL) 
	{
	  error(0, errno, "binnie_open_in: error opening [%s] as bam", filename);
	}
    }
  else if ( !strcasecmp(".sam", filename + filename_len - 4) )
    {
      fp = sam_open(filename, "r", 0);
      if (!fp) 
	{
	  error(0, errno, "binnie_open_in: error opening [%s] as sam", filename);
	}
    }
  else
    {
      error(0, 0, "binnie_open_in: filename [%s] does not end in .bam or .sam", filename);
      return fp;
    }
  
  blog(3, "binnie_open_in: opened fp->fn=[%s]", fp->fn);

  DLOG("binnie_open_in: returning fp=[%p] for filename=[%s]", fp, filename);
  return fp;
}


/*
 * binnie_close
 *
 * Closes the file pointed to by FP.
 *
 */
void binnie_close(samFile *fp)
{
  DLOG("binnie_close: fp=[%p]", fp);
  sam_close(fp);
  DLOG("binnie_close: returning");
}
