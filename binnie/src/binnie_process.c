/*
 * binnie_process.c - processes reads into output bins
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
#include <stdbool.h>
#include "error.h"
#include "progname.h"
#include "xalloc.h"
#include "xstrndup.h"
#include "gl_xlist.h"
#include "gl_avltreehash_list.h"
#include "hash-pjw.h"

/* internationalisation */
#include "gettext.h"

/* binnie includes */
#include "binnie.h"
#include "binnie_log.h"
#include "binnie_process.h"

/* htslib for sam/bam processing */
#include <htslib/sam.h>


/*
 * binnie_process
 * --------------
 * 
 * Processes each read from the original input BAM into output BAMs for 
 * unchanged, bridged, and remap categories depending on the status of 
 * the read in the bridge BAM. 
 *
 * INPUT: pointers to open samFile structures for original and bridge 
 *        BAM files (opened for input and pre-sorted on contig/position) 
 *        and for unchanged, bridged, and remap BAM files (opened for output).
 *
 * OUTPUT: none (side effect is that output files are written to).
 *
 *
 * Procedure in more detail:
 * 
 * 0. Headers are written to output files.
 * 
 * 1. Individual Reads processed into a buffer containing the read data and a result bin. 
 *    (see binnie_read_bin for more details)
 *
 * 2. Check if any of the read's mates are already in the buffer -- if there are, then 
 *    update all reads to increment their mate_count and change the result bin for both 
 *    reads such that if any read is set to Remap, they are all set to Remap or if there 
 *    is any disagreement between Unchanged and Bridged then they are all set to Remap.
 *    (see binnie_read_buffer for more details)
 *
 * 3. When the buffer is full, start writing to output bins, but first perform one final check: 
 *    if all a read's mates have not been added to buffer (or if number of mates is unknown), 
 *    then change its bin to Remap.
 *
 */
bool binnie_process(int buffer_size, int max_buffer_bases, samFile *original_in_fp, samFile *bridge_in_fp, samFile *unchanged_out_fp, samFile *bridged_out_fp, samFile *remap_out_fp)
{
  bool original_done;
  bool bridge_done;
  binnie_read_t *current_bridge_read;
  gl_list_t output_buffer;
  int32_t last_refid;
  int32_t last_pos;
  int32_t buffer_first_pos;
  int32_t buffer_last_pos;
  bam_hdr_t *original_header;
  bam_hdr_t *bridge_header;
  bam_hdr_t *unchanged_header;
  bam_hdr_t *bridged_header;
  bam_hdr_t *remap_header;
  uint32_t read_count;
  bool new_refid;
  int buffer_read_count;
  int buffer_read_count_max;

  DLOG("binnie_process()");

  last_refid = 0;
  last_pos = 0;
  buffer_first_pos = 0;
  buffer_last_pos = 0;
  original_done = false;
  bridge_done = false;
  new_refid = true;
  buffer_read_count = 0;
  buffer_read_count_max = 0;

  /* read BAM/SAM headers */
  blog(3, gettext("reading headers"));
  original_header = sam_hdr_read(original_in_fp);
  DLOG(gettext("binnie_process: original has %d targets"), original_header->n_targets);
  bridge_header = sam_hdr_read(bridge_in_fp);
  DLOG(gettext("binnie_process: bridge has %d targets"), bridge_header->n_targets);


  // TODO we should modify the output header to merge RG entries, etc (and add our PG entry?)
  unchanged_header = original_header;
  remap_header = original_header;
  bridged_header = bridge_header;
  

  /* write BAM/SAM headers to each of the output bins */
  blog(3, gettext("writing headers"));
  sam_hdr_write(unchanged_out_fp, unchanged_header);
  sam_hdr_write(bridged_out_fp, bridged_header);
  sam_hdr_write(remap_out_fp, remap_header);


  /* initialize read buffer */
  DLOG("binnie_process: intializing output buffer");
  output_buffer = gl_list_create_empty(GL_AVLTREEHASH_LIST, 
                                       bbr_equals, 
                                       bbr_hashcode, 
                                       bbr_dispose, 
                                       false);
  
  
  /* Allocate/initialize BAM/SAM read buffers. */
  DLOG(gettext("binnie_process: initializing current bridge read"));
  current_bridge_read = br_init();
  
  /* Read original and bridge in synchrony. */
  DLOG(gettext("binnie_process: beginning read processing loop"));
  read_count = 0;
  do 
    {
      binnie_read_t *original_read;
      binnie_binned_read_t *bbr;
      int ret; 
      bool success;
      int32_t refid;
      int32_t pos;
      int32_t reads_output;

      /* increment read_count */
      read_count++;
      
      DLOG(gettext("binnie_process: processing read [%d]"), read_count);

      DLOG(gettext("binnie_process: initializing original_read"));
      original_read = br_init();


      /* Read from original. */ 
      DLOG(gettext("binnie_process: reading original_read"));
      ret = sam_read1(original_in_fp, original_header, original_read->bam_read);
      if ( ret >= 0 )
        {
          /* We have an original read, compare to current_bridge_read. */
          original_read->bam_read_present = true;

	  /* get refid and pos from original */
	  refid = br_get_refid(original_read);
	  pos = br_get_pos(original_read);
        }
      else if ( ret == -1 )
        {
          /* Reached normal end of file. */
          blog(3, gettext("reached end of original file"));
          original_read->bam_read_present = false;
          original_done = true;
        }
      else 
        {
          /* An error occured reading the original file. */
          errx(BINNIE_EXIT_ERR_READ_ORIG, gettext("binnie_process: error reading from original input file"));
        }

      if (new_refid)
	{
	  if (refid >= 0)
	    blog(1, gettext("processing original reads mapped to reference [%s]"), original_header->target_name[refid]);
	  else 
	    blog(1, gettext("processing original unmapped reads"));
	    
	}
      
      DLOG(gettext("binnie_process: have original read at refid=[%d] pos=[%d]"), refid, pos);

      if ( !(current_bridge_read->bam_read_present) && !bridge_done )
        {
          DLOG(gettext("binnie_process: reading current_bridge_read"));

          /* Read from bridge. */
          ret = sam_read1(bridge_in_fp, bridge_header, current_bridge_read->bam_read);
          if ( ret >= 0 ) 
            {
              /* Have a bridge-mapped read. */
              current_bridge_read->bam_read_present = true;
            }
          else if ( ret == -1 )
            {
              /* Reached normal end of file. */
              blog(3, gettext("reached end of bridge file"));
              current_bridge_read->bam_read_present = false;
              bridge_done = true;
            }
          else 
            {
              /* An error occurred reading the bridge file. */
              errx(BINNIE_EXIT_ERR_READ_BRIDGE, gettext("binnie_process: error reading from bridge input file"));
            }
    
          DLOG(gettext("binnie_process: have bridge-mapped read at refid=[%d] pos=[%d]"), br_get_refid(current_bridge_read), br_get_pos(current_bridge_read));
      
        }

      if (!original_done)
	{
	  DLOG(gettext("binnie_process: checking if original_read equals current_bridge_read. original_read->bam_read_present=[%d] current_bridge_read->bam_read_present=[%d]"), original_read->bam_read_present, current_bridge_read->bam_read_present);
	  if ( current_bridge_read->bam_read_present && br_equals(original_read, current_bridge_read) )
	    {
	      /* have a match for the original read, bin the reads */
	      DLOG(gettext("binnie_process: original_read matches current_bridge_read"));
	      bbr = binnie_read_bin(original_read, current_bridge_read);
	      
	      DLOG(gettext("binnie_process: initializing current_bridge_read"));
	      current_bridge_read = br_init();
	    } 
	  else 
	    {
	      /* original_read doesn't match bridge_read, output the original */
	      DLOG(gettext("binnie_process: original read is not a match for current_bridge_read"));
	      bbr = binnie_read_bin(original_read, NULL);
	    }
	  
	  /* if bbr is NULL, it means that binnie_read_bin wants to discard this read */
	  if (bbr == NULL)
	    {
	      DLOG(gettext("binnie_process: have NULL bbr (binnie_read_bin wants to discard this read), skipping to next iteration of processing loop"));
	      continue;
	    }
	  
	  /* verify refid has not decreased */
	  DLOG(gettext("binnie_process: checking that refid has not decreased.  refid=[%d] last_refid=[%d]"), refid, last_refid);
	  if ( (refid < last_refid) && (refid != -1) && (last_refid != -1) )
	    {
	      errx(BINNIE_EXIT_ERR_BAM_UNSORTED, gettext("binnie_process: sort error -- current refid [%d] was less than the last one [%d]"), refid, last_refid);
	    }
	  
	  /* 
	   * verify refid has not switched from unmapped to mapped 
	   * (all unmapped reads should go at the end) 
	   */
	  DLOG(gettext("binnie_process: checking that refid has not switched from unmapped back to mapped.  refid=[%d] last_refid=[%d]"), refid, last_refid);
	  if ( (last_refid == -1) && (refid != -1) )
	    {
	      errx(BINNIE_EXIT_ERR_BAM_UNSORTED, gettext("binnie_process: sort error -- current refid [%d] was set but last refid was unmapped"), refid);
	    }
	  
	  /* 
	   * if refid has changed, set new_refid flag and 
	   * reset last_pos to 0 (or to -1 if we've changed to unmapped)
	   */
	  DLOG(gettext("binnie_process: checking if refid has changed.  refid=[%d] last_refid=[%d]"), refid, last_refid);
	  if (refid != last_refid)
	    {
	      blog(2, gettext("reference id now [%d]"), refid);
	      DLOG(gettext("have new refid.  refid=[%d] last_refid=[%d]"), refid, last_refid);
	      new_refid = true;
	      if (refid == -1)
		{
		  DLOG(gettext("binnie_process: refid now unmapped.  setting last_pos to -1"));
		  last_pos = -1;
		}
	      else
		{
		  DLOG(gettext("binnie_process: refid still mapped.  resetting last_pos to 0"));
		  last_pos = 0;
		}
	    }
	  else 
	    {
	      new_refid = false;
	    }
	  
	  
	  /* verify that pos has not decreased */
	  DLOG(gettext("binnie_process: checking that pos has not decreased.  pos=[%d] last_pos=[%d]"), pos, last_pos);
	  if ( (pos < last_pos) && (pos != -1) && (last_pos != -1) )
	    {
	      errx(BINNIE_EXIT_ERR_BAM_UNSORTED, gettext("binnie_process: sort error -- current pos [%d] was less than the last one [%d]"), pos, last_pos);
	    }
	  
	  /* 
	   * verify pos has not switched from unmapped to mapped 
	   * (all unmapped reads should go at the end) 
	   */
	  DLOG(gettext("binnie_process: checking that pos has not switched from unmapped to mapped.  pos=[%d] last_pos=[%d]"), pos, last_pos);
	  if ( (last_pos == -1) && (pos != -1) )
	    {
	      errx(BINNIE_EXIT_ERR_BAM_UNSORTED, gettext("binnie_process: sort error -- current pos [%d] was set but last pos was unmapped"), pos);
	    }
	  
	  
	  /* sort order now confirmed - set last_refid and last_pos for next iteration */
	  DLOG(gettext("binnie_process: sort order confirmed. updating last_refid and last_pos from last_refid=[%d] last_pos=[%d]"), last_refid, last_pos);
	  last_refid = refid;
	  last_pos = pos;
	  DLOG(gettext("binnie_process: last_refid and last_pos updated to last_refid=[%d] last_pos=[%d]"), last_refid, last_pos);
	  
	  
	  /* add read to buffer, update mates in buffer, etc */
	  DLOG(gettext("binnie_process: calling binnie_read_buffer"));
	  binnie_read_buffer(bbr, output_buffer);
	  
	  /* update buffer_last_pos */
	  DLOG(gettext("binnie_process: updating buffer_last_pos from buffer_last_pos=[%d]"), buffer_last_pos);
	  buffer_last_pos = bbr->original_pos;
	  DLOG(gettext("binnie_process: updated buffer_last_pos to buffer_last_pos=[%d]"), buffer_last_pos);

	  /* if this is the first read in the buffer, also update buffer_first_pos */
	  if (buffer_read_count == 0)
	    {
	      DLOG(gettext("binnie_process: buffer now contains a single read. updating buffer_first_post from buffer_first_pos=[%d]"), buffer_first_pos);
	      buffer_first_pos = buffer_last_pos;
	      DLOG(gettext("binnie_process: updated buffer_first_post to buffer_first_pos=[%d]"), buffer_first_pos);
	    }
      
	}


      /* update buffer_read_count and set buffer_read_count_max */
      buffer_read_count = gl_list_size(output_buffer);
      if ((buffer_read_count > buffer_read_count_max) && (refid >= 0))
	{
	  buffer_read_count_max = buffer_read_count;
	}


      /* 
       * if original is done or if refid has changed, it's time to flush the buffer completely.  
       * if buffer is "full" (size larger than buffer_size or position range greater than 
       * max_buffer_bases), then we need to flush until the size and/or position are no longer full.  
       */
      DLOG(gettext("binnie_process: beginning buffer output loop.  original_done=[%d] buffer_read_count=[%d] new_refid=[%d] buffer_size=[%d] buffer_last_pos=[%d] buffer_first_pos=[%d] (buffer_last_pos-buffer_first_pos)=[%d] max_buffer_bases=[%d]"), original_done, buffer_read_count, new_refid, buffer_size, buffer_last_pos, buffer_first_pos, (buffer_last_pos-buffer_first_pos), max_buffer_bases);
      reads_output = 0;
      while ( (original_done && (buffer_read_count > 0))
              || (new_refid && (buffer_read_count > 0))
              || ((buffer_size > 0) && (buffer_read_count >= buffer_size))
              || ((max_buffer_bases > 0) && ((buffer_last_pos - buffer_first_pos) >= max_buffer_bases))
              )
        {

          /* get a read from buffer */
	  DLOG("binnie_process: calling gl_list_get_at 0");
          bbr = (binnie_binned_read_t *) gl_list_get_at(output_buffer, 0);
          

          /* write the read to the designated output bin */
          switch (bbr->bin) {
          case BINNIE_UNCHANGED:
	    DLOG(gettext("binnie_process: writing to unchanged output bin."));
            ret = sam_write1(unchanged_out_fp, unchanged_header, bbr->br->bam_read);
	    reads_output++;
            if (ret <= 0)
              {
                errx(BINNIE_EXIT_ERR_WRITE, gettext("binnie_process: could not write to unchanged out file"));
              }
            break;
          case BINNIE_BRIDGED:
	    DLOG(gettext("binnie_process: writing to bridged output bin."));
            ret = sam_write1(bridged_out_fp, bridged_header, bbr->br->bam_read);
	    reads_output++;
            if (ret <= 0)
              {
                errx(BINNIE_EXIT_ERR_WRITE, gettext("binnie_process: could not write to bridged out file"));
              }
            break;
          case BINNIE_REMAP:
	    DLOG(gettext("binnie_process: writing to remap output bin."));
            ret = sam_write1(remap_out_fp, remap_header, bbr->br->bam_read);
	    reads_output++;
            if (ret <= 0)
              {
                errx(BINNIE_EXIT_ERR_WRITE, gettext("binnie_process: could not write to remap out file"));
              }
            break;
          default:
            errx(BINNIE_EXIT_ERR_INVALID_BIN, gettext("binnie_process: invalid bin [%d]"), bbr->bin);
          }
          
          /* remove the read from the buffer (this will call bbr_dispose for us) */
	  DLOG("binnie_process: calling gl_list_remove_at 0");
          success = gl_list_remove_at(output_buffer, 0);
          if (!success)
            {
              errx(BINNIE_EXIT_ERR_BUFFER_REMOVE, gettext("binnie_process: failure removing node from buffer"));
            }

	  /* update read count for next iteration */
          buffer_read_count = gl_list_size(output_buffer);
          if (buffer_read_count > 0)
            {
              /* still have reads in buffer, update buffer_first_pos to position of new first read in buffer */
	      DLOG("binnie_process: calling gl_list_get_at 0");
              bbr = (binnie_binned_read_t *) gl_list_get_at(output_buffer, 0);
              buffer_first_pos = bbr->original_pos;
            }
          else  /* buffer_read_count <= 0 */
            {
              /* no reads in buffer, reset buffer_first_pos and buffer_last_pos to 0 */
              buffer_first_pos = 0;
              buffer_last_pos = 0;
            }

          DLOG(gettext("binnie_process: end of buffer output iteration.  original_done=[%d] buffer_read_count=[%d] new_refid=[%d] buffer_size=[%d] buffer_last_pos=[%d] buffer_first_pos=[%d] (buffer_last_pos-buffer_first_pos)=[%d] max_buffer_bases=[%d]"), original_done, buffer_read_count, new_refid, buffer_size, buffer_last_pos, buffer_first_pos, (buffer_last_pos-buffer_first_pos), max_buffer_bases);

        } /* while original_done ... || new_refid ... || ... > buffer_size || ... > max_buffer_bases */
      DLOG(gettext("binnie_process: finished buffer output loop after outputting [%d] reads."), reads_output);

      DLOG(gettext("binnie_process: done processing read [%d]"), read_count);
    } while ( !original_done );

  DLOG(gettext("binnie_process: finished read processing loop.  original_done=[%d]"), original_done);
  

  DLOG(gettext("binnie_process: checking if bridge is done"));
  if ( !bridge_done || (current_bridge_read->bam_read_present) )
    {
      /* original finished before bridge */
      br_dispose(current_bridge_read);

      /* FATAL ERROR */
      errx(BINNIE_EXIT_ERR_ORIG_TRUNCATED, gettext("binnie_process: original finished but bridge read(s) remain"));
    }


  /* check if buffer still has reads remaining */
  DLOG(gettext("binnie_process: checking if output buffer is empty"));
  buffer_read_count = gl_list_size(output_buffer);
  if (buffer_read_count > 0) 
    {
      /* FATAL ERROR */
      errx(BINNIE_EXIT_ERR_BUFFER_NOT_EMPTY, gettext("output_buffer was not empty at end of binnie_process (%d reads remained)."), buffer_read_count);
    }
  

  blog (3, gettext("freeing the read buffer"));
  gl_list_free (output_buffer);

  blog(1, gettext("finished processing reads. had a maximum of %d reads in buffer (not counting unmapped reads)."), buffer_read_count_max);
  if (buffer_read_count_max >= buffer_size && max_buffer_bases > 0)
    {
      blog(0, gettext("WARNING: buffer was limited by size (%d reads) rather than bases"), buffer_read_count_max);
    }

  DLOG("binnie_process: returning true");
  return true;
} /* binnie_process */


/*
 * binnie_read_bin
 * ---------------
 * 
 * Examines the read and decides in which bin it belongs.
 *
 * INPUT: pointer to binnie_read_t structures for original and bridge-mapped reads, or 
 *        if there is no bridge-mapped read matching this original, bridge_read can be NULL
 * OUTPUT: pointer to binnie_binned_read_t or NULL pointer if this read has no bin
 *
 * SIDE EFFECT: disposes of all the unbinned reads
 *
 * Detail of binning:
 *    --------------------------------------
 *    Original  Bridge    Bin
 *    --------------------------------------
 *    Unmapped  Unmapped  Unchanged
 *    Unmapped  MAPQ >= 0 Bridged
 *    MAPQ == 0 Unmapped  Unchanged
 *    MAPQ == 0 MAPQ == 0 Unchanged*
 *    MAPQ == 0 MAPQ > 0  Remap
 *    MAPQ > 0  Unmapped  Unchanged
 *    MAPQ > 0  MAPQ == 0 Remap
 *    MAPQ > 0  MAPQ > 0  Remap
 *    Deleted   (any)     Remap
 *    Secondary (any)     (dispose)
 *    --------------------------------------
 */
binnie_binned_read_t *binnie_read_bin(binnie_read_t *original_read, binnie_read_t *bridge_read)
{
  binnie_binned_read_t *bbr;
  int32_t original_mapq;
  int32_t bridge_mapq;

  DLOG("binnie_read_bin()");

  if (original_read == NULL) 
    {
      /* original_read must never be NULL: FATAL ERROR */
      errx(BINNIE_EXIT_ERR_NULL, gettext("binnie_read_bin: original_read must not be NULL"));
    }
  
  /* check for deletion of original coodinates */
  if (false) // TODO check for deletion of original coordinates
    {
      /* 
       *    --------------------------------------
       *    Original  Bridge    Bin
       *    --------------------------------------
       *    Deleted   (any)     Remap
       */
      bbr = bbr_init(original_read);
      bbr->bin = BINNIE_REMAP;      
    }
  
  /* check if the original read is a secondary alignment */
  if ( 
      !((original_read->bam_read)->core.flag & BAM_FUNMAP)
      && ((original_read->bam_read)->core.flag & BAM_FSECONDARY)
       )
    {
      /*
       *    --------------------------------------
       *    Original  Bridge    Bin
       *    --------------------------------------
       *    Secondary (any)     (dispose)
       *    --------------------------------------
       */
      bbr = NULL;
      br_dispose(original_read);
      if (bridge_read != NULL)
        {
          br_dispose(bridge_read);
        }
    }
  
  /* check if bridge read is not present */
  if (bridge_read == NULL)
    {
      /* no bridge read (it must be unmapped) */
      /*
       *    --------------------------------------
       *    Original  Bridge    Bin
       *    --------------------------------------
       *    Unmapped  Unmapped  Unchanged
       *    MAPQ == 0 Unmapped  Unchanged
       *    MAPQ > 0  Unmapped  Unchanged
       */
      bbr = bbr_init(original_read);
      bbr->bin = BINNIE_UNCHANGED;
    }
  else
    {
      /* have both an original and a bridge mapped read */
      original_mapq = br_get_mapq(original_read);
      bridge_mapq = br_get_mapq(bridge_read);
      if (original_mapq < 0)
        {
          if (bridge_mapq < 0)
            {
              /* 
               *    --------------------------------------
               *    Original  Bridge    Bin
               *    --------------------------------------
               *    Unmapped  Unmapped  Unchanged
               */
              bbr = bbr_init(original_read);
              bbr->bin = BINNIE_UNCHANGED;
              br_dispose(bridge_read);
            }
          else /* bridge_mapq >= 0 */
            {
              /* 
               *    --------------------------------------
               *    Original  Bridge    Bin
               *    --------------------------------------
               *    Unmapped  MAPQ >= 0 Bridged
               */	      
	      fixup_bridge_from_original(bridge_read, original_read);
              bbr = bbr_init(bridge_read);
              bbr->bin = BINNIE_BRIDGED;
	      bbr->original_refid = br_get_refid(original_read);
	      bbr->original_pos = br_get_pos(original_read);
              br_dispose(original_read);
            }
        }
      else if (original_mapq == 0) 
        {
          if (bridge_mapq < 0)
            {
              /* 
               *    --------------------------------------
               *    Original  Bridge    Bin
               *    --------------------------------------
               *    MAPQ == 0 Unmapped  Unchanged
               */
              bbr = bbr_init(original_read);
              bbr->bin = BINNIE_UNCHANGED;
              br_dispose(bridge_read);
            }
          else if (bridge_mapq == 0)
            {
              /* 
               *    --------------------------------------
               *    Original  Bridge    Bin
               *    --------------------------------------
               *    MAPQ == 0 MAPQ == 0 Unchanged*
               */
              bbr = bbr_init(original_read);
              bbr->bin = BINNIE_UNCHANGED;
              br_dispose(bridge_read);
            }
          else /* bridge_mapq > 0 */
            {
              /* 
               *    --------------------------------------
               *    Original  Bridge    Bin
               *    --------------------------------------
               *    MAPQ == 0 MAPQ > 0  Remap
               */
              bbr = bbr_init(original_read);
              bbr->bin = BINNIE_REMAP;
              br_dispose(bridge_read);
            }
        }
      else /* original_mapq > 0 */ 
        {
          if (bridge_mapq < 0)
            {
              /* 
               *    --------------------------------------
               *    Original  Bridge    Bin
               *    --------------------------------------
               *    MAPQ > 0  Unmapped  Unchanged
               */
              bbr = bbr_init(original_read);
              bbr->bin = BINNIE_UNCHANGED;
              br_dispose(bridge_read);
            }
          else if(bridge_mapq == 0)
            {
              /* 
               *    --------------------------------------
               *    Original  Bridge    Bin
               *    --------------------------------------
               *    MAPQ > 0  MAPQ == 0 Remap
               */
              bbr = bbr_init(original_read);
              bbr->bin = BINNIE_REMAP;
              br_dispose(bridge_read);
            }
          else /* bridge_mapq > 0 */
            {
              /* 
               *    --------------------------------------
               *    Original  Bridge    Bin
               *    --------------------------------------
               *    MAPQ > 0  MAPQ > 0  Remap
               */
              bbr = bbr_init(original_read);
              bbr->bin = BINNIE_REMAP;
              br_dispose(bridge_read);
            }
        }
          
    }

  DLOG("binnie_read_bin: returning bbr in bin [%d]", bbr->bin);
  return bbr;
} /* binnie_read_bin */


/*
 * fixup_bridge_from_original
 * --------------------------
 * 
 * Copies FPAIRED, FREAD1, FREAD2 flags and FI and RG tags from original to bridge. 
 *
 * INPUT: pointers to binnie_read_t structure for bridge and original read
 *
 * SIDE EFFECT: updates bridge_read with data from original_read
 *
 */
void fixup_bridge_from_original (binnie_read_t *bridge_read, binnie_read_t *original_read)
{
  uint8_t *original_tag;
  uint8_t *bridge_tag;

  DLOG(gettext("fixup_bridge_from_original: checking if we need to override bridge_read FPAIRED flag. original_read flag=[%d] bridge_read flag=[%d]"), original_read->bam_read->core.flag, bridge_read->bam_read->core.flag);
  if ( (original_read->bam_read->core.flag & BAM_FPAIRED) 
       && !(bridge_read->bam_read->core.flag & BAM_FPAIRED) )
    {
      bridge_read->bam_read->core.flag = bridge_read->bam_read->core.flag | BAM_FPAIRED;
      blog(4, gettext("set bridge read flag FPAIRED. flag=[%d]"), bridge_read->bam_read->core.flag);
    }

  DLOG(gettext("fixup_bridge_from_original: checking if we need to override bridge_read FREAD1 flag. original_read flag=[%d] bridge_read flag=[%d]"), original_read->bam_read->core.flag, bridge_read->bam_read->core.flag);
  if (original_read->bam_read->core.flag & BAM_FREAD1)
    {
      bridge_read->bam_read->core.flag = bridge_read->bam_read->core.flag | BAM_FREAD1;
      blog(4, gettext("set bridge read flag FREAD1. flag=[%d]"), bridge_read->bam_read->core.flag);
    }	      

  DLOG(gettext("fixup_bridge_from_original: checking if we need to override bridge_read FREAD2 flag. original_read flag=[%d] bridge_read flag=[%d]"), original_read->bam_read->core.flag, bridge_read->bam_read->core.flag);
  if (original_read->bam_read->core.flag & BAM_FREAD2)
    {
      bridge_read->bam_read->core.flag = bridge_read->bam_read->core.flag | BAM_FREAD2;
      blog(4, gettext("set bridge read flag FREAD2. flag=[%d]"), bridge_read->bam_read->core.flag);
    }	 

  DLOG(gettext("fixup_bridge_from_original: checking if we need to override bridge_read FI tag."));
  original_tag = bam_aux_get(original_read->bam_read, "FI");
  if (original_tag > 0)
    {
      /* we have FI tag in original, check for FI tag in bridge */
      bridge_tag = bam_aux_get(bridge_read->bam_read, "FI");
      if (bridge_tag > 0)
	{
	  /* have tag in bridge, delete it */
	  bam_aux_del(bridge_read->bam_read, bridge_tag);
	}
      
      /* add original FI tag to bridge */
      bam_aux_append(bridge_read->bam_read, "FI", 'i', bam_aux_type2size('i'), original_tag);

      blog(4, gettext("set bridge read tag FI=[%d]"), bam_aux2i(bam_aux_get(bridge_read->bam_read, "FI")));
    }

  /* only fixup RG if we are ignoring RG in matching, otherwise we already know they are the same */
  if (ignore_rg)
    {
      DLOG(gettext("fixup_bridge_from_original: checking if we need to override bridge_read RG tag."));
      original_tag = bam_aux_get(original_read->bam_read, "RG");
      if (original_tag)
	{
	  int rg_len;
	  char *rg;
	  
	  /* we have RG tag in original, check for RG tag in bridge */
	  bridge_tag = bam_aux_get(bridge_read->bam_read, "RG");
	  if (bridge_tag)
	    {
	      /* have tag in bridge, delete it */
	      bam_aux_del(bridge_read->bam_read, bridge_tag);
	    }
	  
	  rg = bam_aux2Z(original_tag);
	  rg_len = (int) strlen(rg);
	  DLOG(gettext("have original RG=[%s] len=[%d]"), rg, rg_len);
	  
	  /* add original RG tag to bridge */
	  bam_aux_append(bridge_read->bam_read, "RG", 'Z', rg_len+1, (uint8_t*)rg);
	  
	  blog(4, gettext("set bridge read tag RG=[%s]"), bam_aux2Z(bam_aux_get(bridge_read->bam_read, "RG")));
	}
    }

} /* fixup_bridge_from_original */
	      

/*
 * binnie_read_buffer
 * ------------------
 * 
 * Takes a read to be added and finds all mates (if any) in the buffer that match the 
 * given read and increment the mate_count for each as well as for the new read. 
 * Checks the bin of all buffered reads and this read and if there is any disagreement, 
 * set all bins to remap.
 *
 * INPUT: pointer to binnie_binned_read_t structure and gl_list_t structure for buffer
 *
 * SIDE EFFECT: updates bin and mate_count in the buffer
 *
 */
void binnie_read_buffer (binnie_binned_read_t *bbr, gl_list_t output_buffer)
{
  bool all_bins_agree;
  gl_list_node_t node;
  binnie_binned_read_t *bbri;

  DLOG("binnie_read_buffer()");

  /* search for a node like this one */
  DLOG("binnie_read_buffer: calling gl_list_search");
  node = gl_list_search(output_buffer, bbr);
  if (node == NULL) 
    {
      /* this will be the first node for this template in the buffer, just add it */
      DLOG("binnie_read_buffer: node not found, calling gl_list_add_last");
      gl_list_add_last(output_buffer, bbr);
    }
  else /* node != NULL */
    {
      DLOG("binnie_read_buffer: matching node found");

      /* at least one read for this template is already queued */
      if (bbr->expected_mate_count == 0)
        {
          /* FATAL ERROR: we weren't expecting any mates for this read */
          errx(BINNIE_EXIT_ERR_UNEXPECTED_MATES, gettext("binnie_read_buffer: mate found for read but expected_mate_count was 0.  rg=[%s] qname=[%s]"), br_get_read_group(bbr->br), br_get_qname(bbr->br));
        }

      /* wind to beginning of linked list of buffered reads */
      DLOG("binnie_read_buffer: calling gl_list_node_value");
      bbri = (binnie_binned_read_t *) gl_list_node_value(output_buffer, node);
      while (bbri->prev_mate != NULL)
        {
          bbri = bbri->prev_mate;
        }
      

      /* should now be at first read in buffered linked list */
      if (bbri->prev_mate != NULL)
        {
          /* ERROR: this should have been NULL at this point */
          errx(BINNIE_EXIT_ERR_NOT_NULL, gettext("binnie_read_buffer: should have been at beginning of linked list, but bbri->prev_mate was not NULL (impossible!)"));
        }
      

      /* sweep through linked list to end, processing each read as we go */
      all_bins_agree = true;
      do
        {
          /* increment mate_count to account for the new read for this template */
          bbri->mate_count++;
          
          /* check if the buffered read has an unknown expected mate count */
          if (bbri->expected_mate_count < 0)
            {
              /* check if we know the expected_mate_count for this template */
              if (bbr->expected_mate_count >= 0)
                {
                  /* we do, set the buffered version to what we know */
                  bbri->expected_mate_count = bbr->expected_mate_count;
                }
            }

          /* check that our new read's bin agrees with this one */
          if (bbr->bin != bbri->bin)
            {
              all_bins_agree = false;
            }
          
        } while (bbri->next_mate != NULL);
      

      /* should now be at last read in buffered linked list */ 
      if (bbri->next_mate != NULL)
        {
          /* ERROR: this should have been NULL at this point */
          errx(BINNIE_EXIT_ERR_NOT_NULL, gettext("binnie_read_buffer: should have been at end of linked list, but bbri->next_mate was not NULL (impossible!)"));
        }
      

      /* attach the new read to the end of the linked list and add it to the buffer */
      bbri->next_mate = bbr;
      bbr->prev_mate = bbri;
      DLOG("binnie_read_buffer: calling gl_list_add_last");
      gl_list_add_last(output_buffer, bbr);
      

      /* if all bins did not agree, we need to set them all to remap */
      if (!all_bins_agree)
        {
          /* start at the end with the new read */
          bbri = bbr;

          /* sweep backwards through linked list, resetting all bins to remap */
          do
            {
              bbri->bin = BINNIE_REMAP;
              bbri = bbri->prev_mate;
            } while (bbri != NULL);
        }
      
    } /* node != NULL */

  DLOG("binnie_read_buffer: returning void");
} /* binnie_read_buffer */


/*
 * br_get_refid
 * -----------------------------
 *
 * Gets the target reference id for the binnie_read, or -1 if the read is unmapped.
 *
 * INPUT: pointer to binnie_read_t
 * OUTPUT: int32_t indicating the index of the target in the header
 */
int32_t br_get_refid (const binnie_read_t *br)
{
  int32_t refid;

  DLOG("br_get_refid()");

  refid = -1;

  if ( !((br->bam_read)->core.flag & BAM_FUNMAP) )
    {
      refid = (br->bam_read)->core.tid;
    }
  
  DLOG("br_get_refid: returning refid=[%d]", refid);
  return refid;
}


/*
 * br_get_pos
 * -----------------------------
 *
 * Gets the 0-based position on the target for the binnie_read, or -1 if the read is unmapped.
 *
 * INPUT: pointer to binnie_read_t
 * OUTPUT: int32_t indicating the index of the target in the header
 */
int32_t br_get_pos (const binnie_read_t *br)
{
  int32_t pos;

  DLOG("br_get_pos()");

  pos = -1;
  
  if ( !((br->bam_read)->core.flag & BAM_FUNMAP) )
    {
      pos = (br->bam_read)->core.pos;
    }
  
  DLOG("br_get_pos: returning pos=[%d]", pos);
  return pos;
}


/*
 * br_get_mapq
 * -----------------------------
 *
 * Gets the mapping quality for a binnie_read, or -1 if the read is 
 * unmapped or the mapping quality is not available.
 *
 * INPUT: pointer to binnie_read_t
 * OUTPUT: int32_t indicating mapping quality
 */
int32_t br_get_mapq (const binnie_read_t *br)
{
  int32_t mapq;
  
  DLOG("br_get_mapq()");

  /* set default mapq to -1 */
  mapq = -1;
  if ( !((br->bam_read)->core.flag & BAM_FUNMAP) )
    {
      /* this is a mapped read, get mapping quality */
      mapq = (br->bam_read)->core.qual;
      
      if ( mapq == 255 )
        {
          /* mapq of 255 represents unavailable */
          mapq = -1;
        }
    }
  
  DLOG("br_get_mapq: returning mapq=[%d]", mapq);
  return mapq;
}
 

/*
 * br_get_segment_index
 * -----------------------------
 *
 * Gets the segment index for a binnie_read, from the FI tag if present, 
 * otherwise based on flags (1 for first, 2 for last) or -1 if it can't be determined.
 *
 * INPUT: pointer to binnie_read_t
 * OUTPUT: int32_t indicating segment index 
 */
int32_t br_get_segment_index (const binnie_read_t *br)
{
  int32_t segment_index;
  uint8_t *fi;

  DLOG("br_get_segment_index()");

  fi = bam_aux_get(br->bam_read, "FI");
  if (fi > 0)
    {
      /* we have FI tag, use it */
      segment_index = bam_aux2i(fi);
      DLOG(gettext("br_get_segment_index: tag FI is present.  fi=[%d]"), segment_index);
    }
  else 
    {
      if ((br->bam_read)->core.flag & BAM_FREAD1)
        {
          if ((br->bam_read)->core.flag & BAM_FREAD2)
            {
              /* 
               * FREAD1 && FREAD2 
               *
               * the segment is part of a linear template, but it is neither 
               * the first nor the last segment.
               * this indicates a bad BAM -- we should have had FI in this case
               * FATAL ERROR - bail out
               */
              errx(BINNIE_EXIT_ERR_SEGMENT_INDEX, 
                   gettext("br_get_segment_index: FREAD1 and FREAD2 were set, but FI flag not found for read rg=[%s] qname=[%s]"), 
                   br_get_read_group(br), br_get_qname(br));
            }
          else 
            {
              /*
               * FREAD1 && !FREAD2 
               */
              DLOG(gettext("br_get_segment_index: flag FREAD1 but not FREAD2 is set. this is the first segment in a pair."));
              segment_index = 1;
            }
        }
      else 
        {
          if ((br->bam_read)->core.flag&BAM_FREAD2)
            {
              /* 
               * !FREAD1 && FREAD2 
               */
              DLOG(gettext("br_get_segment_index: flag FREAD2 but not FREAD1 is set. this is the second segment in a pair."));
              segment_index = 2;
            }
          else 
            {
              /* 
               * !FREAD1 && !FREAD2 
               *
               * This may happen for a non-linear template or the index is lost in data processing.
               */
              blog(1, "WARNING: non-linear template or index lost in data processing");
              segment_index = -1;
            }
        }
    }
  
  DLOG("br_get_segment_index: returning segment_index=[%d]", segment_index);

  return segment_index;
}


/*
 * br_get_num_segments
 * --------------------
 *
 * Gets the number of segments for this read's template/qname. 
 * Uses TC tag if present, otherwise based on flags or -1 if it cannot be 
 * determined from the data at this read. 
 *
 * INPUT: pointer to binnie_read_t
 * OUTPUT: int32_t indicating number of segments
 */
int32_t br_get_num_segments (const binnie_read_t *br)
{
  int32_t num_segments;

  DLOG("br_get_num_segments()");

  num_segments = -1;


  uint8_t *tc;
  DLOG("br_get_num_segments: getting tag TC");
  tc = bam_aux_get(br->bam_read, "TC");
  if (tc > 0)
    {
      /* we have TC tag, use it */
      num_segments = bam_aux2i(tc);
      DLOG(gettext("br_get_num_segments: have tag tc=[%d]"), num_segments);
    }
  else /* tc <= 0 */
    {
      /* 
       * If the template has more than 2 segments, the TC tag should be present.
       * since we don't have TC tag, we need to check flags to see whether we have 
       * a single or multiple segments.
       */
      DLOG(gettext("br_get_num_segments: don't have tag TC"));
      if ( !((br->bam_read)->core.flag & BAM_FPAIRED) )
        {
          /*
           * template only has single read
           */
          DLOG(gettext("br_get_num_segments: flag indicates this is not a paired read, only a single segment"));
          num_segments = 1;
        }
      else /* multi-read template */
        {
          /*
           * template having multiple segments in sequencing 
           * check other flags for consistency with this
           */
          DLOG(gettext("br_get_num_segments: flag indicates this read is part of a multiple segment template"));
          if ((br->bam_read)->core.flag & BAM_FREAD1)
            {
              if ((br->bam_read)->core.flag & BAM_FREAD2)
                {
                  /* 
                   * FREAD1 && FREAD2 
                   *
                   * the segment is part of a linear template, but it is neither the first nor the last segment. 
                   */
                  DLOG(gettext("br_get_num_segments: flags FREAD1 and FREAD2 both set.  this read is part of a linear template but is neither the first nor the last segment."));
		  blog(1, gettext("WARNING: unknown number of segments for read rg=[%s] qname=[%s] which is neither the first nor last segment but has no FI tag"), br_get_read_group(br), br_get_qname(br));
                  num_segments = -1;
                }
              else 
                {
                  /* 
                   * FREAD1 && !FREAD2 
                   *
                   * If the template has more than 2 segments, the TC tag should be present. 
                   */
                  DLOG(gettext("br_get_num_segments: flag FREAD1 but not FREAD2 are set, so it must be a paired read."));
                  num_segments = 2;
                }
            }
          else
            {
              if ((br->bam_read)->core.flag & BAM_FREAD2)
                {
                  /*
                   * !FREAD1 && FREAD2
                   * the segment is part of a linear template, but it is neither the first nor the last segment.
                   */
                  DLOG(gettext("br_get_num_segments: flag FREAD2 but not FREAD1 are set, so it must be a paired read."));
                  num_segments = 2;
                }
              else 
                {
                  /*
                   * !FREAD1 && !FREAD2
                   */
                  DLOG(gettext("br_get_num_segments: flags FREAD1 and FREAD2 are both unset. number of segments unknown."));
		  blog(1, gettext("WARNING: unknown number of segments for read rg=[%s] qname=[%s]"), br_get_read_group(br), br_get_qname(br));
                  num_segments = -1;
                }              
            }
        } /* multi-read template */
    } /* tc <= 0 */

  DLOG("br_get_num_segments: rg=[%s] qname=[%s] with flag [%d] has num_segments=[%d]", br_get_read_group(br), br_get_qname(br), (br->bam_read)->core.flag, num_segments);
  return num_segments;
}


/*
 * br_get_read_group
 * -----------------------------
 *
 * Gets the read group ID for a read, from the RG tag if present
 *
 * INPUT: pointer to binnie_read_t
 * OUTPUT: string (char *) pointing to read group id or empty string if no RG tag is present
 */
char *br_get_read_group (const binnie_read_t *br)
{
  char *read_group;
  uint8_t *rg;

  DLOG("br_get_read_group()");

  read_group = "";

  rg = bam_aux_get(br->bam_read, "RG");

  if (rg > 0)
    {
      /* we have RG tag */
      read_group = bam_aux2Z(rg);
      DLOG("br_get_read_group: tag RG is present.  rg=[%s]", read_group);
    }
  else 
    {
      DLOG("br_get_read_group: tag RG is not present.");
    }

  DLOG("br_get_read_group: returning read_group=[%s]", read_group);
  return read_group;
}


/*
 * br_get_qname
 * -----------------------------
 *
 * Gets the qname (read_name) for a read
 *
 * INPUT: pointer to binnie_read_t
 * OUTPUT: string (char *) pointing to qname
 */
char *br_get_qname (const binnie_read_t *br)
{
  char *qname;
  DLOG("br_get_qname()");

  qname = bam_get_qname(br->bam_read); 

  DLOG("br_get_qname: returning qname=[%s]", qname);
  return qname;
}

 

/*
 * br_equals
 * ------------------
 *
 * INPUT: pointers to two binnie_read_t reads to be compared
 * OUTPUT: bool (true if equal, false if not equal)
 *
 */
bool br_equals(const binnie_read_t *br1, const binnie_read_t *br2)
{
  bool equal;
  DLOG("br_equals()");
  
  DLOG("br_equals: ignore_rg=[%d] br1 rg=[%s] uid=[%s] br2 rg=[%s] uid=[%s]", ignore_rg, br_get_read_group(br1), br_get_uid_alloc(br1), br_get_read_group(br2), br_get_uid_alloc(br2));
  equal = ( (ignore_rg || (strcmp(br_get_read_group(br1), br_get_read_group(br2)) == 0))
            && (strcmp(br_get_qname(br1), br_get_qname(br2)) == 0) );
  
  DLOG("br_equals: returning [%d]", equal);
  return equal;
}


/*
 * br_get_uid_alloc
 * ------------------
 *
 * Generates a uid by concatenating read group id and qname 
 * separated by newlines (newlines are not allowed in any of the tags)
 *
 * INPUT: pointer to binnie_read_t 
 * OUTPUT: string (char*) allocated on the heap 
 *
 * SIDE EFFECTS: allocates memory for return pointer. caller must free.
 *
 */
char *br_get_uid_alloc(const binnie_read_t *br)
{
  char *uid;
  char *read_group;
  char *qname;

  DLOG("br_get_uid_alloc()");
  
  /*
  if (ignore_rg)
    {
      read_group = "";
    }
  else
    {
  */
  read_group = br_get_read_group(br);
  /*
    }
  */  

  qname = br_get_qname(br);

  if ( asprintf(&uid, "%s\n%s", read_group, qname) < 0 )
    {
      /* fatal error, cannot continue */
      err(BINNIE_EXIT_ERR_UID, gettext("br_get_uid_alloc: could not format uid"));
    }

  DLOG("br_get_uid_alloc: returning uid=[%s]", uid);
  return uid;
}


/*
 * br_dispose
 * -------------------
 *
 * INPUT: pointer to binnie_read_t read to be disposed of
 *
 */
void br_dispose(const binnie_read_t *br)
{
  DLOG("br_dispose()");

  /* destroy bam1_t struct pointed to by bam_read */
  bam_destroy1(br->bam_read);
  
  /* free the binnie_read_t struct itself */
  free((void *)br);

  DLOG("br_dispose: returning void");
}


/*
 * br_init
 * -------------------
 *
 * OUTPUT: pointer to binnie_read_t with memory allocated (call br_dispose to free)
 *
 */
binnie_read_t *br_init()
{
  binnie_read_t *br;
  
  DLOG("br_init()");
  br = xmalloc(sizeof(binnie_read_t));
  
  br->bam_read_present = false;
  br->bam_read         = bam_init1();

  DLOG("br_init: returning br");
  return br;
}


/*
 * bbr_equals
 * ------------------
 *
 * INPUT: pointers to two binnie_binned_read_t reads to be compared
 * OUTPUT: bool (true if equal, false if not equal)
 *
 */
static bool bbr_equals(const void *elt1, const void *elt2)
{
  const binnie_binned_read_t *bbr1;
  const binnie_binned_read_t *bbr2;
  bool equal;

  DLOG("bbr_equals()");
  bbr1 = elt1;
  bbr2 = elt2;
  
  equal = ( br_equals(bbr1->br, bbr2->br) ); 

  DLOG("bbr_equals: returning [%d]", equal);
  return equal;
}


/*
 * bbr_hashcode
 * --------------------
 *
 * INPUT: pointer to binnie_binned_read_t read to be hashed
 * OUTPUT: size_t hash of the read name
 *
 */
size_t bbr_hashcode(const void *elt)
{
  binnie_binned_read_t *bbr;
  size_t hashcode;
  char *uid;

  DLOG("bbr_hashcode()");
  bbr = elt;

  /* get uid */
  uid = br_get_uid_alloc(bbr->br);
  hashcode = hash_pjw(uid, BINNIE_TABLESIZE);
  DLOG("bbr_hashcode: have hashcode=[%d] for uid=[%s] tablesize=[%z]", hashcode, uid, BINNIE_TABLESIZE);
  free((void *)uid);
  
  DLOG("bbr_hashcode: returning hashcode=[%z]", hashcode);
  return hashcode;
}


/*
 * bbr_dispose
 * -------------------
 *
 * INPUT: pointer to binnie_binned_read_t read to be disposed of
 *
 */
void bbr_dispose(const void *elt)
{
  const binnie_binned_read_t *bbr;

  DLOG("bbr_dispose()");
  bbr = elt;
  
  /* call to br_dispose to delete the br struct and the bam_read struct */
  if (bbr->br)
    br_dispose(bbr->br);

  /* free the binnie_binned_read_t struct itself */
  free((void *)bbr);
  
  DLOG("bbr_dispose: returning void");
}


/*
 * bbr_init
 * -------------------
 *
 * OUTPUT: pointer to binnie_binned_read_t with memory allocated (call bbr_dispose to free)
 *
 */
binnie_binned_read_t *bbr_init(binnie_read_t *br)
{
  binnie_binned_read_t *bbr;
  
  DLOG("bbr_init()");

  bbr = xmalloc(sizeof(binnie_binned_read_t));
  
  bbr->bin                  = BINNIE_UNCHANGED;
  bbr->br                   = br;
  bbr->expected_mate_count  = br_get_num_segments(bbr->br) - 1;
  bbr->mate_count           = 0;
  bbr->original_refid       = br_get_refid(bbr->br);
  bbr->original_pos         = br_get_pos(bbr->br);
  bbr->next_mate            = NULL;
  bbr->prev_mate            = NULL;

  DLOG("bbr_init: returning bbr");
  return bbr;
}


