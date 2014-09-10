/*
 * brindley.c Brindley: standalone co-ordinate liftover.
 *
 * Copyright (c) 2013 Genome Research Ltd. 
 * Author: Nicholas Clarke <nicholas.clarke@sanger.ac.uk>
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
 /*
  * Simple tool to liftover co-ordinates. At the moment, takes in data of the form
  * chr\tposition and a liftover file and outputs the resulting chromosome and position.
 */

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

/* gnulib headers */
#include "progname.h"
#include "xalloc.h"

/* brindley includes */
#include "brindley_log.h"
#include "brindley_coordmap.h"

#define LINE_LENGTH 256

Range* createRange(char* input) {
  char *from_sn = xmalloc(LINE_LENGTH * sizeof(char));
  int from_pos;
  int is2 = sscanf(input, "%s\t%d", from_sn, &from_pos);
  if (is2 == 2) {
    volatile Range* new = xmalloc(sizeof(*new));
    new->start = from_pos-1;
    new->end = from_pos-1;
    char *prefix_sn = xmalloc(LINE_LENGTH * sizeof(char));
    // prefix_sn = strncat(prefix_sn, "chr", 4);
    new->id = strncat(prefix_sn, from_sn, LINE_LENGTH);
    free(from_sn);
    return new;
  } else {
    fprintf(stderr, "%s\n", "Unable to construct range from input.");
    exit(1233);
  }
}

int main(int argc, char *argv[])
{

  /* setup progname */
  set_program_name (argv[0]);

  FILE *in;
  FILE *out;
  char *mapFile;

  if (argc == 3) {
    // Assume in mapfile
    in = fopen(argv[1] ,"r");
    mapFile = argv[2];
    out = stdout;
  } else if (argc == 4) {
    // in mapfile out
    in = fopen(argv[1], "r");
    out = fopen(argv[3], "w");
    mapFile = argv[2];
  } else {
    fprintf(stderr, gettext("Usage: %s [options] <input> <liftover_map> [output]\n"), program_name);
  }

  CoordMap *map = bc_read_file(mapFile);

  if (!in) {
    fprintf(stderr, "%s\nFilename:%s", "Unable to read input file.", argv[1]);
    exit(1234);
  }
  if (!out) {
    fprintf(stderr, "%s\n", "Unable to open output file for writing.");
    exit(1235);
  }

  char line[LINE_LENGTH];
    // Read each line into the thing
  while (fgets(line, LINE_LENGTH-1, in) != NULL) {
    Range* from = createRange(line);
    Range* to = bc_map_range(map, from);
    free(from);
    
    if (to != NULL) {
      fprintf(out, "%s\t%d\n", to->id, to->start+1);
      free(to);
    } else {
      fprintf(out, ".\t.\n");
    }
  }

  bc_free_coordmap(map);

  fclose(out);
  fclose(in);
}

