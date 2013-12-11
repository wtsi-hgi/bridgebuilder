#include "binnie_coordmap.h"
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "xalloc.h"
#include "binnie_log.h"

#define LINE_LENGTH 256

Range* createRange(char* input) {
  char *from_sn = xmalloc(LINE_LENGTH * sizeof(char));
  int from_pos;
  int is2 = sscanf(input, "%s\t%d", from_sn, &from_pos);
  if (is2 == 2) {
    Range* new = xmalloc(sizeof(*new));
    new->start = from_pos-1;
    new->end = from_pos-1;
    new->id = strncat("chr", from_sn, LINE_LENGTH);
    return new;
  } else {
    fprintf(stderr, "%s\n", "Unable to construct range from input.");
    exit(1233);
  }
}

int main(int argc, int *argv[])
{
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
    fprintf(stderr, "Usage: liftover in mapFile [out].\nWhere [out] is not given writes to stdout.\n");
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

  // Read each line into the thing
  char line[LINE_LENGTH];
  while (fgets(line, LINE_LENGTH,in) != NULL) {
    Range* from = createRange(line);
    Range* to = bc_map_range(map, from);

    if (to != NULL) {
     fprintf(out, "%s\t%d\n", to->id, to->start+1);
    }
  }

  fclose(out);
  fclose(in);
}

