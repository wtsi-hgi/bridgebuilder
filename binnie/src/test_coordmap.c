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
    exit(1233);
  }
}

int main()
{
  debug_flag = false;
  DLOG("Test is starting!");
  fflush(stdout);
  const char * filename = "/lustre/scratch113/teams/hgi/users/nc6/svm/hg18ToHg19.coordmap";
  CoordMap *map = bc_read_file(filename);
  
  const char * otherFilename = "/lustre/scratch113/teams/hgi/users/nc6/svm/ibd_ichip.sites.chr-pos";
  const char * outFileName = "/lustre/scratch113/teams/hgi/users/nc6/svm/ibd_ichip.sites.hg19.chr-pos";

  FILE *in = fopen(otherFilename, "r");
  FILE *out = fopen(outFileName, "w");

  if (!in) exit(1234);
  if (!out) exit(1235);

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

