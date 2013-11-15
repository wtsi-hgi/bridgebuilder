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
    volatile Range* new = xmalloc(sizeof(*new));
    new->start = from_pos-1;
    new->end = from_pos-1;
    char *prefix_sn = xmalloc(LINE_LENGTH * sizeof(char));
    prefix_sn = strncat(prefix_sn, "chr", 4);
    new->id = strncat(prefix_sn, from_sn, LINE_LENGTH);
    free(from_sn);
    return new;
  } else {
    exit(1233);
  }
}

void doStuff(CoordMap *map) {
  volatile Range oldR = {100, 200, "chr1"};

  Range *newR = bc_map_range(map, &oldR);
  
  if (newR != NULL) {
    printf("New Pretty Range: %s\t%d\t%d\n", newR->id, newR->start, newR->end);
  } else {
    exit(1);
  }

  free(newR);
}

int main()
{
  DLOG("Test is starting!");
  fflush(stdout);
  const char * filename = "/lustre/scratch113/teams/hgi/users/nc6/svm/hg18ToHg19.coordmap";
  CoordMap *map = bc_read_file(filename);

  const char * otherFilename = "/lustre/scratch113/teams/hgi/users/nc6/svm/ibd_ichip.sites.chr-pos";
  const char * outFileName = "/lustre/scratch113/teams/hgi/users/nc6/svm/ibd_ichip.sites.hg19.chr-pos";

  FILE *infile = fopen(otherFilename, "r");
  FILE *out = fopen(outFileName, "w");

  if (infile == NULL) exit(1234);
  if (out == NULL) exit(1235);


  char line[LINE_LENGTH];
    // Read each line into the thing
  while (fgets(line, LINE_LENGTH-1, infile) != NULL) {
    Range* from = createRange(line);
    Range* to = bc_map_range(map, from);
    free(from);
    
    if (to != NULL) {
      fprintf(out, "%s\t%d\n", to->id, to->start+1);
      free(to);
    }
  }
  

  bc_free_coordmap(map);

  fclose(out);
  fclose(infile);
}

