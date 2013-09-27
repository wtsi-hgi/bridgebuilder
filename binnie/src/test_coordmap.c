#include "binnie_coordmap.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "binnie_log.h"


int main()
{
  debug_flag = true;
  DLOG("Test is starting!");
  fflush(stdout);
  const char * filename = "/media/sf_nc6/scratch/hg18ToHg19.coordmap";

  volatile Range oldR = {107924340, 107934380, "chr12"};

  CoordMap *map = bc_read_file(filename);

  Range *newR = bc_map_range(map, &oldR);

  if (newR != NULL) {
    printf("New Pretty Range: %s\t%d\t%d", newR->id, newR->start, newR->end);
    return 0;
  } else {
    return 1;
  }
}