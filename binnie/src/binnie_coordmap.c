/*
 * binnie_coordmap.c Binnie co-ordinate mapping.
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
#include "config.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "binnie.h"
#include "binnie_log.h"
#include "gl_xlist.h"
#include "gl_avltreehash_list.h"
#include "xalloc.h"
#include "hash-pjw.h"
#include "binnie_coordmap.h"

#define LINE_LENGTH 256
#define LINE_SEP = '\t'
#define MAX_DEPTH 500

typedef struct CoordMap {
  gl_list_t * entries;
} CoordMap;

typedef struct avl_node {
  Range * key;
  Range * data;
  int balance;
  struct avl_node *child[2];
} avl_node;

typedef struct {
  char * key;
  avl_node * data;
} entry;

typedef struct {
  Range* from;
  Range* to;
} RangeMap;

void avl_free(avl_node *thing);

void bc_free_coordmap(CoordMap* coordMap) {
 gl_list_free(*(coordMap->entries));
 free(coordMap->entries);
 free(coordMap);
}

/*
 * entry_equals
 * ------------------
 *
 * INPUT: pointers to two entry structs to be compared
 * OUTPUT: bool (true if equal, false if not equal)
 *
 */
static bool entry_equals(const void *elt1, const void *elt2)
{
  const entry *bbr1;
  const entry *bbr2;

  bbr1 = elt1;
  bbr2 = elt2;
  
  return strcmp(bbr1->key, bbr2->key) == 0;
}

/*
 * entry_hashcode
 * --------------------
 *
 * INPUT: pointer to entry read to be hashed
 * OUTPUT: size_t hash of the read name
 *
 */
size_t entry_hashcode(const void *elt)
{
  size_t hashcode;

  // DLOG("entry_hashcode()");
  const entry *bbr = elt;

  /* get uid */
  char* uid = bbr->key;
  // DLOG("entry_hashcode: calling hash_pjw on uid=[%s]", uid);
  hashcode = hash_pjw(uid, BINNIE_TABLESIZE);
  
  // DLOG("bbr_hashcode: have hashcode=[%zu] for uid=[%s] tablesize=[%zu]", hashcode, uid, BINNIE_TABLESIZE);
  
  // DLOG("bbr_hashcode: returning hashcode=[%zu]", hashcode);
  return hashcode;
}


/*
 * entry_dispose
 * -------------------
 *
 * INPUT: pointer to binnie_binned_read_t read to be disposed of
 *
 */
void entry_dispose(const void *elt)
{
  const entry *bbr;

  DLOG("entry_dispose()");
  bbr = elt;

  /* free the entry struct itself */
  free(bbr->key);
  avl_free(bbr->data);
  free(bbr);

  DLOG("entry_dispose: returning void");
}

Range * clone_range(Range* old) {
  Range* new = xmalloc(sizeof(*new));
  new->start = old->start;
  new->end = old->end;
  new->id = old->id;
  return new;
}

// Allocate a single AVL node.
avl_node* avl_single (Range* key1, Range* data1) {
  avl_node* rn = xmalloc (sizeof *rn);
  Range *key = clone_range(key1);
  Range *data = clone_range(data1);

  if (rn != NULL) {
    rn->key=key;
    rn->data=data;
    rn->balance = 0;
    rn->child[0]=rn->child[1]=NULL;
  }
  return rn;
}

void avl_free(avl_node *tree) {
  if (tree->child[0]!=NULL) {
    avl_free(tree->child[0]);
  }
  if (tree->child[1]!=NULL) {
    avl_free(tree->child[1]);
  }
  free(tree->key->id);
  free(tree->key);
  free(tree->data->id);
  free(tree->data);
  free(tree);

}

int sgn(int x) {
  if (x > 0) {
    return 1;
  } else if (x < 0){
    return -1;
  } else {
    return 1;
  }
}

/*
  Look up a node in the tree. Returns NULL (yuck!) if the tree does not contain the item.
*/
RangeMap *avl_lookup(avl_node *tree, Range* key) {
  int a = tree->key->start < key->start;
  if (a && tree->key->end > key->end) {
    RangeMap *r = xmalloc(sizeof *r);
    r->from = tree->key;
    r->to = tree->data;
    return r;
  } else if (a && tree->key->end > key->start) {
    return NULL;
  } else {
    if (tree->child[a] == NULL) {
      return NULL;
    } else {
     return avl_lookup(tree->child[a], key);
    }
  }
}

// Rotate root in dir and return the new root
avl_node * avl_rotate(avl_node* p, int dir) {
  avl_node *c = p->child[dir];
  p->child[dir] = c->child[!dir];
  c->child[!dir] = p;

  // Update balance!
  int d = (dir == 0) ? -1 : 1;
  if (c->balance * d < 0) {
    p->balance = p->balance + d - c->balance;
  } else {
    p->balance = p->balance  + d;
  }

  if (p->balance * d > 0) {
    c->balance = c->balance + p->balance + d;
  } else {
    c->balance = c->balance + d;
  }

  return c;
}

// Balance the tree after an insert.
// Dir indicates the side of the tree with a positive imbalance - e.g.
// dir=0 suggests a left imbalance, dir=1 suggests a right imbalance.
avl_node * avl_insert_balance(avl_node* p, int dir) {
  int pbal = dir == 0 ? 1 : -1;
  // We need to test whether the sign of c's balance agrees with the dir
  if (sgn(p->child[dir]->balance) != pbal) {
    avl_node *d = avl_rotate(p->child[dir], !dir);
    p->child[dir] = d;
  }
  return avl_rotate(p, dir);
}

/* Insert a node into the tree. We walk down the tree until we insert, storing the path
   which we used to walk down. We then walk back up the same path checking whether we need
   to rebalance.
*/
avl_node* avl_insert(avl_node* tree, Range* key, Range* value) {

  // Iterator
  avl_node *i = tree;
  int idx = 0;
  // Store which nodes we visit
  avl_node *rev[MAX_DEPTH];
  // And which directions we take (corresponds to the direction to the next node.)
  int revd[MAX_DEPTH];

  for (;;) {
    if (i != NULL) {
      // Record we've visited
      int dir = key->start > i->key->start;
      revd[idx] = dir;
      rev[idx] = i;

      // Increment
      idx = idx + 1;
      // and walk
      i = i->child[dir];
    } else {
      idx = idx - 1;
      i = rev[idx];
      break;
    }
  }

  // Insert the new child.
  i->child[revd[idx]] = avl_single(key, value);

  // Walk back up the path.
  while(idx != 0) {
    // Update the balance
    rev[idx]->balance += revd[idx] == 0 ? 1 : -1;

    int bal = rev[idx]->balance;
    if (bal == 0) {
      /* If we get a balance 0 during the walk up, then we have balanced a subtree
      and hence need walk no further.
      */
      break;
    } else if (bal < -1 || bal > 1) {
      // Out of balance, so we need to rebalance.
      rev[idx - 1]->child[revd[idx - 1]] = avl_insert_balance(rev[idx], revd[idx]);
      break;
    }
    idx = idx - 1;
  }

  return rev[idx];
}

CoordMap* bc_read_file(const char *filename) {
  DLOG("bc_read_file()");
  gl_list_t *map = xmalloc(sizeof(gl_list_t));
  *map = gl_list_create_empty(GL_AVLTREEHASH_LIST, 
			      entry_equals, 
			      entry_hashcode, 
			      entry_dispose, 
			      true);
  FILE *fp = fopen(filename, "r");

  if (!fp) exit(1234);

  // Ignore header
  while (fgetc(fp) != '\n') {}
  // Read each line into the thing
  char line[LINE_LENGTH];
  while (fgets(line, LINE_LENGTH, fp) != NULL) {
    // Parse the line
    // Tab separated
    // from_sn from_start      from_end        to_sn   to_start        to_end
    char *from_sn = xmalloc(LINE_LENGTH * sizeof(char));
    char *to_sn = xmalloc(LINE_LENGTH * sizeof(char));
    int from_start, from_end, to_start, to_end;
    sscanf(line, "%s\t%d\t%d\t%s\t%d\t%d", from_sn, &from_start, &from_end, to_sn, &to_start, &to_end);
    DLOG("Input line: %s\t%d\t%d\t%s\t%d\t%d", from_sn, from_start, from_end, to_sn, to_start, to_end);
    Range from_range = { from_start, from_end, from_sn };
    Range to_range = { to_start, to_end, to_sn };
    entry *e_bad = xmalloc(sizeof(entry));
    e_bad->key=from_sn;
    e_bad->data=NULL;

    // Check whether we already have something for this key.
    gl_list_node_t n = gl_list_search(*map, e_bad);
    free(e_bad);
    if (n == NULL) {
      avl_node * newtree = avl_single(&from_range, &to_range);
      entry* e1 = xmalloc(sizeof(entry));
      e1->key=strndup(from_sn, LINE_LENGTH-1);
      e1->data = newtree;
      gl_list_add_last(*map, e1);
    } else {
      entry* e = (entry*) gl_list_node_value(*map, n);
      avl_insert(e->data, &from_range, &to_range);
    }
  }

  fclose(fp);
  CoordMap *cm = xmalloc(sizeof(CoordMap));
  cm->entries = map;
  return cm;
}

Range* bc_map_range(CoordMap* coordMap, Range* oldRef) {
  DLOG("bc_map_range()");
  char * key = oldRef->id;
  //  entry e_bad = {key, NULL};
  entry *e_bad = xmalloc(sizeof(entry));
  e_bad->key=key;
  e_bad->data=NULL;
  gl_list_t* map = coordMap->entries;
  gl_list_node_t n = gl_list_search(*map, e_bad);
  free(e_bad);

  Range *mappedFrom;
  Range *mappedTo;
  if (n == NULL) {
    return NULL;
  } else {
    entry* e = (entry*) gl_list_node_value(*map, n);
    RangeMap* rm = avl_lookup(e->data, oldRef);
    if (rm == NULL) {
      return NULL;
    } else {
      mappedFrom = rm->from;
      mappedTo = rm->to;
      free(rm);
    }
  }

  // Look up offsets in mapped range.
  int toStart = mappedTo->start;
  int toEnd = mappedTo->end;
  if (toStart > toEnd) {
    // If I must program in C, I'm at least going to do silly hacky C things like doing this
    // swap to avoid the 'expense' of another variable.
    toEnd = toStart + toEnd;
    toStart = toEnd - toStart;
    toEnd = toEnd - toStart;
  }
  int start = oldRef->start - mappedFrom->start + toStart;
  int end = start + (oldRef->end - oldRef->start);

  Range* newRef = xmalloc(sizeof *newRef);
  newRef->start = start;
  newRef->end = end;
  newRef->id = mappedTo->id;

  return newRef;
}
