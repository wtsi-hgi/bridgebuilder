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

#include <assert.h>
#define MAX_DEPTH 50

typedef struct avl_node {
  Range * key;
  Range * data;
  int balance;
  struct avl_node *child[2];
} avl_node;

// Allocate a single AVL node.
avl_node* avl_single (Range* key, Range* data) {
  avl_node* n = malloc (sizeof *n);
  if (n != NULL) {
    rn->key=key;
    rn->data=data;
    rn->balance = 0;
    rn->child[0]=rn->child[1]=NULL;
  }

  return n;
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

Range *avl_lookup(avl_node *tree, Range* key) {
  int a = tree->key->start < key->start
  if (a && tree->key->end > key->end) {
    return tree->data;
  } else {
    return avl_lookup(tree->link[a], key);
  }
}

/* Insert a node into the tree. We walk down the tree until we insert, storing the path
   which we used to walk down. We then walk back up the same path checking whether we need
   to rebalance.
*/
void avl_insert(avl_node* tree, Range* key, Range* value) {

  // Iterator
  avl_node *i = tree;
  int idx = 0;
  // Store which nodes we visit
  avl_node *rev[MAX_DEPTH];
  // And which directions we take (corresponds to the direction to the next node.)
  int revd[MAX_DEPTH];

  for (;;) {
    // Record we've visited
    int dir = key->start > tree->key->start;
    revd[idx] = dir;
    rev[idx] = i;

    if (*i != NULL) {
      // Increment
      idx = idx + 1;
      // and walk
      i = i->child[dir];
    } else {
      break;
    }
  }

  // Insert the new child.
  i->child[revd[idx]] = avl_single(key, value);

  // Walk back up the path.
  while(top != 0) {
    // Update the balance
    rev[top]->balance += revd[top] == 0 ? 1 : -1;

    int bal = rev[top]->balance;
    if (bal == 0) {
      /* If we get a balance 0 during the walk up, then we have balanced a subtree
      and hence need walk no further.
      */
      break;
    } else if (bal < -1 || bal > 1) {
      // Out of balance, so we need to rebalance.
      rev[top - 1]->link[revd[top - 1]] = avl_insert_balance(rev[top], redv[top]);
      break;
    }
    top = top - 1;
  }

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

