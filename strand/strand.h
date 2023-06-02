#ifndef STRAND_H
#define STRAND_H

#include "vector.h"
#include "nucleotide.h"
#include <stdbool.h>

typedef struct Strand {
/* public */
  int         (*size)(struct Strand *);
  Nucleotide  (*at)(struct Strand *, int index);
  void        (*push_back)(struct Strand *, Nucleotide);
  void        (*pop_back)(struct Strand *);
  void        (*clear)(struct Strand *);

/* private */
  vector_t *_nucleotides;

} Strand;

/* friend */
bool compare_strand(Strand *, Strand *);

Strand *new_strand();
void free_strand(Strand *);
void print_strand(Strand *);

#endif
