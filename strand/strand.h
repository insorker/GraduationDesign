#ifndef STRAND_H
#define STRAND_H

#include "nucleotide.h"

typedef struct Strand {
/* public */
  int         (*size)(struct Strand *);
  Nucleotide  (*at)(struct Strand *, int index);
  void        (*push_back)(struct Strand *, Nucleotide);
  Nucleotide  (*pop_back)(struct Strand *);
  void        (*clear)(struct Strand *);

/* private */
  int length;
  int _size;
  Nucleotide *items;

  void (*capacity)(struct Strand *, int nlength);
  void (*expand)(struct Strand *);
  void (*shrink)(struct Strand *);

} Strand;

Strand *new_strand();
void free_strand(Strand *);
int  compare_strand(Strand *, Strand *);
void print_strand(Strand *);

#endif
