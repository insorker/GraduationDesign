#ifndef READ_H
#define READ_H

#include "strand.h"

typedef struct Read {
/* public extends */
  int         (*size)(struct Read *);
  Nucleotide  (*at)(struct Read *, int index);
  void        (*push_back)(struct Read *, Nucleotide);
  Nucleotide  (*pop_back)(struct Read *);
  void        (*clear)(struct Read *);

/* private */
  Strand *super;

} Read;

Read *new_read();
void free_read(Read *);
void print_read(Read *);

#endif
