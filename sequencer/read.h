#ifndef READ_H
#define READ_H

#include "strand.h"
#include "vector.h"

typedef enum ReadErrorType {
  READ_ERROR_NONE,
  READ_ERROR_SUB,
  READ_ERROR_DEL,
  READ_ERROR_INS
} ReadErrorType;

typedef struct Read {
/* super */
  Strand *super;

/* public extends */
  int           (*size)(struct Read *);
  Nucleotide    (*at)(struct Read *, int index);
  void          (*push_back)(struct Read *, Nucleotide);
  void          (*pop_back)(struct Read *);
  void          (*clear)(struct Read *);

/* public */
  vector_t *error;
  vector_t *anchor;

} Read;

Read *new_read();
Read *copy_read(Read *);
void free_read(Read *);
void print_read(Read *);
void print_read_error(Read *);
void print_read_anchor(Read *);

#endif
