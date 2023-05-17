#ifndef TR_READ_H
#define TR_READ_H

#include "read.h"

typedef enum TRReadState {
  TR_READ_STATE_CREDIBLE,
  TR_READ_STATE_VARIANT,
  TR_READ_STATE_OMMITED
} TRReadState;

typedef struct TRRead {
/* super */
  Read *super;

/* public extends */
  int         (*size)(struct TRRead *);
  Nucleotide  (*at)(struct TRRead *, int index);
  void        (*push_back)(struct TRRead *, Nucleotide);
  Nucleotide  (*pop_back)(struct TRRead *);
  void        (*clear)(struct TRRead *);

/* public */
  TRReadState state;

} TRRead;

TRRead *new_tr_read();
void free_tr_read(TRRead *);
void copy_tr_read(TRRead *, Read *);
void print_tr_read(TRRead *);

#endif
