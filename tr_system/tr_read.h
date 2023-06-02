#ifndef TR_READ_H
#define TR_READ_H

#include "read.h"

typedef enum TRReadState {
  TR_READ_STATE_ACTIVE,
  TR_READ_STATE_INACTIVE,
  TR_READ_STATE_VARIANT,
  TR_READ_STATE_OMITTED
} TRReadState;

typedef struct TRRead {
/* super */
  Read *super;

/* public extends */
  int         (*size)(struct TRRead *);
  Nucleotide  (*at)(struct TRRead *, int index);
  void        (*push_back)(struct TRRead *, Nucleotide);
  void        (*pop_back)(struct TRRead *);
  void        (*clear)(struct TRRead *);

/* public */
  TRReadState state;

} TRRead;

TRRead *new_tr_read();
void copy_tr_read(TRRead *, Read *);
void free_tr_read(TRRead *);

void print_tr_read(TRRead *);
void print_tr_read_state(TRRead *);

#endif
