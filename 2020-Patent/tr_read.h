#ifndef TR_READ_H
#define TR_READ_H

#include "read.h"

#define TR_READ_DELAY 10
#define TR_READ_SEARCH 20

typedef enum TRReadState {
  TR_READ_STATE_ACTIVE,
  TR_READ_STATE_VARIANT,
  TR_READ_STATE_INACTIVE,
  TR_READ_STATE_OMITTED
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
  int delay;

} TRRead;

TRRead *new_tr_read();
void free_tr_read(TRRead *);
void copy_tr_read(TRRead *, Read *);
void print_tr_read(TRRead *);

#endif
