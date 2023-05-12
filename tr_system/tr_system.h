#ifndef TR_SYSTEM_H
#define TR_SYSTEM_H

#include "strand.h"
#include "read.h"
#include "tr_read.h"
#include "tr_window.h"

typedef struct TRSystem {
/* public */
  void (*set_reads)(struct TRSystem *, Read *reads[], int nreads);
  Strand *(*get_consensus)(struct TRSystem *);

/* private */
  int nreads;       // number of reads
  TRRead **reads;   // list of reads
  int trw_width;    // look-ahead window width
  TRWindow *trw;    // look-ahead window
  
  void        (*set_window)(struct TRSystem *);

} TRSystem;

TRSystem *new_tr_system(int trw_width);
void free_tr_system(TRSystem *);

#endif
