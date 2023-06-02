#ifndef TR_WINDOW_H
#define TR_WINDOW_H

#include "nucleotide.h"
#include "tr_read.h"

typedef struct TRWindow {
/* public */
  int nrow;
  int ncol;

  Strand *(*get_consensus)(
      struct TRWindow *trw,
      TRRead *reads[],
      int *pcmp);

} TRWindow;

TRWindow *new_tr_window(int nrow, int ncol);
void free_tr_window(TRWindow *);

#endif
