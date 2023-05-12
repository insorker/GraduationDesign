#include "tr_window.h"
#include <malloc.h>

Strand *tr_window_get_consensus(
    TRWindow *trw,
    TRRead *reads[],
    int *reads_pcmp);

Nucleotide *tr_window_new_consensus(TRWindow *);

TRWindow *new_tr_window(int nrow, int ncol) {
  TRWindow *trw = (TRWindow *)malloc(sizeof(TRWindow));

/* public */
  trw->nrow = nrow;
  trw->ncol = ncol;

  trw->get_consensus = tr_window_get_consensus;
  
  return trw;
}

void free_tr_window(TRWindow *trw) {
  free(trw);
}

Strand *tr_window_get_consensus(
    TRWindow *trw,
    TRRead *reads[],
    int *reads_pcmp)
{
  Strand *consensus = new_strand();

  for (int i = 0; i < trw->ncol; i++) {
    double weight[NUCLEOTIDE_SIZE] = {};
    Nucleotide weight_max = 0;

    for (int j = 0; j < trw->nrow; j++) {
      int p = reads_pcmp[j] + i;

      if (reads[j]->state == TR_READ_STATE_CREDIBLE) {
        if (reads[j]->size(reads[j]) <= p) {
          weight[NUCLEOTIDE_N]++;
        }
        else {
          Nucleotide n = reads[j]->at(reads[j], p);

          weight[n]++;
        }
      }
      else {
        continue;
      }
    }

    for (int j = 0; j < NUCLEOTIDE_SIZE; j++) {
      if (weight[j] > weight_max) {
        consensus->push_back(consensus, j);
      }
    }
  }

  return consensus;
}
