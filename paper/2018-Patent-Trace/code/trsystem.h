#ifndef TRSYSTEM_H
#define TRSYSTEM_H

#include "generator.h"
#include "strand.h"

struct TRWindow {
  int         nrow;       // number of rows
  int         ncol;       // number of cols
  Nucleotide  weight[NUCLEOTIDE_NUM]; // weight of ACGTN
  Nucleotide *consensus;  // plurality consensus
};
typedef struct TRWindow TRWindow;

struct TRSystem {
  int       nreads;     // number of reads
  Strand  **reads;      // reads produced by the polynucleotide sequencer
  int      *pos_cmp;    // position of comparison
  int      *pos_lah;    // position of look ahead window
  TRWindow *win_cmp;    // position of comparison
  TRWindow *win_lah;    // look ahead window
  Strand   *consensus;  // plurality consensus strand
};
typedef struct TRSystem TRSystem;

TRWindow *trwin_new(const int nrow, const int ncol);
void trwin_free(TRWindow *win);
void trwin_weight_reset();
Nucleotide trwin_weight_cmp(Nucleotide weight[]);
void trwin_consensus(
    TRWindow *win,
    Strand *reads[],
    int *pos_cmp,
    Nucleotide (*weight_cmp)(Nucleotide [])
    );

TRSystem *trs_new(TRConfig trc, Strand *reads[]);
void trs_free(TRSystem *trs);
void trs_run(TRSystem *trs);
void trs_print(TRSystem *trs);

#endif
