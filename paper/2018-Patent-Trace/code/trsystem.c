#include "trsystem.h"
#include <malloc.h>
#include <string.h>
#define TRS_RUN_TEST 1

TRWindow *trwin_new(const int nrow, const int ncol) {
  TRWindow *win =
    (TRWindow *)malloc(sizeof(TRWindow));

  win->nrow = nrow;
  win->ncol = ncol;
  trwin_weight_reset(win);
  win->consensus =
    (Nucleotide *)malloc(sizeof(Nucleotide) * ncol);

  return win;
}

void trwin_free(TRWindow *win) {
  free(win->consensus);
  free(win);
}

void trwin_weight_reset(TRWindow *win) {
  memset(win->weight, 0, sizeof(win->weight));
}

/* deprecated, better using explicit rules */
Nucleotide trwin_weight_cmp(Nucleotide weight[]) {
  Nucleotide res = N;
  int cnt = 0;

  for (int i = 0; i < NUCLEOTIDE_NUM; i++) {
    if (weight[i] > cnt) {
      res = i;
      cnt = weight[i];
    }
  }

  return res;
}

void trwin_consensus(
    TRWindow *win,
    Strand *reads[],
    int *pos_cmp,
    Nucleotide (*weight_cmp)(Nucleotide [])
    )
{
  for (int i = 0; i < win->ncol; i++) {
    trwin_weight_reset(win);
    for (int j = 0; j < win->nrow; j++) {
      int p = pos_cmp[j] + i;

      if (reads[j]->state == OMITTED) continue;
      if (reads[j]->len <= p) {
        win->weight[N]++;
        continue;
      }

      win->weight[reads[j]->seq[p]]++;
    }
    win->consensus[i] = weight_cmp(win->weight);
  }
}

TRSystem *trs_new(TRConfig trc, Strand *reads[]) {
  TRSystem *trs = (TRSystem *)malloc(sizeof(TRSystem));

  trs->nreads = trc.num;
  trs->reads = reads;
  trs->pos_cmp = (int *)malloc(sizeof(int) * trc.num);
  memset(trs->pos_cmp, 0, sizeof(int) * trc.num);
  trs->pos_lah = (int *)malloc(sizeof(int) * trc.num);
  for (int i = 0; i < trc.num; i++) {
    trs->pos_lah[i] = 1;
  }
  trs->win_cmp = trwin_new(trc.num, 1);
  trs->win_lah = trwin_new(trc.num, trc.width);
  trs->consensus = strand_new(CONSENSUS);
  
  return trs;
}

void trs_free(TRSystem *trs) {
  free(trs->pos_cmp);
  free(trs->pos_lah);
  trwin_free(trs->win_cmp);
  trwin_free(trs->win_lah);
  strand_free(trs->consensus);
  free(trs);
}

void trs_run(TRSystem *trs) {
  Nucleotide baseCallConsensus = A;

  do {
    // position of comparison
    trwin_consensus(trs->win_cmp, trs->reads, trs->pos_cmp,
        trwin_weight_cmp);
    baseCallConsensus = trs->win_cmp->consensus[0];
    if (baseCallConsensus == N) {
      break;
    }
    else {
      strand_append(trs->consensus, baseCallConsensus);
    }

    for (int i = 0; i < trs->nreads; i++) {
      Strand *read = trs->reads[i];

      if (read->state == OMITTED) continue;
      if (read->len <= trs->pos_cmp[i]) {
        read->state = OMITTED;
        continue;
      }
      if (read->seq[trs->pos_cmp[i]] != baseCallConsensus) {
        read->state = VARIANT;
#if TRS_RUN_TEST
      printf("Variant read %2d at %3d\n", i, trs->pos_cmp[i]);
#endif
        continue;
      }
    }

    // look ahead window
    trwin_consensus(trs->win_lah, trs->reads, trs->pos_lah,
        trwin_weight_cmp);
// #if TRS_RUN_TEST
//     nucleotide_print(trs->win_lah->consensus, trs->win_lah->ncol);
// #endif

    for (int i = 0; i < trs->nreads; i++) {
      Strand *read = trs->reads[i];

      if (read->state == OMITTED) continue;
      else if (read->state == VARIANT) {
        Nucleotide *rn;   // read nucleotide
        Nucleotide *wn;   // win nucleotide
        int rl;           // read nucleotide length
        int wl;           // win nucleotide length

        wn = trs->win_lah->consensus;
        wl = trs->win_lah->ncol;

        // substitution
        rn = read->seq + trs->pos_lah[i];
        rl = read->len - trs->pos_lah[i] < wl ?
          read->len - trs->pos_lah[i] : wl;
// #if TRS_RUN_TEST
//     nucleotide_print(rn, rl);
// #endif
        if (nucleotide_cmp(rn, wn, rl, wl) == 0) {
          trs->pos_cmp[i]++;
          trs->pos_lah[i]++;
          read->state = CREDIBLE;
#if TRS_RUN_TEST
          printf("%2d: sub\n", i);
#endif
          continue;
        }

        // deletion
        rn = read->seq + trs->pos_cmp[i];
        rl = read->len - trs->pos_cmp[i] < wl ?
          read->len - trs->pos_cmp[i] : wl;
        if (nucleotide_cmp(rn, wn, rl, wl) == 0) {
          read->state = CREDIBLE;
#if TRS_RUN_TEST
          printf("%2d: del\n", i);
#endif
          continue;
        }

        // insertion
        if (trs->pos_lah[i] < read->len &&
            read->seq[trs->pos_lah[i]] == baseCallConsensus)
        {
          rn = read->seq + trs->pos_lah[i] + 1;
          rl = read->len - trs->pos_lah[i] + 1 ?
            read->len - trs->pos_lah[i] + 1 : wl;
          if (nucleotide_cmp(rn, wn, rl, wl) == 0) {
            trs->pos_cmp[i] += 2;
            trs->pos_lah[i] += 2;
            read->state = CREDIBLE;
#if TRS_RUN_TEST
          printf("%2d: ins\n", i);
#endif
            continue;
          }
        }

#if TRS_RUN_TEST
          printf("%2d: omitted\n", i);
#endif

        read->state = OMITTED;
      }
      else {
        trs->pos_cmp[i]++;
        trs->pos_lah[i]++;
        continue;
      }
    }
    // strand_print(trs->consensus);
  } while (baseCallConsensus != N);
}

void trs_print(TRSystem *trs) {
  strand_print(trs->consensus);
}
