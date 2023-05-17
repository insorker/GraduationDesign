#include "tr_system.h"
#include <malloc.h>
#include <assert.h>

/* public */
void    tr_system_set_reads(TRSystem *, Read *reads[], int nreads);
Strand *tr_system_get_consensus(TRSystem *);

/* private */
void        tr_system_set_window(TRSystem *);


TRSystem *new_tr_system(int trw_width) {
  TRSystem *trs = (TRSystem *)malloc(sizeof(TRSystem));

/* public */
  trs->set_reads = tr_system_set_reads;
  trs->get_consensus = tr_system_get_consensus;

/* private */
  trs->nreads = 0;
  trs->reads = NULL;
  trs->trw_width = trw_width;
  trs->trw = NULL;

  trs->set_window = tr_system_set_window;

  return trs;
}

void free_tr_system(TRSystem *trs) {
  for (int i = 0; i < trs->nreads; i++) {
    free_tr_read(trs->reads[i]);
  }
  free(trs->reads);
  free(trs->trw);
  free(trs);
}

void tr_system_set_reads(TRSystem *trs, Read *reads[], int nreads) {
  for (int i = 0; i < trs->nreads; i++) {
    free_tr_read(trs->reads[i]);
  }
  free(trs->reads);

  trs->nreads = nreads;
  trs->reads = (TRRead **)malloc(nreads * sizeof(TRRead *));

  for (int i = 0; i < nreads; i++) {
    trs->reads[i] = new_tr_read();
    copy_tr_read(trs->reads[i], reads[i]);
  }

  trs->set_window(trs);
}

Strand *tr_system_get_consensus(TRSystem *trs) {
  assert(trs->nreads != 0);
  assert(trs->reads != NULL);

  Strand *consensus = new_strand();
  int *pcmp = (int *)malloc(trs->nreads * sizeof(int));
  int nreads_left = trs->nreads;

  for (int i = 0; i < trs->nreads; i++) {
    trs->reads[i]->state = TR_READ_STATE_CREDIBLE;
  }
  for (int i = 0; i < trs->nreads; i++) {
    pcmp[i] = 0;
  }

  while (nreads_left) {
    /* initialize */
    Nucleotide single_concensus = NUCLEOTIDE_N;

    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state != TR_READ_STATE_OMMITED) {
        read->state = TR_READ_STATE_CREDIBLE;
      }
      if (pcmp[i] >= read->size(read)) {
        read->state = TR_READ_STATE_OMMITED;
        nreads_left--;
      }
    }
    if (nreads_left == 0) {
      break;
    }

    /* core algorithm */

    // get single concensus
    double weight[NUCLEOTIDE_SIZE] = {};
    double weight_max = 0;

    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state != TR_READ_STATE_OMMITED) {
        Nucleotide n = read->at(read, pcmp[i]);

        weight[n]++;
      }
    }
    
    for (int i = 0; i < NUCLEOTIDE_SIZE; i++) {
      if (weight[i] > weight_max) {
        single_concensus = i;
        weight_max = weight[i];
      }
    }

    consensus->push_back(consensus, single_concensus);

    // set reads' state
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state != TR_READ_STATE_OMMITED) {
        Nucleotide n = read->at(read, pcmp[i]);

        if (n != single_concensus) {
          read->state = TR_READ_STATE_VARIANT;
        }
      }
    }

    // get look-ahead window consensus
    Strand *win_consensus = trs->trw->get_consensus(trs->trw, trs->reads, pcmp);
    // get comparison consensus
    Strand *cmp_consensus = new_strand();

    cmp_consensus->push_back(cmp_consensus, single_concensus);
    for (int i = 0; i < win_consensus->size(win_consensus) - 1; i++) {
      cmp_consensus->push_back(
          cmp_consensus,
          win_consensus->at(
            win_consensus, i));
    }

    // move pcmp
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state == TR_READ_STATE_CREDIBLE) {
        pcmp[i]++;
      }
      else if (read->state == TR_READ_STATE_VARIANT) {
        // get read slice start from look-ahead window
        Strand *read_slice_win = new_strand();
        // get read slice start from comparison
        Strand *read_slice_cmp = new_strand();


        for (int j = 0; j < trs->trw_width; j++) {
          if (pcmp[i] + 1 + j < read->size(read)) {
            read_slice_win->push_back(read_slice_win, read->at(read, pcmp[i] + 1 + j));
            read_slice_cmp->push_back(read_slice_cmp, read->at(read, pcmp[i] + j));
          }
          else if (pcmp[i] + j < read->size(read)) {
            read_slice_win->push_back(read_slice_win, NUCLEOTIDE_N);
            read_slice_cmp->push_back(read_slice_cmp, read->at(read, pcmp[i] + j));
          }
          else {
            read_slice_win->push_back(read_slice_win, NUCLEOTIDE_N);
            read_slice_cmp->push_back(read_slice_cmp, NUCLEOTIDE_N);
          }
        }

        // substitution
        if (compare_strand(win_consensus, read_slice_win) == 0) {
          pcmp[i]++;
        }
        // deletion
        else if (compare_strand(win_consensus, read_slice_cmp) == 0) {
          // do nothing
        }
        // insertion
        else if (compare_strand(cmp_consensus, read_slice_win) == 0){
          pcmp[i] += 2;
        }
        else {
          read->state = TR_READ_STATE_OMMITED;
          nreads_left--;
        }

        free_strand(read_slice_win);
        free_strand(read_slice_cmp);
      }
    }
    
    free_strand(win_consensus);
    free_strand(cmp_consensus);
  }

  free(pcmp);

  return consensus;
}

void tr_system_set_window(TRSystem *trs) {
  free_tr_window(trs->trw);
  trs->trw = new_tr_window(trs->nreads, trs->trw_width);
}
