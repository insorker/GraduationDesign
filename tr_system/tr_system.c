#include "tr_system.h"
#include <malloc.h>
#include <assert.h>

/* public */
void    tr_system_set_reads(TRSystem *, Read *reads[], int nreads);
Strand *tr_system_get_consensus(TRSystem *);

/* private */
void    tr_system_set_window(TRSystem *);


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

#if TR_NORMAL == 1
Strand *tr_system_get_consensus(TRSystem *trs) {
  assert(trs->nreads != 0);
  assert(trs->reads != NULL);

  Strand *consensus = new_strand();
  int *pcmp = (int *)malloc(trs->nreads * sizeof(int));

  for (int i = 0; i < trs->nreads; i++) {
    trs->reads[i]->state = TR_READ_STATE_ACTIVE;
    pcmp[i] = 0;
  }

  while (1) {
    /* initialize */
    int nreads_active = trs->nreads;
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state == TR_READ_STATE_OMITTED) {
        nreads_active--;
      }
      else if (pcmp[i] >= read->size(read)) {
        read->state = TR_READ_STATE_OMITTED;
        nreads_active--;
      }
      else if (read->state != TR_READ_STATE_OMITTED) {
        read->state = TR_READ_STATE_ACTIVE;
      }
    }
    if (nreads_active == 0) {
      break;
    }

    /* core algorithm */

    // get single concensus
    Nucleotide single_consensus = NUCLEOTIDE_N;
    double weight[NUCLEOTIDE_SIZE] = {};
    double weight_max = 0;

    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state != TR_READ_STATE_OMITTED) {
        Nucleotide n = read->at(read, pcmp[i]);

        weight[n]++;
      }
    }
    
    for (int i = 0; i < NUCLEOTIDE_SIZE; i++) {
      if (weight[i] > weight_max) {
        single_consensus = i;
        weight_max = weight[i];
      }
    }

    consensus->push_back(consensus, single_consensus);

    // set variant reads' state
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state != TR_READ_STATE_OMITTED) {
        Nucleotide n = read->at(read, pcmp[i]);

        if (n != single_consensus) {
          read->state = TR_READ_STATE_VARIANT;
        }
      }
    }

    // get look-ahead window consensus
    Strand *win_consensus = trs->trw->get_consensus(trs->trw, trs->reads, pcmp);
    // get comparison consensus
    Strand *cmp_consensus = new_strand();

    cmp_consensus->push_back(cmp_consensus, single_consensus);
    for (int i = 0; i < win_consensus->size(win_consensus) - 1; i++) {
      cmp_consensus->push_back(
          cmp_consensus,
          win_consensus->at(win_consensus, i));
    }

    // move pcmp
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state == TR_READ_STATE_ACTIVE) {
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
        if (compare_strand(win_consensus, read_slice_win) == true) {
          pcmp[i]++;
        }
        // deletion
        else if (compare_strand(win_consensus, read_slice_cmp) == true) {
          // do nothing
        }
        // insertion
        else if (compare_strand(cmp_consensus, read_slice_win) == true){
          pcmp[i] += 2;
        }
        else {
          read->state = TR_READ_STATE_OMITTED;

#if TR_VERBOSE
          printf("%d\n", pcmp[i]);
          print_tr_read(read);
          printf("- "),
            print_strand(consensus);
          printf("- win-c"),
            print_strand(win_consensus);
          printf("- cmp-c"),
            print_strand(win_consensus);
          printf("- win-r"),
            print_strand(read_slice_win);
          printf("- cmp-r"),
            print_strand(read_slice_cmp);
#endif
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
#endif

#if TR_INDETERMINATED == 1 && TR_ALIGNMENT == 0
#define TR_READ_DELAY   10
#define TR_READ_SEARCH  20
Strand *tr_system_get_consensus(TRSystem *trs) {
  assert(trs->nreads != 0);
  assert(trs->reads != NULL);

  Strand *consensus = new_strand();
  int *pcmp = (int *)malloc(trs->nreads * sizeof(int));
  int *delay = (int *)malloc(trs->nreads * sizeof(int));

  for (int i = 0; i < trs->nreads; i++) {
    trs->reads[i]->state = TR_READ_STATE_ACTIVE;
    pcmp[i] = 0;
    delay[i] = TR_READ_DELAY;
  }

  while (1) {
    /* initialize */
    int nreads_active = trs->nreads;
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state == TR_READ_STATE_OMITTED) {
        nreads_active--;
      }
      else if (pcmp[i] >= read->size(read)) {
        read->state = TR_READ_STATE_OMITTED;
        nreads_active--;
      }
      else if (read->state == TR_READ_STATE_INACTIVE) {
        nreads_active--;

        if (--(delay[i]) == 0) {
          int pcmp_next = 0;

          for (int j = pcmp[i] + TR_READ_DELAY; j <= pcmp[i] + TR_READ_SEARCH; j++) {
            if (j + 1 < read->size(read)) {
              int fault = 0, min_fault = TR_READ_DELAY + 1;

              for (int k = 0; k < TR_READ_DELAY; k++) {
                if (consensus->at(consensus, consensus->size(consensus) - 1 - k)
                    != read->at(read, j - k))
                {
                  fault++;
                }
              }
              if (fault < min_fault && fault <= TR_READ_DELAY / 5) {
                min_fault = fault;
                pcmp_next = j;
              }
            }
          }
          if (pcmp_next != 0) {
#if TR_VERBOSE
            printf("\nset active %d\n", i);
            printf("  pcmp:%d ", pcmp[i]);
#endif
            pcmp[i] = pcmp_next + 1;
            read->state = TR_READ_STATE_ACTIVE;
            nreads_active++;
#if TR_VERBOSE
            printf("%d\n", pcmp[i]);
            print_strand(consensus);
            print_tr_read(read);
#endif
          }
          else {
            read->state = TR_READ_STATE_OMITTED;
          }

          delay[i] = TR_READ_DELAY;
        }
      }
      else if (read->state != TR_READ_STATE_OMITTED) {
        read->state = TR_READ_STATE_ACTIVE;
      }
    }
    if (nreads_active == 0) {
      break;
    }

    /* core algorithm */

    // get single concensus
    Nucleotide single_consensus = NUCLEOTIDE_N;
    double weight[NUCLEOTIDE_SIZE] = {};
    double weight_max = 0;

    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state != TR_READ_STATE_OMITTED
          && read->state != TR_READ_STATE_INACTIVE)
      {
        Nucleotide n = read->at(read, pcmp[i]);

        weight[n]++;
      }
    }
    
    for (int i = 0; i < NUCLEOTIDE_SIZE; i++) {
      if (weight[i] > weight_max) {
        single_consensus = i;
        weight_max = weight[i];
      }
    }

    consensus->push_back(consensus, single_consensus);

    // set variant reads' state
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state != TR_READ_STATE_OMITTED
          && read->state != TR_READ_STATE_INACTIVE) {
        Nucleotide n = read->at(read, pcmp[i]);

        if (n != single_consensus) {
          read->state = TR_READ_STATE_VARIANT;
        }
      }
    }

    // get look-ahead window consensus
    Strand *win_consensus = trs->trw->get_consensus(trs->trw, trs->reads, pcmp);
    // get comparison consensus
    Strand *cmp_consensus = new_strand();

    cmp_consensus->push_back(cmp_consensus, single_consensus);
    for (int i = 0; i < win_consensus->size(win_consensus) - 1; i++) {
      cmp_consensus->push_back(
          cmp_consensus,
          win_consensus->at(win_consensus, i));
    }

    // move pcmp
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state == TR_READ_STATE_ACTIVE) {
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
        if (compare_strand(win_consensus, read_slice_win) == true) {
          pcmp[i]++;
        }
        // deletion
        else if (compare_strand(win_consensus, read_slice_cmp) == true) {
          // do nothing
        }
        // insertion
        else if (compare_strand(cmp_consensus, read_slice_win) == true){
          pcmp[i] += 2;
        }
        else {
          read->state = TR_READ_STATE_INACTIVE;

#if TR_VERBOSE
          printf("%d\n", pcmp[i]);
          print_tr_read(read);
          printf("- "),
            print_strand(consensus);
          printf("- win-c"),
            print_strand(win_consensus);
          printf("- cmp-c"),
            print_strand(win_consensus);
          printf("- win-r"),
            print_strand(read_slice_win);
          printf("- cmp-r"),
            print_strand(read_slice_cmp);
#endif
        }

        free_strand(read_slice_win);
        free_strand(read_slice_cmp);
      }
    }
    
    free_strand(win_consensus);
    free_strand(cmp_consensus);
  }

  free(pcmp);
  free(delay);

  return consensus;
}
#endif

#if TR_INDETERMINATED == 0 && TR_ALIGNMENT == 1
Strand *tr_system_get_partial_consensus(TRSystem *trs) {
  assert(trs->nreads != 0);
  assert(trs->reads != NULL);

  Strand *consensus = new_strand();
  int *pcmp = (int *)malloc(trs->nreads * sizeof(int));

  for (int i = 0; i < trs->nreads; i++) {
    trs->reads[i]->state = TR_READ_STATE_ACTIVE;
    pcmp[i] = 0;
  }

  while (1) {
    /* initialize */
    int nreads_active = trs->nreads;
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state == TR_READ_STATE_OMITTED) {
        nreads_active--;
      }
      else if (pcmp[i] >= read->size(read)) {
        read->state = TR_READ_STATE_OMITTED;
        nreads_active--;
      }
      else if (read->state != TR_READ_STATE_OMITTED) {
        read->state = TR_READ_STATE_ACTIVE;
      }
    }
    if (nreads_active == 0) {
      break;
    }

    /* core algorithm */

    // get single concensus
    Nucleotide single_consensus = NUCLEOTIDE_N;
    double weight[NUCLEOTIDE_SIZE] = {};
    double weight_max = 0;

    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state != TR_READ_STATE_OMITTED) {
        Nucleotide n = read->at(read, pcmp[i]);

        weight[n]++;
      }
    }
    
    for (int i = 0; i < NUCLEOTIDE_SIZE; i++) {
      if (weight[i] > weight_max) {
        single_consensus = i;
        weight_max = weight[i];
      }
    }

    consensus->push_back(consensus, single_consensus);

    // set variant reads' state
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state != TR_READ_STATE_OMITTED) {
        Nucleotide n = read->at(read, pcmp[i]);

        if (n != single_consensus) {
          read->state = TR_READ_STATE_VARIANT;
        }
      }
    }

    // get look-ahead window consensus
    Strand *win_consensus = trs->trw->get_consensus(trs->trw, trs->reads, pcmp);
    // get comparison consensus
    Strand *cmp_consensus = new_strand();

    cmp_consensus->push_back(cmp_consensus, single_consensus);
    for (int i = 0; i < win_consensus->size(win_consensus) - 1; i++) {
      cmp_consensus->push_back(
          cmp_consensus,
          win_consensus->at(win_consensus, i));
    }

    // move pcmp
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state == TR_READ_STATE_ACTIVE) {
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
        if (compare_strand(win_consensus, read_slice_win) == true) {
          pcmp[i]++;
        }
        // deletion
        else if (compare_strand(win_consensus, read_slice_cmp) == true) {
          // do nothing
        }
        // insertion
        else if (compare_strand(cmp_consensus, read_slice_win) == true){
          pcmp[i] += 2;
        }
        else {
          read->state = TR_READ_STATE_OMITTED;

#if TR_VERBOSE
          printf("%d\n", pcmp[i]);
          print_tr_read(read);
          printf("- "),
            print_strand(consensus);
          printf("- win-c"),
            print_strand(win_consensus);
          printf("- cmp-c"),
            print_strand(win_consensus);
          printf("- win-r"),
            print_strand(read_slice_win);
          printf("- cmp-r"),
            print_strand(read_slice_cmp);
#endif
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
Strand *tr_system_get_consensus(TRSystem *trs) {
  assert(trs->nreads != 0);
  assert(trs->reads != NULL);

  Strand *consensus = new_strand();
  int nreads = trs->nreads;
  int nanchor = trs->reads[0]->super->anchor->size(trs->reads[0]->super->anchor);
  int *panchor = (int *)malloc(sizeof(int) * nreads);
  TRRead **tr_reads = trs->reads;
  TRRead **reads = (TRRead **)malloc(sizeof(TRRead *) * trs->nreads);
  Read ***reads_slice = (Read ***)malloc(sizeof(Read **) * (nanchor + 1));

  for (int i = 0; i < nreads; i++) {
    reads[i] = new_tr_read();
    panchor[i] = -1;

    for (int j = 0; j < tr_reads[i]->size(tr_reads[i]); j++) {
      reads[i]->push_back(reads[i], tr_reads[i]->at(tr_reads[i], j));
    }

    vector_t *error = tr_reads[i]->super->error;
    for (int j = 0; j < error->size(error); j++) {
      reads[i]->super->error->push_back(
          reads[i]->super->error,
          error->at(error, j));
    }

    vector_t *anchor = tr_reads[i]->super->anchor;
    for (int j = 0; j < anchor->size(anchor); j++) {
      reads[i]->super->anchor->push_back(
          reads[i]->super->anchor,
          anchor->at(anchor, j));
    }
  }

  for (int i = 0; i <= nanchor; i++) {
    reads_slice[i] = (Read **)malloc(sizeof(Read *) * trs->nreads);

    for (int j = 0; j < nreads; j++) {
      reads_slice[i][j] = new_read();
      vector_t *anchor = reads[j]->super->anchor;
      int st = panchor[j] == -1 ? 0 : *(int *)anchor->at(anchor, panchor[j]);
      int ed = ++panchor[j] < nanchor ? *(int *)anchor->at(anchor, panchor[j]) : reads[j]->size(reads[j]);
      ed += trs->trw_width;

      // printf("%d %d\n", st, ed);
      // for (int i = 0; i < anchor->size(anchor); i++) {
      //   printf("%d ", *(int *)anchor->at(anchor, i));
      // }
      // printf("\n");
      for (int k = st; k < ed; k++) {
        if (k < reads[j]->size(reads[j])) {
          reads_slice[i][j]->push_back(
              reads_slice[i][j],
              reads[j]->at(reads[j], k));
        }
        else {
          reads_slice[i][j]->push_back(
              reads_slice[i][j],
              NUCLEOTIDE_N);
        }
      }
    }

    trs->set_reads(trs, reads_slice[i], nreads);
    Strand *partial_consensus = tr_system_get_partial_consensus(trs);
    for (int j = 0; j < partial_consensus->size(partial_consensus) - trs->trw_width; j++) {
      consensus->push_back(
          consensus,
          partial_consensus->at(partial_consensus, j));
    }
    // print_strand(partial_consensus);

    free_strand(partial_consensus);
    for (int j = 0; j < nreads; j++) {
      free_read(reads_slice[i][j]);
    }
    free(reads_slice[i]);
  }

  free(panchor);
  for (int i = 0; i < nreads; i++) {
    free_tr_read(reads[i]);
  }
  free(reads);
  free(reads_slice);

  return consensus;
}
#endif

#if TR_INDETERMINATED == 1 && TR_ALIGNMENT == 1
#define TR_READ_DELAY   10
#define TR_READ_SEARCH  20

Strand *tr_system_get_partial_consensus(TRSystem *trs) {
  assert(trs->nreads != 0);
  assert(trs->reads != NULL);

  Strand *consensus = new_strand();
  int *pcmp = (int *)malloc(trs->nreads * sizeof(int));
  int *delay = (int *)malloc(trs->nreads * sizeof(int));

  for (int i = 0; i < trs->nreads; i++) {
    trs->reads[i]->state = TR_READ_STATE_ACTIVE;
    pcmp[i] = 0;
    delay[i] = TR_READ_DELAY;
  }

  while (1) {
    /* initialize */
    int nreads_active = trs->nreads;
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state == TR_READ_STATE_OMITTED) {
        nreads_active--;
      }
      else if (pcmp[i] >= read->size(read)) {
        read->state = TR_READ_STATE_OMITTED;
        nreads_active--;
      }
      else if (read->state == TR_READ_STATE_INACTIVE) {
        nreads_active--;

        if (--(delay[i]) == 0) {
          int pcmp_next = 0;

          for (int j = pcmp[i] + TR_READ_DELAY; j <= pcmp[i] + TR_READ_SEARCH; j++) {
            if (j + 1 < read->size(read)) {
              int fault = 0, min_fault = TR_READ_DELAY + 1;

              for (int k = 0; k < TR_READ_DELAY; k++) {
                if (consensus->at(consensus, consensus->size(consensus) - 1 - k)
                    != read->at(read, j - k))
                {
                  fault++;
                }
              }
              if (fault < min_fault && fault <= TR_READ_DELAY / 5) {
                min_fault = fault;
                pcmp_next = j;
              }
            }
          }
          if (pcmp_next != 0) {
#if TR_VERBOSE
            printf("\nset active %d\n", i);
            printf("  pcmp:%d ", pcmp[i]);
#endif
            pcmp[i] = pcmp_next + 1;
            read->state = TR_READ_STATE_ACTIVE;
            nreads_active++;
#if TR_VERBOSE
            printf("%d\n", pcmp[i]);
            print_strand(consensus);
            print_tr_read(read);
#endif
          }
          else {
            read->state = TR_READ_STATE_OMITTED;
          }

          delay[i] = TR_READ_DELAY;
        }
      }
      else if (read->state != TR_READ_STATE_OMITTED) {
        read->state = TR_READ_STATE_ACTIVE;
      }
    }
    if (nreads_active == 0) {
      break;
    }

    /* core algorithm */

    // get single concensus
    Nucleotide single_consensus = NUCLEOTIDE_N;
    double weight[NUCLEOTIDE_SIZE] = {};
    double weight_max = 0;

    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state != TR_READ_STATE_OMITTED
          && read->state != TR_READ_STATE_INACTIVE)
      {
        Nucleotide n = read->at(read, pcmp[i]);

        weight[n]++;
      }
    }
    
    for (int i = 0; i < NUCLEOTIDE_SIZE; i++) {
      if (weight[i] > weight_max) {
        single_consensus = i;
        weight_max = weight[i];
      }
    }

    consensus->push_back(consensus, single_consensus);

    // set variant reads' state
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state != TR_READ_STATE_OMITTED
          && read->state != TR_READ_STATE_INACTIVE) {
        Nucleotide n = read->at(read, pcmp[i]);

        if (n != single_consensus) {
          read->state = TR_READ_STATE_VARIANT;
        }
      }
    }

    // get look-ahead window consensus
    Strand *win_consensus = trs->trw->get_consensus(trs->trw, trs->reads, pcmp);
    // get comparison consensus
    Strand *cmp_consensus = new_strand();

    cmp_consensus->push_back(cmp_consensus, single_consensus);
    for (int i = 0; i < win_consensus->size(win_consensus) - 1; i++) {
      cmp_consensus->push_back(
          cmp_consensus,
          win_consensus->at(win_consensus, i));
    }

    // move pcmp
    for (int i = 0; i < trs->nreads; i++) {
      TRRead *read = trs->reads[i];

      if (read->state == TR_READ_STATE_ACTIVE) {
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
        if (compare_strand(win_consensus, read_slice_win) == true) {
          pcmp[i]++;
        }
        // deletion
        else if (compare_strand(win_consensus, read_slice_cmp) == true) {
          // do nothing
        }
        // insertion
        else if (compare_strand(cmp_consensus, read_slice_win) == true){
          pcmp[i] += 2;
        }
        else {
          read->state = TR_READ_STATE_INACTIVE;

#if TR_VERBOSE
          printf("%d\n", pcmp[i]);
          print_tr_read(read);
          printf("- "),
            print_strand(consensus);
          printf("- win-c"),
            print_strand(win_consensus);
          printf("- cmp-c"),
            print_strand(win_consensus);
          printf("- win-r"),
            print_strand(read_slice_win);
          printf("- cmp-r"),
            print_strand(read_slice_cmp);
#endif
        }

        free_strand(read_slice_win);
        free_strand(read_slice_cmp);
      }
    }
    
    free_strand(win_consensus);
    free_strand(cmp_consensus);
  }

  free(pcmp);
  free(delay);

  return consensus;
}
Strand *tr_system_get_consensus(TRSystem *trs) {
  assert(trs->nreads != 0);
  assert(trs->reads != NULL);

  Strand *consensus = new_strand();
  int nreads = trs->nreads;
  int nanchor = trs->reads[0]->super->anchor->size(trs->reads[0]->super->anchor);
  int *panchor = (int *)malloc(sizeof(int) * nreads);
  TRRead **tr_reads = trs->reads;
  TRRead **reads = (TRRead **)malloc(sizeof(TRRead *) * trs->nreads);
  Read ***reads_slice = (Read ***)malloc(sizeof(Read **) * (nanchor + 1));

  for (int i = 0; i < nreads; i++) {
    reads[i] = new_tr_read();
    panchor[i] = -1;

    for (int j = 0; j < tr_reads[i]->size(tr_reads[i]); j++) {
      reads[i]->push_back(reads[i], tr_reads[i]->at(tr_reads[i], j));
    }

    vector_t *error = tr_reads[i]->super->error;
    for (int j = 0; j < error->size(error); j++) {
      reads[i]->super->error->push_back(
          reads[i]->super->error,
          error->at(error, j));
    }

    vector_t *anchor = tr_reads[i]->super->anchor;
    for (int j = 0; j < anchor->size(anchor); j++) {
      reads[i]->super->anchor->push_back(
          reads[i]->super->anchor,
          anchor->at(anchor, j));
    }
  }

  for (int i = 0; i <= nanchor; i++) {
    reads_slice[i] = (Read **)malloc(sizeof(Read *) * trs->nreads);

    for (int j = 0; j < nreads; j++) {
      reads_slice[i][j] = new_read();
      vector_t *anchor = reads[j]->super->anchor;
      int st = panchor[j] == -1 ? 0 : *(int *)anchor->at(anchor, panchor[j]);
      int ed = ++panchor[j] < nanchor ? *(int *)anchor->at(anchor, panchor[j]) : reads[j]->size(reads[j]);
      ed += trs->trw_width;

      // printf("%d %d\n", st, ed);
      // for (int i = 0; i < anchor->size(anchor); i++) {
      //   printf("%d ", *(int *)anchor->at(anchor, i));
      // }
      // printf("\n");
      for (int k = st; k < ed; k++) {
        if (k < reads[j]->size(reads[j])) {
          reads_slice[i][j]->push_back(
              reads_slice[i][j],
              reads[j]->at(reads[j], k));
        }
        else {
          reads_slice[i][j]->push_back(
              reads_slice[i][j],
              NUCLEOTIDE_N);
        }
      }
    }

    trs->set_reads(trs, reads_slice[i], nreads);
    Strand *partial_consensus = tr_system_get_partial_consensus(trs);
    for (int j = 0; j < partial_consensus->size(partial_consensus) - trs->trw_width; j++) {
      consensus->push_back(
          consensus,
          partial_consensus->at(partial_consensus, j));
    }
    // print_strand(partial_consensus);

    free_strand(partial_consensus);
    for (int j = 0; j < nreads; j++) {
      free_read(reads_slice[i][j]);
    }
    free(reads_slice[i]);
  }

  free(panchor);
  for (int i = 0; i < nreads; i++) {
    free_tr_read(reads[i]);
  }
  free(reads);
  free(reads_slice);

  return consensus;
}

#endif

void tr_system_set_window(TRSystem *trs) {
  free_tr_window(trs->trw);
  trs->trw = new_tr_window(trs->nreads, trs->trw_width);
}
